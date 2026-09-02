//! The BAM file writer.
//!
//! Ported from `htsjdk.samtools.BAMFileWriter.writeHeader` and `writeAlignment`, on top of
//! `BlockCompressedOutputStream`.
//!
//! A BAM file is the SAM text header and the sequence dictionary, both written *inside* the
//! BGZF stream, followed by the records and the empty terminator block. Two details in the
//! framing are easy to get subtly wrong, and neither shows up when reading the result back:
//!
//! - the header text length is written **without** a null terminator, while each sequence name
//!   length is written **with** one, so the same `writeString` helper is called two different
//!   ways four lines apart;
//! - the sequence dictionary is redundant with the text header, and htsjdk writes both, so a
//!   writer that emits only one produces a shorter, still-readable file.

use std::io::{self, Write};

use htsjdk_bgzf::{BgzfWriter, Deflater, DEFAULT_COMPRESSION_LEVEL};

use crate::bin::BIN_GENOMIC_SPAN;
use crate::header::SamHeader;
use crate::index::{BamIndexer, Chunk};
use crate::record::{BamRecord, EncodeError};

/// `BAMFileConstants.BAM_MAGIC`.
pub const BAM_MAGIC: [u8; 4] = *b"BAM\x01";

/// Writes a BAM file: header, dictionary, records, terminator.
pub struct BamWriter<W: Write> {
    bgzf: BgzfWriter<W>,
    /// Reference lengths, kept for the too-large-reference bin rule.
    reference_lengths: Vec<i32>,
    /// Present when an index is being built alongside the file.
    indexer: Option<BamIndexer>,
}

/// `BAMFileWriter.writeHeader(BinaryCodec, SAMFileHeader)`: the BAM header binary content, written
/// to any sink: `BAM\1` magic, the length-prefixed header text, then the binary sequence dictionary.
/// Shared by [`BamWriter::new`] and [`write_bam_header_block`] so the two cannot drift.
pub(crate) fn write_header_binary<W: Write>(
    w: &mut W,
    header: &SamHeader,
    keep_existing_version_number: bool,
) -> io::Result<()> {
    let text = if keep_existing_version_number {
        header.encode()
    } else {
        header.encode_replacing_version()
    };

    w.write_all(&BAM_MAGIC)?;

    // `writeString(headerText, true, false)`: length prefix, no null terminator. The
    // length counts UTF-16 units, as everywhere else in htsjdk.
    let text_bytes: Vec<u8> = text.encode_utf16().map(|u| (u & 0xFF) as u8).collect();
    w.write_all(&(text_bytes.len() as i32).to_le_bytes())?;
    w.write_all(&text_bytes)?;

    // The dictionary again, in binary. Redundant with the text, and written anyway.
    w.write_all(&(header.sequences.len() as i32).to_le_bytes())?;
    for seq in &header.sequences {
        let name: Vec<u8> = seq.name.encode_utf16().map(|u| (u & 0xFF) as u8).collect();
        // `writeString(name, true, true)`: the length here DOES include the terminator.
        w.write_all(&((name.len() + 1) as i32).to_le_bytes())?;
        w.write_all(&name)?;
        w.write_all(&[0])?;
        w.write_all(&seq.length.to_le_bytes())?;
    }
    Ok(())
}

/// `BAMFileWriter.writeHeader(OutputStream, SAMFileHeader)`: the header binary content written to a
/// fresh `BlockCompressedOutputStream` and `flush()`ed, so the result is complete BGZF block(s) with
/// **no** EOF terminator (htsjdk never closes that stream). This is the leading segment of a
/// block-copy reheader; the copied data blocks and the terminator follow it.
pub fn write_bam_header_block(header: &SamHeader) -> io::Result<Vec<u8>> {
    let mut bgzf = BgzfWriter::new(Vec::new());
    // The reheader path is the one that passes true: it is copying an existing file's blocks, and
    // its whole purpose is to leave everything but the header alone.
    write_header_binary(&mut bgzf, header, true)?;
    bgzf.into_inner_without_terminator()
}

/// `SAMFileWriterImpl.setHeader`'s `header.setSortOrder(this.sortOrder)`.
///
/// The writer's own sort order comes from `SAMFileWriterFactory.makeBAMWriter`, which reads it off
/// the header with `getSortOrder()` and hands it straight back. So a header carrying a recognised
/// `SO` keeps it, and one carrying none gets `unsorted` -- which is why a BAM written from a header
/// with no `@HD SO` still has one.
///
/// `getSortOrder()` answers `unsorted` for a value it does not recognise as well, so an `SO` of
/// something else would be replaced rather than kept. No golden reaches that, and the branch is
/// written to match the reference rather than to be convenient.
fn stamp_sort_order(header: &SamHeader) -> SamHeader {
    const RECOGNISED: [&str; 5] = [
        "unsorted",
        "queryname",
        "coordinate",
        "duplicate",
        "unknown",
    ];
    let mut stamped = header.clone();
    let order = match header.attributes.get("SO") {
        Some(value) if RECOGNISED.contains(&value) => value.to_string(),
        _ => "unsorted".to_string(),
    };
    stamped.set_sort_order(&order);
    stamped
}

impl<W: Write> BamWriter<W> {
    /// `BAMFileWriter.writeHeader`: magic, header text, then the binary dictionary.
    ///
    /// The BGZF underneath takes htsjdk's own defaults: the JDK deflater at level five. A caller
    /// running under GATK writes the same records through [`Self::with_compression`], because
    /// GATK replaces the static deflater factory and `GATKConfig` sets the level to two.
    pub fn new(inner: W, header: &SamHeader) -> io::Result<Self> {
        Self::with_compression(inner, header, DEFAULT_COMPRESSION_LEVEL, Deflater::Jdk)
    }

    /// The same writer with the BGZF compression named, which is the only thing that differs
    /// between a file htsjdk writes and the same file written by GATK.
    pub fn with_compression(
        inner: W,
        header: &SamHeader,
        level: u32,
        deflater: Deflater,
    ) -> io::Result<Self> {
        let mut bgzf = BgzfWriter::with_deflater(inner, level, deflater);
        // `SAMFileWriterImpl.setHeader` stamps the sort order onto the header before writing it,
        // defaulting to `unsorted`, and only then encodes -- with the version REPLACED, which is
        // the half `write_bam_header_block` does not do. See the `bam-header-version` golden.
        let stamped = stamp_sort_order(header);
        write_header_binary(&mut bgzf, &stamped, false)?;

        Ok(BamWriter {
            bgzf,
            reference_lengths: header.sequences.iter().map(|s| s.length).collect(),
            indexer: None,
        })
    }

    /// Builds a BAI index alongside the file, as `SAMFileWriterFactory.setCreateIndex(true)`
    /// does.
    ///
    /// This must be enabled before the first record: the index records the virtual file
    /// pointer around every record, and there is no way to recover the ones already written.
    pub fn with_index(mut self) -> Self {
        self.indexer = Some(BamIndexer::new(self.reference_lengths.clone()));
        self
    }

    /// Whether this reference is too long for the BAI bin field.
    ///
    /// `BAMRecordCodec.warnIfReferenceIsTooLargeForBinField` forces the bin to 0 for these,
    /// after warning once. It lives here rather than on the record because it needs the
    /// sequence dictionary, which the record does not carry.
    fn reference_too_large_for_bin(&self, reference_index: i32) -> bool {
        reference_index >= 0
            && self
                .reference_lengths
                .get(reference_index as usize)
                .is_some_and(|&len| len > BIN_GENOMIC_SPAN)
    }

    /// `BAMFileWriter.writeAlignment`.
    pub fn write(&mut self, record: &BamRecord) -> Result<(), WriteError> {
        let forced_bin = self.reference_too_large_for_bin(record.reference_index);
        let bytes = if forced_bin {
            record.encode_with_bin(0)
        } else {
            record.encode()
        }
        .map_err(WriteError::Encode)?;

        // htsjdk takes the pointer *before* encoding and again after, so the chunk spans
        // exactly this record's bytes. Taking it after the write for the start, or including
        // the next record, shifts every chunk in the index.
        let start_offset = self.bgzf.file_pointer();
        self.bgzf.write_all(&bytes).map_err(WriteError::Io)?;
        let stop_offset = self.bgzf.file_pointer();

        if let Some(indexer) = &mut self.indexer {
            let index_bin = if forced_bin {
                0
            } else if record.alignment_start != crate::bin::NO_ALIGNMENT_START {
                crate::bin::compute_indexing_bin(record.alignment_start, record.alignment_end())
                    .unwrap_or(0)
            } else {
                0
            };
            indexer.process(
                record.reference_index,
                record.alignment_start,
                record.alignment_end(),
                index_bin,
                record.read_unmapped(),
                Chunk {
                    start: start_offset,
                    end: stop_offset,
                },
            );
        }
        Ok(())
    }

    /// Closes the BGZF stream, emitting the empty terminator block.
    pub fn finish(self) -> io::Result<W> {
        self.bgzf.into_inner()
    }

    /// Closes the stream and returns the file alongside its BAI index.
    ///
    /// Panics unless [`Self::with_index`] was called.
    pub fn finish_with_index(self) -> io::Result<(W, Vec<u8>)> {
        let indexer = self.indexer.expect("with_index was not enabled");
        let index = indexer.finish();
        Ok((self.bgzf.into_inner()?, index))
    }
}

/// Why a record could not be written.
#[derive(Debug)]
pub enum WriteError {
    Encode(EncodeError),
    Io(io::Error),
}

impl std::fmt::Display for WriteError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            WriteError::Encode(e) => write!(f, "cannot encode record: {e:?}"),
            WriteError::Io(e) => write!(f, "{e}"),
        }
    }
}

impl std::error::Error for WriteError {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::header::{SamHeader, SequenceRecord};

    fn header() -> SamHeader {
        let mut h = SamHeader::new();
        h.set_sort_order("coordinate");
        h.sequences.push(SequenceRecord::new("chr1", 250_000_000));
        h
    }

    /// The uncompressed prefix, checked field by field before any record is written.
    #[test]
    fn the_file_opens_with_the_magic_then_the_header_text() {
        let h = header();
        let w = BamWriter::new(Vec::new(), &h).unwrap();
        let bytes = w.finish().unwrap();
        let plain = htsjdk_bgzf::decompress_all(&bytes).unwrap();

        assert_eq!(&plain[0..4], b"BAM\x01");
        let text_len = i32::from_le_bytes(plain[4..8].try_into().unwrap()) as usize;
        let text = &plain[8..8 + text_len];
        assert_eq!(text, h.encode().as_bytes());
        assert_ne!(
            text.last(),
            Some(&0),
            "the header text is NOT null terminated; only sequence names are"
        );

        let mut p = 8 + text_len;
        assert_eq!(
            i32::from_le_bytes(plain[p..p + 4].try_into().unwrap()),
            1,
            "the dictionary is written again in binary, redundantly with the text"
        );
        p += 4;
        let name_len = i32::from_le_bytes(plain[p..p + 4].try_into().unwrap()) as usize;
        assert_eq!(
            name_len, 5,
            "sequence name length DOES include its terminator"
        );
        p += 4;
        assert_eq!(&plain[p..p + 4], b"chr1");
        assert_eq!(plain[p + 4], 0);
        p += name_len;
        assert_eq!(
            i32::from_le_bytes(plain[p..p + 4].try_into().unwrap()),
            250_000_000
        );
    }

    #[test]
    fn an_empty_file_still_ends_with_the_terminator_block() {
        let w = BamWriter::new(Vec::new(), &header()).unwrap();
        let bytes = w.finish().unwrap();
        assert!(
            bytes.ends_with(&htsjdk_bgzf::EMPTY_GZIP_BLOCK),
            "a BAM without the EOF block is a truncated BAM to every reader that checks"
        );
    }

    /// A reference longer than the bin scheme can address forces the bin to 0. The record
    /// itself cannot know this, so the writer applies it.
    #[test]
    fn an_over_long_reference_forces_the_bin_to_zero() {
        let mut h = SamHeader::new();
        h.sequences
            .push(SequenceRecord::new("big", BIN_GENOMIC_SPAN + 1));
        let mut rec = BamRecord {
            read_name: "r".into(),
            reference_index: 0,
            alignment_start: 100,
            cigar: crate::cigar::Cigar::new(vec![crate::cigar::CigarElement {
                length: 4,
                op: crate::cigar::Op::M,
            }]),
            read_bases: b"ACGT".to_vec(),
            base_qualities: vec![30; 4],
            ..Default::default()
        };
        rec.mapping_quality = 60;

        let mut w = BamWriter::new(Vec::new(), &h).unwrap();
        w.write(&rec).unwrap();
        let plain = htsjdk_bgzf::decompress_all(&w.finish().unwrap()).unwrap();

        // Find the record: it starts after magic + text + dictionary.
        let text_len = i32::from_le_bytes(plain[4..8].try_into().unwrap()) as usize;
        let mut p = 8 + text_len + 4;
        let name_len = i32::from_le_bytes(plain[p..p + 4].try_into().unwrap()) as usize;
        p += 4 + name_len + 4;

        let bin = u16::from_le_bytes(plain[p + 14..p + 16].try_into().unwrap());
        assert_eq!(bin, 0, "a reference past BIN_GENOMIC_SPAN gets bin 0");

        // And on a normal-length reference the same record keeps its computed bin.
        let mut small = SamHeader::new();
        small
            .sequences
            .push(SequenceRecord::new("small", BIN_GENOMIC_SPAN));
        let mut w2 = BamWriter::new(Vec::new(), &small).unwrap();
        w2.write(&rec).unwrap();
        let plain2 = htsjdk_bgzf::decompress_all(&w2.finish().unwrap()).unwrap();
        let text_len2 = i32::from_le_bytes(plain2[4..8].try_into().unwrap()) as usize;
        let mut q = 8 + text_len2 + 4;
        let name_len2 = i32::from_le_bytes(plain2[q..q + 4].try_into().unwrap()) as usize;
        q += 4 + name_len2 + 4;
        assert_ne!(
            u16::from_le_bytes(plain2[q + 14..q + 16].try_into().unwrap()),
            0
        );
    }

    /// The compression is the only difference between a file htsjdk writes and GATK's own.
    ///
    /// The two deflaters' agreement with the reference is measured where each has an oracle:
    /// `gkl-deflate`'s suites and `gatk-rs`'s goldens. What is asserted here is that this writer
    /// hands the choice through and changes nothing else about the file -- the header, the records
    /// and the terminator are the format's.
    #[test]
    fn the_compression_reaches_the_blocks_and_nothing_else() {
        let header = header();
        let record = BamRecord {
            read_name: "r1".to_string(),
            reference_index: 0,
            alignment_start: 100,
            ..BamRecord::default()
        };

        let mut jdk = BamWriter::new(Vec::new(), &header).unwrap();
        jdk.write(&record).unwrap();
        let jdk = jdk.finish().unwrap();

        // Level TWO with GKL is the pair a real `gatk` run writes, and it is igzip's, so a build
        // without ISA-L refuses rather than answering with zlib's bytes.
        let mut gatk = BamWriter::with_compression(Vec::new(), &header, 5, Deflater::Gkl).unwrap();
        gatk.write(&record).unwrap();
        let gatk = gatk.finish().unwrap();

        // The two files decompress to the same bytes and are not the same file.
        let (a, b) = (
            htsjdk_bgzf::decompress_all(&jdk).unwrap(),
            htsjdk_bgzf::decompress_all(&gatk).unwrap(),
        );
        assert_eq!(a, b, "the records are the records whatever compressed them");
        assert!(jdk.starts_with(&[0x1f, 0x8b]) && gatk.starts_with(&[0x1f, 0x8b]));
    }
}
