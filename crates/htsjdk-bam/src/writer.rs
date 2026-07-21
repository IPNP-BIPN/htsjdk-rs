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

use htsjdk_bgzf::BgzfWriter;

use crate::bin::BIN_GENOMIC_SPAN;
use crate::header::SamHeader;
use crate::record::{BamRecord, EncodeError};

/// `BAMFileConstants.BAM_MAGIC`.
pub const BAM_MAGIC: [u8; 4] = *b"BAM\x01";

/// Writes a BAM file: header, dictionary, records, terminator.
pub struct BamWriter<W: Write> {
    bgzf: BgzfWriter<W>,
    /// Reference lengths, kept for the too-large-reference bin rule.
    reference_lengths: Vec<i32>,
}

impl<W: Write> BamWriter<W> {
    /// `BAMFileWriter.writeHeader`: magic, header text, then the binary dictionary.
    pub fn new(inner: W, header: &SamHeader) -> io::Result<Self> {
        let mut bgzf = BgzfWriter::new(inner);
        let text = header.encode();

        bgzf.write_all(&BAM_MAGIC)?;

        // `writeString(headerText, true, false)`: length prefix, no null terminator. The
        // length counts UTF-16 units, as everywhere else in htsjdk.
        let text_bytes: Vec<u8> = text.encode_utf16().map(|u| (u & 0xFF) as u8).collect();
        bgzf.write_all(&(text_bytes.len() as i32).to_le_bytes())?;
        bgzf.write_all(&text_bytes)?;

        // The dictionary again, in binary. Redundant with the text, and written anyway.
        bgzf.write_all(&(header.sequences.len() as i32).to_le_bytes())?;
        for seq in &header.sequences {
            let name: Vec<u8> = seq.name.encode_utf16().map(|u| (u & 0xFF) as u8).collect();
            // `writeString(name, true, true)`: the length here DOES include the terminator.
            bgzf.write_all(&((name.len() + 1) as i32).to_le_bytes())?;
            bgzf.write_all(&name)?;
            bgzf.write_all(&[0])?;
            bgzf.write_all(&seq.length.to_le_bytes())?;
        }

        Ok(BamWriter {
            bgzf,
            reference_lengths: header.sequences.iter().map(|s| s.length).collect(),
        })
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
        let bytes = if self.reference_too_large_for_bin(record.reference_index) {
            record.encode_with_bin(0)
        } else {
            record.encode()
        }
        .map_err(WriteError::Encode)?;
        self.bgzf.write_all(&bytes).map_err(WriteError::Io)
    }

    /// Closes the BGZF stream, emitting the empty terminator block.
    pub fn finish(self) -> io::Result<W> {
        self.bgzf.into_inner()
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
}
