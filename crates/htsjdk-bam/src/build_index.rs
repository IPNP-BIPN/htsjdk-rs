//! `BuildBamIndex` / `BAMIndexer` fed from a BAM read off disk.
//!
//! Ports the read side of `htsjdk.samtools.BuildBamIndex`: read a coordinate-sorted BAM and build its
//! `.bai`. The index itself is [`BamIndexer`], already proven byte-identical to htsjdk's `BAMIndexer`
//! when driven by [`crate::writer::BamWriter`] (the `bai_conformance` goldens). The only new thing
//! here is feeding it from a **read** rather than a write: for each record we take the BGZF virtual
//! file pointer before and after its bytes, exactly as `BAMFileWriter` does, so the chunk spans the
//! record's own bytes and the `.bai` matches. Because a `.bai` is a property of the BAM and not of
//! how it was produced, `build_bam_index(bam)` reproduces the index `BamWriter::with_index` would
//! have written for that same BAM.

use std::io::Read;

use htsjdk_bgzf::{BgzfError, BgzfReader};

use crate::index::{BamIndexer, Chunk};
use crate::record::{BamRecord, DecodeError};
use crate::writer::BAM_MAGIC;

/// Why an index could not be built.
#[derive(Debug)]
pub enum BuildIndexError {
    NotABam,
    Bgzf(BgzfError),
    Truncated,
    Decode(DecodeError),
    /// `BAMIndexer.processAlignment` refused the record, which only a read past
    /// `BIN_GENOMIC_SPAN` makes it do: `htsjdk.samtools.SAMException` with this message, caused by
    /// `IllegalStateException: Read position too high for BAI bin indexing.` `BuildBamIndex` has
    /// already created the `.bai` and leaves it empty.
    Index(String),
}

impl From<BgzfError> for BuildIndexError {
    fn from(e: BgzfError) -> Self {
        BuildIndexError::Bgzf(e)
    }
}

impl From<DecodeError> for BuildIndexError {
    fn from(e: DecodeError) -> Self {
        BuildIndexError::Decode(e)
    }
}

fn read_exact_or_eof<R: Read>(r: &mut R, buf: &mut [u8]) -> Result<bool, BuildIndexError> {
    let mut filled = 0;
    while filled < buf.len() {
        match r.read(&mut buf[filled..]) {
            Ok(0) => {
                if filled == 0 {
                    return Ok(false); // clean EOF at a record boundary
                }
                return Err(BuildIndexError::Truncated);
            }
            Ok(n) => filled += n,
            Err(_) => return Err(BuildIndexError::Truncated),
        }
    }
    Ok(true)
}

fn read_i32<R: Read>(r: &mut R) -> Result<i32, BuildIndexError> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b)
        .map_err(|_| BuildIndexError::Truncated)?;
    Ok(i32::from_le_bytes(b))
}

/// Consume the BAM header binary content, returning the per-reference names and lengths (`l_ref`).
/// Leaves the reader positioned at the first record, so `virtual_pos()` is the first record's
/// virtual offset.
fn consume_header<R: Read>(
    reader: &mut BgzfReader<R>,
) -> Result<(Vec<String>, Vec<i32>), BuildIndexError> {
    let mut magic = [0u8; 4];
    reader
        .read_exact(&mut magic)
        .map_err(|_| BuildIndexError::Truncated)?;
    if magic != BAM_MAGIC {
        return Err(BuildIndexError::NotABam);
    }
    let l_text = read_i32(reader)?;
    std::io::copy(
        &mut reader.by_ref().take(l_text as u64),
        &mut std::io::sink(),
    )
    .map_err(|_| BuildIndexError::Truncated)?;
    let n_ref = read_i32(reader)?;
    let mut names = Vec::with_capacity(n_ref.max(0) as usize);
    let mut lengths = Vec::with_capacity(n_ref.max(0) as usize);
    for _ in 0..n_ref {
        let l_name = read_i32(reader)?;
        let mut name = Vec::new();
        reader
            .by_ref()
            .take(l_name as u64)
            .read_to_end(&mut name)
            .map_err(|_| BuildIndexError::Truncated)?;
        // `l_name` counts the NUL terminator.
        name.pop();
        names.push(String::from_utf8_lossy(&name).into_owned());
        lengths.push(read_i32(reader)?);
    }
    Ok((names, lengths))
}

/// `BuildBamIndex`: read a BAM (BGZF-framed) and return its `.bai` bytes.
pub fn build_bam_index(bam: &[u8]) -> Result<Vec<u8>, BuildIndexError> {
    let mut reader = BgzfReader::new(bam);
    let (reference_names, reference_lengths) = consume_header(&mut reader)?;
    let mut indexer = BamIndexer::new(reference_lengths);

    loop {
        // htsjdk takes the pointer before the record's bytes (the block_size int included) and again
        // after, so the chunk spans exactly this record.
        let start_offset = reader.virtual_pos();
        let mut size_bytes = [0u8; 4];
        if !read_exact_or_eof(&mut reader, &mut size_bytes)? {
            break;
        }
        let block_size = i32::from_le_bytes(size_bytes);
        let mut payload = vec![0u8; block_size as usize];
        reader
            .read_exact(&mut payload)
            .map_err(|_| BuildIndexError::Truncated)?;
        let stop_offset = reader.virtual_pos();

        // Decode the record from its full bytes (the length prefix followed by the payload).
        let mut full = Vec::with_capacity(4 + payload.len());
        full.extend_from_slice(&size_bytes);
        full.extend_from_slice(&payload);
        let (record, _) = BamRecord::decode(&full)?.ok_or(BuildIndexError::Truncated)?;

        indexer
            .process(
                record.reference_index,
                record.alignment_start,
                record.alignment_end(),
                record.read_unmapped(),
                Chunk {
                    start: start_offset,
                    end: stop_offset,
                },
            )
            .map_err(|_| {
                BuildIndexError::Index(format!(
                    "Exception creating BAM index for record {}",
                    crate::text::sam_string(&record, &reference_names)
                ))
            })?;
    }

    Ok(indexer.finish())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cigar::{Cigar, CigarElement, Op};
    use crate::header::{SamHeader, SequenceRecord};
    use crate::record::BamRecord;
    use crate::writer::BamWriter;

    fn bam_with(records: &[BamRecord]) -> Vec<u8> {
        let mut h = SamHeader::new();
        h.set_sort_order("coordinate");
        h.sequences.push(SequenceRecord::new("chr1", 100_000));
        let mut w = BamWriter::new(Vec::new(), &h).unwrap();
        for r in records {
            w.write(r).unwrap();
        }
        w.finish().unwrap()
    }

    fn mapped(name: &str, start: i32) -> BamRecord {
        BamRecord {
            read_name: name.to_string(),
            reference_index: 0,
            alignment_start: start,
            mapping_quality: 60,
            cigar: Cigar::new(vec![CigarElement {
                length: 10,
                op: Op::M,
            }]),
            flags: 0,
            ..BamRecord::default()
        }
    }

    #[test]
    fn produces_a_bai_with_the_right_magic_and_reference_count() {
        let bam = bam_with(&[mapped("a", 10), mapped("b", 20)]);
        let bai = build_bam_index(&bam).unwrap();
        assert_eq!(&bai[..4], b"BAI\x01");
        // n_ref = 1.
        assert_eq!(i32::from_le_bytes(bai[4..8].try_into().unwrap()), 1);
    }

    #[test]
    fn a_non_bam_input_is_rejected() {
        let mut w = htsjdk_bgzf::BgzfWriter::new(Vec::new());
        std::io::Write::write_all(&mut w, b"XXXXnot bam").unwrap();
        let framed = w.into_inner().unwrap();
        assert!(matches!(
            build_bam_index(&framed),
            Err(BuildIndexError::NotABam)
        ));
    }

    fn on_a_long_reference(starts: &[i32]) -> Vec<u8> {
        let mut h = SamHeader::new();
        h.set_sort_order("coordinate");
        h.sequences.push(SequenceRecord::new("big", 600_000_000));
        let mut w = BamWriter::new(Vec::new(), &h).unwrap();
        for (i, &start) in starts.iter().enumerate() {
            let mut r = mapped(&format!("r{i}"), start);
            r.read_bases = b"ACGTACGTAC".to_vec();
            r.base_qualities = vec![40; 10];
            w.write(&r).unwrap();
        }
        w.finish().unwrap()
    }

    /// `BuildBamIndex` on a read past 2^29 used to panic in the linear index. htsjdk 5.0.0 throws
    /// this `SAMException`, measured with the GATK 4.7.0.0 jar.
    #[test]
    fn a_read_past_the_bin_span_is_refused_not_a_panic() {
        let bam = on_a_long_reference(&[100, 550_000_000]);
        match build_bam_index(&bam) {
            Err(BuildIndexError::Index(message)) => assert_eq!(
                message,
                "Exception creating BAM index for record \
                 r1\t0\tbig\t550000000\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII"
            ),
            other => panic!("{other:?}"),
        }
    }

    /// The bin field of a read on an over-long reference is 0; the index bins it by position.
    #[test]
    fn a_read_on_an_over_long_reference_is_binned_by_its_position() {
        let bai = build_bam_index(&on_a_long_reference(&[100])).unwrap();
        let index = crate::index::read_bai(&bai).unwrap();
        assert_eq!(index.references[0].bins[0].bin_number, 4681);
    }
}
