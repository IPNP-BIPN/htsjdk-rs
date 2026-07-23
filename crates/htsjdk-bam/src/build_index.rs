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

use crate::bin::{compute_indexing_bin, BIN_GENOMIC_SPAN, NO_ALIGNMENT_START};
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

/// Consume the BAM header binary content, returning the per-reference lengths (`l_ref`). Leaves the
/// reader positioned at the first record, so `virtual_pos()` is the first record's virtual offset.
fn consume_header<R: Read>(reader: &mut BgzfReader<R>) -> Result<Vec<i32>, BuildIndexError> {
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
    let mut lengths = Vec::with_capacity(n_ref.max(0) as usize);
    for _ in 0..n_ref {
        let l_name = read_i32(reader)?;
        std::io::copy(
            &mut reader.by_ref().take(l_name as u64),
            &mut std::io::sink(),
        )
        .map_err(|_| BuildIndexError::Truncated)?;
        lengths.push(read_i32(reader)?);
    }
    Ok(lengths)
}

/// `BuildBamIndex`: read a BAM (BGZF-framed) and return its `.bai` bytes.
pub fn build_bam_index(bam: &[u8]) -> Result<Vec<u8>, BuildIndexError> {
    let mut reader = BgzfReader::new(bam);
    let reference_lengths = consume_header(&mut reader)?;
    let mut indexer = BamIndexer::new(reference_lengths.clone());

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

        let index_bin = index_bin_for(&record, &reference_lengths);
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

    Ok(indexer.finish())
}

/// The indexing bin `BAMFileWriter` would record: 0 for a reference too large for the 6-level scheme
/// or for an unplaced read, else `computeIndexingBin`.
fn index_bin_for(record: &BamRecord, reference_lengths: &[i32]) -> i32 {
    let too_large = record.reference_index >= 0
        && reference_lengths
            .get(record.reference_index as usize)
            .is_some_and(|&len| len > BIN_GENOMIC_SPAN);
    if too_large || record.alignment_start == NO_ALIGNMENT_START {
        0
    } else {
        compute_indexing_bin(record.alignment_start, record.alignment_end()).unwrap_or(0)
    }
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
}
