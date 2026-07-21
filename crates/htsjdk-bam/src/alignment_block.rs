//! Alignment blocks: the contiguous read-to-reference runs of a CIGAR.
//!
//! Ported from `htsjdk.samtools.SAMUtils.getAlignmentBlocks` and `htsjdk.samtools.AlignmentBlock`
//! at htsjdk 4.2.0.
//!
//! Both coordinates are **1-based**, matching the Java, because callers index the reference with
//! `getReferenceStart() - 1` and the read with `getReadStart() - 1` and any port that quietly
//! rebased them to 0 would put every off-by-one in the caller instead of here.
//!
//! One property that looks like an implementation detail and is not: **consecutive `M` operators
//! produce separate blocks**. The loop appends a block per element and never merges, so `10M10M`
//! gives two blocks covering the same span as `20M` gives one. Downstream code that indexes by
//! position within a block therefore sees different offsets for the same read, which is exactly
//! how Picard's `BAD_CYCLES` divergence surfaces.

use crate::cigar::{Cigar, Op};

/// One contiguous run of aligned bases.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AlignmentBlock {
    /// 1-based offset of the first base of this block in the read.
    pub read_start: i32,
    /// 1-based position of the first base of this block on the reference.
    pub reference_start: i32,
    pub length: i32,
}

/// `SAMUtils.getAlignmentBlocks(cigar, alignmentStart, ...)`.
///
/// `H` and `P` advance neither coordinate; `S` and `I` advance only the read; `N` and `D` advance
/// only the reference. `M`, `=` and `X` emit a block.
pub fn alignment_blocks(cigar: &Cigar, alignment_start: i32) -> Vec<AlignmentBlock> {
    let mut blocks = Vec::new();
    let mut read_base = 1;
    let mut ref_base = alignment_start;

    for e in &cigar.elements {
        let len = e.length as i32;
        match e.op {
            Op::H | Op::P => {}
            Op::S => read_base += len,
            Op::N | Op::D => ref_base += len,
            Op::I => read_base += len,
            Op::M | Op::Eq | Op::X => {
                blocks.push(AlignmentBlock {
                    read_start: read_base,
                    reference_start: ref_base,
                    length: len,
                });
                read_base += len;
                ref_base += len;
            }
        }
    }
    blocks
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cigar::CigarElement;

    fn cigar(spec: &[(u32, Op)]) -> Cigar {
        Cigar::new(
            spec.iter()
                .map(|&(length, op)| CigarElement { length, op })
                .collect(),
        )
    }

    #[test]
    fn a_simple_match_is_one_block() {
        let b = alignment_blocks(&cigar(&[(20, Op::M)]), 100);
        assert_eq!(
            b,
            [AlignmentBlock {
                read_start: 1,
                reference_start: 100,
                length: 20
            }]
        );
    }

    /// A deletion advances the reference and not the read, so the second block starts at read
    /// base 11 and reference base 115.
    #[test]
    fn a_deletion_splits_the_read_into_two_blocks() {
        let b = alignment_blocks(&cigar(&[(10, Op::M), (5, Op::D), (10, Op::M)]), 100);
        assert_eq!(
            b,
            [
                AlignmentBlock {
                    read_start: 1,
                    reference_start: 100,
                    length: 10
                },
                AlignmentBlock {
                    read_start: 11,
                    reference_start: 115,
                    length: 10
                },
            ]
        );
    }

    /// An insertion advances the read and not the reference.
    #[test]
    fn an_insertion_advances_only_the_read() {
        let b = alignment_blocks(&cigar(&[(10, Op::M), (3, Op::I), (10, Op::M)]), 100);
        assert_eq!(b[1].read_start, 14);
        assert_eq!(b[1].reference_start, 110);
    }

    #[test]
    fn soft_clips_advance_the_read_and_hard_clips_advance_nothing() {
        let b = alignment_blocks(&cigar(&[(5, Op::S), (10, Op::M)]), 100);
        assert_eq!(b[0].read_start, 6);
        assert_eq!(b[0].reference_start, 100);

        let b = alignment_blocks(&cigar(&[(5, Op::H), (10, Op::M)]), 100);
        assert_eq!(b[0].read_start, 1);
        assert_eq!(b[0].reference_start, 100);
    }

    /// The loop never merges, so adjacent match operators stay separate blocks. Any code that
    /// works in per-block offsets sees a different picture for `10M10M` than for `20M`, on reads
    /// that align identically.
    #[test]
    fn adjacent_match_operators_are_not_merged() {
        let split = alignment_blocks(&cigar(&[(10, Op::M), (10, Op::M)]), 100);
        let joined = alignment_blocks(&cigar(&[(20, Op::M)]), 100);
        assert_eq!(split.len(), 2);
        assert_eq!(joined.len(), 1);
        assert_eq!(split[1].read_start, 11);
        assert_eq!(split[1].reference_start, 110);
    }

    #[test]
    fn eq_and_x_emit_blocks_like_m() {
        let b = alignment_blocks(&cigar(&[(5, Op::Eq), (2, Op::X), (5, Op::Eq)]), 1);
        assert_eq!(b.len(), 3);
        assert_eq!(b[2].read_start, 8);
        assert_eq!(b[2].reference_start, 8);
    }

    #[test]
    fn an_unmapped_cigar_has_no_blocks() {
        assert!(alignment_blocks(&cigar(&[]), 0).is_empty());
    }
}
