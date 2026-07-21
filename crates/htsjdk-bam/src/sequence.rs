//! IUPAC base comparison.
//!
//! Ported from `htsjdk.samtools.util.SequenceUtil.basesEqual`, `isNoCall` and the static `bases`
//! table at htsjdk 4.2.0.
//!
//! The table maps each byte to a 4-bit mask of the nucleotides it stands for, and `basesEqual`
//! compares the **masks**, not the bytes. Three consequences, none of them obvious from the name:
//!
//!  1. Comparison is **case-insensitive**, because the static block copies every entry from `A`
//!     to `Z` into its lowercase counterpart. `a` equals `A`.
//!  2. `.` is given the same mask as `N`, so `.` equals `N` and equals `n`.
//!  3. Every byte *not* in the table keeps mask 0, and 0 equals 0. So two bytes that are not
//!     bases at all - `*` against `#`, say - compare **equal**, while either against a real base
//!     compares unequal. A reimplementation that compared uppercased bytes would agree on every
//!     valid input and disagree here.
//!
//! It is set equality and not set overlap: `M` is `A|C` and `A` is `A`, and those masks differ,
//! so `basesEqual(A, M)` is false. `basesEqualAmbiguous` is the overlap test, and Picard's
//! alignment metrics do not use it.

const A_MASK: u8 = 1;
const C_MASK: u8 = 2;
const G_MASK: u8 = 4;
const T_MASK: u8 = 8;

/// `SequenceUtil.BASES_ARRAY_LENGTH`.
const BASES_ARRAY_LENGTH: usize = 127;

/// The static `bases` table, built exactly as htsjdk's static block builds it.
fn bases_table() -> [u8; BASES_ARRAY_LENGTH] {
    let mut t = [0u8; BASES_ARRAY_LENGTH];
    t[b'A' as usize] = A_MASK;
    t[b'C' as usize] = C_MASK;
    t[b'G' as usize] = G_MASK;
    t[b'T' as usize] = T_MASK;
    t[b'M' as usize] = A_MASK | C_MASK;
    t[b'R' as usize] = A_MASK | G_MASK;
    t[b'W' as usize] = A_MASK | T_MASK;
    t[b'S' as usize] = C_MASK | G_MASK;
    t[b'Y' as usize] = C_MASK | T_MASK;
    t[b'K' as usize] = G_MASK | T_MASK;
    t[b'V' as usize] = A_MASK | C_MASK | G_MASK;
    t[b'H' as usize] = A_MASK | C_MASK | T_MASK;
    t[b'D' as usize] = A_MASK | G_MASK | T_MASK;
    t[b'B' as usize] = C_MASK | G_MASK | T_MASK;
    t[b'N' as usize] = A_MASK | C_MASK | G_MASK | T_MASK;
    // The lowercase copy runs over the whole A..Z range, so it also copies the zeros of the
    // letters that are not IUPAC codes. Then '.' is set, *after* the loop, so it has no
    // lowercase twin - which does not matter, since '.' has no case.
    let mut i = b'A';
    while i <= b'Z' {
        t[(i + 32) as usize] = t[i as usize];
        i += 1;
    }
    t[b'.' as usize] = A_MASK | C_MASK | G_MASK | T_MASK;
    t
}

static BASES: std::sync::LazyLock<[u8; BASES_ARRAY_LENGTH]> = std::sync::LazyLock::new(bases_table);

/// `SequenceUtil.basesEqual(lhs, rhs)`.
///
/// The bounds check is htsjdk's: a byte outside the table's range is **unequal to everything**,
/// including to another out-of-range byte. That is the one case where two unknown bytes do not
/// compare equal, and it turns on the table's length rather than on anything about bases.
pub fn bases_equal(lhs: u8, rhs: u8) -> bool {
    if lhs as usize >= BASES_ARRAY_LENGTH || rhs as usize >= BASES_ARRAY_LENGTH {
        return false;
    }
    BASES[lhs as usize] == BASES[rhs as usize]
}

/// `SequenceUtil.isNoCall(base)`.
///
/// Note this is a plain byte test and not a table lookup, so it is *not* the same predicate as
/// "has the N mask". `.` and `N` and `n` are no-calls; `B`, whose mask is also ambiguous, is not.
pub fn is_no_call(base: u8) -> bool {
    base == b'N' || base == b'n' || base == b'.'
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn comparison_ignores_case() {
        assert!(bases_equal(b'a', b'A'));
        assert!(bases_equal(b'g', b'G'));
        assert!(!bases_equal(b'a', b'g'));
    }

    #[test]
    fn dot_is_the_same_as_n() {
        assert!(bases_equal(b'.', b'N'));
        assert!(bases_equal(b'.', b'n'));
        assert!(!bases_equal(b'.', b'A'));
    }

    /// The mask test is equality, not overlap, so an ambiguity code does not match the bases it
    /// stands for.
    #[test]
    fn an_ambiguity_code_does_not_match_its_members() {
        assert!(!bases_equal(b'M', b'A'), "M is A|C, which is not A");
        assert!(!bases_equal(b'N', b'A'));
        assert!(bases_equal(b'M', b'm'));
    }

    /// Two bytes that are not bases share mask 0 and therefore compare equal. This is the case a
    /// byte-comparing reimplementation gets wrong in the direction of being *more* correct.
    #[test]
    fn two_non_bases_compare_equal_to_each_other() {
        assert!(bases_equal(b'*', b'#'));
        assert!(bases_equal(b'Z', b'Q'), "neither letter is an IUPAC code");
        assert!(!bases_equal(b'*', b'A'));
    }

    /// ...unless one of them is outside the table, where the bounds check rejects first.
    #[test]
    fn a_byte_past_the_table_is_equal_to_nothing() {
        assert!(!bases_equal(200, 200), "the same byte, still unequal");
        assert!(!bases_equal(200, b'A'));
        assert!(bases_equal(126, 125), "both inside the table, both mask 0");
    }

    /// `isNoCall` is a byte test, not a mask test, so it disagrees with the table.
    #[test]
    fn is_no_call_is_not_the_ambiguity_mask() {
        assert!(is_no_call(b'N') && is_no_call(b'n') && is_no_call(b'.'));
        assert!(!is_no_call(b'B'), "ambiguous, but not a no-call");
        assert!(
            bases_equal(b'.', b'N') && is_no_call(b'.') && is_no_call(b'N'),
            "here the two predicates happen to agree"
        );
    }
}

/// `SequenceUtil.isBisulfiteConverted(read, reference, negativeStrand)`.
///
/// On the positive strand a reference `C` read as `T` is a conversion, not a mismatch; on the
/// negative strand a reference `G` read as `A` is. Both tests go through [`bases_equal`], so
/// they fold case and accept `.` for `N` along with everything else that predicate accepts.
pub fn is_bisulfite_converted(read: u8, reference: u8, negative_strand: bool) -> bool {
    if negative_strand {
        bases_equal(reference, b'G') && bases_equal(read, b'A')
    } else {
        bases_equal(reference, b'C') && bases_equal(read, b'T')
    }
}

/// `SequenceUtil.bisulfiteBasesEqual(negativeStrand, read, reference)`.
pub fn bisulfite_bases_equal(negative_strand: bool, read: u8, reference: u8) -> bool {
    bases_equal(read, reference) || is_bisulfite_converted(read, reference, negative_strand)
}

/// `SequenceUtil.countMismatches(read, referenceBases, referenceOffset, bisulfiteSequence, matchAmbiguousRef)`,
/// restricted to `matchAmbiguousRef = false`, which is what every Picard caller passes.
///
/// Three things decide the answer and none is about arithmetic:
///
///   * the walk is over **alignment blocks**, so inserted and clipped read bases are not compared
///     at all and deleted reference bases are stepped over;
///   * the comparison is [`bases_equal`], so `N` in the read against `N` in the reference is a
///     **match**, and any byte outside the IUPAC table matches any other such byte;
///   * `referenceOffset` is subtracted from the reference block start, so the caller can pass a
///     slice of the contig rather than the whole thing. Picard's GC-bias collector passes 0 and
///     the whole contig.
///
/// htsjdk wraps the whole loop in `try { } catch (Exception e)` and rethrows as `SAMException`,
/// which turns an out-of-range reference into an error rather than a silent wrong answer. Here
/// that is the slice index panicking, which has the same effect and the same cause.
pub fn count_mismatches(
    read_bases: &[u8],
    blocks: &[crate::alignment_block::AlignmentBlock],
    reference_bases: &[u8],
    reference_offset: i32,
    negative_strand: bool,
    bisulfite_sequence: bool,
) -> i32 {
    let mut mismatches = 0;
    for block in blocks {
        let read_block_start = (block.read_start - 1) as usize;
        let reference_block_start = (block.reference_start - 1 - reference_offset) as usize;
        for i in 0..block.length as usize {
            let read = read_bases[read_block_start + i];
            let reference = reference_bases[reference_block_start + i];
            let matches = if bisulfite_sequence {
                bisulfite_bases_equal(negative_strand, read, reference)
            } else {
                bases_equal(read, reference)
            };
            if !matches {
                mismatches += 1;
            }
        }
    }
    mismatches
}

/// `SequenceUtil.countInsertedBases(cigar)`: the summed length of `I` elements.
pub fn count_inserted_bases(cigar: &crate::cigar::Cigar) -> i32 {
    cigar
        .elements
        .iter()
        .filter(|e| e.op == crate::cigar::Op::I)
        .map(|e| e.length as i32)
        .sum()
}

/// `SequenceUtil.countDeletedBases(cigar)`: the summed length of `D` elements.
///
/// `N` is **not** counted. A spliced read's skipped reference is not a deletion here, which
/// matters for RNA-seq alignments where `N` runs are long and common.
pub fn count_deleted_bases(cigar: &crate::cigar::Cigar) -> i32 {
    cigar
        .elements
        .iter()
        .filter(|e| e.op == crate::cigar::Op::D)
        .map(|e| e.length as i32)
        .sum()
}

#[cfg(test)]
mod counter_tests {
    use super::*;
    use crate::alignment_block::alignment_blocks;
    use crate::cigar::{Cigar, CigarElement, Op};

    fn cigar(spec: &[(u32, Op)]) -> Cigar {
        Cigar::new(
            spec.iter()
                .map(|&(length, op)| CigarElement { length, op })
                .collect(),
        )
    }

    #[test]
    fn mismatches_are_counted_only_inside_alignment_blocks() {
        let c = cigar(&[(4, Op::M), (2, Op::I), (4, Op::M)]);
        let blocks = alignment_blocks(&c, 1);
        // Read bases 4..6 are the insertion and are never compared, however wrong they look.
        let read = b"AAAAGGAAAA";
        let reference = b"AAAAAAAA";
        assert_eq!(
            count_mismatches(read, &blocks, reference, 0, false, false),
            0,
            "the inserted G's are not compared to anything"
        );
    }

    /// A deletion steps the reference forward without consuming read bases.
    #[test]
    fn a_deletion_steps_the_reference() {
        let c = cigar(&[(4, Op::M), (3, Op::D), (4, Op::M)]);
        let blocks = alignment_blocks(&c, 1);
        let read = b"AAAAAAAA";
        let reference = b"AAAAGGGAAAA";
        assert_eq!(
            count_mismatches(read, &blocks, reference, 0, false, false),
            0
        );
    }

    /// `N` against `N` is a match, because the comparison is mask equality.
    #[test]
    fn n_matches_n() {
        let c = cigar(&[(4, Op::M)]);
        let blocks = alignment_blocks(&c, 1);
        assert_eq!(
            count_mismatches(b"ANNA", &blocks, b"ANNA", 0, false, false),
            0
        );
        assert_eq!(
            count_mismatches(b"ANNA", &blocks, b"AAAA", 0, false, false),
            2,
            "N against A is a mismatch"
        );
    }

    #[test]
    fn the_reference_offset_shifts_the_reference_index() {
        let c = cigar(&[(4, Op::M)]);
        let blocks = alignment_blocks(&c, 101);
        let reference = b"GGGGAAAA";
        // Without the offset the read would be compared against reference[100..], out of range.
        assert_eq!(
            count_mismatches(b"AAAA", &blocks, reference, 96, false, false),
            0
        );
    }

    /// Bisulfite: C→T on the forward strand and G→A on the reverse are conversions, not
    /// mismatches, and each applies on one strand only.
    #[test]
    fn bisulfite_conversions_are_not_mismatches_on_their_own_strand() {
        let c = cigar(&[(4, Op::M)]);
        let blocks = alignment_blocks(&c, 1);
        assert_eq!(
            count_mismatches(b"TTTT", &blocks, b"CCCC", 0, false, true),
            0,
            "C read as T on the forward strand"
        );
        assert_eq!(
            count_mismatches(b"TTTT", &blocks, b"CCCC", 0, true, true),
            4,
            "the same pair on the reverse strand is a mismatch"
        );
        assert_eq!(
            count_mismatches(b"AAAA", &blocks, b"GGGG", 0, true, true),
            0,
            "G read as A on the reverse strand"
        );
    }

    #[test]
    fn without_bisulfite_a_conversion_is_a_mismatch() {
        let c = cigar(&[(4, Op::M)]);
        let blocks = alignment_blocks(&c, 1);
        assert_eq!(
            count_mismatches(b"TTTT", &blocks, b"CCCC", 0, false, false),
            4
        );
    }

    /// `N` in the CIGAR is a skip, not a deletion, and the counters keep them apart.
    #[test]
    fn skips_are_not_deletions() {
        let c = cigar(&[(4, Op::M), (100, Op::N), (4, Op::M), (3, Op::D), (2, Op::I)]);
        assert_eq!(count_deleted_bases(&c), 3, "the N run is not counted");
        assert_eq!(count_inserted_bases(&c), 2);
    }

    #[test]
    fn a_cigar_with_no_indels_counts_zero() {
        let c = cigar(&[(10, Op::M), (5, Op::S)]);
        assert_eq!(count_inserted_bases(&c), 0);
        assert_eq!(count_deleted_bases(&c), 0);
    }
}
