//! `htsjdk.samtools.DuplicateScoringStrategy`: how a duplicate set picks its survivor.
//!
//! The score is a `short`, and every step of it is written to stay far from overflow rather than to
//! be the obvious arithmetic. Reproducing it means reproducing that care, because the places it
//! saturates are the places two records tie and the tie-break decides which read a BAM keeps.
//!
//! * `SUM_OF_BASE_QUALITIES` sums qualities **at or above 15** and clamps the sum to
//!   `Short.MAX_VALUE / 2`, because two long high-quality reads would otherwise exceed a `short`
//!   when their two scores are added.
//! * `TOTAL_MAPPED_REFERENCE_LENGTH` takes the cigar's reference length, clamped the same way, and
//!   adds the mate's under `assumeMateCigar`.
//! * `RANDOM` hashes the **read name** with `Murmur3(1)`, masks to 14 bits, and subtracts
//!   `Short.MIN_VALUE / 4`. Both ends of a template hash the same name, so a pair scores together.
//!
//! And then, whatever the strategy, a vendor-failed record is discounted by `Short.MIN_VALUE / 2`,
//! once per end, which is why every branch above is capped at a quarter or a half of the range.
//!
//! The arithmetic is Java's: `short` addition wraps, and the clamps are what keep it from doing so.
//! The port uses `i16` with wrapping operations for the same reason, so a corpus that reaches the
//! cap agrees rather than panicking in debug and diverging in release.

use crate::cigar::Cigar;
use crate::murmur3::Murmur3;
use crate::record::BamRecord;

const READ_PAIRED: u16 = 0x1;
const READ_UNMAPPED: u16 = 0x4;
const MATE_UNMAPPED: u16 = 0x8;
const VENDOR_FAILED: u16 = 0x200;

const SHORT_MAX: i32 = 32767;
const SHORT_MIN: i32 = -32768;

/// `DuplicateScoringStrategy.ScoringStrategy`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScoringStrategy {
    SumOfBaseQualities,
    TotalMappedReferenceLength,
    Random,
}

/// `getSumOfBaseQualities`: the sum of the qualities that are **15 or more**.
fn sum_of_base_qualities(record: &BamRecord) -> i32 {
    record
        .base_qualities
        .iter()
        .filter(|&&q| q >= 15)
        .map(|&q| q as i32)
        .sum()
}

fn reference_length(cigar: &Cigar) -> i32 {
    cigar
        .elements
        .iter()
        .filter(|element| element.op.consumes_reference_bases())
        .map(|element| element.length as i32)
        .sum()
}

/// `computeDuplicateScore(record, strategy, assumeMateCigar)`.
///
/// `mate_cigar` is what `SAMUtils.getMateCigar` would return: `None` where the record carries no
/// `MC`, which under `assume_mate_cigar` is the case htsjdk throws on and a caller must not reach.
pub fn compute_duplicate_score(
    record: &BamRecord,
    strategy: ScoringStrategy,
    assume_mate_cigar: bool,
    mate_cigar: Option<&Cigar>,
) -> i16 {
    let mut score: i16 = 0;
    match strategy {
        ScoringStrategy::SumOfBaseQualities => {
            score = score.wrapping_add(sum_of_base_qualities(record).min(SHORT_MAX / 2) as i16);
        }
        ScoringStrategy::TotalMappedReferenceLength => {
            if record.flags & READ_UNMAPPED == 0 {
                score = reference_length(&record.cigar).min(SHORT_MAX / 2) as i16;
            }
            if assume_mate_cigar
                && record.flags & READ_PAIRED != 0
                && record.flags & MATE_UNMAPPED == 0
            {
                let mate = mate_cigar.expect("assume_mate_cigar with no MC is htsjdk's throw");
                score = score.wrapping_add(reference_length(mate).min(SHORT_MAX / 2) as i16);
            }
        }
        ScoringStrategy::Random => {
            // 14 bits of the name's hash, then shifted up so the result is non-negative: a number
            // between 0 and Short.MAX_VALUE / 2, which leaves room for the pair's second score and
            // for the vendor discount below.
            let hashed =
                Murmur3::new(1).hash_unencoded_chars(&record.read_name) & 0b11_1111_1111_1111;
            score = score.wrapping_add(hashed as i16);
            score = score.wrapping_sub((SHORT_MIN / 4) as i16);
        }
    }
    if record.flags & VENDOR_FAILED != 0 {
        score = score.wrapping_add((SHORT_MIN / 2) as i16);
    }
    score
}

/// `SAMUtils.getCanonicalRecordName`: the read name, with `/1` or `/2` appended for a paired read.
fn canonical_record_name(record: &BamRecord) -> String {
    const FIRST_OF_PAIR: u16 = 0x40;
    if record.flags & READ_PAIRED == 0 {
        return record.read_name.clone();
    }
    if record.flags & FIRST_OF_PAIR != 0 {
        format!("{}/1", record.read_name)
    } else {
        format!("{}/2", record.read_name)
    }
}

/// `String.compareTo`, which answers a **difference** and not a sign.
///
/// It returns the difference of the first differing UTF-16 code units, and the difference of the
/// lengths when one string is a prefix of the other. That matters here because
/// `DuplicateScoringStrategy.compare` returns it verbatim: comparing `read2` with `read5` answers
/// -3, not -1, and a port that normalized to a sign would disagree with the reference on every
/// tie-break while sorting identically.
fn java_string_compare(a: &str, b: &str) -> i32 {
    let (a_units, b_units): (Vec<u16>, Vec<u16>) =
        (a.encode_utf16().collect(), b.encode_utf16().collect());
    for (x, y) in a_units.iter().zip(b_units.iter()) {
        if x != y {
            return *x as i32 - *y as i32;
        }
    }
    a_units.len() as i32 - b_units.len() as i32
}

/// `DuplicateScoringStrategy.compare`: negative when `first` is the better record.
///
/// Paired beats unpaired before any score is computed, the scores are compared **descending** (the
/// higher score sorts first), and the name breaks a tie so the answer does not depend on input
/// order. Every branch returns the reference's own number rather than a sign: the score branch
/// returns a difference of shorts, and the tie-break returns `String.compareTo`.
pub fn compare(
    first: &BamRecord,
    second: &BamRecord,
    strategy: ScoringStrategy,
    assume_mate_cigar: bool,
    first_mate_cigar: Option<&Cigar>,
    second_mate_cigar: Option<&Cigar>,
) -> i32 {
    let first_paired = first.flags & READ_PAIRED != 0;
    let second_paired = second.flags & READ_PAIRED != 0;
    if first_paired != second_paired {
        return if first_paired { -1 } else { 1 };
    }
    let cmp = compute_duplicate_score(second, strategy, assume_mate_cigar, second_mate_cigar)
        as i32
        - compute_duplicate_score(first, strategy, assume_mate_cigar, first_mate_cigar) as i32;
    if cmp != 0 {
        return cmp;
    }
    java_string_compare(
        &canonical_record_name(first),
        &canonical_record_name(second),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(name: &str, quals: Vec<u8>, flags: u16) -> BamRecord {
        BamRecord {
            read_name: name.to_string(),
            base_qualities: quals,
            flags,
            ..Default::default()
        }
    }

    #[test]
    fn qualities_below_fifteen_do_not_count() {
        let low = record("a", vec![14, 14, 14], 0);
        let high = record("a", vec![15, 15, 15], 0);
        assert_eq!(
            compute_duplicate_score(&low, ScoringStrategy::SumOfBaseQualities, false, None),
            0
        );
        assert_eq!(
            compute_duplicate_score(&high, ScoringStrategy::SumOfBaseQualities, false, None),
            45
        );
    }

    #[test]
    fn the_sum_is_capped_at_a_quarter_of_a_short() {
        let long = record("a", vec![40; 5000], 0);
        assert_eq!(
            compute_duplicate_score(&long, ScoringStrategy::SumOfBaseQualities, false, None),
            (SHORT_MAX / 2) as i16
        );
    }

    #[test]
    fn a_vendor_failure_is_discounted_by_half_the_range() {
        let clean = record("a", vec![30, 30], 0);
        let failed = record("a", vec![30, 30], VENDOR_FAILED);
        let clean_score =
            compute_duplicate_score(&clean, ScoringStrategy::SumOfBaseQualities, false, None);
        let failed_score =
            compute_duplicate_score(&failed, ScoringStrategy::SumOfBaseQualities, false, None);
        assert_eq!(
            failed_score,
            clean_score.wrapping_add((SHORT_MIN / 2) as i16)
        );
        assert!(failed_score < clean_score);
    }

    #[test]
    fn both_ends_of_a_template_get_the_same_random_score() {
        let first = record("read1", vec![30], READ_PAIRED | 0x40);
        let second = record("read1", vec![10], READ_PAIRED | 0x80);
        assert_eq!(
            compute_duplicate_score(&first, ScoringStrategy::Random, false, None),
            compute_duplicate_score(&second, ScoringStrategy::Random, false, None)
        );
    }

    #[test]
    fn paired_beats_unpaired_before_any_score() {
        let paired = record("a", vec![], READ_PAIRED);
        let unpaired = record("z", vec![40; 100], 0);
        assert_eq!(
            compare(
                &paired,
                &unpaired,
                ScoringStrategy::SumOfBaseQualities,
                false,
                None,
                None
            ),
            -1
        );
    }

    #[test]
    fn a_tie_is_broken_by_the_canonical_name_and_keeps_its_magnitude() {
        let a = record("read2", vec![30], READ_PAIRED | 0x40);
        let b = record("read5", vec![30], READ_PAIRED | 0x40);
        // String.compareTo answers the difference of the code units, and compare returns it
        // verbatim: '2' - '5' is -3, not -1.
        assert_eq!(
            compare(
                &a,
                &b,
                ScoringStrategy::SumOfBaseQualities,
                false,
                None,
                None
            ),
            -3
        );
    }
}
