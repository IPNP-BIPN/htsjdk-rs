//! Coordinate ordering of records.
//!
//! Ports `htsjdk.samtools.SAMRecordCoordinateComparator` at tag 4.2.0: the order a `SO:coordinate`
//! file is written in, which SortSam, MarkDuplicates, and most of the alignment pipeline depend on.
//!
//! The primary key is the file-order position: reference index, then alignment start, with an
//! **unmapped read (reference index -1) sorting last**. Records at the same position are then broken
//! by strand (forward before reverse) and, within the same strand, by read name, flags, mapping
//! quality, mate reference index, mate alignment start, and inferred insert size, in that order.
//!
//! Risk R11/R5 in the plan flagged this comparator as *suspected non-total*. Reading it settles the
//! question: it is a **total order**. Two records that agree on every field above compare equal, so
//! their relative order is decided by the sort's stability, not by the comparator. Byte-identity
//! therefore needs a **stable, in-memory sort** (no spill reordering), which is why the oracle
//! contract pins `MAX_RECORDS_IN_RAM` for coordinate-sorting tools. See decision 0021.

use std::cmp::Ordering;

use crate::record::BamRecord;

const READ_NEGATIVE_STRAND: u16 = 0x10;
/// `SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX`.
const NO_ALIGNMENT_REFERENCE_INDEX: i32 = -1;

/// `SAMRecordCoordinateComparator.fileOrderCompare`: reference index then alignment start, with an
/// unmapped read (reference index -1) after every mapped one.
pub fn file_order_compare(a: &BamRecord, b: &BamRecord) -> Ordering {
    let r1 = a.reference_index;
    let r2 = b.reference_index;
    if r1 == NO_ALIGNMENT_REFERENCE_INDEX {
        return if r2 == NO_ALIGNMENT_REFERENCE_INDEX {
            Ordering::Equal
        } else {
            Ordering::Greater
        };
    }
    if r2 == NO_ALIGNMENT_REFERENCE_INDEX {
        return Ordering::Less;
    }
    r1.cmp(&r2)
        .then_with(|| a.alignment_start.cmp(&b.alignment_start))
}

/// `SAMRecordCoordinateComparator.compare`: the full coordinate order with all tie-breaks.
pub fn compare(a: &BamRecord, b: &BamRecord) -> Ordering {
    let cmp = file_order_compare(a, b);
    if cmp != Ordering::Equal {
        return cmp;
    }

    let neg1 = a.flags & READ_NEGATIVE_STRAND != 0;
    let neg2 = b.flags & READ_NEGATIVE_STRAND != 0;
    if neg1 == neg2 {
        // compareInts on each field in turn; the read name first, as htsjdk does.
        a.read_name
            .cmp(&b.read_name)
            .then_with(|| a.flags.cmp(&b.flags))
            .then_with(|| a.mapping_quality.cmp(&b.mapping_quality))
            .then_with(|| a.mate_reference_index.cmp(&b.mate_reference_index))
            .then_with(|| a.mate_alignment_start.cmp(&b.mate_alignment_start))
            .then_with(|| a.inferred_insert_size.cmp(&b.inferred_insert_size))
    } else {
        // Forward (false) before reverse (true).
        if neg1 {
            Ordering::Greater
        } else {
            Ordering::Less
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rec(ref_index: i32, start: i32, flags: u16, name: &str) -> BamRecord {
        BamRecord {
            reference_index: ref_index,
            alignment_start: start,
            flags,
            read_name: name.to_string(),
            ..Default::default()
        }
    }

    #[test]
    fn reference_index_then_start_is_the_primary_key() {
        assert_eq!(
            compare(&rec(0, 100, 0, "a"), &rec(1, 5, 0, "a")),
            Ordering::Less
        );
        assert_eq!(
            compare(&rec(0, 200, 0, "a"), &rec(0, 100, 0, "a")),
            Ordering::Greater
        );
    }

    #[test]
    fn an_unmapped_read_sorts_after_every_mapped_read() {
        let unmapped = rec(NO_ALIGNMENT_REFERENCE_INDEX, 0, 0x4, "z");
        let mapped = rec(5, 999_999, 0, "a");
        assert_eq!(compare(&unmapped, &mapped), Ordering::Greater);
        assert_eq!(compare(&mapped, &unmapped), Ordering::Less);
        // Two unmapped reads fall through to the tie-breaks (both have refIndex -1).
        let u2 = rec(NO_ALIGNMENT_REFERENCE_INDEX, 0, 0x4, "a");
        assert_eq!(compare(&unmapped, &u2), Ordering::Greater); // "z" > "a"
    }

    #[test]
    fn at_the_same_position_forward_sorts_before_reverse() {
        let fwd = rec(0, 100, 0, "z");
        let rev = rec(0, 100, READ_NEGATIVE_STRAND, "a");
        // Strand decides before the read name, so the forward "z" still precedes the reverse "a".
        assert_eq!(compare(&fwd, &rev), Ordering::Less);
    }

    #[test]
    fn same_position_and_strand_falls_through_to_the_name() {
        assert_eq!(
            compare(&rec(0, 100, 0, "a"), &rec(0, 100, 0, "b")),
            Ordering::Less
        );
    }

    #[test]
    fn fully_equal_records_tie() {
        assert_eq!(
            compare(&rec(0, 100, 0, "a"), &rec(0, 100, 0, "a")),
            Ordering::Equal
        );
    }

    /// A small consistency check on the strand branch: within one position, all forward reads
    /// precede all reverse reads, and a stable sort keeps equal records in input order.
    #[test]
    fn sorting_partitions_a_position_into_forward_then_reverse() {
        let mut recs = [
            rec(0, 100, READ_NEGATIVE_STRAND, "a"),
            rec(0, 100, 0, "b"),
            rec(0, 100, READ_NEGATIVE_STRAND, "c"),
            rec(0, 100, 0, "a"),
        ];
        recs.sort_by(compare);
        let seen: Vec<(&str, bool)> = recs
            .iter()
            .map(|r| (r.read_name.as_str(), r.flags & READ_NEGATIVE_STRAND != 0))
            .collect();
        assert_eq!(seen, [("a", false), ("b", false), ("a", true), ("c", true)]);
    }
}
