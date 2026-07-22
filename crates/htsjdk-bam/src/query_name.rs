//! Queryname ordering of records.
//!
//! Ports `htsjdk.samtools.SAMRecordQueryNameComparator` at tag 4.2.0: the order a
//! `SO:queryname` file is written in, and the order tools like SortSam and FastqToSam impose. The
//! primary key is the read name (`String.compareTo`), and equal names are then broken, in order,
//! by: paired before unpaired, first-of-pair before second, forward before reverse, primary before
//! secondary, non-supplementary before supplementary, and finally the `HI` (hit index) tag.
//!
//! `compareReadNames` is plain `String.compareTo`, i.e. lexicographic by UTF-16 code unit. Read
//! names are ASCII in every file this targets, where that equals byte order, so the port compares
//! the `&str`s directly and says so rather than pretending to a generality it does not exercise.

use std::cmp::Ordering;

use crate::record::BamRecord;
use crate::tag::{Tag, TagValue};

const READ_PAIRED: u16 = 0x1;
const READ_NEGATIVE_STRAND: u16 = 0x10;
const FIRST_OF_PAIR: u16 = 0x40;
const SECOND_OF_PAIR: u16 = 0x80;
const SECONDARY_ALIGNMENT: u16 = 0x100;
const SUPPLEMENTARY_ALIGNMENT: u16 = 0x800;

fn has(rec: &BamRecord, bit: u16) -> bool {
    rec.flags & bit != 0
}

/// `SAMRecord.getIntegerAttribute(SAMTag.HI)`.
fn hit_index(rec: &BamRecord) -> Option<i64> {
    match rec.tags.get(Tag::new(b"HI")) {
        Some(TagValue::Int(v)) => Some(*v),
        _ => None,
    }
}

/// `SAMRecordQueryNameComparator.compareReadNames`: `readName1.compareTo(readName2)`.
pub fn compare_read_names(name1: &str, name2: &str) -> Ordering {
    name1.cmp(name2)
}

/// `SAMRecordQueryNameComparator.compare`: the full queryname order with all tie-breaks.
pub fn compare(a: &BamRecord, b: &BamRecord) -> Ordering {
    // fileOrderCompare: the read name decides unless it ties.
    let cmp = compare_read_names(&a.read_name, &b.read_name);
    if cmp != Ordering::Equal {
        return cmp;
    }

    let a_paired = has(a, READ_PAIRED);
    let b_paired = has(b, READ_PAIRED);
    if a_paired || b_paired {
        if !a_paired {
            return Ordering::Greater; // an unpaired read sorts after a paired one
        }
        if !b_paired {
            return Ordering::Less;
        }
        if has(a, FIRST_OF_PAIR) && has(b, SECOND_OF_PAIR) {
            return Ordering::Less;
        }
        if has(a, SECOND_OF_PAIR) && has(b, FIRST_OF_PAIR) {
            return Ordering::Greater;
        }
    }

    if has(a, READ_NEGATIVE_STRAND) != has(b, READ_NEGATIVE_STRAND) {
        // Forward (false) before reverse (true).
        return if has(a, READ_NEGATIVE_STRAND) {
            Ordering::Greater
        } else {
            Ordering::Less
        };
    }
    if has(a, SECONDARY_ALIGNMENT) != has(b, SECONDARY_ALIGNMENT) {
        // Primary before secondary.
        return if has(b, SECONDARY_ALIGNMENT) {
            Ordering::Less
        } else {
            Ordering::Greater
        };
    }
    if has(a, SUPPLEMENTARY_ALIGNMENT) != has(b, SUPPLEMENTARY_ALIGNMENT) {
        // Non-supplementary before supplementary.
        return if has(b, SUPPLEMENTARY_ALIGNMENT) {
            Ordering::Less
        } else {
            Ordering::Greater
        };
    }

    // The HI tag last: a record with one sorts after a record without.
    match (hit_index(a), hit_index(b)) {
        (Some(h1), Some(h2)) => h1.cmp(&h2),
        (Some(_), None) => Ordering::Greater,
        (None, Some(_)) => Ordering::Less,
        (None, None) => Ordering::Equal,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rec(name: &str, flags: u16) -> BamRecord {
        BamRecord {
            read_name: name.to_string(),
            flags,
            ..Default::default()
        }
    }

    #[test]
    fn read_name_is_the_primary_key_and_is_plain_string_order() {
        // String.compareTo: "r10" < "r2" because '1' < '2' at the second character.
        assert_eq!(compare(&rec("r10", 0), &rec("r2", 0)), Ordering::Less);
        assert_eq!(compare(&rec("b", 0), &rec("a", 0)), Ordering::Greater);
        // A prefix sorts before the longer name.
        assert_eq!(compare(&rec("r", 0), &rec("r1", 0)), Ordering::Less);
    }

    #[test]
    fn a_paired_read_sorts_before_an_unpaired_one_of_the_same_name() {
        let paired = rec("x", READ_PAIRED | FIRST_OF_PAIR);
        let unpaired = rec("x", 0);
        assert_eq!(compare(&paired, &unpaired), Ordering::Less);
        assert_eq!(compare(&unpaired, &paired), Ordering::Greater);
    }

    #[test]
    fn first_of_pair_sorts_before_second_of_pair() {
        let first = rec("x", READ_PAIRED | FIRST_OF_PAIR);
        let second = rec("x", READ_PAIRED | SECOND_OF_PAIR);
        assert_eq!(compare(&first, &second), Ordering::Less);
        assert_eq!(compare(&second, &first), Ordering::Greater);
    }

    #[test]
    fn forward_sorts_before_reverse() {
        let fwd = rec("x", 0);
        let rev = rec("x", READ_NEGATIVE_STRAND);
        assert_eq!(compare(&fwd, &rev), Ordering::Less);
    }

    #[test]
    fn primary_before_secondary_and_before_supplementary() {
        assert_eq!(
            compare(&rec("x", 0), &rec("x", SECONDARY_ALIGNMENT)),
            Ordering::Less
        );
        assert_eq!(
            compare(&rec("x", 0), &rec("x", SUPPLEMENTARY_ALIGNMENT)),
            Ordering::Less
        );
    }

    #[test]
    fn the_hit_index_tag_breaks_a_full_tie() {
        let mut a = rec("x", 0);
        let mut b = rec("x", 0);
        a.tags.insert(Tag::new(b"HI"), TagValue::Int(1));
        b.tags.insert(Tag::new(b"HI"), TagValue::Int(2));
        assert_eq!(compare(&a, &b), Ordering::Less);
        // A record with HI sorts after one without.
        let without = rec("x", 0);
        assert_eq!(compare(&a, &without), Ordering::Greater);
        assert_eq!(compare(&without, &a), Ordering::Less);
    }

    #[test]
    fn two_identical_records_tie() {
        assert_eq!(compare(&rec("x", 0), &rec("x", 0)), Ordering::Equal);
    }
}
