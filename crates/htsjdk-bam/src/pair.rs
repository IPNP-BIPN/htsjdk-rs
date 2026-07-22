//! Mate pairing information.
//!
//! Ports `htsjdk.samtools.SamPairUtil.setMateInfo` and `computeInsertSize` at tag 4.2.0: given the
//! two ends of a template, stamp each with its mate's reference, position, strand, unmapped flag,
//! mapping quality (`MQ`) and, optionally, mate cigar (`MC`), and set the inferred insert size.
//! This is what FixMateInformation and MergeBamAlignment rely on to make a pair self-consistent.
//!
//! Three cases, each reproduced exactly: both ends mapped, both unmapped, and one of each. The
//! insert size is only nonzero when both ends are mapped to the same reference; its sign is `+` on
//! the first end and `-` on the second, and its magnitude counts from 5' end to 5' end with a `±1`
//! adjustment so that abutting reads have insert size `±1`, not `0`.

use crate::record::BamRecord;
use crate::tag::{Tag, TagValue};

const READ_UNMAPPED: u16 = 0x4;
const MATE_UNMAPPED: u16 = 0x8;
const READ_NEGATIVE_STRAND: u16 = 0x10;
const MATE_NEGATIVE_STRAND: u16 = 0x20;

const NO_ALIGNMENT_REFERENCE_INDEX: i32 = -1;
const NO_ALIGNMENT_START: i32 = 0;

fn is_unmapped(rec: &BamRecord) -> bool {
    rec.flags & READ_UNMAPPED != 0
}

fn is_negative_strand(rec: &BamRecord) -> bool {
    rec.flags & READ_NEGATIVE_STRAND != 0
}

fn set_flag(rec: &mut BamRecord, bit: u16, value: bool) {
    if value {
        rec.flags |= bit;
    } else {
        rec.flags &= !bit;
    }
}

fn mq(mapping_quality: u8) -> TagValue {
    TagValue::Int(mapping_quality as i64)
}

fn mc(cigar_text: String) -> TagValue {
    TagValue::Str(cigar_text)
}

/// `SamPairUtil.computeInsertSize`.
pub fn compute_insert_size(first_end: &BamRecord, second_end: &BamRecord) -> i32 {
    if is_unmapped(first_end) || is_unmapped(second_end) {
        return 0;
    }
    if first_end.reference_index != second_end.reference_index {
        return 0;
    }
    let first_5prime = if is_negative_strand(first_end) {
        first_end.alignment_end()
    } else {
        first_end.alignment_start
    };
    let second_5prime = if is_negative_strand(second_end) {
        second_end.alignment_end()
    } else {
        second_end.alignment_start
    };
    let adjustment = if second_5prime >= first_5prime { 1 } else { -1 };
    second_5prime - first_5prime + adjustment
}

/// `SamPairUtil.setMateInfo(rec1, rec2, setMateCigar)`.
pub fn set_mate_info(rec1: &mut BamRecord, rec2: &mut BamRecord, set_mate_cigar: bool) {
    let mc_tag = Tag::new(b"MC");
    let mq_tag = Tag::new(b"MQ");

    if !is_unmapped(rec1) && !is_unmapped(rec2) {
        // Both mapped: cross-copy coordinates, strand, mapping quality, and (optionally) cigar.
        let (r1_ref, r1_start, r1_neg, r1_mapq, r1_cigar) = snapshot(rec1);
        let (r2_ref, r2_start, r2_neg, r2_mapq, r2_cigar) = snapshot(rec2);

        rec1.mate_reference_index = r2_ref;
        rec1.mate_alignment_start = r2_start;
        set_flag(rec1, MATE_NEGATIVE_STRAND, r2_neg);
        set_flag(rec1, MATE_UNMAPPED, false);
        rec1.tags.insert(mq_tag, mq(r2_mapq));

        rec2.mate_reference_index = r1_ref;
        rec2.mate_alignment_start = r1_start;
        set_flag(rec2, MATE_NEGATIVE_STRAND, r1_neg);
        set_flag(rec2, MATE_UNMAPPED, false);
        rec2.tags.insert(mq_tag, mq(r1_mapq));

        if set_mate_cigar {
            rec1.tags.insert(mc_tag, mc(r2_cigar));
            rec2.tags.insert(mc_tag, mc(r1_cigar));
        } else {
            rec1.tags.remove(mc_tag);
            rec2.tags.remove(mc_tag);
        }

        let insert_size = compute_insert_size(rec1, rec2);
        rec1.inferred_insert_size = insert_size;
        rec2.inferred_insert_size = -insert_size;
    } else if is_unmapped(rec1) && is_unmapped(rec2) {
        // Both unmapped: clear coordinates, keep strand, and cross-set the mate-unmapped flag.
        clear_both_unmapped(rec1, is_negative_strand(rec2), mq_tag, mc_tag);
        clear_both_unmapped(rec2, is_negative_strand(rec1), mq_tag, mc_tag);
    } else {
        // Exactly one mapped: copy its coordinates onto the unmapped mate.
        if is_unmapped(rec1) {
            set_one_mapped(rec2, rec1, set_mate_cigar, mq_tag, mc_tag);
        } else {
            set_one_mapped(rec1, rec2, set_mate_cigar, mq_tag, mc_tag);
        }
    }
}

fn snapshot(rec: &BamRecord) -> (i32, i32, bool, u8, String) {
    (
        rec.reference_index,
        rec.alignment_start,
        is_negative_strand(rec),
        rec.mapping_quality,
        rec.cigar.to_text(),
    )
}

fn clear_both_unmapped(rec: &mut BamRecord, mate_negative: bool, mq_tag: Tag, mc_tag: Tag) {
    rec.reference_index = NO_ALIGNMENT_REFERENCE_INDEX;
    rec.alignment_start = NO_ALIGNMENT_START;
    rec.mate_reference_index = NO_ALIGNMENT_REFERENCE_INDEX;
    rec.mate_alignment_start = NO_ALIGNMENT_START;
    set_flag(rec, MATE_NEGATIVE_STRAND, mate_negative);
    set_flag(rec, MATE_UNMAPPED, true);
    rec.tags.remove(mq_tag);
    rec.tags.remove(mc_tag);
    rec.inferred_insert_size = 0;
}

fn set_one_mapped(
    mapped: &mut BamRecord,
    unmapped: &mut BamRecord,
    set_mate_cigar: bool,
    mq_tag: Tag,
    mc_tag: Tag,
) {
    let m_ref = mapped.reference_index;
    let m_start = mapped.alignment_start;
    let m_neg = is_negative_strand(mapped);
    let m_mapq = mapped.mapping_quality;
    let m_cigar = mapped.cigar.to_text();

    // The unmapped read takes the mapped read's coordinates.
    unmapped.reference_index = m_ref;
    unmapped.alignment_start = m_start;
    let u_neg = is_negative_strand(unmapped);

    mapped.mate_reference_index = m_ref;
    mapped.mate_alignment_start = m_start;
    set_flag(mapped, MATE_NEGATIVE_STRAND, u_neg);
    set_flag(mapped, MATE_UNMAPPED, true);
    mapped.tags.remove(mq_tag);
    mapped.tags.remove(mc_tag);
    mapped.inferred_insert_size = 0;

    unmapped.mate_reference_index = m_ref;
    unmapped.mate_alignment_start = m_start;
    set_flag(unmapped, MATE_NEGATIVE_STRAND, m_neg);
    set_flag(unmapped, MATE_UNMAPPED, false);
    unmapped.tags.insert(mq_tag, mq(m_mapq));
    if set_mate_cigar {
        unmapped.tags.insert(mc_tag, mc(m_cigar));
    } else {
        unmapped.tags.remove(mc_tag);
    }
    unmapped.inferred_insert_size = 0;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cigar::{Cigar, CigarElement, Op};

    fn rec(flags: u16, ref_index: i32, start: i32, mapq: u8, cigar_len: u32) -> BamRecord {
        BamRecord {
            flags,
            reference_index: ref_index,
            alignment_start: start,
            mapping_quality: mapq,
            cigar: Cigar::new(vec![CigarElement {
                length: cigar_len,
                op: Op::M,
            }]),
            ..Default::default()
        }
    }

    #[test]
    fn insert_size_counts_5prime_to_5prime_with_a_unit_adjustment() {
        // Forward read at 100 (36M) and reverse read at 200 (36M): 5' positions 100 and 235.
        let fwd = rec(0, 0, 100, 60, 36);
        let rev = rec(READ_NEGATIVE_STRAND, 0, 200, 60, 36);
        // second5 (235) - first5 (100) + 1 = 136.
        assert_eq!(compute_insert_size(&fwd, &rev), 136);
        // Symmetric magnitude the other way: first5=235, second5=100, adjustment -1 -> -136.
        assert_eq!(compute_insert_size(&rev, &fwd), -136);
    }

    #[test]
    fn insert_size_is_zero_across_references_or_when_unmapped() {
        assert_eq!(
            compute_insert_size(&rec(0, 0, 100, 60, 36), &rec(0, 1, 100, 60, 36)),
            0
        );
        assert_eq!(
            compute_insert_size(&rec(READ_UNMAPPED, 0, 0, 0, 36), &rec(0, 0, 100, 60, 36)),
            0
        );
    }

    #[test]
    fn both_mapped_cross_sets_mate_fields_and_insert_size() {
        let mut r1 = rec(0x1 | 0x40, 0, 100, 60, 36); // paired, first of pair, forward
        let mut r2 = rec(0x1 | 0x80 | READ_NEGATIVE_STRAND, 0, 200, 50, 36); // reverse
        set_mate_info(&mut r1, &mut r2, true);

        assert_eq!(r1.mate_reference_index, 0);
        assert_eq!(r1.mate_alignment_start, 200);
        assert!(r1.flags & MATE_NEGATIVE_STRAND != 0); // mate (r2) is reverse
        assert!(r1.flags & MATE_UNMAPPED == 0);
        assert_eq!(r1.tags.get(Tag::new(b"MQ")), Some(&TagValue::Int(50)));
        assert_eq!(
            r1.tags.get(Tag::new(b"MC")),
            Some(&TagValue::Str("36M".into()))
        );
        assert_eq!(r1.inferred_insert_size, 136);
        assert_eq!(r2.inferred_insert_size, -136);
    }

    #[test]
    fn set_mate_cigar_false_removes_the_mc_tag() {
        let mut r1 = rec(0x1, 0, 100, 60, 36);
        let mut r2 = rec(0x1, 0, 200, 60, 36);
        r1.tags
            .insert(Tag::new(b"MC"), TagValue::Str("stale".into()));
        set_mate_info(&mut r1, &mut r2, false);
        assert_eq!(r1.tags.get(Tag::new(b"MC")), None);
    }

    #[test]
    fn one_mapped_copies_coordinates_to_the_unmapped_mate() {
        let mut mapped = rec(0x1 | 0x40, 2, 500, 60, 36);
        let mut unmapped = rec(0x1 | 0x80 | READ_UNMAPPED, -1, 0, 0, 0);
        set_mate_info(&mut mapped, &mut unmapped, true);

        // The unmapped read inherits the mapped read's coordinates.
        assert_eq!(unmapped.reference_index, 2);
        assert_eq!(unmapped.alignment_start, 500);
        assert!(unmapped.flags & MATE_UNMAPPED == 0);
        assert_eq!(unmapped.tags.get(Tag::new(b"MQ")), Some(&TagValue::Int(60)));
        // The mapped read's mate is the unmapped one.
        assert!(mapped.flags & MATE_UNMAPPED != 0);
        assert_eq!(mapped.tags.get(Tag::new(b"MQ")), None);
        assert_eq!(mapped.inferred_insert_size, 0);
    }
}
