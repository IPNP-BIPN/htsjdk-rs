//! Conformance for the CRAM record's three flag words and its mate chain, against
//! `htsjdk.samtools.cram.structure.CRAMCompressionRecord`.
//!
//! Goldens from `tools/cram-conformance/CramRecordFlagsDump.java` in the pinned oracle.
//!
//! The rows that justify the suite:
//!
//! ```text
//! mask    cram  511  255
//! chain   triple                  0  1   0  0   200  250
//! chain   triple                  1  33  1  0   300  0
//! chain   triple                  2  17  0  0   100  -250
//! chain   pair-mate-no-reference  0  1   0  -1  0    0
//! ```
//!
//! The two narrow words are masked to a byte. A chain becomes a ring, and the template length is
//! computed once and negated once, so the middle of a triple keeps the zero it was built with. And
//! a mate on no reference loses its position.

use std::io::Read;

use htsjdk_cram::record_flags::{
    compute_insert_size, restore_mate_info, Flags, MateRecord, NO_ALIGNMENT_REFERENCE_INDEX,
};

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/cram_record_flags.txt.gz");
    let file = std::fs::File::open(&path).expect("corpus");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("corpus is gzip");
    text
}

fn rows<'a>(corpus: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    let prefix = format!("{kind}\t");
    corpus
        .lines()
        .filter_map(|line| line.strip_prefix(&prefix))
        .map(|rest| rest.split('\t').collect())
        .collect()
}

/// Every predicate a record answers from its three words.
#[test]
fn the_predicates_are_the_reference_predicates() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "flags") {
        let flags = Flags {
            bam: row[0].parse().expect("bam"),
            cram: row[1].parse().expect("cram"),
            mate: row[2].parse().expect("mate"),
        };
        let answered = [
            ("detached", flags.is_detached()),
            ("mateDownstream", flags.has_mate_downstream()),
            ("forceQuality", flags.is_force_preserve_quality_scores()),
            ("unknownBases", flags.is_unknown_bases()),
            ("paired", flags.is_read_paired()),
            ("unmapped", flags.is_segment_unmapped()),
            ("first", flags.is_first_segment()),
            ("last", flags.is_last_segment()),
            ("secondary", flags.is_secondary_alignment()),
        ]
        .iter()
        .map(|(name, value)| format!("{name}={value}"))
        .collect::<Vec<_>>()
        .join(",");
        assert_eq!(answered, row[3], "flags {}/{}/{}", row[0], row[1], row[2]);
        compared += 1;
    }
    assert_eq!(compared, 28, "flag words compared");
}

/// The two narrow words are masked to a byte on the way out rather than refused on the way in.
#[test]
fn the_narrow_words_are_masked_to_a_byte() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "mask") {
        let (field, stored, returned) = (
            row[0],
            row[1].parse::<i32>().expect("stored"),
            row[2].parse::<i32>().expect("returned"),
        );
        let flags = Flags {
            bam: 0,
            cram: stored,
            mate: stored,
        };
        let actual = match field {
            "cram" => flags.cram_flags(),
            "mate" => flags.mate_flags(),
            other => panic!("{other}"),
        };
        assert_eq!(actual, returned, "mask {field} {stored}");
        compared += 1;
    }
    assert_eq!(compared, 6, "masked values compared");
}

/// A chain restored into a ring, record by record.
#[test]
fn a_chain_becomes_the_ring_the_reference_made() {
    let corpus = corpus();
    let mut compared = 0;
    let mut labels: Vec<&str> = Vec::new();
    for row in rows(&corpus, "chain") {
        if !labels.contains(&row[0]) {
            labels.push(row[0]);
        }
        compared += 1;
    }

    for label in &labels {
        let expected: Vec<Vec<&str>> = rows(&corpus, "chain")
            .into_iter()
            .filter(|row| row[0] == *label)
            .collect();
        let mut records = chain(label);
        assert_eq!(records.len(), expected.len(), "{label} length");
        restore_mate_info(&mut records);

        for (index, record) in records.iter().enumerate() {
            let row = &expected[index];
            assert_eq!(record.flags.bam.to_string(), row[2], "{label}[{index}] bam");
            assert_eq!(
                record.flags.mate_flags().to_string(),
                row[3],
                "{label}[{index}] mate flags"
            );
            assert_eq!(
                record.mate_reference_index.to_string(),
                row[4],
                "{label}[{index}] mate reference"
            );
            assert_eq!(
                record.mate_alignment_start.to_string(),
                row[5],
                "{label}[{index}] mate start"
            );
            assert_eq!(
                record.template_size.to_string(),
                row[6],
                "{label}[{index}] template size"
            );
        }
    }
    assert_eq!(compared, 18, "chain rows compared");
    assert_eq!(labels.len(), 9, "chains compared");
}

/// The corpus's chains, built the way the dump built them. Every record is 50 bases with no read
/// features, so its alignment end is its start plus 49.
fn chain(label: &str) -> Vec<MateRecord> {
    match label {
        "pair-forward" => vec![mapped(0, 100, false), mapped(0, 200, false)],
        "pair-reverse-second" => vec![mapped(0, 100, false), mapped(0, 200, true)],
        "pair-same-start" => vec![mapped(0, 100, false), mapped(0, 100, false)],
        "pair-second-before" => vec![mapped(0, 200, false), mapped(0, 100, false)],
        "pair-other-reference" => vec![mapped(0, 100, false), mapped(1, 200, false)],
        "pair-unmapped-mate" => vec![mapped(0, 100, false), unmapped()],
        "pair-mate-no-reference" => vec![
            mapped(0, 100, false),
            mapped(NO_ALIGNMENT_REFERENCE_INDEX, 200, false),
        ],
        "triple" => vec![
            mapped(0, 100, false),
            mapped(0, 200, false),
            mapped(0, 300, true),
        ],
        "single" => vec![mapped(0, 100, false)],
        other => panic!("{other}"),
    }
}

fn mapped(reference_index: i32, alignment_start: i32, negative_strand: bool) -> MateRecord {
    MateRecord {
        flags: Flags {
            bam: 0x1 | if negative_strand { 0x10 } else { 0 },
            cram: 0,
            mate: 0,
        },
        reference_index,
        alignment_start,
        alignment_end: alignment_start + 49,
        mate_reference_index: -1,
        mate_alignment_start: 0,
        template_size: 0,
        records_to_next_fragment: -1,
    }
}

fn unmapped() -> MateRecord {
    MateRecord {
        flags: Flags {
            bam: 0x1 | 0x4,
            cram: 0,
            mate: 0,
        },
        reference_index: 0,
        alignment_start: 200,
        alignment_end: 249,
        mate_reference_index: -1,
        mate_alignment_start: 0,
        template_size: 0,
        records_to_next_fragment: -1,
    }
}

/// The template length, stated on its own: zero where either end is unmapped or the two are on
/// different references, and otherwise never zero, because of the sign that is added to it.
#[test]
fn the_template_length_is_never_zero_between_two_mapped_ends() {
    let corpus = corpus();
    for row in rows(&corpus, "insert") {
        let (first, last) = (
            row[1].parse::<i32>().expect("first"),
            row[2].parse::<i32>().expect("last"),
        );
        assert_eq!(first, -last, "{} is negated on the last", row[0]);
    }

    // Two ends at the same position are one apart, not zero apart.
    let same = compute_insert_size(&mapped(0, 100, false), &mapped(0, 100, false));
    assert_eq!(same, 1);
    // And the sign follows which end comes first.
    assert_eq!(
        compute_insert_size(&mapped(0, 200, false), &mapped(0, 100, false)),
        -101
    );
    // A negative-strand end measures from its alignment end.
    assert_eq!(
        compute_insert_size(&mapped(0, 100, false), &mapped(0, 200, true)),
        150
    );
    // And an unmapped end makes it zero however the two are placed.
    assert_eq!(compute_insert_size(&mapped(0, 100, false), &unmapped()), 0);
}

/// Detaching sets one bit, clears another, and forgets the distance to the next fragment.
#[test]
fn detaching_is_the_reference_state() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "detach") {
        let (before, before_next, after, after_next) = (
            row[0].parse::<i32>().expect("before"),
            row[1].parse::<i32>().expect("before next"),
            row[2].parse::<i32>().expect("after"),
            row[3].parse::<i32>().expect("after next"),
        );
        let mut record = mapped(0, 100, false);
        record.flags.cram = before;
        record.records_to_next_fragment = before_next;
        assert_eq!(record.flags.cram_flags(), before);
        record.set_to_detached_state();
        assert_eq!(record.flags.cram_flags(), after, "detach from {before}");
        assert_eq!(record.records_to_next_fragment, after_next);
        compared += 1;
    }
    assert_eq!(compared, 3, "detachments compared");
}
