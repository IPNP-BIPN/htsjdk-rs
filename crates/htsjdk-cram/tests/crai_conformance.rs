//! Conformance for the CRAM index, against `htsjdk.samtools.cram.CRAIEntry` and `CRAIIndex`.
//!
//! Goldens from `tools/cram-conformance/CramCraiDump.java` in the pinned oracle.
//!
//! The rows that justify the suite:
//!
//! ```text
//! sort       unmapped-last  1:100:10:0:0:1;-1:1:1:0:0:1;0:100:10:0:0:1  0:100:10:0:0:1;1:100:10:0:0:1;-1:1:1:0:0:1
//! intersect  0:100:0:0:0:1  0:100:0:0:0:1  false
//! intersect  0:100:10:0:0:1  0:100:0:0:0:1  false
//! find       whole-sequence-by-span  0  100  0  0:100:10:0:0:1;0:200:10:0:1:1;0:300:10:0:2:1
//! ```
//!
//! Unmapped sorts last whatever its start says. Two identical entries of span zero do not
//! intersect, and neither does a zero-span entry with one that contains it, because the test is a
//! midpoint comparison rather than an overlap. And a span below one is a wildcard, not an empty
//! query.

use std::io::Read;

use htsjdk_cram::crai::{find, leftmost, write_index, CraiEntry, CraiError};

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/cram_crai.txt.gz");
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

/// `sequence:start:span:container:offset:size`, the notation the dump writes a list in.
fn entry(specification: &str) -> CraiEntry {
    let parts: Vec<&str> = specification.split(':').collect();
    CraiEntry::new(
        parts[0].parse().expect("sequence"),
        parts[1].parse().expect("start"),
        parts[2].parse().expect("span"),
        parts[3].parse().expect("container"),
        parts[4].parse().expect("slice offset"),
        parts[5].parse().expect("slice size"),
    )
    .expect("built")
}

fn entries(column: &str) -> Vec<CraiEntry> {
    if column == "-" {
        return Vec::new();
    }
    column.split(';').map(entry).collect()
}

fn show(entries: &[CraiEntry]) -> String {
    if entries.is_empty() {
        return "-".to_string();
    }
    entries
        .iter()
        .map(|entry| {
            format!(
                "{}:{}:{}:{}:{}:{}",
                entry.sequence_id,
                entry.alignment_start,
                entry.alignment_span,
                entry.container_start_byte_offset,
                entry.slice_byte_offset_from_compression_header_start,
                entry.slice_byte_size
            )
        })
        .collect::<Vec<_>>()
        .join(";")
}

/// The line an entry makes, with its tabs shown as spaces the way the dump shows them.
#[test]
fn an_entry_serializes_to_the_reference_line() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "entry") {
        let entry = CraiEntry::new(
            row[0].parse().expect("sequence"),
            row[1].parse().expect("start"),
            row[2].parse().expect("span"),
            row[3].parse().expect("container"),
            row[4].parse().expect("slice offset"),
            row[5].parse().expect("slice size"),
        )
        .expect("built");
        assert_eq!(
            entry.serialize().replace('\t', " "),
            row[6],
            "entry {}",
            row[0]
        );
        compared += 1;
    }
    assert_eq!(compared, 5, "entries serialized");

    // Negative starts, spans and offsets all go through: the constructor checks one thing only.
    let negative = CraiEntry::new(-1, -1, -1, 500, 10, 20).expect("built");
    assert_eq!(
        negative.serialize().replace('\t', " "),
        "-1 -1 -1 500 10 20"
    );
}

/// And the entry a line makes.
#[test]
fn a_line_parses_to_the_reference_entry() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "parse") {
        // The dump shows the line with spaces, because its own rows are tab separated.
        let line = row[0].replace(' ', "\t");
        let entry = CraiEntry::parse(&line).expect("parsed");
        assert_eq!(entry.serialize().replace('\t', " "), row[1], "parse {line}");
        compared += 1;
    }
    assert_eq!(compared, 3, "lines parsed");
}

/// Sorting, which is what writing an index does to it.
#[test]
fn the_order_is_the_reference_order() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "sort") {
        let (label, before, after) = (row[0], entries(row[1]), row[2]);
        let mut sorted = before.clone();
        sorted.sort();
        assert_eq!(show(&sorted), after, "sort {label}");

        // And the bytes an index of them makes, which is the sort plus a newline each. The dump
        // shows its newlines as `|` and its tabs as spaces.
        let index = rows(&corpus, "index")
            .into_iter()
            .find(|index| index[0] == label)
            .expect("an index row");
        let written = String::from_utf8(write_index(&before)).expect("utf8");
        assert_eq!(
            written.replace('\n', "|").replace('\t', " ").trim_end(),
            index[1],
            "index {label}"
        );
        compared += 1;
    }
    assert_eq!(compared, 5, "orders compared");
}

/// The overlap test, at every boundary the corpus has.
#[test]
fn the_overlap_test_is_the_reference_arithmetic() {
    let corpus = corpus();
    let mut compared = 0;
    let mut zero_span_rows = 0;
    for row in rows(&corpus, "intersect") {
        let (first, second, expected) = (entry(row[0]), entry(row[1]), row[2]);
        assert_eq!(
            first.intersects(&second).to_string(),
            expected,
            "intersect {} {}",
            row[0],
            row[1]
        );
        if first.alignment_span == 0 || second.alignment_span == 0 {
            // A zero span never intersects, which an overlap test would not say.
            assert_eq!(expected, "false");
            zero_span_rows += 1;
        }
        compared += 1;
    }
    assert_eq!(compared, 9, "pairs compared");
    assert_eq!(zero_span_rows, 2, "of them with a zero span");

    // Touching is not overlapping: 100..110 and 110..120 do not intersect, and 100..110 and
    // 109..110 do.
    assert!(!entry("0:100:10:0:0:1").intersects(&entry("0:110:10:0:0:1")));
    assert!(entry("0:100:10:0:0:1").intersects(&entry("0:109:1:0:0:1")));

    // An unmapped entry never intersects, not even with itself.
    assert!(!entry("-1:100:10:0:0:1").intersects(&entry("-1:100:10:0:0:1")));
}

/// Querying a list, including the two wildcards.
#[test]
fn a_query_finds_what_the_reference_found() {
    let corpus = corpus();
    let list = entries("0:100:10:0:0:1;0:200:10:0:1:1;0:300:10:0:2:1;1:100:10:0:3:1;-1:1:1:0:4:1");
    let mut compared = 0;
    for row in rows(&corpus, "find") {
        let (label, sequence_id, start, span, expected) = (
            row[0],
            row[1].parse::<i32>().expect("sequence"),
            row[2].parse::<i32>().expect("start"),
            row[3].parse::<i32>().expect("span"),
            row[4],
        );
        assert_eq!(
            show(&find(&list, sequence_id, start, span)),
            expected,
            "find {label}"
        );
        compared += 1;
    }
    assert_eq!(compared, 6, "queries compared");

    // A start below one and a span below one are both wildcards, so they find the same three
    // entries a query over the whole sequence would.
    assert_eq!(find(&list, 0, 0, 10).len(), 3);
    assert_eq!(find(&list, 0, 100, 0).len(), 3);
    // And an unmapped query finds the unmapped entry, because the overlap test that would refuse
    // it is never reached.
    assert_eq!(find(&list, -1, 0, 0).len(), 1);

    for row in rows(&corpus, "leftmost") {
        let expected = row[1];
        let found = if row[0] == "empty" {
            leftmost(&[])
        } else {
            leftmost(&list)
        };
        let shown = found.map(|entry| show(&[entry])).unwrap_or("-".to_string());
        assert_eq!(shown, expected, "leftmost {}", row[0]);
    }
}

/// What building or parsing an entry refuses.
#[test]
fn the_failures_are_the_reference_failures() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "err") {
        let (what, detail, class, message) = (row[0], row[1], row[2], row[3]);
        let error = match what {
            "multi-ref" => CraiEntry::new(-2, 1, 1, 0, 0, 1).expect_err("refused"),
            _ => CraiEntry::parse(&detail.replace(' ', "\t")).expect_err("refused"),
        };
        assert_eq!(error.java_exception(), class, "err {what}");
        assert_eq!(error.message(), message, "err {what}");
        compared += 1;
    }
    assert_eq!(compared, 5, "refusals compared");

    // The line constructor does not make the multi-reference check the other one makes, so a line
    // naming -2 parses where a call naming it does not.
    assert!(matches!(
        CraiEntry::new(-2, 1, 1, 0, 0, 1),
        Err(CraiError::MultiReference)
    ));
    assert_eq!(
        CraiEntry::parse("-2\t1\t1\t0\t0\t1")
            .expect("parsed")
            .sequence_id,
        -2
    );
}
