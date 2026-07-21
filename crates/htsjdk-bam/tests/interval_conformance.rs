//! Conformance for `.interval_list` bodies against htsjdk's `IntervalList`.
//!
//! Goldens from `tools/interval-conformance/IntervalListDump.java` in the pinned oracle.
//!
//! The property under test is the **order**. `IntervalList.sorted()` uses
//! `IntervalCoordinateComparator`, which keys on the contig's index in the sequence dictionary,
//! while `Interval.compareTo` keys on the contig *name*. The corpus's dictionary is
//! `chr1 chr2 chr10 chrX`, chosen so those two disagree: a port that sorted naturally puts
//! `chr10` before `chr2` and produces a valid file in the wrong order.

use std::io::Read;

use htsjdk_bam::interval::{Interval, IntervalList};

fn corpus() -> Vec<(String, String)> {
    let p =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/interval_list.txt.gz");
    let f = std::fs::File::open(&p).expect("corpus");
    let mut s = String::new();
    flate2::read::GzDecoder::new(f)
        .read_to_string(&mut s)
        .expect("corpus is gzip");
    s.lines()
        .filter(|l| !l.is_empty())
        .map(|l| match l.split_once('\t') {
            Some((name, value)) => (name.to_string(), unescape(value)),
            // The `empty` case has no payload at all.
            None => (l.to_string(), String::new()),
        })
        .collect()
}

fn unescape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut chars = s.chars();
    while let Some(c) = chars.next() {
        if c != '\\' {
            out.push(c);
            continue;
        }
        match chars.next() {
            Some('n') => out.push('\n'),
            Some('t') => out.push('\t'),
            Some('\\') => out.push('\\'),
            Some(other) => {
                out.push('\\');
                out.push(other);
            }
            None => out.push('\\'),
        }
    }
    out
}

fn dictionary() -> Vec<String> {
    ["chr1", "chr2", "chr10", "chrX"]
        .iter()
        .map(|s| s.to_string())
        .collect()
}

fn iv(c: &str, s: i32, e: i32) -> Interval {
    Interval::new(c, s, e)
}

fn named(c: &str, s: i32, e: i32, neg: bool, n: Option<&str>) -> Interval {
    Interval::with_strand_and_name(c, s, e, neg, n)
}

fn list(intervals: Vec<Interval>) -> IntervalList {
    let mut l = IntervalList::new(dictionary());
    l.intervals = intervals;
    l
}

fn cases() -> Vec<(String, String)> {
    let mut out = Vec::new();
    let mut push = |name: &str, l: &IntervalList| out.push((name.to_string(), l.write_body()));

    push("one_interval", &list(vec![iv("chr1", 100, 200)]));

    let unsorted = list(vec![
        iv("chr10", 1, 10),
        iv("chrX", 1, 10),
        iv("chr2", 1, 10),
        iv("chr1", 1, 10),
    ]);
    push("unsorted_as_given", &unsorted);
    push("unsorted_sorted", &unsorted.sorted());

    let named_list = list(vec![
        named("chr1", 1, 10, false, Some("plus_named")),
        named("chr1", 20, 30, true, Some("minus_named")),
        named("chr1", 40, 50, false, None),
        named("chr1", 60, 70, true, None),
    ]);
    push("strands_and_names", &named_list.sorted());

    let tail = list(vec![
        named("chr1", 1, 10, true, Some("zzz")),
        named("chr1", 1, 10, false, Some("zzz")),
        named("chr1", 1, 10, false, Some("aaa")),
        named("chr1", 1, 10, false, None),
    ]);
    push("comparator_tail", &tail.sorted());

    let overlapping = list(vec![
        named("chr1", 1, 10, false, Some("first")),
        named("chr1", 5, 20, false, Some("second")),
        named("chr1", 21, 30, false, Some("abutting")),
        named("chr1", 32, 40, false, Some("separated")),
    ]);
    // `uniqued()` with no argument is `uniqued(true)`, so the no-argument form concatenates.
    // The two golden cases are identical for exactly that reason.
    push("uniqued", &overlapping.uniqued(true));
    push("uniqued_concatenated", &overlapping.uniqued(true));
    push("sorted_not_uniqued", &overlapping.sorted());

    let to_pad = list(vec![
        named("chr1", 5, 10, false, Some("near_start")),
        named("chr2", 500, 600, false, Some("middle")),
    ]);
    push("padded", &to_pad.padded(100, 50, |_| 100_000));

    push(
        "single_base",
        &list(vec![named("chr1", 42, 42, false, Some("point"))]),
    );
    push("empty", &list(vec![]));

    out
}

#[test]
fn every_interval_list_body_matches_htsjdks() {
    let golden = corpus();
    let ours = cases();
    assert_eq!(golden.len(), ours.len(), "case lists differ in length");

    let mut mismatches = Vec::new();
    for (i, (name, text)) in ours.iter().enumerate() {
        let (gname, gtext) = &golden[i];
        assert_eq!(gname, name, "case {i}: the lists have drifted");
        if gtext != text {
            mismatches.push(format!("{name}\n  htsjdk: {gtext:?}\n  ours  : {text:?}"));
        }
    }
    assert!(
        mismatches.is_empty(),
        "{} of {} cases diverge:\n{}",
        mismatches.len(),
        ours.len(),
        mismatches.join("\n")
    );
}

/// The finding, asserted on the goldens: sorting follows the dictionary, not the alphabet.
#[test]
fn the_sort_follows_the_dictionary_and_not_the_alphabet() {
    let golden = corpus();
    let (_, sorted) = golden
        .iter()
        .find(|(n, _)| n == "unsorted_sorted")
        .expect("unsorted_sorted");
    let contigs: Vec<&str> = sorted
        .lines()
        .filter_map(|l| l.split('\t').next())
        .collect();
    assert_eq!(contigs, ["chr1", "chr2", "chr10", "chrX"]);
    let mut alphabetical = contigs.clone();
    alphabetical.sort();
    assert_ne!(
        contigs, alphabetical,
        "the corpus must distinguish the two orderings, or it proves nothing"
    );
}

/// `uniqued()` takes the concatenating branch, so the two golden cases are byte-identical.
/// If htsjdk ever changed the no-argument default, this would fail rather than pass quietly.
#[test]
fn the_no_argument_uniqued_concatenates_names() {
    let golden = corpus();
    let get = |n: &str| {
        golden
            .iter()
            .find(|(name, _)| name == n)
            .map(|(_, v)| v.clone())
            .expect(n)
    };
    assert_eq!(get("uniqued"), get("uniqued_concatenated"));
    assert!(get("uniqued").contains("first|second|abutting"));
}
