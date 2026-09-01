//! `QueryInterval` and the chunk arithmetic against the reference's own answers.
//!
//! Both sides print the objects as htsjdk's `toString` does, so a mismatch reads as
//! `0:100-200` rather than as three numbers, and the test parses that form back.
//!
//! While the suite is `golden-pending` the dump is named by `QUERY_DUMP` (decision 0008).

use std::path::Path;

use htsjdk_bam::index::Chunk;
use htsjdk_bam::query::{
    chunks_are_adjacent, chunks_overlap, compare_chunks, display, optimize_chunk_list,
    optimize_intervals, QueryInterval,
};

/// `0:100-200`.
fn interval(text: &str) -> QueryInterval {
    let (reference, rest) = text.split_once(':').expect("ref:start-end");
    let (start, end) = rest.split_once('-').expect("start-end");
    QueryInterval::new(
        reference.parse().expect("a reference index"),
        start.parse().expect("a start"),
        end.parse().expect("an end"),
    )
}

fn intervals(text: &str) -> Vec<QueryInterval> {
    if text == "[]" {
        return Vec::new();
    }
    text.split(',').map(interval).collect()
}

/// `10:0-20:40`, which is `blockAddress:blockOffset` twice.
fn chunk(text: &str) -> Chunk {
    let parts: Vec<&str> = text.split(['-', ':']).collect();
    assert_eq!(
        parts.len(),
        4,
        "blockAddress:blockOffset-blockAddress:blockOffset"
    );
    let vfp = |block: &str, offset: &str| -> u64 {
        (block.parse::<u64>().expect("a block address") << 16)
            | offset.parse::<u64>().expect("an offset")
    };
    Chunk {
        start: vfp(parts[0], parts[1]),
        end: vfp(parts[2], parts[3]),
    }
}

fn chunks(text: &str) -> Vec<Chunk> {
    if text == "[]" {
        return Vec::new();
    }
    // A chunk carries two colons and one dash, so the list separator is the comma alone.
    text.split(',').map(chunk).collect()
}

fn show_chunk(c: &Chunk) -> String {
    format!(
        "{}:{}-{}:{}",
        c.start >> 16,
        c.start & 0xFFFF,
        c.end >> 16,
        c.end & 0xFFFF
    )
}

fn show_chunks(list: &[Chunk]) -> String {
    if list.is_empty() {
        return "[]".to_string();
    }
    list.iter().map(show_chunk).collect::<Vec<_>>().join(",")
}

fn show_intervals(list: &[QueryInterval]) -> String {
    if list.is_empty() {
        return "[]".to_string();
    }
    list.iter().map(display).collect::<Vec<_>>().join(",")
}

#[test]
fn every_answer_matches_the_reference() {
    let golden = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/query.txt.gz");
    let dump = match std::env::var("QUERY_DUMP") {
        Ok(path) => std::fs::read_to_string(path).expect("the dump named by QUERY_DUMP"),
        Err(_) if golden.exists() => {
            panic!("the golden landed: read it here instead of skipping, and drop this branch")
        }
        Err(_) => {
            println!(
                "skipped: the query golden is still pending. Run the suite and point QUERY_DUMP at \
                 tools/conformance/pending/query.QueryDump.txt"
            );
            return;
        }
    };

    let mut rows = 0;
    for line in dump.lines() {
        let fields: Vec<&str> = line.trim().split('\t').collect();
        match fields.as_slice() {
            ["compare", a, b, expected] => {
                let ours = interval(a).compare_to(&interval(b));
                assert_eq!(ours.to_string(), *expected, "compare {a} {b}");
            }
            ["overlaps", a, b, expected] => {
                let ours = interval(a).overlaps(&interval(b));
                assert_eq!(ours.to_string(), *expected, "overlaps {a} {b}");
            }
            ["abuts", a, b, expected] => {
                let ours = interval(a).ends_at_start_of(&interval(b));
                assert_eq!(ours.to_string(), *expected, "abuts {a} {b}");
            }
            ["optimize", input, expected] => {
                let ours = optimize_intervals(&intervals(input));
                assert_eq!(show_intervals(&ours), *expected, "optimize {input}");
            }
            ["chunkcmp", a, b, expected] => {
                let ours = match compare_chunks(&chunk(a), &chunk(b)) {
                    std::cmp::Ordering::Less => -1,
                    std::cmp::Ordering::Equal => 0,
                    std::cmp::Ordering::Greater => 1,
                };
                assert_eq!(ours.to_string(), *expected, "chunkcmp {a} {b}");
            }
            ["chunkovl", a, b, expected] => {
                let ours = chunks_overlap(&chunk(a), &chunk(b));
                assert_eq!(ours.to_string(), *expected, "chunkovl {a} {b}");
            }
            ["chunkadj", a, b, expected] => {
                let ours = chunks_are_adjacent(&chunk(a), &chunk(b));
                assert_eq!(ours.to_string(), *expected, "chunkadj {a} {b}");
            }
            ["chunkopt", minimum, input, expected] => {
                let ours =
                    optimize_chunk_list(&chunks(input), minimum.parse().expect("a virtual offset"));
                assert_eq!(show_chunks(&ours), *expected, "chunkopt {minimum} {input}");
            }
            _ => panic!("unrecognized dump line: {line}"),
        }
        rows += 1;
    }
    assert!(rows > 100, "the dump has every pair, not a sample");
    println!("rows={rows}");
}
