//! `QueryInterval` and the chunk arithmetic against the reference's own answers.
//!
//! Both sides print the objects as htsjdk's `toString` does, so a mismatch reads as
//! `0:100-200` rather than as three numbers, and the test parses that form back.
//!
//! The golden is committed and re-derived by the `query` suite on every run; the dump can
//! still be overridden with an environment variable while a harness change is being checked.

use std::io::Read;
use std::path::Path;

use htsjdk_bam::cigar::{Cigar, CigarElement, Op};
use htsjdk_bam::index::Chunk;
use htsjdk_bam::iterator_filter::{
    compare_interval_to_record, FilteringIteratorState, IntervalComparison, MultipleIntervalsFilter,
};
use htsjdk_bam::query::{
    chunks_are_adjacent, chunks_overlap, compare_chunks, display, optimize_chunk_list,
    optimize_intervals, region_to_bins, QueryInterval,
};
use htsjdk_bam::record::BamRecord;

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

const READ_UNMAPPED: u16 = 0x4;

/// Record `i` of the filter corpus, as `QueryDump.filterRecord` builds it.
fn filter_record(i: usize) -> BamRecord {
    let mut record = BamRecord {
        read_name: format!("read{i}"),
        reference_index: (i % 2) as i32,
        alignment_start: 100 * (i as i32 + 1),
        read_bases: b"ACGTACGTAC".to_vec(),
        base_qualities: vec![30; 10],
        ..Default::default()
    };
    if i % 3 == 2 {
        record.flags = READ_UNMAPPED;
    } else {
        record.cigar = Cigar::new(vec![CigarElement {
            length: 10,
            op: Op::M,
        }]);
        record.mapping_quality = 60;
    }
    record
}

fn comparison_name(comparison: IntervalComparison) -> &'static str {
    match comparison {
        IntervalComparison::Before => "BEFORE",
        IntervalComparison::After => "AFTER",
        IntervalComparison::Overlapping => "OVERLAPPING",
        IntervalComparison::Contained => "CONTAINED",
    }
}

fn state_name(state: FilteringIteratorState) -> &'static str {
    match state {
        FilteringIteratorState::MatchesFilter => "MATCHES_FILTER",
        FilteringIteratorState::StopIteration => "STOP_ITERATION",
        FilteringIteratorState::ContinueIteration => "CONTINUE_ITERATION",
    }
}

#[test]
fn every_answer_matches_the_reference() {
    // The golden was produced by the pinned container on real x86-64 and is re-derived on every
    // run; `QUERY_DUMP` still overrides it, which is how a harness change is checked before CI
    // sees it.
    let dump = match std::env::var("QUERY_DUMP") {
        Ok(path) => std::fs::read_to_string(path).expect("the dump named by QUERY_DUMP"),
        Err(_) => {
            let golden = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/query.txt.gz");
            let file = std::fs::File::open(&golden).expect("the committed golden");
            let mut text = String::new();
            flate2::read::GzDecoder::new(file)
                .read_to_string(&mut text)
                .expect("the golden decompresses");
            text
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
            ["cmprec", interval_text, index, expected] => {
                let ours = compare_interval_to_record(
                    &interval(interval_text),
                    &filter_record(index.parse().expect("a record index")),
                );
                assert_eq!(
                    comparison_name(ours),
                    *expected,
                    "cmprec {interval_text} {index}"
                );
            }
            ["filter", contained, index, expected] => {
                // The filter is stateful, so the run is replayed from the top for each row rather
                // than kept between rows: the dump records the answer at that point in the walk.
                let intervals = vec![
                    QueryInterval::new(0, 100, 200),
                    QueryInterval::new(0, 500, 900),
                    QueryInterval::new(1, 100, 100000),
                ];
                let upto: usize = index.parse().expect("a record index");
                let mut filter =
                    MultipleIntervalsFilter::new(intervals, contained.parse().expect("a flag"));
                let mut state = FilteringIteratorState::ContinueIteration;
                for i in 0..=upto {
                    state = filter.compare_to_filter(&filter_record(i));
                }
                assert_eq!(state_name(state), *expected, "filter {contained} {index}");
            }
            ["bins", start, end, expected] => {
                let ours = match region_to_bins(
                    start.parse().expect("a start"),
                    end.parse().expect("an end"),
                ) {
                    None => "null".to_string(),
                    Some(bins) if bins.is_empty() => "[]".to_string(),
                    Some(bins) => bins
                        .iter()
                        .map(|b| b.to_string())
                        .collect::<Vec<_>>()
                        .join(","),
                };
                assert_eq!(ours, *expected, "bins {start} {end}");
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
