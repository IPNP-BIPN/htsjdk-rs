//! The record filters against the reference's own answers, in both forms of `filterOut`.
//!
//! The corpus is built by index on both sides rather than read from a fixture, so the Rust half is
//! the same twelve records the Java half made: name `read<i>`, unmapped when `i` is odd, and an
//! `RG` of `rg1`, `rg2` or nothing as `i % 3` decides.
//!
//! While the suite is `golden-pending` the dump is named by `FILTER_DUMP`; the committed corpus may
//! only come from the pinned container on real x86-64 (decision 0008).

use std::collections::HashSet;
use std::path::Path;

use htsjdk_bam::filter::{AlignedFilter, ReadNameFilter, SamRecordFilter, TagFilter};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};

const READ_UNMAPPED: u16 = 0x4;

/// Record `i` of the harness's corpus.
fn record(i: usize) -> BamRecord {
    let mut record = BamRecord {
        read_name: format!("read{i}"),
        read_bases: b"ACGT".to_vec(),
        base_qualities: vec![30; 4],
        ..Default::default()
    };
    if i % 2 == 1 {
        record.flags = READ_UNMAPPED;
        record.reference_index = -1;
        record.alignment_start = 0;
    } else {
        record.reference_index = 0;
        record.alignment_start = 100 + i as i32;
        record.mapping_quality = 60;
    }
    match i % 3 {
        0 => {
            record
                .tags
                .insert(Tag::new(b"RG"), TagValue::Str("rg1".to_string()));
        }
        1 => {
            record
                .tags
                .insert(Tag::new(b"RG"), TagValue::Str("rg2".to_string()));
        }
        _ => {}
    }
    record
}

fn filter_for(label: &str) -> Box<dyn SamRecordFilter> {
    let names: HashSet<String> = ["read0", "read3", "read4", "read11"]
        .iter()
        .map(|s| s.to_string())
        .collect();
    let values = vec![TagValue::Str("rg1".to_string())];
    match label {
        "aligned_include" => Box::new(AlignedFilter {
            include_aligned: true,
        }),
        "aligned_exclude" => Box::new(AlignedFilter {
            include_aligned: false,
        }),
        "name_include" => Box::new(ReadNameFilter {
            names,
            include_reads: true,
        }),
        "name_exclude" => Box::new(ReadNameFilter {
            names,
            include_reads: false,
        }),
        "tag_include" => Box::new(TagFilter {
            tag: Tag::new(b"RG"),
            values,
            include_reads: true,
        }),
        "tag_exclude" => Box::new(TagFilter {
            tag: Tag::new(b"RG"),
            values,
            include_reads: false,
        }),
        other => panic!("unknown filter {other}"),
    }
}

#[test]
fn every_decision_matches_the_reference() {
    let golden = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/filter.txt.gz");
    let dump = match std::env::var("FILTER_DUMP") {
        Ok(path) => std::fs::read_to_string(path).expect("the dump named by FILTER_DUMP"),
        Err(_) if golden.exists() => {
            panic!("the golden landed: read it here instead of skipping, and drop this branch")
        }
        Err(_) => {
            println!(
                "skipped: the filter golden is still pending. Run the suite and point FILTER_DUMP \
                 at tools/conformance/pending/filter.FilterDump.txt"
            );
            return;
        }
    };

    let (mut singles, mut pairs) = (0, 0);
    for line in dump.lines() {
        let fields: Vec<&str> = line.trim().split('\t').collect();
        match fields.as_slice() {
            ["single", label, index, expected] => {
                let i: usize = index.parse().expect("a record index");
                let ours = filter_for(label).filter_out(&record(i));
                assert_eq!(ours.to_string(), *expected, "single {label} {i}");
                singles += 1;
            }
            ["pair", label, first, second, expected] => {
                let (a, b): (usize, usize) = (
                    first.parse().expect("an index"),
                    second.parse().expect("an index"),
                );
                let ours = filter_for(label).filter_out_pair(&record(a), &record(b));
                assert_eq!(ours.to_string(), *expected, "pair {label} {a},{b}");
                pairs += 1;
            }
            _ => panic!("unrecognized dump line: {line}"),
        }
    }
    assert!(singles > 0 && pairs > 0, "both forms ran");
    println!("single={singles} pair={pairs}");
}
