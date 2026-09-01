//! `DuplicateScoringStrategy` against the reference's own scores and comparisons.
//!
//! The corpus is rebuilt by index rather than read from a fixture, so this file and
//! `tools/duplicate-scoring-conformance/DuplicateScoringDump.java` have to agree about what record
//! `i` is; each says so in the same words.
//!
//! The golden is committed and re-derived by the `duplicate-scoring` suite on every run; the dump can
//! still be overridden with an environment variable while a harness change is being checked.

use std::io::Read;
use std::path::Path;

use htsjdk_bam::cigar::{Cigar, CigarElement, Op};
use htsjdk_bam::duplicate_scoring::{compare, compute_duplicate_score, ScoringStrategy};
use htsjdk_bam::record::BamRecord;

const READ_PAIRED: u16 = 0x1;
const READ_UNMAPPED: u16 = 0x4;
const MATE_UNMAPPED: u16 = 0x8;
const FIRST_OF_PAIR: u16 = 0x40;
const SECOND_OF_PAIR: u16 = 0x80;
const VENDOR_FAILED: u16 = 0x200;

fn record(i: usize) -> BamRecord {
    let length = if i == 7 { 5000 } else { 8 };
    let bases: Vec<u8> = (0..length).map(|b| b"ACGT"[b % 4]).collect();
    let quals: Vec<u8> = (0..length)
        .map(|b| {
            if i.is_multiple_of(2) && b.is_multiple_of(2) {
                10
            } else {
                20 + (i % 5) as u8
            }
        })
        .collect();
    let mut record = BamRecord {
        read_name: format!("read{i}"),
        read_bases: bases,
        base_qualities: quals,
        ..Default::default()
    };
    if i % 4 == 3 {
        record.flags |= READ_UNMAPPED;
        record.reference_index = -1;
        record.alignment_start = 0;
    } else {
        record.reference_index = 0;
        record.alignment_start = 100 + i as i32;
        record.mapping_quality = 60;
        record.cigar = Cigar::new(vec![CigarElement {
            length: length as u32,
            op: Op::M,
        }]);
    }
    if i % 3 != 2 {
        record.flags |= READ_PAIRED | MATE_UNMAPPED;
        record.flags |= if i.is_multiple_of(3) {
            FIRST_OF_PAIR
        } else {
            SECOND_OF_PAIR
        };
    }
    if i % 5 == 4 {
        record.flags |= VENDOR_FAILED;
    }
    record
}

fn strategy_for(name: &str) -> ScoringStrategy {
    match name {
        "SUM_OF_BASE_QUALITIES" => ScoringStrategy::SumOfBaseQualities,
        "TOTAL_MAPPED_REFERENCE_LENGTH" => ScoringStrategy::TotalMappedReferenceLength,
        "RANDOM" => ScoringStrategy::Random,
        other => panic!("unknown strategy {other}"),
    }
}

#[test]
fn every_score_and_comparison_matches_the_reference() {
    // The golden was produced by the pinned container on real x86-64 and is re-derived on every
    // run; `DUPLICATE_SCORING_DUMP` still overrides it, which is how a local run checks a change to the
    // harness before CI does.
    let dump = match std::env::var("DUPLICATE_SCORING_DUMP") {
        Ok(path) => {
            std::fs::read_to_string(path).expect("the dump named by DUPLICATE_SCORING_DUMP")
        }
        Err(_) => {
            let golden =
                Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/duplicate_scoring.txt.gz");
            let file = std::fs::File::open(&golden).expect("the committed golden");
            let mut text = String::new();
            flate2::read::GzDecoder::new(file)
                .read_to_string(&mut text)
                .expect("the golden decompresses");
            text
        }
    };

    let (mut scores, mut comparisons) = (0, 0);
    for line in dump.lines() {
        let fields: Vec<&str> = line.trim().split('\t').collect();
        match fields.as_slice() {
            ["score", strategy, index, expected] => {
                let i: usize = index.parse().expect("a record index");
                let ours = compute_duplicate_score(&record(i), strategy_for(strategy), false, None);
                assert_eq!(ours.to_string(), *expected, "score {strategy} {i}");
                scores += 1;
            }
            ["compare", strategy, first, second, expected] => {
                let (a, b): (usize, usize) = (
                    first.parse().expect("an index"),
                    second.parse().expect("an index"),
                );
                let ours = compare(
                    &record(a),
                    &record(b),
                    strategy_for(strategy),
                    false,
                    None,
                    None,
                );
                // htsjdk's compare returns a difference of scores, not a sign; the port returns the
                // same number, so this is an equality and not a comparison of signs.
                assert_eq!(ours.to_string(), *expected, "compare {strategy} {a},{b}");
                comparisons += 1;
            }
            _ => panic!("unrecognized dump line: {line}"),
        }
    }
    assert!(scores > 0 && comparisons > 0, "both families ran");
    println!("score={scores} compare={comparisons}");
}
