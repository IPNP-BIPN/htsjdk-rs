//! Conformance for the cigar rebuilt from read features, against `getCigarForReadFeatures`.
//!
//! Goldens from `tools/cram-conformance/CramRecordCigarDump.java` in the pinned oracle.
//!
//! The rows that justify the suite:
//!
//! ```text
//! case       deletion-at-one  8  D@1 len=2         2D8M
//! case       bases-only       8  b@1 bases=ACGT    8M
//! case       zero-read-length 0  X@4 base=T ref=A  4M
//! roundtrip  x-operator  8X  ACGTACGT  1  8M  changed
//! ```
//!
//! The cigar is not stored anywhere: the matches come back as the gaps between feature positions.
//! A deletion at the first position leaves the whole read after it, a `Bases` feature carrying
//! read bases is dropped on the floor, a read length of 0 takes the accumulated length instead,
//! and a record written with `8X` comes back as `8M` because a substitution and a match are the
//! same operator here.

use std::io::Read;

use htsjdk_bam::cigar::{Cigar, CigarElement, Op};
use htsjdk_cram::read_features::{
    cigar_for_read_features, create_read_features, Qualities, ReadFeature, NO_CODE,
};

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/cram_record_cigar.txt.gz");
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

fn cigar(text: &str) -> Cigar {
    let mut elements = Vec::new();
    let mut length = 0u32;
    for byte in text.bytes() {
        if byte.is_ascii_digit() {
            length = length * 10 + u32::from(byte - b'0');
        } else {
            let op = match byte {
                b'M' => Op::M,
                b'I' => Op::I,
                b'D' => Op::D,
                b'N' => Op::N,
                b'S' => Op::S,
                b'H' => Op::H,
                b'P' => Op::P,
                b'=' => Op::Eq,
                b'X' => Op::X,
                other => panic!("cigar operator {}", other as char),
            };
            elements.push(CigarElement { length, op });
            length = 0;
        }
    }
    Cigar::new(elements)
}

/// The dump's `describe`, read back: `X@4 base=T ref=A`, `D@5 len=2`, `S@1 bases=TTT`.
fn parse_features(text: &str) -> Vec<ReadFeature> {
    if text == "-" {
        return Vec::new();
    }
    text.split(',')
        .map(|entry| {
            let mut parts = entry.split(' ');
            let head = parts.next().expect("operator");
            let (operator, position) = head.split_at(1);
            let position: i32 = position.trim_start_matches('@').parse().expect("position");
            let field = |name: &str| -> Option<String> {
                entry
                    .split(' ')
                    .find_map(|part| part.strip_prefix(&format!("{name}=")))
                    .map(|value| value.to_string())
            };
            let length = || -> i32 { field("len").expect("len").parse().expect("length") };
            let bytes = |name: &str| -> Vec<u8> { field(name).expect(name).into_bytes() };
            let base = |name: &str| -> u8 { field(name).expect(name).as_bytes()[0] };
            let quality = || -> i8 { field("quality").expect("quality").parse().expect("quality") };

            match operator {
                "X" => ReadFeature::Substitution {
                    position,
                    base: base("base"),
                    reference_base: base("ref"),
                    code: NO_CODE,
                },
                "B" => ReadFeature::ReadBase {
                    position,
                    base: base("base"),
                    quality: quality(),
                },
                "i" => ReadFeature::InsertBase {
                    position,
                    base: base("base"),
                },
                "I" => ReadFeature::Insertion {
                    position,
                    sequence: bytes("bases"),
                },
                "S" => ReadFeature::SoftClip {
                    position,
                    sequence: bytes("bases"),
                },
                "b" => ReadFeature::Bases {
                    position,
                    bases: bytes("bases"),
                },
                "q" => ReadFeature::Scores {
                    position,
                    scores: bytes("scores"),
                },
                "Q" => ReadFeature::BaseQualityScore {
                    position,
                    quality: quality(),
                },
                "D" => ReadFeature::Deletion {
                    position,
                    length: length(),
                },
                "N" => ReadFeature::RefSkip {
                    position,
                    length: length(),
                },
                "P" => ReadFeature::Padding {
                    position,
                    length: length(),
                },
                "H" => ReadFeature::HardClip {
                    position,
                    length: length(),
                },
                other => panic!("operator {other}"),
            }
        })
        .collect()
}

/// Every hand-built feature list rebuilds into the cigar the reference rebuilt.
#[test]
fn every_feature_list_rebuilds_the_reference_cigar() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "case") {
        let (label, read_length, features, expected) = (row[0], row[1], row[2], row[3]);
        let read_length: i32 = read_length.parse().expect("read length");
        let mine = cigar_for_read_features(&parse_features(features), read_length);
        assert_eq!(mine.to_text(), expected, "{label}");
        compared += 1;
    }
    assert_eq!(compared, 26, "feature lists compared");
}

/// The matches are the gaps: nothing in a list of substitutions says a base matched, and the
/// whole read still comes back as one M.
#[test]
fn the_matches_are_the_gaps_between_the_features() {
    let corpus = corpus();
    let of = |label: &str| -> String {
        rows(&corpus, "case")
            .into_iter()
            .find(|row| row[0] == label)
            .unwrap_or_else(|| panic!("{label}"))[3]
            .to_string()
    };
    for label in [
        "empty",
        "one-substitution",
        "substitution-at-one",
        "substitution-at-end",
        "two-substitutions",
        "adjacent-substitutions",
        "read-base-is-also-m",
    ] {
        assert_eq!(of(label), "8M", "{label}");
    }
}

/// A feature that consumes no read bases winds the read cursor back, which is what keeps the
/// trailing match at the right length.
#[test]
fn a_reference_only_operator_winds_the_read_cursor_back() {
    let corpus = corpus();
    let of = |label: &str| -> String {
        rows(&corpus, "case")
            .into_iter()
            .find(|row| row[0] == label)
            .unwrap_or_else(|| panic!("{label}"))[3]
            .to_string()
    };
    assert_eq!(of("deletion"), "4M2D4M");
    assert_eq!(of("ref-skip"), "4M2N4M");
    assert_eq!(of("padding"), "4M2P4M");
    // At the first position there is nothing before it, and the whole read follows.
    assert_eq!(of("deletion-at-one"), "2D8M");
}

/// The three features the switch does not name contribute nothing, including the one that
/// carries read bases.
#[test]
fn the_features_the_switch_ignores_contribute_nothing() {
    let corpus = corpus();
    let of = |label: &str| -> String {
        rows(&corpus, "case")
            .into_iter()
            .find(|row| row[0] == label)
            .unwrap_or_else(|| panic!("{label}"))[3]
            .to_string()
    };
    assert_eq!(of("base-quality-score-only"), "8M");
    assert_eq!(of("scores-only"), "8M");
    assert_eq!(of("bases-only"), "8M", "and Bases carries read bases");
    assert_eq!(of("bases-then-substitution"), "8M");
}

/// A record's own cigar, through its features and back, and the only shapes that change.
#[test]
fn the_round_trip_loses_only_the_x_and_the_equals() {
    let corpus = corpus();
    let reference = rows(&corpus, "reference")[0][0].as_bytes().to_vec();
    let mut compared = 0;
    let mut changed = Vec::new();
    for row in rows(&corpus, "roundtrip") {
        let (label, text, bases, start, expected, stability) =
            (row[0], row[1], row[2], row[3], row[4], row[5]);
        let start: i32 = start.parse().expect("start");
        let qualities = vec![40u8; bases.len()];
        let features = create_read_features(
            &cigar(text),
            start,
            bases.as_bytes(),
            Qualities::Present(&qualities),
            &reference,
        )
        .unwrap_or_else(|e| panic!("{label}: {}", e.message()));

        let read_length = cigar(text).read_length() as i32;
        let rebuilt = cigar_for_read_features(&features, read_length);
        assert_eq!(rebuilt.to_text(), expected, "{label}");
        assert_eq!(
            if rebuilt.to_text() == text {
                "same"
            } else {
                "changed"
            },
            stability,
            "{label}: stability"
        );
        if stability == "changed" {
            changed.push(label);
        }
        compared += 1;
    }
    assert_eq!(compared, 13, "round trips compared");
    assert_eq!(
        changed,
        ["x-operator", "eq-operator"],
        "only the two operators the features cannot carry"
    );
}
