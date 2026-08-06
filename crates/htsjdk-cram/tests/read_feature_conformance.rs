//! Conformance for the CRAM read features, against `CRAMRecordReadFeatures`.
//!
//! Goldens from `tools/cram-conformance/CramReadFeatureDump.java` in the pinned oracle.
//!
//! The rows that justify the suite:
//!
//! ```text
//! case     x-over-matching-bases  8X      ACGTACGT  IIIIIIII  1   0
//! feature  insertion-of-three     0  i  3  base=T
//! feature  soft-clip-of-three     0  S  1  bases=TTT
//! case     n-past-the-end         4M      NNNN      IIII      22  3
//! err      empty-quals-non-acgtn  4M  ACGM  <empty array>  1  ArrayIndexOutOfBoundsException  Index 3 out of bounds for length 0
//! ```
//!
//! An `X` over bases that match emits nothing, so the cigar's own claim is not consulted. An
//! insertion of three bases is three features and a soft clip of three is one, decided five lines
//! apart in the same loop. Four `N`s placed at the end of the reference produce three features,
//! because the fourth is compared against the `N` the reference reads past its end. And the
//! missing-quality test is an identity test, so an empty array that is not the singleton is
//! indexed like any other.

use std::io::Read;

use htsjdk_bam::cigar::{Cigar, CigarElement, Op};
use htsjdk_cram::read_features::{
    alignment_end, create_read_features, Qualities, ReadFeature, ReadFeatureError,
    MISSING_QUALITY_SCORE,
};

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/cram_read_features.txt.gz");
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

fn reference(corpus: &str) -> Vec<u8> {
    rows(corpus, "reference")[0][0].as_bytes().to_vec()
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

/// `SEQ="*"` is the empty array htsjdk fills with `N`s.
fn bases(text: &str) -> Vec<u8> {
    if text == "*" {
        Vec::new()
    } else {
        text.as_bytes().to_vec()
    }
}

/// `QUAL="*"` leaves the NULL_QUALS singleton; the dump's `<empty array>` is a different array
/// with the same contents, which is the whole point of the distinction.
fn qualities(text: &str) -> Option<Vec<u8>> {
    match text {
        "*" => None,
        "<empty array>" => Some(Vec::new()),
        _ => Some(text.bytes().map(|b| b - 33).collect()),
    }
}

fn quals_of(parsed: &Option<Vec<u8>>) -> Qualities<'_> {
    match parsed {
        None => Qualities::Missing,
        Some(values) => Qualities::Present(values),
    }
}

/// The dump's `payload`: everything the feature carries beyond its position.
fn payload(feature: &ReadFeature) -> String {
    match feature {
        ReadFeature::Substitution {
            base,
            reference_base,
            code,
            ..
        } => format!(
            "base={} ref={} code={code}",
            *base as char, *reference_base as char
        ),
        ReadFeature::ReadBase { base, quality, .. } => {
            format!("base={} quality={quality}", *base as char)
        }
        ReadFeature::InsertBase { base, .. } => format!("base={}", *base as char),
        ReadFeature::SoftClip { sequence, .. } | ReadFeature::Insertion { sequence, .. } => {
            format!("bases={}", String::from_utf8_lossy(sequence))
        }
        ReadFeature::Bases { bases, .. } => format!("bases={}", String::from_utf8_lossy(bases)),
        ReadFeature::Scores { scores, .. } => format!("scores={}", String::from_utf8_lossy(scores)),
        ReadFeature::BaseQualityScore { quality, .. } => format!("quality={quality}"),
        ReadFeature::Deletion { length, .. }
        | ReadFeature::RefSkip { length, .. }
        | ReadFeature::Padding { length, .. }
        | ReadFeature::HardClip { length, .. } => format!("length={length}"),
    }
}

/// The twelve operator letters, taken from the reference's own classes.
#[test]
fn the_operators_are_the_reference_operators() {
    let corpus = corpus();
    let letter = |class: &str| -> u8 {
        rows(&corpus, "op")
            .into_iter()
            .find(|row| row[0] == class)
            .unwrap_or_else(|| panic!("{class}"))[1]
            .as_bytes()[0]
    };
    let position = 1;
    let cases: Vec<(&str, ReadFeature)> = vec![
        (
            "BaseQualityScore",
            ReadFeature::BaseQualityScore {
                position,
                quality: 0,
            },
        ),
        (
            "Bases",
            ReadFeature::Bases {
                position,
                bases: vec![],
            },
        ),
        (
            "Deletion",
            ReadFeature::Deletion {
                position,
                length: 1,
            },
        ),
        (
            "HardClip",
            ReadFeature::HardClip {
                position,
                length: 1,
            },
        ),
        (
            "InsertBase",
            ReadFeature::InsertBase {
                position,
                base: b'A',
            },
        ),
        (
            "Insertion",
            ReadFeature::Insertion {
                position,
                sequence: vec![],
            },
        ),
        (
            "Padding",
            ReadFeature::Padding {
                position,
                length: 1,
            },
        ),
        (
            "ReadBase",
            ReadFeature::ReadBase {
                position,
                base: b'A',
                quality: 0,
            },
        ),
        (
            "RefSkip",
            ReadFeature::RefSkip {
                position,
                length: 1,
            },
        ),
        (
            "Scores",
            ReadFeature::Scores {
                position,
                scores: vec![],
            },
        ),
        (
            "SoftClip",
            ReadFeature::SoftClip {
                position,
                sequence: vec![],
            },
        ),
        (
            "Substitution",
            ReadFeature::Substitution {
                position,
                base: b'A',
                reference_base: b'C',
                code: -1,
            },
        ),
    ];
    assert_eq!(cases.len(), 12, "every class in the package");
    for (class, feature) in cases {
        assert_eq!(feature.operator(), letter(class), "{class}");
    }

    let consts = rows(&corpus, "consts");
    assert_eq!(consts[0][0], MISSING_QUALITY_SCORE.to_string());
}

/// Every case the reference walked produces the same features here, in the same order.
#[test]
fn every_case_produces_the_features_the_reference_produced() {
    let corpus = corpus();
    let reference = reference(&corpus);
    let features_by_case: Vec<Vec<&str>> = rows(&corpus, "feature");
    let ends: Vec<Vec<&str>> = rows(&corpus, "alignend");

    let mut compared = 0;
    for row in rows(&corpus, "case") {
        let (label, cigar_text, base_text, qual_text, start, count) =
            (row[0], row[1], row[2], row[3], row[4], row[5]);
        let start: i32 = start.parse().expect("start");
        let read_bases = bases(base_text);
        let parsed = qualities(qual_text);

        let mine = create_read_features(
            &cigar(cigar_text),
            start,
            &read_bases,
            quals_of(&parsed),
            &reference,
        )
        .unwrap_or_else(|e| panic!("{label}: {}", e.message()));

        assert_eq!(mine.len().to_string(), count, "{label}: feature count");

        let expected: Vec<&Vec<&str>> = features_by_case
            .iter()
            .filter(|feature| feature[0] == label)
            .collect();
        assert_eq!(expected.len(), mine.len(), "{label}: rows");
        for (index, feature) in mine.iter().enumerate() {
            let row = expected[index];
            assert_eq!(row[1], index.to_string(), "{label}: index");
            assert_eq!(
                (feature.operator() as char).to_string(),
                row[2],
                "{label}/{index}: operator"
            );
            assert_eq!(
                feature.position().to_string(),
                row[3],
                "{label}/{index}: position"
            );
            assert_eq!(payload(feature), row[4], "{label}/{index}: payload");
        }

        let end = ends
            .iter()
            .find(|row| row[0] == label)
            .unwrap_or_else(|| panic!("{label}: no alignment end"));
        let read_length: i32 = end[2].parse().expect("read length");
        assert_eq!(
            alignment_end(&mine, start, read_length).to_string(),
            end[3],
            "{label}: alignment end"
        );
        compared += 1;
    }
    assert_eq!(compared, 20, "cases compared");
}

/// Nothing is stored for a base that matches, and the cigar's claim is not consulted.
#[test]
fn a_match_is_stored_as_nothing_whatever_the_cigar_says() {
    let corpus = corpus();
    let count = |label: &str| -> i32 {
        rows(&corpus, "case")
            .into_iter()
            .find(|row| row[0] == label)
            .unwrap_or_else(|| panic!("{label}"))[5]
            .parse()
            .expect("count")
    };
    assert_eq!(count("perfect-match"), 0);
    assert_eq!(count("x-over-matching-bases"), 0, "an X that matches");
    assert_eq!(count("eq-over-mismatching-bases"), 6, "an = that does not");
}

/// An insertion of n bases is n features and a soft clip of n is one, decided five lines apart.
#[test]
fn an_insertion_is_split_and_a_soft_clip_is_not() {
    let corpus = corpus();
    let of = |label: &str| -> Vec<Vec<&str>> {
        rows(&corpus, "feature")
            .into_iter()
            .filter(|row| row[0] == label)
            .collect()
    };
    let inserted: Vec<Vec<&str>> = of("insertion-of-three")
        .into_iter()
        .filter(|row| row[2] == "i")
        .collect();
    assert_eq!(inserted.len(), 3, "one feature per inserted base");
    assert_eq!(
        inserted.iter().map(|row| row[3]).collect::<Vec<_>>(),
        ["3", "4", "5"],
        "consecutive one-based positions"
    );

    let clipped = of("soft-clip-of-three");
    assert_eq!(clipped[0][2], "S");
    assert_eq!(clipped[0][4], "bases=TTT", "one feature for all three");
}

/// The two ways a mismatch is stored, and the quality score the second one carries twice.
#[test]
fn a_mismatch_splits_on_the_alphabet_and_not_on_the_cigar() {
    let corpus = corpus();
    let of = |label: &str| -> Vec<Vec<&str>> {
        rows(&corpus, "feature")
            .into_iter()
            .filter(|row| row[0] == label)
            .collect()
    };
    // The read base is outside ACGTN.
    let read_base = of("read-base-not-acgtn");
    assert_eq!(read_base[0][2], "B");
    assert_eq!(read_base[0][4], "base=M quality=40");
    // The reference base is, which is the same branch.
    let reference_base = of("reference-not-acgtn");
    assert_eq!(reference_base.len(), 4);
    assert!(reference_base.iter().all(|row| row[2] == "B"));
    // And with no qualities at all it is the missing score rather than a lookup.
    let missing = of("null-quals-non-acgtn");
    assert_eq!(
        missing[0][4],
        format!("base=M quality={MISSING_QUALITY_SCORE}")
    );
}

/// Both ways the reference walks off its own array, with the exception it raises.
#[test]
fn the_two_out_of_bounds_paths_are_the_reference_errors() {
    let corpus = corpus();
    let reference = reference(&corpus);
    let mut compared = 0;
    for row in rows(&corpus, "err") {
        let (label, cigar_text, base_text, qual_text, start, class, message) =
            (row[0], row[1], row[2], row[3], row[4], row[5], row[6]);
        let start: i32 = start.parse().expect("start");
        let read_bases = bases(base_text);
        let parsed = qualities(qual_text);

        let error = create_read_features(
            &cigar(cigar_text),
            start,
            &read_bases,
            quals_of(&parsed),
            &reference,
        )
        .expect_err(label);
        assert_eq!(class, "ArrayIndexOutOfBoundsException", "{label}");
        assert!(
            matches!(error, ReadFeatureError::IndexOutOfBounds { .. }),
            "{label}: {error:?}"
        );
        assert_eq!(error.message(), message, "{label}");
        compared += 1;
    }
    assert_eq!(compared, 2, "out-of-bounds paths compared");
}
