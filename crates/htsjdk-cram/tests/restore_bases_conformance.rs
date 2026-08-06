//! Conformance for the bases restored from read features, against `restoreReadBases`.
//!
//! Goldens from `tools/cram-conformance/CramRestoreBasesDump.java` in the pinned oracle.
//!
//! The rows that justify the suite:
//!
//! ```text
//! case  bases-over-an-insertion  1  8  false  I@1 bases=GGG,b@1 bases=TTT  TTTACGTT
//! case  past-the-end             20 8  false  -                           CGGCANNN
//! lookup 0    N
//! lookup 93   =
//! err   high-inserted-base  1  8  false  i@3 base=233  ArrayIndexOutOfBoundsException  Index -23 out of bounds for length 127
//! ```
//!
//! The bases come from three sources: the features say what differs, the reference supplies the
//! rest, and the matrix turns a code back into a base. `Bases` and `ReadBase` are applied in a
//! second pass and overwrite what the first one wrote. The trailing fill stops at the end of the
//! reference and leaves the array's zeros, which the lookup then turns into `N`. And that lookup
//! is a 127-byte table indexed by a **signed** byte, so `]` is an `=` and `0xE9` is index -23.

use std::io::Read;

use htsjdk_cram::read_features::{
    normalize_base, restore_read_bases, to_bam_read_bases, ReadFeature, ReadFeatureError, NO_CODE,
};
use htsjdk_cram::substitution_matrix::SubstitutionMatrix;

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/cram_restore_bases.txt.gz");
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

fn reference(corpus: &str) -> (Vec<u8>, i32) {
    let row = &rows(corpus, "reference")[0];
    (
        row[0].as_bytes().to_vec(),
        row[1].parse().expect("region start"),
    )
}

fn matrix(corpus: &str) -> SubstitutionMatrix {
    let hex = rows(corpus, "matrix")[0][0];
    let bytes: Vec<u8> = (0..hex.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&hex[i..i + 2], 16).expect("hex"))
        .collect();
    SubstitutionMatrix::from_encoded(bytes.try_into().expect("five bytes"))
}

/// The dump's `escape`, applied to the restored bases.
fn escape(bases: &[u8]) -> String {
    if bases.is_empty() {
        return "-".to_string();
    }
    bases
        .iter()
        .map(|b| {
            if (0x20..=0x7E).contains(b) {
                (*b as char).to_string()
            } else {
                format!("\\u{:04X}", *b as u16)
            }
        })
        .collect()
}

/// The dump's `describe`, read back.
fn parse_features(text: &str) -> Vec<ReadFeature> {
    if text == "-" {
        return Vec::new();
    }
    text.split(',')
        .map(|entry| {
            let head = entry.split(' ').next().expect("operator");
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
            match operator {
                "X" => ReadFeature::Substitution {
                    position,
                    base: b'N',
                    reference_base: b'N',
                    code: field("code").expect("code").parse().expect("code"),
                },
                "B" => ReadFeature::ReadBase {
                    position,
                    base: field("base").expect("base").as_bytes()[0],
                    quality: field("quality").expect("quality").parse().expect("quality"),
                },
                // An inserted base is written as a number, because it is not always printable.
                "i" => ReadFeature::InsertBase {
                    position,
                    base: field("base").expect("base").parse::<u16>().expect("base") as u8,
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
                    quality: field("quality").expect("quality").parse().expect("quality"),
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

/// Every case the reference restored comes back the same way here.
#[test]
fn every_case_restores_the_bases_the_reference_restored() {
    let corpus = corpus();
    let (reference, region_start) = reference(&corpus);
    let matrix = matrix(&corpus);
    let mut compared = 0;
    for row in rows(&corpus, "case") {
        let (label, start, length, unknown, features, expected) =
            (row[0], row[1], row[2], row[3], row[4], row[5]);
        let mine = restore_read_bases(
            &parse_features(features),
            unknown == "true",
            start.parse().expect("start"),
            length.parse().expect("read length"),
            &reference,
            region_start,
            &matrix,
        )
        .unwrap_or_else(|e| panic!("{label}: {}", e.message()));
        assert_eq!(escape(&mine), expected, "{label}");
        compared += 1;
    }
    assert_eq!(compared, 24, "cases compared");
}

/// The second pass overwrites what the first one wrote, whatever wrote it.
#[test]
fn read_base_and_bases_are_applied_last_and_win() {
    let corpus = corpus();
    let of = |label: &str| -> String {
        rows(&corpus, "case")
            .into_iter()
            .find(|row| row[0] == label)
            .unwrap_or_else(|| panic!("{label}"))[5]
            .to_string()
    };
    // Over the reference fill.
    assert_eq!(of("read-base"), "ACGMTGCA");
    assert_eq!(of("bases"), "TTTTTGCA");
    // Over an insertion that already wrote those positions.
    assert_eq!(of("bases-over-an-insertion"), "TTTACGTT");
    // And over a substitution at the same position.
    assert_eq!(of("read-base-over-a-substitution"), "ACGMTGCA");
    assert_eq!(
        of("substitution-code-0"),
        "ACGATGCA",
        "which the substitution alone would have made an A"
    );
}

/// The trailing fill stops at the end of the reference, and the zeros it leaves become N.
#[test]
fn past_the_end_of_the_reference_the_array_keeps_its_zeros() {
    let corpus = corpus();
    let of = |label: &str| -> String {
        rows(&corpus, "case")
            .into_iter()
            .find(|row| row[0] == label)
            .unwrap_or_else(|| panic!("{label}"))[5]
            .to_string()
    };
    assert_eq!(of("past-the-end"), "CGGCANNN");
    assert_eq!(of("entirely-past-the-end"), "NNNN");
    assert_eq!(of("past-the-end-with-a-feature"), "TCGGCANN");
    // Nothing wrote an N: a zero byte becomes one on the way through the lookup.
    let mut zeros = vec![0u8; 3];
    to_bam_read_bases(&mut zeros).expect("in range");
    assert_eq!(zeros, b"NNN");
}

/// The lookup every restored base goes through, at every boundary that matters.
#[test]
fn the_lookup_is_a_127_byte_table_indexed_by_a_signed_byte() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "lookup") {
        let value: u16 = row[0].parse().expect("value");
        let expected = row[1];
        let mut byte = [value as u8];
        match to_bam_read_bases(&mut byte) {
            Ok(()) => assert_eq!(escape(&byte), expected, "{value}"),
            Err(error) => {
                assert_eq!(expected, "ArrayIndexOutOfBoundsException", "{value}");
                assert!(
                    matches!(
                        error,
                        ReadFeatureError::IndexOutOfBounds { length: 127, .. }
                    ),
                    "{value}: {error:?}"
                );
            }
        }
        compared += 1;
    }
    assert_eq!(compared, 18, "lookup entries compared");

    // The two that are not obvious: a NUL is an N, and `]` is an `=` because the table is built
    // by adding 32 to every BAM read base.
    let mut bracket = *b"]";
    to_bam_read_bases(&mut bracket).expect("in range");
    assert_eq!(bracket, [b'=']);
}

/// A substitution is resolved against a normalized reference base, so an IUPAC code on the
/// reference resolves as though it were N.
#[test]
fn a_substitution_resolves_against_a_normalized_reference_base() {
    let corpus = corpus();
    let matrix = matrix(&corpus);
    for row in rows(&corpus, "matrixrow") {
        let reference_base = row[0].as_bytes()[0];
        for entry in row[1].split(',') {
            let (code, base) = entry.split_once('=').expect("code=base");
            let code: u8 = code.parse().expect("code");
            assert_eq!(
                matrix.base(reference_base, code).expect("a base"),
                base.as_bytes()[0],
                "{}[{code}]",
                reference_base as char
            );
        }
    }
    assert_eq!(normalize_base(b'M'), b'N');
    assert_eq!(normalize_base(b'a'), b'A');

    let restored = rows(&corpus, "case")
        .into_iter()
        .find(|row| row[0] == "substitution-against-iupac")
        .expect("row")[5]
        .to_string();
    assert_eq!(restored, "ARWS", "the M resolved as an N and gave an A");
}

/// Everything that is refused, with the reference's own exception and message.
#[test]
fn the_failures_are_the_reference_failures() {
    let corpus = corpus();
    let (reference, region_start) = reference(&corpus);
    let matrix = matrix(&corpus);
    let mut compared = 0;
    for row in rows(&corpus, "err") {
        let (label, start, length, unknown, features, class, message) =
            (row[0], row[1], row[2], row[3], row[4], row[5], row[6]);
        let error = restore_read_bases(
            &parse_features(features),
            unknown == "true",
            start.parse().expect("start"),
            length.parse().expect("read length"),
            &reference,
            region_start,
            &matrix,
        )
        .expect_err(label);
        assert_eq!(error.message(), message, "{label}");
        match class {
            "ArrayIndexOutOfBoundsException" => assert!(
                matches!(error, ReadFeatureError::IndexOutOfBounds { .. }),
                "{label}: {error:?}"
            ),
            "IllegalArgumentException" => assert!(
                matches!(error, ReadFeatureError::Matrix(_)),
                "{label}: {error:?}"
            ),
            other => panic!("{label}: {other}"),
        }
        compared += 1;
    }
    assert_eq!(compared, 6, "failures compared");
}

/// The two shortcuts at the top return the empty sequence rather than a run of Ns.
#[test]
fn unknown_bases_and_a_zero_read_length_return_nothing() {
    let corpus = corpus();
    let (reference, region_start) = reference(&corpus);
    let matrix = matrix(&corpus);
    for (unknown, length) in [(true, 8), (false, 0)] {
        let restored = restore_read_bases(
            &[ReadFeature::Substitution {
                position: 1,
                base: b'N',
                reference_base: b'N',
                code: NO_CODE,
            }],
            unknown,
            1,
            length,
            &reference,
            region_start,
            &matrix,
        )
        .expect("no work at all");
        assert!(
            restored.is_empty(),
            "and the features are not even looked at"
        );
    }
}
