//! Conformance for the three experimental codecs, against `GolombIntegerCodec`,
//! `GolombRiceIntegerCodec` and `GolombLongCodec`.
//!
//! Goldens from `tools/cram-conformance/CramGolombDump.java` in the pinned oracle, which drives
//! each codec through the slice-blocks streams the reader and writer really use and takes the raw
//! core block back out.
//!
//! The rows that justify the suite:
//!
//! ```text
//! golomb  0  4  -1   60    3
//! rice    0  8  255  7f80  255
//! rice    0  8  256  8000  256
//! err  golomb  offset=0 m=1  0  IllegalArgumentException  M parameter must be at least 2.
//! ```
//!
//! Golomb does not round-trip a value whose offset sum is negative and does not say so.
//! Golomb-Rice's parameter is a power, not a divisor: built with 8 it divides by 256. And the
//! divisor Golomb refuses is one Golomb-Rice takes without a word.

use std::io::Read;

use htsjdk_cram::bit_stream::{BitInputStream, BitOutputStream};
use htsjdk_cram::golomb::{
    parse_params, serialize_params, GolombError, GolombIntegerCodec, GolombLongCodec,
    GolombRiceIntegerCodec,
};

fn corpus() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/cram_golomb.txt.gz");
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

fn hex(bytes: &[u8]) -> String {
    if bytes.is_empty() {
        return "-".to_string();
    }
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

/// Golomb: a unary quotient, then a remainder at one of two widths.
#[test]
fn golomb_writes_and_reads_what_the_reference_did() {
    let corpus = corpus();
    let mut compared = 0;
    let mut corrupted = 0;
    for row in rows(&corpus, "golomb") {
        let (offset, m, value, expected, back) = (
            row[0].parse::<i32>().expect("offset"),
            row[1].parse::<i32>().expect("m"),
            row[2].parse::<i32>().expect("value"),
            row[3],
            row[4].parse::<i32>().expect("read back"),
        );
        let codec = GolombIntegerCodec::new(offset, m).expect("built");
        let mut out = BitOutputStream::new();
        codec.write(&mut out, value).expect("written");
        let bytes = out.into_bytes();
        assert_eq!(hex(&bytes), expected, "golomb {offset}/{m}/{value}");

        let mut input = BitInputStream::new(&bytes);
        assert_eq!(
            codec.read(&mut input).expect("read"),
            back,
            "golomb {offset}/{m}/{value} read back"
        );
        if back != value {
            corrupted += 1;
        }
        compared += 1;
    }
    assert_eq!(compared, 48, "Golomb values compared");

    // The rows where the round trip does not hold are exactly the ones whose offset sum is below
    // zero. Nothing in the codec reports them.
    assert_eq!(corrupted, 2, "values the codec silently changed");
}

/// Golomb-Rice, whose parameter is a power rather than a divisor.
#[test]
fn golomb_rice_writes_and_reads_what_the_reference_did() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "rice") {
        let (offset, log2m, value, expected, back) = (
            row[0].parse::<i32>().expect("offset"),
            row[1].parse::<i32>().expect("log2m"),
            row[2].parse::<i32>().expect("value"),
            row[3],
            row[4].parse::<i32>().expect("read back"),
        );
        let codec = GolombRiceIntegerCodec::new(offset, log2m);
        let mut out = BitOutputStream::new();
        codec.write(&mut out, value).expect("written");
        let bytes = out.into_bytes();
        assert_eq!(hex(&bytes), expected, "rice {offset}/{log2m}/{value}");

        let mut input = BitInputStream::new(&bytes);
        assert_eq!(
            codec.read(&mut input).expect("read"),
            back,
            "rice {offset}/{log2m}/{value} read back"
        );
        compared += 1;
    }
    assert_eq!(compared, 36, "Golomb-Rice values compared");

    // Built with 8 it divides by 256, so 255 fits behind a single zero and 256 does not. That is
    // the whole of the parameter's meaning, and the reference calls it m.
    let codec = GolombRiceIntegerCodec::new(0, 8);
    let mut out = BitOutputStream::new();
    codec.write(&mut out, 255).expect("written");
    assert_eq!(hex(&out.into_bytes()), "7f80");
}

/// Golomb-Long: the same arithmetic on a long.
#[test]
fn golomb_long_writes_and_reads_what_the_reference_did() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "long") {
        let (offset, m, value, expected, back) = (
            row[0].parse::<i64>().expect("offset"),
            row[1].parse::<i32>().expect("m"),
            row[2].parse::<i64>().expect("value"),
            row[3],
            row[4].parse::<i64>().expect("read back"),
        );
        let codec = GolombLongCodec::new(offset, m).expect("built");
        let mut out = BitOutputStream::new();
        codec.write(&mut out, value).expect("written");
        let bytes = out.into_bytes();
        assert_eq!(hex(&bytes), expected, "long {offset}/{m}/{value}");

        let mut input = BitInputStream::new(&bytes);
        assert_eq!(
            codec.read(&mut input).expect("read"),
            back,
            "long {offset}/{m}/{value} read back"
        );
        compared += 1;
    }
    assert_eq!(compared, 14, "Golomb-Long values compared");
}

/// Values in a row, where the unary prefixes butt against each other. For a power-of-two divisor
/// the three codecs agree byte for byte, which is the only place they can be compared directly.
#[test]
fn a_sequence_packs_the_way_the_reference_packed_it() {
    let corpus = corpus();
    let mut compared = 0;
    let mut encodings = Vec::new();
    for row in rows(&corpus, "seq") {
        let (name, params, values, expected) = (row[0], row[1], row[2], row[3]);
        let numbers: Vec<i32> = params
            .split(' ')
            .filter_map(|part| part.split_once('='))
            .map(|(_, value)| value.parse().expect("param"))
            .collect();
        let values: Vec<i64> = values
            .split(',')
            .map(|value| value.parse().expect("value"))
            .collect();

        let mut out = BitOutputStream::new();
        for value in &values {
            match name {
                "golomb" => GolombIntegerCodec::new(numbers[0], numbers[1])
                    .expect("built")
                    .write(&mut out, *value as i32),
                "rice" => GolombRiceIntegerCodec::new(numbers[0], numbers[1])
                    .write(&mut out, *value as i32),
                "long" => GolombLongCodec::new(i64::from(numbers[0]), numbers[1])
                    .expect("built")
                    .write(&mut out, *value),
                other => panic!("{other}"),
            }
            .expect("written");
        }
        let bytes = hex(&out.into_bytes());
        assert_eq!(bytes, expected, "seq {name} {params}");
        encodings.push(bytes);
        compared += 1;
    }
    assert_eq!(compared, 3, "sequences compared");
    assert_eq!(encodings[0], encodings[1], "golomb and rice agree at m = 4");
    assert_eq!(encodings[1], encodings[2], "and so does the long");
}

/// The encoding parameters, which are two ITF8s whatever the codec.
#[test]
fn the_encoding_parameters_are_the_reference_bytes() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "ser") {
        let (offset, m, expected, reparsed) = (
            row[1].parse::<i32>().expect("offset"),
            row[2].parse::<i32>().expect("m"),
            row[3],
            row[4],
        );
        let bytes = serialize_params(offset, m);
        assert_eq!(hex(&bytes), expected, "ser {} {offset}/{m}", row[0]);
        let (offset, m) = parse_params(&bytes).expect("parsed");
        assert_eq!(hex(&serialize_params(offset, m)), reparsed);
        compared += 1;
    }
    assert_eq!(compared, 4, "parameter sets compared");
}

/// What each refuses, with the reference's own wording, and the one that refuses nothing.
#[test]
fn the_failures_are_the_reference_failures() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "err") {
        let (name, params, class, message) = (row[0], row[1], row[3], row[4]);
        let error = if params == "read-length" {
            match name {
                "golomb" => GolombIntegerCodec::new(0, 4)
                    .expect("built")
                    .read_with_length(4)
                    .expect_err("refused"),
                "rice" => GolombRiceIntegerCodec::new(0, 2)
                    .read_with_length(4)
                    .expect_err("refused"),
                "long" => GolombLongCodec::new(0, 4)
                    .expect("built")
                    .read_with_length(4)
                    .map(|_| 0)
                    .expect_err("refused"),
                other => panic!("{other}"),
            }
        } else {
            let numbers: Vec<i32> = params
                .split(' ')
                .filter_map(|part| part.split_once('='))
                .map(|(_, value)| value.parse().expect("param"))
                .collect();
            match name {
                "golomb" => GolombIntegerCodec::new(numbers[0], numbers[1]).expect_err("refused"),
                "long" => {
                    GolombLongCodec::new(i64::from(numbers[0]), numbers[1]).expect_err("refused")
                }
                other => panic!("{other}"),
            }
        };
        assert_eq!(error.java_exception(), class, "err {name} {params}");
        assert_eq!(error.message(), message, "err {name} {params}");
        compared += 1;
    }
    assert_eq!(compared, 7, "refusals compared");

    // Golomb-Rice takes the divisor the other two refuse, and the corpus carries the row it wrote
    // with it rather than a refusal.
    assert!(matches!(
        GolombIntegerCodec::new(0, 1),
        Err(GolombError::MTooSmall)
    ));
    let accepted = rows(&corpus, "rice")
        .into_iter()
        .filter(|row| row[1] == "1")
        .count();
    assert!(accepted > 0, "Golomb-Rice wrote with the refused parameter");
}
