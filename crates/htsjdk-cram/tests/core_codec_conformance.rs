//! Conformance for the three core integer codecs, against `BetaIntegerCodec`,
//! `GammaIntegerCodec` and `SubexponentialIntegerCodec`.
//!
//! Goldens from `tools/cram-conformance/CramCoreCodecDump.java` in the pinned oracle, which drives
//! the codecs through the slice-blocks streams the reader and writer really use, and takes the raw
//! core block back out.
//!
//! The rows that justify the suite:
//!
//! ```text
//! gamma   0  4     20    4
//! subexp  0  2  3  60    3
//! subexp  0  2  4  80    4
//! err  gamma  offset=0  0  IllegalArgumentException  Gamma codec handles only positive values.  Value 0 + Offset 0 <= 0
//! ```
//!
//! Gamma writes `length - 1` zeros and then the value with its top bit, so 4 is `001 00`.
//! Subexponential changes shape at `2^k`: 3 goes in two bits behind one zero, 4 behind a unary
//! prefix and without its top bit. And the two spaces in Gamma's message are the reference's.

use std::io::Read;

use htsjdk_cram::bit_stream::{BitInputStream, BitOutputStream};
use htsjdk_cram::core_codecs::{
    BetaIntegerCodec, CodecError, GammaIntegerCodec, SubexponentialIntegerCodec,
};

fn corpus() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/cram_core_codecs.txt.gz");
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

fn write_one(
    write: impl FnOnce(&mut BitOutputStream) -> Result<(), CodecError>,
) -> Result<Vec<u8>, CodecError> {
    let mut out = BitOutputStream::new();
    write(&mut out)?;
    Ok(out.into_bytes())
}

/// Beta: a fixed width, and every value the reference wrote.
#[test]
fn beta_writes_and_reads_what_the_reference_did() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "beta") {
        let (offset, bits, value, expected, back) = (
            row[0].parse::<i32>().expect("offset"),
            row[1].parse::<i32>().expect("bits"),
            row[2].parse::<i32>().expect("value"),
            row[3],
            row[4].parse::<i32>().expect("read back"),
        );
        let codec = BetaIntegerCodec::new(offset, bits);
        let bytes = write_one(|out| codec.write(out, value)).expect("written");
        assert_eq!(hex(&bytes), expected, "beta {offset}/{bits}/{value}");
        let mut input = BitInputStream::new(&bytes);
        assert_eq!(codec.read(&mut input).expect("back"), back);
        compared += 1;
    }
    assert_eq!(compared, 13, "beta values compared");
}

/// Gamma: a length prefix of zeros, and a bit length that comes from a floating-point log.
#[test]
fn gamma_writes_and_reads_what_the_reference_did() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "gamma") {
        let (offset, value, expected, back) = (
            row[0].parse::<i32>().expect("offset"),
            row[1].parse::<i32>().expect("value"),
            row[2],
            row[3].parse::<i32>().expect("read back"),
        );
        let codec = GammaIntegerCodec::new(offset);
        let bytes = write_one(|out| codec.write(out, value)).expect("written");
        assert_eq!(hex(&bytes), expected, "gamma {offset}/{value}");
        let mut input = BitInputStream::new(&bytes);
        assert_eq!(codec.read(&mut input).expect("back"), back);
        compared += 1;
    }
    assert_eq!(compared, 34, "gamma values compared");
}

/// Subexponential: two regimes, and the split between them at `2^k`.
#[test]
fn subexponential_writes_and_reads_what_the_reference_did() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "subexp") {
        let (offset, k, value, expected, back) = (
            row[0].parse::<i32>().expect("offset"),
            row[1].parse::<i32>().expect("k"),
            row[2].parse::<i32>().expect("value"),
            row[3],
            row[4].parse::<i32>().expect("read back"),
        );
        let codec = SubexponentialIntegerCodec::new(offset, k);
        let bytes = write_one(|out| codec.write(out, value)).expect("written");
        assert_eq!(hex(&bytes), expected, "subexp {offset}/{k}/{value}");
        let mut input = BitInputStream::new(&bytes);
        assert_eq!(codec.read(&mut input).expect("back"), back);
        compared += 1;
    }
    assert_eq!(compared, 48, "subexponential values compared");
}

/// The bit length comes from a floating-point log, and every power of two is where that could
/// go wrong. The corpus walks them up to `2^31 - 1`.
#[test]
fn the_bit_length_is_right_at_every_power_of_two() {
    let corpus = corpus();
    let gamma: Vec<(i32, &str)> = rows(&corpus, "gamma")
        .into_iter()
        .filter(|row| row[0] == "0")
        .map(|row| (row[1].parse().expect("value"), row[2]))
        .collect();

    // A power of two is the first value needing another bit, so its encoding is one byte-pattern
    // longer than its predecessor's.
    for (value, encoded) in &gamma {
        let bits = 2 * (32 - (*value as u32).leading_zeros()) - 1;
        let bytes = bits.div_ceil(8) as usize;
        assert_eq!(
            encoded.len() / 2,
            bytes,
            "gamma {value} takes {bits} bits, so {bytes} bytes"
        );
    }
    assert!(
        gamma.iter().any(|(value, _)| *value == i32::MAX),
        "and the corpus reaches the top of the range"
    );
}

/// Values in a row pack against each other with no alignment between them.
#[test]
fn a_sequence_packs_with_no_alignment() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "seq") {
        let (name, params, values, expected) = (row[0], row[1], row[2], row[3]);
        let values: Vec<i32> = values
            .split(',')
            .map(|v| v.parse().expect("value"))
            .collect();
        let bytes = write_one(|out| {
            for value in &values {
                match name {
                    "beta" => BetaIntegerCodec::new(0, 3).write(out, *value)?,
                    "gamma" => GammaIntegerCodec::new(0).write(out, *value)?,
                    "subexp" => SubexponentialIntegerCodec::new(0, 1).write(out, *value)?,
                    other => panic!("{other}"),
                }
            }
            Ok(())
        })
        .expect("written");
        assert_eq!(hex(&bytes), expected, "{name} {params}");
        compared += 1;
    }
    assert_eq!(compared, 3, "sequences compared");
}

/// Everything each codec refuses, with the reference's own message.
#[test]
fn the_failures_are_the_reference_failures() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "err") {
        let (name, params, value, class, message) = (row[0], row[1], row[2], row[3], row[4]);
        let value: i32 = value.parse().expect("value");
        let numbers: Vec<i32> = params
            .split(' ')
            .filter_map(|part| part.split_once('='))
            .map(|(_, v)| v.parse().expect("param"))
            .collect();
        let error = write_one(|out| match name {
            "beta" => BetaIntegerCodec::new(numbers[0], numbers[1]).write(out, value),
            "gamma" => GammaIntegerCodec::new(numbers[0]).write(out, value),
            "subexp" => SubexponentialIntegerCodec::new(numbers[0], numbers[1]).write(out, value),
            other => panic!("{other}"),
        })
        .expect_err("refused");
        assert_eq!(error.java_exception(), class, "{name} {params} {value}");
        assert_eq!(error.message(), message, "{name} {params} {value}");
        compared += 1;
    }
    assert_eq!(compared, 7, "failures compared");

    // Beta is the only one of the three with an upper bound at all.
    let above: Vec<&str> = rows(&corpus, "err")
        .into_iter()
        .filter(|row| row[4].contains("greater than or equal to limit"))
        .map(|row| row[0])
        .collect();
    assert_eq!(above, ["beta", "beta"]);
}
