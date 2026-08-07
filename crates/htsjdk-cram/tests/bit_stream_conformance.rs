//! Conformance for the core bit stream, against `DefaultBitOutputStream` and
//! `DefaultBitInputStream`.
//!
//! Goldens from `tools/cram-conformance/CramBitStreamDump.java` in the pinned oracle.
//!
//! The rows that justify the suite:
//!
//! ```text
//! write  a-byte-too-wide-for-its-width  byte 0xFF in 3 bits  e0    3
//! write  long-12-bits                   long 0xABC in 12 bits  abc0  4
//! read   read-the-padding               80  readBits 1,7  true,0
//! err    zero-bits-into-a-partial-buffer  byte 0x0A in 4 bits then byte 0x01 in 0 bits  ArrayIndexOutOfBoundsException  Index 8 out of bounds for length 8
//! ```
//!
//! Bits go in most significant first, so a 12-bit write leaves its low four in the buffer and the
//! flush pads them with zeros on the right. The padding is indistinguishable from data: reading one
//! bit and then seven from `80` gives a `true` and a zero, and nothing says which of the seven were
//! written. And a write of zero bits against a partial buffer indexes a table of eight masks with
//! eight.
//!
//! Every row carries what it was given, refusals included, so nothing here is rebuilt from a label.

use std::io::Read;

use htsjdk_cram::bit_stream::{BitError, BitInputStream, BitOutputStream};

fn corpus() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/cram_bit_stream.txt.gz");
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

fn unhex(text: &str) -> Vec<u8> {
    if text == "-" {
        return Vec::new();
    }
    (0..text.len())
        .step_by(2)
        .map(|at| u8::from_str_radix(&text[at..at + 2], 16).expect("hex"))
        .collect()
}

/// Replay the write the row describes. The description is the dump's own, so the operation comes
/// from the golden rather than from a table this test keeps.
fn replay_write(description: &str) -> Result<(Vec<u8>, i32), BitError> {
    let mut out = BitOutputStream::new();
    for step in description.split(" then ") {
        let step = step.trim();
        if step == "nothing" {
            continue;
        }
        let words: Vec<&str> = step.split_whitespace().collect();
        match words.as_slice() {
            // "byte 0x01 in 1 bits" and the second half of "byte 0x0A in 4 then 0x05 in 4",
            // which the dump writes without repeating the word "byte".
            ["byte", value, "in", bits, ..] => {
                out.write_byte_bits(parse_hex(value) as u8, bits.parse().expect("bits"))?
            }
            ["long", value, "in", bits, ..] => {
                out.write_long_bits(parse_hex(value) as i64, bits.parse().expect("bits"))?
            }
            ["int", value, "in", bits, ..] => {
                out.write_int_bits(parse_hex(value) as i32, bits.parse().expect("bits"))?
            }
            [value, "in", bits] => {
                out.write_byte_bits(parse_hex(value) as u8, bits.parse().expect("bits"))?
            }
            [bit, "repeated", repeat] => out.write_bits(
                bit.parse().expect("a boolean"),
                repeat.parse().expect("a repeat count"),
            )?,
            other => panic!("{other:?}"),
        }
    }
    let buffered = out.buffered_bits();
    Ok((out.into_bytes(), buffered))
}

fn parse_hex(text: &str) -> i128 {
    let text = text.strip_prefix("0x").unwrap_or(text);
    i128::from_str_radix(text, 16).expect("a hex value")
}

/// Every write the reference made, and the bits it left in the buffer.
#[test]
fn every_write_lands_the_reference_bytes() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "write") {
        let (label, description, expected, buffered) = (row[0], row[1], row[2], row[3]);
        let (bytes, left) = replay_write(description).unwrap_or_else(|error| {
            panic!("{label}: {}", error.message());
        });
        assert_eq!(hex(&bytes), expected, "write {label}: {description}");
        // The flush pads the partial byte, so what is left in the buffer is what the padding
        // covered. The dump reports it, which is the only way to see it at all.
        assert_eq!(left.to_string(), buffered, "write {label} buffered bits");
        compared += 1;
    }
    assert_eq!(compared, 32, "writes compared");
}

/// Every read the reference made, over bytes the row carries.
#[test]
fn every_read_returns_the_reference_values() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "read") {
        let (label, input, operation, expected) = (row[0], unhex(row[1]), row[2], row[3]);
        let (kind, widths) = operation.split_once(' ').expect("an operation");
        let widths: Vec<i32> = widths
            .split(',')
            .map(|width| width.parse().expect("a width"))
            .collect();

        let mut stream = BitInputStream::new(&input);
        let mut values = Vec::new();
        for width in widths {
            values.push(match kind {
                // A width of one is read through readBit, which returns a boolean, and the dump
                // prints it as one.
                "readBits" if width == 1 => stream
                    .read_bit()
                    .unwrap_or_else(|error| panic!("{label}: {}", error.message()))
                    .to_string(),
                "readBits" => stream
                    .read_bits(width)
                    .unwrap_or_else(|error| panic!("{label}: {}", error.message()))
                    .to_string(),
                "readLongBits" => stream
                    .read_long_bits(width)
                    .unwrap_or_else(|error| panic!("{label}: {}", error.message()))
                    .to_string(),
                other => panic!("{other}"),
            });
        }
        assert_eq!(values.join(","), expected, "read {label}");
        compared += 1;
    }
    assert_eq!(compared, 11, "reads compared");
}

/// Everything the stream refuses, with the class and the message the reference gave.
#[test]
fn the_failures_are_the_reference_failures() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "err") {
        let (label, given, class, message) = (row[0], row[1], row[2], row[3]);
        let error = if let Some((input, operation)) = given.split_once(" read") {
            // A read: the bytes, then the widths asked for.
            let input = unhex(input);
            let (kind, widths) = operation.split_once(' ').expect("an operation");
            let mut stream = BitInputStream::new(&input);
            let mut failure = None;
            for width in widths.split(',') {
                let width: i32 = width.parse().expect("a width");
                let result = match kind {
                    "Bits" if width == 1 => stream.read_bit().map(i64::from),
                    "Bits" => stream.read_bits(width).map(i64::from),
                    "LongBits" => stream.read_long_bits(width),
                    other => panic!("{other}"),
                };
                if let Err(error) = result {
                    failure = Some(error);
                    break;
                }
            }
            failure.expect("refused")
        } else {
            replay_write(given).expect_err("refused")
        };
        assert_eq!(error.java_exception(), class, "err {label}");
        assert_eq!(error.message(), message, "err {label}");
        compared += 1;
    }
    assert_eq!(compared, 11, "refusals compared");

    // The one refusal that is neither a width nor an end of stream: a write of zero bits against a
    // buffer that already holds some.
    let mut out = BitOutputStream::new();
    out.write_byte_bits(0x0A, 4).expect("four bits");
    assert_eq!(
        out.write_byte_bits(1, 0),
        Err(BitError::MaskIndexOutOfBounds)
    );
}

/// The padding a flush adds is data as far as anything downstream can tell.
#[test]
fn the_padding_is_indistinguishable_from_what_was_written() {
    // One bit of true, flushed, is the same byte as eight bits of 0x80.
    let mut one = BitOutputStream::new();
    one.write_bits(true, 1).expect("one bit");
    let mut eight = BitOutputStream::new();
    eight.write_byte_bits(0x80, 8).expect("eight bits");
    assert_eq!(one.into_bytes(), eight.into_bytes());

    // And the corpus says so too: reading one bit and then seven from that byte gives a true and a
    // zero, with nothing to say the seven were never written.
    let corpus = corpus();
    let padding = rows(&corpus, "read")
        .into_iter()
        .find(|row| row[0] == "read-the-padding")
        .expect("the padding row");
    assert_eq!(padding[1], "80");
    assert_eq!(padding[3], "true,0");
}
