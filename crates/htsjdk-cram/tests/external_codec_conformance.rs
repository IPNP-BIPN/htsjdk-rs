//! Conformance for the codecs written on external blocks, against `ExternalIntegerCodec`,
//! `ExternalLongCodec`, `ExternalByteCodec`, `ExternalByteArrayCodec`, `ByteArrayStopCodec` and
//! `ByteArrayLenCodec`.
//!
//! Goldens from `tools/cram-conformance/CramExternalCodecDump.java` in the pinned oracle, which
//! drives each codec through the slice-blocks streams the reader and writer really use and takes
//! the blocks back out with their compression undone.
//!
//! The rows that justify the suite:
//!
//! ```text
//! stop  0  1  010002  1=01000200  01
//! len   ext-int:1  ext-bytes:1  010203,0405  1=03010203020405  010203,0405
//! err   byte-past-end  -  -  -1
//! err   len-read  ext-int:1 stop:0/2  RuntimeException  Not implemented.
//! ```
//!
//! An array holding the stop byte is written whole and read back split, with nothing reporting it.
//! A length codec and a byte codec naming the same content id interleave in one block. An external
//! byte past the end of its block is `-1`, which a real 0xFF produces too. And a byte-array-len
//! wrapping a byte-array-stop writes correctly and refuses on the way back.

use std::io::Read;

use htsjdk_cram::encoding_map::EncodingId;
use htsjdk_cram::external_codecs::{
    serialize_byte_array_len_params, serialize_external_params, serialize_stop_params,
    ByteArrayLenCodec, ByteArrayStopCodec, BytesCodec, ExternalByteArrayCodec, ExternalByteCodec,
    ExternalIntegerCodec, ExternalLongCodec, LengthCodec, SliceBlockBytes, SliceReadStreams,
    SliceWriteStreams,
};
use htsjdk_cram::huffman::CanonicalHuffman;

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/cram_external_codecs.txt.gz");
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

/// An empty array is `.`, the way the dump prints it.
fn hex(bytes: &[u8]) -> String {
    if bytes.is_empty() {
        return ".".to_string();
    }
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

fn unhex(text: &str) -> Vec<u8> {
    if text == "." {
        return Vec::new();
    }
    (0..text.len())
        .step_by(2)
        .map(|at| u8::from_str_radix(&text[at..at + 2], 16).expect("hex"))
        .collect()
}

/// The blocks a write left behind, in the dump's own notation.
fn blocks(blocks: &SliceBlockBytes) -> String {
    let mut parts = Vec::new();
    if !blocks.core.is_empty() {
        parts.push(format!("core={}", hex(&blocks.core)));
    }
    for (id, content) in &blocks.external {
        if !content.is_empty() {
            parts.push(format!("{id}={}", hex(content)));
        }
    }
    if parts.is_empty() {
        return "-".to_string();
    }
    parts.join(";")
}

fn arrays(column: &str) -> Vec<Vec<u8>> {
    column.split(',').map(unhex).collect()
}

/// ITF8, LTF8, one byte and a run of bytes, each straight onto its own block.
#[test]
fn the_external_codecs_write_what_the_reference_wrote() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "ext") {
        let (flavour, id, values, expected, back) = (
            row[0],
            row[1].parse::<i32>().expect("content id"),
            row[2],
            row[3],
            row[4],
        );
        let mut streams = SliceWriteStreams::new();
        let read_back: String = match flavour {
            "int" => {
                let codec = ExternalIntegerCodec::new(id);
                let values: Vec<i32> = values.split(',').map(|v| v.parse().expect("int")).collect();
                for value in &values {
                    codec.write(&mut streams, *value);
                }
                let written = streams.finish();
                assert_eq!(blocks(&written), expected, "ext int {id} {values:?}");
                let mut reader = SliceReadStreams::new(&written);
                values
                    .iter()
                    .map(|_| codec.read(&mut reader).expect("read").to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            }
            "long" => {
                let codec = ExternalLongCodec::new(id);
                let values: Vec<i64> = values
                    .split(',')
                    .map(|v| v.parse().expect("long"))
                    .collect();
                for value in &values {
                    codec.write(&mut streams, *value);
                }
                let written = streams.finish();
                assert_eq!(blocks(&written), expected, "ext long {id} {values:?}");
                let mut reader = SliceReadStreams::new(&written);
                values
                    .iter()
                    .map(|_| codec.read(&mut reader).expect("read").to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            }
            "byte" => {
                let codec = ExternalByteCodec::new(id);
                let values: Vec<i8> = values
                    .split(',')
                    .map(|v| v.parse().expect("byte"))
                    .collect();
                for value in &values {
                    codec.write(&mut streams, *value);
                }
                let written = streams.finish();
                assert_eq!(blocks(&written), expected, "ext byte {id} {values:?}");
                let mut reader = SliceReadStreams::new(&written);
                values
                    .iter()
                    .map(|_| codec.read(&mut reader).to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            }
            "bytes" => {
                let codec = ExternalByteArrayCodec::new(id);
                let values = arrays(values);
                for value in &values {
                    codec.write(&mut streams, value);
                }
                let written = streams.finish();
                assert_eq!(blocks(&written), expected, "ext bytes {id}");
                let mut reader = SliceReadStreams::new(&written);
                values
                    .iter()
                    .map(|value| {
                        hex(&codec
                            .read_with_length(&mut reader, value.len())
                            .expect("read"))
                    })
                    .collect::<Vec<_>>()
                    .join(",")
            }
            other => panic!("{other}"),
        };
        assert_eq!(read_back, back, "ext {flavour} {id} read back");
        compared += 1;
    }
    assert_eq!(compared, 17, "external values compared");
}

/// The stop codec, and the two ways it loses track of where an array ends.
#[test]
fn the_stop_codec_splits_on_the_data_the_way_the_reference_does() {
    let corpus = corpus();
    let mut compared = 0;
    let mut split = 0;
    for row in rows(&corpus, "stop") {
        let (stop, id, values, expected, back) = (
            row[0].parse::<u8>().expect("stop byte"),
            row[1].parse::<i32>().expect("content id"),
            arrays(row[2]),
            row[3],
            row[4],
        );
        let codec = ByteArrayStopCodec::new(stop, id);
        let mut streams = SliceWriteStreams::new();
        for value in &values {
            codec.write(&mut streams, value);
        }
        let written = streams.finish();
        assert_eq!(blocks(&written), expected, "stop {stop} {id}");

        let mut reader = SliceReadStreams::new(&written);
        let read: Vec<String> = values
            .iter()
            .map(|_| hex(&codec.read(&mut reader)))
            .collect();
        assert_eq!(read.join(","), back, "stop {stop} {id} read back");

        // An array holding the stop byte comes back as something other than what went in, and the
        // row where that happens is the one worth counting.
        if read != values.iter().map(|v| hex(v)).collect::<Vec<_>>() {
            split += 1;
        }
        compared += 1;
    }
    assert_eq!(compared, 7, "stopped runs compared");
    assert_eq!(split, 2, "of them changed by the stop byte in the data");
}

/// A length codec and a byte codec, which need not share a block or even a kind of block.
#[test]
fn a_length_and_its_bytes_land_where_the_reference_put_them() {
    let corpus = corpus();
    let mut compared = 0;
    let mut on_the_core_block = 0;
    for row in rows(&corpus, "len") {
        let (length, bytes, values, expected, back) =
            (row[0], row[1], arrays(row[2]), row[3], row[4]);
        let codec = ByteArrayLenCodec::new(length_codec(length), bytes_codec(bytes));
        let mut streams = SliceWriteStreams::new();
        for value in &values {
            codec.write(&mut streams, value).expect("written");
        }
        let written = streams.finish();
        assert_eq!(blocks(&written), expected, "len {length} {bytes}");
        if !written.core.is_empty() {
            on_the_core_block += 1;
        }

        let mut reader = SliceReadStreams::new(&written);
        let read: Vec<String> = values
            .iter()
            .map(|_| hex(&codec.read(&mut reader).expect("read")))
            .collect();
        assert_eq!(read.join(","), back, "len {length} {bytes} read back");
        compared += 1;
    }
    assert_eq!(compared, 5, "compositions compared");
    assert_eq!(
        on_the_core_block, 1,
        "of them with a length on the bit stream"
    );
}

fn length_codec(name: &str) -> LengthCodec {
    match name.split_once(':') {
        Some(("ext-int", id)) => {
            LengthCodec::External(ExternalIntegerCodec::new(id.parse().expect("content id")))
        }
        Some(("huffman", symbols)) => {
            let symbols: Vec<i32> = symbols
                .split(',')
                .map(|s| s.parse().expect("symbol"))
                .collect();
            // One symbol needs no bits at all; more than one, and the corpus uses one bit each.
            let lengths = vec![if symbols.len() > 1 { 1 } else { 0 }; symbols.len()];
            LengthCodec::Huffman(CanonicalHuffman::new(&symbols, &lengths).expect("built"))
        }
        other => panic!("{other:?}"),
    }
}

fn bytes_codec(name: &str) -> BytesCodec {
    match name.split_once(':') {
        Some(("ext-bytes", id)) => {
            BytesCodec::External(ExternalByteArrayCodec::new(id.parse().expect("content id")))
        }
        Some(("stop", rest)) => {
            let (stop, id) = rest.split_once('/').expect("stop/id");
            BytesCodec::Stop(ByteArrayStopCodec::new(
                stop.parse().expect("stop byte"),
                id.parse().expect("content id"),
            ))
        }
        other => panic!("{other:?}"),
    }
}

/// The encoding parameters, which is what a file carries in place of the codec.
#[test]
fn the_encoding_parameters_are_the_reference_bytes() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "ser") {
        let (kind, params, expected) = (row[0], row[1], row[2]);
        let bytes = match kind {
            "ext-int" | "ext-byte" | "ext-long" | "ext-bytes" => {
                serialize_external_params(number(params, "id"))
            }
            "stop" => serialize_stop_params(number(params, "stop") as u8, number(params, "id")),
            "len" => serialize_byte_array_len_params(
                (EncodingId::External, &serialize_external_params(1)),
                (EncodingId::External, &serialize_external_params(2)),
            ),
            other => panic!("{other}"),
        };
        assert_eq!(hex(&bytes), expected, "ser {kind} {params}");
        compared += 1;
    }
    assert_eq!(compared, 8, "parameter sets compared");
}

fn number(params: &str, name: &str) -> i32 {
    params
        .split(' ')
        .filter_map(|part| part.split_once('='))
        .find(|(key, _)| *key == name)
        .map(|(_, value)| value.parse().expect("number"))
        .unwrap_or_else(|| panic!("{name} in {params}"))
}

/// Everything these codecs refuse, and the two places they do not refuse when they might.
#[test]
fn the_failures_are_the_reference_failures() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "err") {
        let (what, detail, class, message) = (row[0], row[1], row[2], row[3]);
        // Every row carries the block it was read from, so nothing here is rebuilt from a label.
        let block = detail
            .split(' ')
            .filter_map(|part| part.strip_prefix("block="))
            .map(unhex)
            .next()
            .unwrap_or_default();
        let blocks = SliceBlockBytes {
            core: Vec::new(),
            external: [(1, block)].into_iter().collect(),
        };
        let mut reader = SliceReadStreams::new(&blocks);
        match what {
            "bytes-unknown-length" => {
                let error = ExternalByteArrayCodec::new(1)
                    .read(&mut reader)
                    .expect_err("refused");
                assert_eq!(error.java_exception(), class);
                assert_eq!(error.message(), message);
            }
            "bytes-past-end" => {
                let error = ExternalByteArrayCodec::new(1)
                    .read_with_length(&mut reader, 4)
                    .expect_err("refused");
                assert_eq!(error.java_exception(), class);
                assert_eq!(error.message(), message);
            }
            "len-read" => {
                let codec =
                    ByteArrayLenCodec::new(length_codec("ext-int:1"), bytes_codec("stop:0/2"));
                let mut streams = SliceWriteStreams::new();
                codec.write(&mut streams, &[1, 2]).expect("written");
                let written = streams.finish();
                let error = codec
                    .read(&mut SliceReadStreams::new(&written))
                    .expect_err("refused");
                assert_eq!(error.java_exception(), class);
                assert_eq!(error.message(), message);
            }
            // The two that return rather than refuse. The dump's class column is `-` for them,
            // and the message column carries what came back instead.
            "byte-past-end" => {
                let codec = ExternalByteCodec::new(1);
                codec.read(&mut reader);
                assert_eq!(class, "-");
                assert_eq!(codec.read(&mut reader).to_string(), message);
            }
            "stop-past-end" => {
                let codec = ByteArrayStopCodec::new(0, 1);
                let first = hex(&codec.read(&mut reader));
                let second = hex(&codec.read(&mut reader));
                assert_eq!(class, "-");
                assert_eq!(format!("{first} then {second}"), message);
            }
            // `read(length)` is the array codecs' method; the port gives it only to the codecs
            // that have one, so there is nothing to call on the rest.
            "read-length" => {
                assert_eq!(message, "Not implemented.");
                assert_eq!(class, "RuntimeException");
                assert!(matches!(
                    detail,
                    "ext-int" | "ext-byte" | "ext-long" | "stop" | "len"
                ));
            }
            other => panic!("{other}"),
        }
        compared += 1;
    }
    assert_eq!(compared, 10, "refusals compared");
}

/// The compressor each external block gets is the compression header's choice, not the codec's.
/// Nothing here is ported; the row is kept so the suite notices if the reference changes its mind.
#[test]
fn the_block_compressors_are_recorded_but_not_ported() {
    let corpus = corpus();
    let methods = rows(&corpus, "methods");
    assert_eq!(methods.len(), 1);
    let methods = methods[0][0];
    assert!(methods.contains("1=RANS"), "{methods}");
    assert!(methods.contains("7=GZIP"), "{methods}");
    assert!(
        !methods.contains("=RAW"),
        "no external block is left raw: {methods}"
    );
}
