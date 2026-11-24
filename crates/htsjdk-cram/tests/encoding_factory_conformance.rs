//! Conformance for `EncodingFactory`: which encoding a data series type and an identifier resolve
//! to, fall-through included.
//!
//! Goldens from `tools/cram-conformance/CramEncodingFactoryDump.java` in the pinned oracle, which
//! asks for every one of the four types against every one of the ten identifiers and records the
//! class that came back, or the refusal.
//!
//! The rows that justify the suite:
//!
//! ```text
//! make  INT   BYTE_ARRAY_LEN   010101010102  ByteArrayLenEncoding   LenEncoding: Content ID: 1 ByteEncoding: Content ID: 2
//! make  LONG  BYTE_ARRAY_STOP  0001          ByteArrayStopEncoding  Content ID: 1 StopByte: 0
//! err   BYTE  BYTE_ARRAY_LEN   010101010102  IllegalArgumentException  Encoding not found: value type=BYTE, encoding id=BYTE_ARRAY_LEN
//! ```
//!
//! An INT named with a byte array identifier gets a byte array encoding, because only the BYTE arm
//! of the reference's switch ends in a break. BYTE never falls through; INT and LONG always do.

use std::io::Read;

use htsjdk_cram::encoding_factory::{create_encoding, DataSeriesType, FactoryError};
use htsjdk_cram::encoding_map::EncodingId;

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/cram_encoding_factory.txt.gz");
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

fn unhex(text: &str) -> Vec<u8> {
    if text == "-" {
        return Vec::new();
    }
    (0..text.len())
        .step_by(2)
        .map(|at| u8::from_str_radix(&text[at..at + 2], 16).expect("hex"))
        .collect()
}

fn value_type(name: &str) -> DataSeriesType {
    match name {
        "BYTE" => DataSeriesType::Byte,
        "INT" => DataSeriesType::Int,
        "LONG" => DataSeriesType::Long,
        "BYTE_ARRAY" => DataSeriesType::ByteArray,
        other => panic!("{other}"),
    }
}

fn encoding_id(name: &str) -> EncodingId {
    match name {
        "NULL" => EncodingId::Null,
        "EXTERNAL" => EncodingId::External,
        "GOLOMB" => EncodingId::Golomb,
        "HUFFMAN" => EncodingId::Huffman,
        "BYTE_ARRAY_LEN" => EncodingId::ByteArrayLen,
        "BYTE_ARRAY_STOP" => EncodingId::ByteArrayStop,
        "BETA" => EncodingId::Beta,
        "SUBEXPONENTIAL" => EncodingId::Subexponential,
        "GOLOMB_RICE" => EncodingId::GolombRice,
        "GAMMA" => EncodingId::Gamma,
        other => panic!("{other}"),
    }
}

/// Every pair the factory resolves, and the class and description it resolves to.
#[test]
fn every_pair_resolves_the_way_the_reference_resolved_it() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "make") {
        let (type_name, id_name, params, class, described) =
            (row[0], row[1], unhex(row[2]), row[3], row[4]);
        let encoding = create_encoding(value_type(type_name), encoding_id(id_name), &params)
            .unwrap_or_else(|error| panic!("{type_name}/{id_name}: {}", error.message()));
        assert_eq!(encoding.java_class(), class, "{type_name}/{id_name}");
        assert_eq!(encoding.describe(), described, "{type_name}/{id_name}");
        compared += 1;
    }
    assert_eq!(compared, 27, "pairs resolved");
}

/// Every pair the factory refuses, with the message that names both halves of what was asked for.
#[test]
fn every_refusal_is_the_reference_refusal() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "err") {
        let (type_name, id_name, params, class, message) =
            (row[0], row[1], unhex(row[2]), row[3], row[4]);
        // The one refusal htsjdk 6.0.0 lifts (#1854) and this port takes: see
        // docs/htsjdk-upstream.md. 5.0.0 could not read such a slice at all, so no output moves.
        if (type_name, id_name) == ("BYTE_ARRAY", "HUFFMAN") {
            let encoding = create_encoding(value_type(type_name), encoding_id(id_name), &params)
                .expect("htsjdk 6.0.0 reads it");
            assert_eq!(encoding.java_class(), "CanonicalHuffmanByteArrayEncoding");
            compared += 1;
            continue;
        }
        let error = create_encoding(value_type(type_name), encoding_id(id_name), &params)
            .expect_err("refused");
        assert_eq!(error.java_exception(), class, "{type_name}/{id_name}");
        assert_eq!(error.message(), message, "{type_name}/{id_name}");
        compared += 1;
    }
    assert_eq!(compared, 22, "pairs refused");
}

/// The fall-through, stated as itself rather than left implicit in the row counts.
#[test]
fn only_the_byte_arm_stops_at_its_own_identifiers() {
    let stop_params = [0x00u8, 0x01];
    let len_params = [0x01u8, 0x01, 0x01, 0x01, 0x01, 0x02];

    // BYTE is the one arm with a break, so a byte array identifier is refused there.
    for (id, params) in [
        (EncodingId::ByteArrayStop, &stop_params[..]),
        (EncodingId::ByteArrayLen, &len_params[..]),
    ] {
        assert!(matches!(
            create_encoding(DataSeriesType::Byte, id, params),
            Err(FactoryError::NotFound { .. })
        ));
        // And accepted for both of the arms that fall through.
        for value_type in [DataSeriesType::Int, DataSeriesType::Long] {
            let encoding = create_encoding(value_type, id, params).expect("fell through");
            assert!(
                encoding.java_class().starts_with("ByteArray"),
                "{} got {}",
                value_type.name(),
                encoding.java_class()
            );
        }
    }

    // An INT with an external identifier stops in its own arm rather than falling into LONG's.
    let external =
        create_encoding(DataSeriesType::Int, EncodingId::External, &[0x01]).expect("its own arm");
    assert_eq!(external.java_class(), "ExternalIntegerEncoding");

    // A LONG with GOLOMB gets the long codec, and an INT with GOLOMB the integer one, so the
    // fall-through is ordered rather than a free-for-all.
    assert_eq!(
        create_encoding(DataSeriesType::Long, EncodingId::Golomb, &[0x00, 0x04])
            .expect("long arm")
            .java_class(),
        "GolombLongEncoding"
    );
    assert_eq!(
        create_encoding(DataSeriesType::Int, EncodingId::Golomb, &[0x00, 0x04])
            .expect("int arm")
            .java_class(),
        "GolombIntegerEncoding"
    );
}

/// NULL matches nothing anywhere, which is what makes the refusal reachable from every type.
#[test]
fn the_null_identifier_is_refused_by_every_type() {
    let corpus = corpus();
    let refused: Vec<&str> = rows(&corpus, "err")
        .into_iter()
        .filter(|row| row[1] == "NULL")
        .map(|row| row[0])
        .collect();
    assert_eq!(refused, ["BYTE", "INT", "LONG", "BYTE_ARRAY"]);

    for value_type in [
        DataSeriesType::Byte,
        DataSeriesType::Int,
        DataSeriesType::Long,
        DataSeriesType::ByteArray,
    ] {
        assert_eq!(
            create_encoding(value_type, EncodingId::Null, &[])
                .expect_err("refused")
                .message(),
            format!(
                "Encoding not found: value type={}, encoding id=NULL",
                value_type.name()
            )
        );
    }
}

/// htsjdk 6.0.0 #1854: a QS series encoded HUFFMAN. The record reader builds QS as a BYTE_ARRAY
/// for every slice, so 5.0.0 refused the whole slice whether or not a record kept its qualities.
/// A LONG series named HUFFMAN still falls through to the refusal.
#[test]
fn a_byte_array_series_can_be_huffman() {
    use htsjdk_cram::codecs::read_byte_array;
    use htsjdk_cram::external_codecs::{SliceReadStreams, SliceWriteStreams};
    use htsjdk_cram::huffman::{serialize_byte_params, CanonicalHuffman};

    // One symbol of length 0: the constant series, which writes no bits.
    let constant = create_encoding(
        DataSeriesType::ByteArray,
        EncodingId::Huffman,
        &serialize_byte_params(&[30], &[0]),
    )
    .expect("read");
    let blocks = SliceWriteStreams::new().finish();
    let mut streams = SliceReadStreams::new(&blocks);
    assert_eq!(
        read_byte_array(&constant, &mut streams, Some(4)).unwrap(),
        [30u8; 4]
    );

    // Two symbols of one bit each, read off the core bit stream a symbol at a time.
    let (symbols, bit_lengths) = ([20i8, 40], [1, 1]);
    let two = create_encoding(
        DataSeriesType::ByteArray,
        EncodingId::Huffman,
        &serialize_byte_params(&symbols, &bit_lengths),
    )
    .expect("read");
    let huffman = CanonicalHuffman::new(&symbols, &bit_lengths).unwrap();
    let mut write = SliceWriteStreams::new();
    for symbol in [40i8, 20, 20, 40, 40] {
        huffman.write(write.core(), symbol).unwrap();
    }
    let blocks = write.finish();
    let mut streams = SliceReadStreams::new(&blocks);
    assert_eq!(
        read_byte_array(&two, &mut streams, Some(5)).unwrap(),
        [40u8, 20, 20, 40, 40]
    );

    assert!(matches!(
        create_encoding(
            DataSeriesType::Long,
            EncodingId::Huffman,
            &[0x01, 0x2a, 0x01, 0x00]
        ),
        Err(FactoryError::NotFound { .. })
    ));
}
