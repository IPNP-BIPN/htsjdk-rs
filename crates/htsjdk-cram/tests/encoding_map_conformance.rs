//! Conformance for the data series encoding map, against
//! `htsjdk.samtools.cram.structure.CompressionHeaderEncodingMap`.
//!
//! Goldens from `tools/cram-conformance/CramEncodingMapDump.java` in the pinned oracle.
//!
//! The rows that justify the suite:
//!
//! ```text
//! written  BF,CF,RI,RL,AP,RG,RN,NF,MF,NS,NP,TS,TL,MQ,FN,FP,FC,BA,QS,BS,IN,DL,RS,SC,PD,HC
//! ignored  TC,TN
//! refuse   encoding-id-negative  java.lang.ArrayIndexOutOfBoundsException  Index -1 out of ...
//! ```
//!
//! Twenty-six of the thirty-two data series are written, in the enum's ordinal order rather than
//! the order the constructor populates them in. Two more are read and then dropped. And an encoding
//! id byte of 255 arrives as index -1, because the read is signed and the index is unchecked.

use std::io::Read;

use htsjdk_cram::encoding_map::{
    DataSeries, EncodingId, EncodingMap, EncodingMapError, DATA_SERIES, ENCODING_ID_COUNT, NOT_READ,
};

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/cram_encoding_map.txt.gz");
    let file = std::fs::File::open(&path).expect("corpus");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("corpus is gzip");
    text
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
        .map(|i| u8::from_str_radix(&text[i..i + 2], 16).expect("hex"))
        .collect()
}

fn rows<'a>(corpus: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    let prefix = format!("{kind}\t");
    corpus
        .lines()
        .filter_map(|line| line.strip_prefix(&prefix))
        .map(|rest| rest.split('\t').collect())
        .collect()
}

/// The two counts, measured rather than assumed.
#[test]
fn there_are_thirty_two_data_series_and_ten_encodings() {
    let corpus = corpus();
    let sizes = rows(&corpus, "sizes");
    assert_eq!(sizes.len(), 1);
    assert_eq!(sizes[0][0], DATA_SERIES.len().to_string());
    assert_eq!(sizes[0][1], ENCODING_ID_COUNT.to_string());
}

/// Every data series, its type, and the content id this implementation assigns.
///
/// The content ids are not in the specification, so a reader must discover them from the map. They
/// are pinned here because a port that writes different ones writes a different file.
#[test]
fn every_data_series_has_the_ordinal_type_and_content_id_the_reference_gives_it() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "series") {
        let ordinal: usize = row[0].parse().expect("ordinal");
        let series = DataSeries(ordinal);
        assert_eq!(series.canonical_name(), row[1], "ordinal {ordinal}");
        assert_eq!(series.series_type().name(), row[2], "{}: type", row[1]);
        assert_eq!(
            series.content_id().to_string(),
            row[3],
            "{}: content id",
            row[1]
        );
        assert_eq!(
            DataSeries::by_canonical_name(&row[1].as_bytes().try_into().expect("two bytes")),
            Some(series),
            "{}: found by name",
            row[1]
        );
        compared += 1;
    }
    assert_eq!(compared, 32, "data series compared");
}

/// The ten encodings, and which of them may live in an external block.
#[test]
fn every_encoding_id_is_the_one_the_reference_numbers_it() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "encid") {
        let id: i32 = row[0].parse().expect("id");
        let encoding = EncodingId::from_id(id).expect("a known id");
        assert_eq!(encoding.name(), row[1], "id {id}");
        assert_eq!(
            encoding.is_external().to_string(),
            row[2],
            "{}: external",
            row[1]
        );
        compared += 1;
    }
    assert_eq!(compared, 10, "encodings compared");
    assert_eq!(EncodingId::from_id(ENCODING_ID_COUNT), None);
}

/// Every map the reference wrote parses, and writing it back gives the same bytes.
#[test]
fn every_map_round_trips_to_the_bytes_the_reference_wrote() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "map") {
        let (label, byte_size, map_size) = (row[0], row[1], row[2]);
        let bytes = unhex(row[3]);
        assert_eq!(bytes.len().to_string(), byte_size, "{label}: declared size");

        let map = EncodingMap::read(&bytes).expect("the reference's own map parses");
        assert_eq!(map.len().to_string(), map_size, "{label}: entries");
        assert_eq!(hex(&map.write()), row[3], "{label}: written back");

        let prefixed = map.write_prefixed();
        let (again, consumed) =
            EncodingMap::read_prefixed(&prefixed).expect("the prefixed form parses");
        assert_eq!(consumed, prefixed.len(), "{label}: prefixed length");
        assert_eq!(again, map, "{label}: prefixed round trip");
        compared += 1;
    }
    assert_eq!(compared, 2, "maps compared");
}

/// Every entry, field by field, so a wrong byte says which series it belongs to.
#[test]
fn every_entry_carries_the_encoding_and_parameters_the_reference_recorded() {
    let corpus = corpus();
    let maps: std::collections::HashMap<&str, Vec<u8>> = rows(&corpus, "map")
        .into_iter()
        .map(|row| (row[0], unhex(row[3])))
        .collect();

    let mut compared = 0;
    for row in rows(&corpus, "entry") {
        let (label, index) = (row[0], row[1].parse::<usize>().expect("index"));
        let map = EncodingMap::read(&maps[label]).expect("parses");
        let series = map.series()[index];
        assert_eq!(series.canonical_name(), row[2], "{label}/{index}: series");

        let descriptor = map.get(series).expect("a descriptor for every series");
        assert_eq!(
            (descriptor.id as i32).to_string(),
            row[3],
            "{label}/{index}: encoding id"
        );
        assert_eq!(
            descriptor.parameters.len().to_string(),
            row[4],
            "{label}/{index}: parameter length"
        );
        assert_eq!(
            hex(&descriptor.parameters),
            row[5],
            "{label}/{index}: parameters"
        );
        compared += 1;
    }
    assert_eq!(compared, 52, "entries compared");
}

/// The series htsjdk writes, in the order it writes them, which is the enum's ordinal order and
/// not the alphabetical order its constructor populates in.
#[test]
fn twenty_six_series_are_written_in_ordinal_order() {
    let corpus = corpus();
    let written = rows(&corpus, "written");
    assert_eq!(written.len(), 1);
    let names: Vec<&str> = written[0][0].split(',').collect();
    assert_eq!(names.len(), 26, "series written");

    // Ordinal order: each name's ordinal is strictly greater than the last.
    let mut previous = None;
    for name in &names {
        let series = DataSeries::by_canonical_name(&name.as_bytes().try_into().expect("two bytes"))
            .expect("a known series");
        if let Some(last) = previous {
            assert!(series > last, "{name} follows its predecessor by ordinal");
        }
        previous = Some(series);
    }

    // And the six that are missing are the ones the module names.
    let missing: Vec<&str> = DATA_SERIES
        .iter()
        .map(|(name, _, _)| *name)
        .filter(|name| !names.contains(name))
        .collect();
    assert_eq!(missing, vec!["TC", "TN", "BB", "QQ", "TM", "TV"]);

    // The alphabetical order the constructor uses is a different order, which is the point.
    let mut alphabetical = names.clone();
    alphabetical.sort_unstable();
    assert_ne!(alphabetical, names, "ordinal order is not alphabetical");
}

/// The two series that are parsed and then dropped, so a map can hold fewer entries than its own
/// count declared.
#[test]
fn tc_and_tn_are_read_and_dropped() {
    let corpus = corpus();
    let ignored = rows(&corpus, "ignored");
    assert_eq!(ignored.len(), 1);
    assert_eq!(ignored[0][0].split(',').collect::<Vec<_>>(), NOT_READ);

    let mut map = Vec::new();
    map.push(3u8);
    for name in [b"BF", b"TC", b"TN"] {
        map.extend_from_slice(name);
        map.push(EncodingId::External as u8);
        map.push(1);
        map.push(1);
    }
    let read = EncodingMap::read(&map).expect("parses");
    assert_eq!(read.len(), 1, "three entries declared, one kept");
}

/// The refusals, and the fact that two of them are not CRAM errors at all.
#[test]
fn an_unknown_series_is_a_cram_error_and_an_unknown_encoding_is_an_array_index() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "refuse") {
        let (case, class, message) = (row[0], row[1], row[2]);
        let ours = match case {
            "unknown-canonical-name" => EncodingMapError::UnknownDataSeries(*b"ZZ"),
            "name-with-a-control-byte" => EncodingMapError::UnknownDataSeries([b'Z', 1]),
            "encoding-id-past-the-end" => EncodingMapError::EncodingIdOutOfBounds(10),
            "encoding-id-negative" => EncodingMapError::EncodingIdOutOfBounds(-1),
            other => panic!("{other}: no such case"),
        };
        // The dump escapes anything outside printable ASCII, which the control-byte name needs.
        let escaped: String = ours
            .message()
            .chars()
            .map(|c| {
                if (' '..='~').contains(&c) {
                    c.to_string()
                } else {
                    format!("\\u{:04x}", c as u32)
                }
            })
            .collect();
        assert_eq!(escaped, message, "{case}: message");
        match case {
            "unknown-canonical-name" | "name-with-a-control-byte" => {
                assert_eq!(class, "htsjdk.samtools.cram.CRAMException")
            }
            _ => assert_eq!(class, "java.lang.ArrayIndexOutOfBoundsException"),
        }
        compared += 1;
    }
    assert_eq!(compared, 4, "refusals compared");
}
