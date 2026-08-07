//! Conformance for codec negotiation, against `CompressionHeaderFactory` and the default
//! `CompressionHeaderEncodingMap`.
//!
//! Goldens from `tools/cram-conformance/CramNegotiationDump.java` in the pinned oracle.
//!
//! The rows that justify the suite:
//!
//! ```text
//! compressor  all-the-same  1000 bytes  29  29  36  RANS
//! tag         XBZ  61626300,646566676800  BYTE_ARRAY_STOP  09e058425a
//! tag         XBB  630100000001,63020000000102  BYTE_ARRAY_LEN  0104e05842420104e0584242
//! dictionary  same-tags-opposite-order  XAc,XBc;XBc,XAc  .;XAcXBc  indexes=1,1
//! ```
//!
//! A tie between GZIP and rANS goes to rANS. A `Z` of two sizes becomes a stop-byte array whose
//! stop byte is a tab. A `B` of two sizes under the threshold becomes a length-prefixed array whose
//! two halves are both external on the same content id. And two records whose tags differ only in
//! order share one dictionary entry.
//!
//! The gzip lengths are the reference's: htsjdk compresses with the JDK's `Deflater`, whose output
//! length is its zlib's business and not this crate's. Every row carries all three lengths, so the
//! rule is checked without one.

use std::io::Read;

use htsjdk_cram::encoding_map::EncodingId;
use htsjdk_cram::negotiation::{
    best_compressor, tag_encoding, tag_id_as_int, tag_id_dictionary, unused_byte, Compressor,
};

fn corpus() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/cram_negotiation.txt.gz");
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

/// The compressor, chosen by the three lengths the reference measured.
#[test]
fn the_compressor_is_the_reference_choice() {
    let corpus = corpus();
    let mut compared = 0;
    let mut ties = 0;
    for row in rows(&corpus, "compressor") {
        let (label, gzip, rans0, rans1, chosen) = (
            row[0],
            row[2].parse::<usize>().expect("gzip length"),
            row[3].parse::<usize>().expect("rans0 length"),
            row[4].parse::<usize>().expect("rans1 length"),
            row[5],
        );
        assert_eq!(
            best_compressor(gzip, rans0, rans1).name(),
            chosen,
            "compressor {label}"
        );
        if gzip == rans0 || gzip == rans1 {
            // A tie between GZIP and either rANS goes to rANS, which is the order of the
            // comparisons rather than of the compressions.
            assert_eq!(chosen, "RANS", "compressor {label} is a tie");
            ties += 1;
        }
        compared += 1;
    }
    assert_eq!(compared, 7, "compressor choices compared");
    assert_eq!(ties, 1, "of them a tie");

    // And the rule stated on its own, at the boundary the corpus does not have.
    assert_eq!(best_compressor(10, 10, 10), Compressor::Rans);
    assert_eq!(best_compressor(9, 10, 10), Compressor::Gzip);
    assert_eq!(best_compressor(10, 11, 10), Compressor::Rans);
}

/// A tag's encoding, from its type and the sizes of its values.
#[test]
fn a_tag_gets_the_reference_encoding() {
    let corpus = corpus();
    let mut compared = 0;
    let mut stop_encodings = 0;
    for row in rows(&corpus, "tag") {
        let (name, values, id, parameters) = (row[0], row[1], row[2], row[3]);
        let tag_type = name.as_bytes()[2];
        let tag_id = tag_id_as_int(&[name.as_bytes()[0], name.as_bytes()[1], tag_type]);

        let values: Vec<Vec<u8>> = values.split(',').map(unhex).collect();
        let sizes: Vec<usize> = values.iter().map(|value| value.len()).collect();
        let data: Vec<u8> = values.concat();

        let encoding = tag_encoding(tag_type, tag_id, &sizes, &data);
        assert_eq!(encoding_name(encoding.id), id, "tag {name}");
        assert_eq!(hex(&encoding.parameters), parameters, "tag {name}");
        if encoding.id == EncodingId::ByteArrayStop {
            stop_encodings += 1;
        }
        compared += 1;
    }
    assert_eq!(compared, 11, "tag encodings compared");
    assert_eq!(stop_encodings, 2, "of them stop-byte arrays");
}

fn encoding_name(id: EncodingId) -> &'static str {
    id.name()
}

/// The stop byte of a `Z` is a tab, chosen rather than searched for, and the one a long `B` gets is
/// searched for.
#[test]
fn the_stop_byte_is_chosen_for_one_type_and_found_for_the_other() {
    let corpus = corpus();

    // A Z of two sizes: the parameters begin with the stop byte, and it is 0x09.
    let z = rows(&corpus, "tag")
        .into_iter()
        .find(|row| row[0] == "XBZ")
        .expect("the two-size Z");
    assert_eq!(&unhex(z[3])[..1], b"\t");

    // A B whose values are both over the threshold: the stop byte is one its data never uses.
    let b = rows(&corpus, "tag")
        .into_iter()
        .find(|row| row[0] == "XCB")
        .expect("the long B");
    let data: Vec<u8> = b[1].split(',').map(unhex).collect::<Vec<_>>().concat();
    assert_eq!(i32::from(unhex(b[3])[0]), unused_byte(&data));

    // And the search returns -1 only when every byte value appears.
    let all: Vec<u8> = (0..=255u8).collect();
    assert_eq!(unused_byte(&all), -1);
    assert_eq!(unused_byte(&[]), 0);
}

/// The dictionary, and the index each record is given.
#[test]
fn the_dictionary_is_the_reference_dictionary() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "dictionary") {
        let (label, input, expected, indexes) = (row[0], row[1], row[2], row[3]);

        let records: Vec<Vec<[u8; 3]>> = input
            .split(';')
            .map(|record| {
                record
                    .split(',')
                    .map(|name| {
                        let bytes = name.as_bytes();
                        [bytes[0], bytes[1], bytes[2]]
                    })
                    .collect()
            })
            .collect();

        let (groups, given) = tag_id_dictionary(&records);
        let shown = groups
            .iter()
            .map(|group| {
                if group.is_empty() {
                    ".".to_string()
                } else {
                    group
                        .iter()
                        .map(|id| String::from_utf8_lossy(id).to_string())
                        .collect()
                }
            })
            .collect::<Vec<_>>()
            .join(";");
        assert_eq!(shown, expected, "dictionary {label}");
        assert_eq!(
            format!(
                "indexes={}",
                given
                    .iter()
                    .map(|index| index.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            ),
            indexes,
            "dictionary {label} indexes"
        );
        compared += 1;
    }
    assert_eq!(compared, 4, "dictionaries compared");

    // Group 0 is the empty one every container carries, so a record's index is never zero unless
    // it has no tags at all.
    let (groups, indexes) = tag_id_dictionary(&[Vec::new(), vec![*b"XAc"]]);
    assert!(groups[0].is_empty());
    assert_eq!(indexes, [0, 1]);
}

/// The fixed table: every data series the default strategy names, and the encoding it names it
/// with. Nothing here depends on the records.
#[test]
fn the_data_series_table_is_fixed() {
    let corpus = corpus();
    let series = rows(&corpus, "series");
    assert_eq!(series.len(), 32, "data series named");

    // Three series are stop-byte arrays and the rest external. The three are the ones whose
    // values have no fixed width and no length of their own: the read name, an insertion and a
    // soft clip.
    let external = series.iter().filter(|row| row[2] == "EXTERNAL").count();
    let stop: Vec<&str> = series
        .iter()
        .filter(|row| row[2] == "BYTE_ARRAY_STOP")
        .map(|row| row[0])
        .collect();
    assert_eq!(stop, ["RN", "IN", "SC"]);

    // And six of the thirty-two are not in the map at all. Three are the series htsjdk does not
    // read, two are the ones it does not write, and a reader that expects every series to be
    // named finds nothing for them rather than a default.
    let absent: Vec<&str> = series
        .iter()
        .filter(|row| row[2] == "-")
        .map(|row| row[0])
        .collect();
    assert_eq!(absent, ["TC", "TN", "BB", "QQ", "TM", "TV"]);
    assert_eq!(external + stop.len() + absent.len(), series.len());

    // And the content ids are handed out in the table's own order, starting at one.
    let first: Vec<&str> = series.iter().take(6).map(|row| row[3]).collect();
    assert_eq!(first, ["01", "02", "03", "04", "05", "06"]);
}
