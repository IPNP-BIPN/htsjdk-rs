//! Conformance for the tag encoding map, against `htsjdk.samtools.cram.structure.CompressionHeader`
//! and `ReadTag`.
//!
//! Goldens from `tools/cram-conformance/CramTagEncodingMapDump.java` in the pinned oracle.
//!
//! The rows that justify the suite:
//!
//! ```text
//! order   tagged               MDZ,NMc,XXf
//! order   reverse-order        MDZ,NMc,XXf
//! order   same-name-two-types  XXZ,XXi
//! order   large-integer        NMi
//! guard   32                   tagID 32 overlaps with data series content ID
//! ```
//!
//! Two files whose records introduce the same three tags in opposite orders produce the same map,
//! because the order is the packed key's. One name at two types is two entries. An integer of 1 to
//! 4 is written as `NMc` and one of 100000 as `NMi`, from the same Java type. And the collision
//! guard fires only over a range no printable tag can reach.

use std::io::Read;

use htsjdk_cram::tag_encoding_map::{
    int_to_name_type_3, int_to_name_type_4, name_type_to_int, overlap_message,
    overlaps_data_series_content_id, TagEncodingMap,
};

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/cram_tag_encoding_map.txt.gz");
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

/// The dump escapes anything outside printable ASCII, and a tag of three spaces round trips
/// through a form that has none of them escaped but is still worth comparing that way.
fn escape(text: &[u8]) -> String {
    text.iter()
        .map(|byte| {
            if (0x20..=0x7e).contains(byte) {
                (*byte as char).to_string()
            } else {
                format!("\\u{:04x}", *byte as u32)
            }
        })
        .collect()
}

/// The key is the tag: two name bytes and the type, packed into twenty-four bits.
#[test]
fn every_tag_packs_into_the_integer_the_reference_packs_it_into() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "id") {
        let name: [u8; 2] = row[0].as_bytes().try_into().expect("two bytes");
        let tag_type = row[1].as_bytes()[0];
        let key = name_type_to_int(name, tag_type);
        assert_eq!(key.to_string(), row[2], "{}{}", row[0], row[1]);
        assert_eq!(escape(&int_to_name_type_3(key)), row[3], "three-byte form");
        assert_eq!(escape(&int_to_name_type_4(key)), row[4], "four-byte form");
        compared += 1;
    }
    assert_eq!(compared, 9, "tag ids compared");
}

/// Every map the reference wrote parses, and writing it back gives the same bytes.
#[test]
fn every_map_round_trips_to_the_bytes_the_reference_wrote() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "tmap") {
        let (label, byte_size, map_size) = (row[0], row[1], row[2]);
        let bytes = unhex(row[3]);
        assert_eq!(bytes.len().to_string(), byte_size, "{label}: declared size");

        let map = TagEncodingMap::read(&bytes).expect("the reference's own map parses");
        assert_eq!(map.len().to_string(), map_size, "{label}: entries");
        assert_eq!(hex(&map.write()), row[3], "{label}: written back");

        let prefixed = map.write_prefixed();
        let (again, consumed) =
            TagEncodingMap::read_prefixed(&prefixed).expect("the prefixed form parses");
        assert_eq!(consumed, prefixed.len(), "{label}: prefixed length");
        assert_eq!(again, map, "{label}: prefixed round trip");
        compared += 1;
    }
    assert_eq!(compared, 5, "maps compared");
}

/// Every entry, field by field.
#[test]
fn every_entry_carries_the_key_encoding_and_parameters_the_reference_recorded() {
    let corpus = corpus();
    let maps: std::collections::HashMap<&str, Vec<u8>> = rows(&corpus, "tmap")
        .into_iter()
        .map(|row| (row[0], unhex(row[3])))
        .collect();

    let mut compared = 0;
    for row in rows(&corpus, "tentry") {
        let (label, index) = (row[0], row[1].parse::<usize>().expect("index"));
        let map = TagEncodingMap::read(&maps[label]).expect("parses");
        let key = map.tag_ids()[index];
        assert_eq!(key.to_string(), row[2], "{label}/{index}: key");
        assert_eq!(
            escape(&int_to_name_type_3(key)),
            row[3],
            "{label}/{index}: name and type"
        );

        let descriptor = map.get(key).expect("a descriptor for every key");
        assert_eq!(
            (descriptor.id as i32).to_string(),
            row[4],
            "{label}/{index}: encoding id"
        );
        assert_eq!(
            descriptor.parameters.len().to_string(),
            row[5],
            "{label}/{index}: parameter length"
        );
        assert_eq!(
            hex(&descriptor.parameters),
            row[6],
            "{label}/{index}: parameters"
        );
        compared += 1;
    }
    assert_eq!(compared, 9, "entries compared");
}

/// The order the tags come out in is the packed key's, not the order the records introduced them.
///
/// Two files that differ only in that order produce byte-identical maps, which is the assertion.
#[test]
fn the_write_order_is_the_keys_and_not_the_records() {
    let corpus = corpus();
    let orders: std::collections::HashMap<&str, &str> = rows(&corpus, "order")
        .into_iter()
        .map(|row| (row[0], row[1]))
        .collect();
    let maps: std::collections::HashMap<&str, &str> = rows(&corpus, "tmap")
        .into_iter()
        .map(|row| (row[0], row[3]))
        .collect();

    assert_eq!(orders["tagged"], "MDZ,NMc,XXf");
    assert_eq!(
        orders["reverse-order"], orders["tagged"],
        "the introduction order is not in the file"
    );
    assert_eq!(
        maps["reverse-order"], maps["tagged"],
        "and the two files' maps are byte-identical"
    );

    // One name at two types is two entries, ordered by their type character.
    assert_eq!(orders["same-name-two-types"], "XXZ,XXi");
    // An integer narrows to the smallest type that holds it, so the same Java type gives two keys.
    assert_eq!(orders["large-integer"], "NMi");
    assert_eq!(orders["untagged"], "-");

    // And every order row is what our own reader gives back.
    let mut compared = 0;
    for (label, expected) in &orders {
        let map = TagEncodingMap::read(&unhex(maps[label])).expect("parses");
        let mine: Vec<String> = map
            .tag_ids()
            .iter()
            .map(|key| escape(&int_to_name_type_3(*key)))
            .collect();
        let joined = if mine.is_empty() {
            "-".to_string()
        } else {
            mine.join(",")
        };
        assert_eq!(&joined, expected, "{label}");
        compared += 1;
    }
    assert_eq!(compared, 5, "orders compared");
}

/// The guard over a range no printable tag can enter.
#[test]
fn the_collision_guard_fires_only_where_no_tag_can_reach() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "guard") {
        let tag_id: i32 = row[0].parse().expect("tag id");
        let mine = if overlaps_data_series_content_id(tag_id) {
            overlap_message(tag_id)
        } else {
            "accepted".to_string()
        };
        assert_eq!(mine, row[1], "tag id {tag_id}");
        compared += 1;
    }
    assert_eq!(compared, 4, "guard answers compared");

    // The smallest printable tag is far above the range the guard covers.
    assert!(!overlaps_data_series_content_id(name_type_to_int(
        *b"  ", b' '
    )));
}
