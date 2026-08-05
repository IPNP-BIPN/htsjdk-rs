//! Conformance for the compression header's preservation map, against
//! `htsjdk.samtools.cram.structure.CompressionHeader`.
//!
//! Goldens from `tools/cram-conformance/CramPreservationMapDump.java` in the pinned oracle.
//!
//! The rows that justify the suite:
//!
//! ```text
//! pmap     four-unmapped  21  5  05524e01415000525201534d1b1b1b1b1b54440100
//! boolean  2   false
//! refuse   unknown-key             java.lang.RuntimeException      Unknown preservation map key: ZZ
//! refuse   no-tag-dictionary       htsjdk.samtools.cram.CRAMException  substitution matrix and ...
//! tdgroup  tagged  1  MDZ,NMc,XXf
//! ```
//!
//! The map begins with a hardcoded 5, a boolean is `== 1` so a 2 is false and nothing complains,
//! an unknown key is a plain `RuntimeException` while a missing mandatory key is a `CRAMException`
//! naming both, and the tag dictionary's first group is empty even when every record carries tags.

use std::io::Read;

use htsjdk_cram::preservation_map::{
    parse_dictionary, PreservationMap, PreservationMapError, BASES_SIZE, WRITTEN_MAP_SIZE,
};

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/cram_preservation_map.txt.gz");
    let file = std::fs::File::open(&path).expect("corpus");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("corpus is gzip");
    text
}

fn hex(bytes: &[u8]) -> String {
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

/// The substitution matrix's width, measured rather than assumed.
#[test]
fn the_substitution_matrix_is_five_bytes() {
    let corpus = corpus();
    let sizes = rows(&corpus, "sizes");
    assert_eq!(sizes.len(), 1);
    assert_eq!(sizes[0][0], BASES_SIZE.to_string());
}

/// Every map the reference wrote parses, and writing it back gives the same bytes.
///
/// This is the whole claim: the count is a constant, the order is htsjdk's, and reproducing both
/// is what byte-identity means here.
#[test]
fn every_map_round_trips_to_the_bytes_the_reference_wrote() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "pmap") {
        let (label, byte_size, map_size) = (row[0], row[1], row[2]);
        let bytes = unhex(row[3]);
        assert_eq!(bytes.len().to_string(), byte_size, "{label}: declared size");
        assert_eq!(
            map_size,
            WRITTEN_MAP_SIZE.to_string(),
            "{label}: the count the writer wrote"
        );

        let map = PreservationMap::read(&bytes).expect("the reference's own map parses");
        assert_eq!(hex(&map.write()), row[3], "{label}: written back");
        // The prefixed form is the map plus its own ITF8 length, which is what sits in the header.
        let prefixed = map.write_prefixed();
        let (again, consumed) =
            PreservationMap::read_prefixed(&prefixed).expect("the prefixed form parses");
        assert_eq!(consumed, prefixed.len(), "{label}: prefixed length");
        assert_eq!(again, map, "{label}: prefixed round trip");
        compared += 1;
    }
    assert_eq!(compared, 5, "maps compared");
}

/// The three flags, as the reference's own reader reported them.
#[test]
fn the_three_flags_read_as_the_reference_read_them() {
    let corpus = corpus();
    let maps: std::collections::HashMap<&str, Vec<u8>> = rows(&corpus, "pmap")
        .into_iter()
        .map(|row| (row[0], unhex(row[3])))
        .collect();

    let mut compared = 0;
    for row in rows(&corpus, "flags") {
        let label = row[0];
        let map = PreservationMap::read(maps.get(label).expect("a map for every flags row"))
            .expect("parses");
        assert_eq!(map.preserve_read_names.to_string(), row[1], "{label}: RN");
        assert_eq!(map.ap_delta.to_string(), row[2], "{label}: AP");
        assert_eq!(map.reference_required.to_string(), row[3], "{label}: RR");
        compared += 1;
    }
    assert_eq!(compared, 5, "flag rows compared");
}

/// The substitution matrix and the tag dictionary, as the map carries them.
#[test]
fn the_matrix_and_the_dictionary_are_read_where_the_golden_found_them() {
    let corpus = corpus();
    let maps: std::collections::HashMap<&str, Vec<u8>> = rows(&corpus, "pmap")
        .into_iter()
        .map(|row| (row[0], unhex(row[3])))
        .collect();

    let mut compared = 0;
    for row in rows(&corpus, "sm") {
        let map = PreservationMap::read(&maps[row[0]]).expect("parses");
        assert_eq!(hex(&map.substitution_matrix), row[1], "{}: SM", row[0]);
        compared += 1;
    }
    for row in rows(&corpus, "td") {
        let map = PreservationMap::read(&maps[row[0]]).expect("parses");
        let bytes = unhex(row[2]);
        assert_eq!(bytes.len().to_string(), row[1], "{}: TD length", row[0]);
        assert_eq!(
            parse_dictionary(&bytes),
            map.tag_id_dictionary,
            "{}",
            row[0]
        );
        compared += 1;
    }
    assert_eq!(compared, 10, "matrix and dictionary rows compared");
}

/// The dictionary's groups, and the empty one that is always first.
#[test]
fn the_first_dictionary_group_is_empty_even_when_every_record_has_tags() {
    let corpus = corpus();
    let maps: std::collections::HashMap<&str, Vec<u8>> = rows(&corpus, "pmap")
        .into_iter()
        .map(|row| (row[0], unhex(row[3])))
        .collect();

    let mut compared = 0;
    let mut saw_populated_group = false;
    for row in rows(&corpus, "tdgroup") {
        let (label, index) = (row[0], row[1].parse::<usize>().expect("index"));
        let map = PreservationMap::read(&maps[label]).expect("parses");
        let group = &map.tag_id_dictionary[index];
        let mine: Vec<String> = group
            .iter()
            .map(|id| String::from_utf8_lossy(id).to_string())
            .collect();
        let expected = if mine.is_empty() {
            "-".to_string()
        } else {
            mine.join(",")
        };
        assert_eq!(expected, row[2], "{label}/{index}");
        if index == 0 {
            assert!(
                group.is_empty(),
                "{label}: group 0 is the record with no tags"
            );
        } else {
            saw_populated_group = true;
        }
        compared += 1;
    }
    assert_eq!(compared, 6, "dictionary groups compared");
    assert!(
        saw_populated_group,
        "a file whose records carry tags is in the corpus"
    );
}

/// A boolean is `== 1`: everything else is false, and nothing raises.
#[test]
fn a_boolean_is_one_rather_than_non_zero() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "boolean") {
        let value: u8 = row[0].parse().expect("byte");
        let mut map = Vec::new();
        map.push(WRITTEN_MAP_SIZE as u8);
        map.extend_from_slice(b"RN");
        map.push(value);
        map.extend_from_slice(b"AP");
        map.push(1);
        map.extend_from_slice(b"RR");
        map.push(1);
        map.extend_from_slice(b"SM");
        map.extend_from_slice(&[0u8; BASES_SIZE]);
        map.extend_from_slice(b"TD");
        map.push(4);
        map.extend_from_slice(b"NMi\x00");

        let read = PreservationMap::read(&map).expect("parses");
        assert_eq!(read.preserve_read_names.to_string(), row[1], "RN = {value}");
        compared += 1;
    }
    assert_eq!(compared, 5, "boolean values compared");
}

/// The two refusals, and the fact that they are different kinds of failure.
///
/// An unknown key is a plain `RuntimeException`; a missing mandatory key is a `CRAMException` whose
/// message names both keys whichever one is absent.
#[test]
fn the_refusals_carry_the_messages_the_reference_carries() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "refuse") {
        let (case, class, message) = (row[0], row[1], row[2]);
        let ours = match case {
            "unknown-key" => PreservationMapError::UnknownKey(*b"ZZ"),
            "no-substitution-matrix" | "no-tag-dictionary" => {
                PreservationMapError::MissingMatrixOrDictionary
            }
            other => panic!("{other}: no such case"),
        };
        assert_eq!(ours.message(), message, "{case}: message");
        // The class distinction is the finding: one is a CRAMException and one is not.
        match case {
            "unknown-key" => assert_eq!(class, "java.lang.RuntimeException"),
            _ => assert_eq!(class, "htsjdk.samtools.cram.CRAMException"),
        }
        compared += 1;
    }
    assert_eq!(compared, 3, "refusals compared");
}

/// The compression header does not depend on the reads: four files differing in record count and
/// read length produce the same block, and only the tags move it.
#[test]
fn the_header_is_the_same_for_files_that_differ_only_in_their_reads() {
    let corpus = corpus();
    let blocks = rows(&corpus, "block");
    assert_eq!(blocks.len(), 5);

    let untagged: Vec<&Vec<&str>> = blocks.iter().filter(|row| row[0] != "tagged").collect();
    assert_eq!(untagged.len(), 4);
    for row in &untagged[1..] {
        assert_eq!(row[1], untagged[0][1], "{}: length", row[0]);
        assert_eq!(row[2], untagged[0][2], "{}: digest", row[0]);
    }

    let tagged = blocks
        .iter()
        .find(|row| row[0] == "tagged")
        .expect("tagged");
    assert_ne!(tagged[2], untagged[0][2], "tags move the header");
}
