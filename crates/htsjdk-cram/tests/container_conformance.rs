//! Conformance for the CRAM file definition and container header, against `CramHeader` and
//! `ContainerHeader`.
//!
//! Goldens from `tools/cram-conformance/CramContainerDump.java` in the pinned oracle.
//!
//! The rows that justify the suite are the ones no reading of the specification produces:
//!
//! ```text
//! id   far-too-long-to-fit-in-twenty  6661722d746f6f2d6c6f6e672d746f2d6669742d
//! hdr  four-unmapped  2  15  -1  4542278  0  0  0  0  1  -  1339669765
//! ```
//!
//! The first is a 29-character id written as 20 bytes with no error. The second is the container
//! that ends every CRAM, whose `alignmentStart` is 4542278, which is `0x454F46`, which is `EOF`.

use std::io::Read;

use htsjdk_cram::container::{container_headers, ContainerHeader, FileDefinition};

fn corpus() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/cram_container.txt.gz");
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

/// The `id` rows: padding and truncation, which the writer does in silence.
#[test]
fn every_file_id_is_padded_and_truncated_as_the_reference_does() {
    let corpus = corpus();
    let mut compared = 0;
    for line in corpus.lines() {
        let Some(rest) = line.strip_prefix("id\t") else {
            continue;
        };
        let (label, expected) = rest.split_once('\t').expect("a label and bytes");
        let id = if label == "<empty>" { "" } else { label };
        assert_eq!(
            hex(&FileDefinition::new(3, 0, id).id),
            expected,
            "id {label:?}"
        );
        compared += 1;
    }
    assert_eq!(compared, 4, "id rows compared");
}

/// The `bytes` rows carry each container header exactly as it sits in the file, so parsing them
/// and rendering the result reproduces the `hdr` row. That is the whole layout under test: a field
/// read at the wrong width shifts every field after it.
#[test]
fn every_container_header_parses_as_the_reference_parses_it() {
    let corpus = corpus();
    let mut compared = 0;

    for line in corpus.lines() {
        let Some(rest) = line.strip_prefix("bytes\t") else {
            continue;
        };
        let fields: Vec<&str> = rest.split('\t').collect();
        let (label, index, bytes) = (fields[0], fields[1], unhex(fields[2]));

        let header = ContainerHeader::read(&bytes, 3).expect("the golden's bytes parse");
        assert_eq!(
            header.byte_length,
            bytes.len(),
            "{label}/{index}: the header consumed {} of {} bytes",
            header.byte_length,
            bytes.len()
        );

        let landmarks = if header.landmarks.is_empty() {
            "-".to_string()
        } else {
            header
                .landmarks
                .iter()
                .map(i32::to_string)
                .collect::<Vec<_>>()
                .join(",")
        };
        let mine = format!(
            "hdr\t{label}\t{index}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{landmarks}\t{}",
            header.blocks_byte_size,
            header.reference_context_id,
            header.alignment_start,
            header.alignment_span,
            header.record_count,
            header.global_record_counter,
            header.base_count,
            header.block_count,
            header.checksum
        );
        let golden = corpus
            .lines()
            .find(|l| l.starts_with(&format!("hdr\t{label}\t{index}\t")))
            .unwrap_or_else(|| panic!("{label}/{index}: no hdr row"));
        assert_eq!(mine, golden, "{label}/{index}");
        compared += 1;
    }

    assert_eq!(compared, 11, "container headers compared");
}

/// The `file` rows: walking a whole file must find the same containers, which is the layout and the
/// skipping arithmetic together.
#[test]
fn walking_a_file_finds_the_containers_the_reference_found() {
    let corpus = corpus();
    let mut compared = 0;

    for line in corpus.lines() {
        let Some(rest) = line.strip_prefix("file\t") else {
            continue;
        };
        let fields: Vec<&str> = rest.split('\t').collect();
        let (label, expected_count) = (fields[0], fields[2].parse::<usize>().expect("a count"));

        // Rebuild the file from its definition and the container headers and blocks the golden
        // recorded, which is enough to walk: the blocks themselves are opaque here.
        let def_row = corpus
            .lines()
            .find(|l| l.starts_with(&format!("def\t{label}\t")))
            .expect("a def row");
        let def_fields: Vec<&str> = def_row.split('\t').collect();
        let version: Vec<&str> = def_fields[3].split('.').collect();
        let definition = FileDefinition {
            major: version[0].parse().expect("a major"),
            minor: version[1].parse().expect("a minor"),
            id: unhex(def_fields[4]).try_into().expect("twenty bytes"),
        };

        let mut cram = definition.write();
        for index in 0..expected_count {
            let bytes_row = corpus
                .lines()
                .find(|l| l.starts_with(&format!("bytes\t{label}\t{index}\t")))
                .expect("a bytes row");
            let header_bytes = unhex(bytes_row.split('\t').nth(3).expect("the bytes"));
            let header = ContainerHeader::read(&header_bytes, definition.major).expect("parses");
            cram.extend_from_slice(&header_bytes);
            cram.extend(std::iter::repeat_n(
                0u8,
                header.blocks_byte_size.max(0) as usize,
            ));
        }

        let (parsed_definition, headers) = container_headers(&cram).expect("the file walks");
        assert_eq!(parsed_definition, definition, "{label}: file definition");
        assert_eq!(headers.len(), expected_count, "{label}: container count");
        assert!(
            headers.last().expect("at least one").is_eof(),
            "{label}: the last container is the EOF one"
        );
        compared += 1;
    }

    assert_eq!(compared, 4, "files walked");
}
