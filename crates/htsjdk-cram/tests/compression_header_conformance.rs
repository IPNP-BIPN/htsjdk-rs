//! Conformance for the compression header as a whole, against
//! `htsjdk.samtools.cram.structure.CompressionHeader`.
//!
//! Goldens from `tools/cram-conformance/CramCompressionHeaderDump.java` in the pinned oracle,
//! which writes a header, reads it back, and writes it again.
//!
//! The rows that justify the suite:
//!
//! ```text
//! header  2.1  true  true  true  OQZ;XAZ  0001008...  (174 bytes)
//! header  3.0  true  true  true  OQZ;XAZ  0001008...  (178 bytes)
//! back    3.0  true  true  true  OQZ;XAZ  same        -
//! err     wrong-block-type  MAPPED_SLICE_HEADER  RuntimeIOException  Compression header block expected, found: MAPPED_SLICE
//! ```
//!
//! The same header is four bytes longer in version 3, and those four are a CRC-32 over everything
//! before them. A header read and written again is byte-identical. And a block of the wrong kind
//! is refused with a message naming what was found.

use std::io::Read;

use htsjdk_cram::compression_header::{crc32, CompressionHeader};
use htsjdk_cram::preservation_map::PreservationMap;

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/cram_compression_header.txt.gz");
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
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

fn unhex(text: &str) -> Vec<u8> {
    (0..text.len())
        .step_by(2)
        .map(|at| u8::from_str_radix(&text[at..at + 2], 16).expect("hex"))
        .collect()
}

fn major(version: &str) -> u8 {
    version
        .split('.')
        .next()
        .expect("major")
        .parse()
        .expect("major")
}

/// Every header the reference wrote reads back into the same three maps, and writing them again
/// gives the same bytes.
#[test]
fn a_header_read_and_written_again_is_the_reference_block() {
    let corpus = corpus();
    let headers = rows(&corpus, "header");
    let backs = rows(&corpus, "back");
    assert_eq!(headers.len(), 6, "headers written");
    assert_eq!(backs.len(), headers.len(), "and read back");

    for (header, back) in headers.iter().zip(&backs) {
        let (version, block) = (header[0], unhex(header[5]));
        let parsed = CompressionHeader::read_block(&block, major(version))
            .unwrap_or_else(|error| panic!("{version}: {}", error.message()));

        // The three flags, as the reference reported them after reading.
        assert_eq!(
            parsed.preservation.preserve_read_names.to_string(),
            back[1],
            "{version} RN"
        );
        assert_eq!(
            parsed.preservation.ap_delta.to_string(),
            back[2],
            "{version} AP"
        );
        assert_eq!(
            parsed.preservation.reference_required.to_string(),
            back[3],
            "{version} RR"
        );
        assert_eq!(tags(&parsed.preservation), back[4], "{version} TD");

        // And the block written again. The reference recorded "same"; anything else would be its
        // hex, and the port has to produce whichever it recorded.
        let again = parsed.write_block(major(version));
        let expected = if back[5] == "same" {
            header[5]
        } else {
            back[5]
        };
        assert_eq!(hex(&again), expected, "{version} rewritten");
    }
}

fn tags(preservation: &PreservationMap) -> String {
    if preservation.tag_id_dictionary.is_empty() {
        return "-".to_string();
    }
    preservation
        .tag_id_dictionary
        .iter()
        .map(|group| {
            if group.is_empty() {
                ".".to_string()
            } else {
                group
                    .iter()
                    .map(|tag| String::from_utf8_lossy(tag).to_string())
                    .collect()
            }
        })
        .collect::<Vec<_>>()
        .join(";")
}

/// The version changes the block and nothing inside it: four bytes of checksum, and the content
/// identical up to them.
#[test]
fn the_checksum_is_the_only_difference_between_the_two_versions() {
    let corpus = corpus();
    let headers = rows(&corpus, "header");
    let (v21, v30): (Vec<_>, Vec<_>) = headers.iter().partition(|row| row[0] == "2.1");
    assert_eq!(v21.len(), 2);

    for (old, new) in v21.iter().zip(&v30) {
        let (old, new) = (unhex(old[5]), unhex(new[5]));
        assert_eq!(new.len(), old.len() + 4, "version 3 adds four bytes");
        assert_eq!(new[..old.len()], old[..], "and changes nothing before them");

        // Those four are a CRC-32 over everything before them, little-endian.
        let checksum = u32::from_le_bytes(new[old.len()..].try_into().expect("four bytes"));
        assert_eq!(
            checksum,
            crc32(&old),
            "the checksum covers header and content"
        );
    }

    // The checksum differs between two headers that differ only in their flags, which is the point
    // of carrying it.
    let checksums: Vec<u32> = v30
        .iter()
        .map(|row| {
            let block = unhex(row[5]);
            u32::from_le_bytes(block[block.len() - 4..].try_into().expect("four bytes"))
        })
        .collect();
    assert_ne!(checksums[0], checksums[1]);
}

/// The block's own fields, taken out of the content the reference reported.
#[test]
fn the_block_is_raw_and_says_it_carries_a_compression_header() {
    let corpus = corpus();
    let sections = rows(&corpus, "section");
    let content = sections
        .iter()
        .find(|row| row[0] == "content")
        .map(|row| unhex(row[1]))
        .expect("content");
    assert_eq!(
        sections
            .iter()
            .find(|row| row[0] == "content-type")
            .expect("type")[1],
        "COMPRESSION_HEADER"
    );
    assert_eq!(
        sections
            .iter()
            .find(|row| row[0] == "compression")
            .expect("method")[1],
        "RAW"
    );

    // The content parses into the three maps, and writing them gives it back byte for byte.
    let parsed = CompressionHeader::read_content(&content).expect("three maps");
    assert_eq!(hex(&parsed.write_content()), hex(&content));

    // And the block the port writes around that content is what the reference's header row shows.
    let header = rows(&corpus, "header")
        .into_iter()
        .find(|row| row[0] == "3.0" && row[1] == "true" && row[2] == "true")
        .expect("the same header");
    assert_eq!(hex(&parsed.write_block(3)), header[5]);
}

/// What reading a compression header refuses.
#[test]
fn the_failures_are_the_reference_failures() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "err") {
        let (what, detail, class, message) = (row[0], row[1], row[2], row[3]);
        let error = match what {
            // The detail is the content that was wrapped in a block, so the row carries its own
            // input rather than a label to rebuild one from.
            "no-substitution-matrix" | "no-tag-dictionary" | "unknown-key" => {
                CompressionHeader::read_content(&unhex(detail)).expect_err("refused")
            }
            "wrong-block-type" => {
                // A raw slice header block: method 0, content type 2, no content id, sizes of 3.
                let mut block = vec![0u8, 2, 0, 3, 3, 1, 2, 3];
                let checksum = crc32(&block);
                block.extend_from_slice(&checksum.to_le_bytes());
                CompressionHeader::read_block(&block, 3).expect_err("refused")
            }
            other => panic!("{other}"),
        };
        assert_eq!(error.java_exception(), class, "err {what}");
        assert_eq!(error.message(), message, "err {what}");
        compared += 1;
    }
    assert_eq!(compared, 4, "refusals compared");
}
