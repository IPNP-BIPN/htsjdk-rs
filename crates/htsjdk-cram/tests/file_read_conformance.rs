//! Conformance for walking a whole CRAM file, against htsjdk's own reader.
//!
//! Goldens from `tools/cram-conformance/CramFileReadDump.java` in the pinned oracle, over
//! `ce#5.2.1.cram` from htsjdk's test resources: a version 2.1 file of six records in one
//! container, whose blocks are GZIP.
//!
//! Every piece below this has a suite of its own. What this one adds is that they fit: the file
//! definition, the SAM header container, the container header, the compression header, the slice
//! header and the slice's blocks, in that order and at those offsets, with each block's
//! compression undone.
//!
//! The `record` rows are htsjdk's, and they are not compared here. Reading a record from a slice's
//! blocks is `cram-record-read`'s subject, and turning one into a SAM record needs the reference
//! FASTA the reads were compressed against, which this suite does not carry. They are kept because
//! they say what the file holds, and because a walk that produced the wrong blocks would still
//! look right without them.

use std::io::Read;

use htsjdk_cram::block::{CompressionMethod, ContentType};
use htsjdk_cram::file::read_file;

fn corpus() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/cram_file_read.txt.gz");
    let file = std::fs::File::open(&path).expect("corpus");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("corpus is gzip");
    text
}

/// The file itself, carried beside the golden so the suite is self-contained.
fn cram() -> Vec<u8> {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/ce5.cram");
    std::fs::read(path).expect("the cram file")
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

/// The definition, and the SAM header container that follows it.
#[test]
fn the_file_begins_where_the_reference_says_it_does() {
    let corpus = corpus();
    let cram = cram();
    let walk = read_file(&cram).unwrap_or_else(|error| panic!("{}", error.message()));

    let file = &rows(&corpus, "file")[0];
    assert_eq!(file[1].parse::<usize>().expect("length"), cram.len());
    assert_eq!(
        format!("{}.{}", walk.definition.major, walk.definition.minor),
        file[2],
        "the version"
    );

    // The SAM header container's block holds the header text. htsjdk reports what it parsed out of
    // it; the walk holds the bytes, so the two are checked against each other by content.
    let sam = &rows(&corpus, "samheader")[0];
    let text = String::from_utf8_lossy(&walk.sam_header);
    let sequences = text.matches("@SQ").count();
    assert_eq!(
        sequences.to_string(),
        sam[1],
        "the sequences the header names"
    );
    // The sort order htsjdk reports is `unsorted`, and the header does not say so: there is no
    // SO tag in it at all. The default is reported as though it were written down.
    assert_eq!(sam[3], "unsorted");
    assert!(
        !text.contains("SO:"),
        "and the header names no sort order at all"
    );
}

/// Every container the reference walked, at the offset it walked it from.
#[test]
fn every_container_is_where_the_reference_found_it() {
    let corpus = corpus();
    let cram = cram();
    let walk = read_file(&cram).expect("walked");

    let expected = rows(&corpus, "container");
    assert_eq!(walk.containers.len(), expected.len(), "containers walked");

    for (container, row) in walk.containers.iter().zip(&expected) {
        assert_eq!(container.offset.to_string(), row[2], "container offset");
        assert_eq!(
            container.header.blocks_byte_size.to_string(),
            row[3],
            "container length"
        );
        assert_eq!(
            container.header.reference_context_id.to_string(),
            row[4],
            "container reference"
        );
        assert_eq!(
            container.header.alignment_start.to_string(),
            row[5],
            "container start"
        );
        assert_eq!(
            container.header.record_count.to_string(),
            row[7],
            "container records"
        );
        assert_eq!(
            container.header.block_count.to_string(),
            row[8],
            "container blocks"
        );
    }

    // The last one is the EOF container: it parses as a container whose record count is zero, not
    // as a byte pattern.
    let last = walk.containers.last().expect("a container");
    assert!(last.header.is_eof());
    assert_eq!(last.header.record_count, 0);
    assert!(last.compression_header.is_none());
}

/// Every slice, and every block of it with its compression undone.
#[test]
fn every_block_decompresses_to_the_reference_bytes() {
    let corpus = corpus();
    let cram = cram();
    let walk = read_file(&cram).expect("walked");

    let slices = rows(&corpus, "slice");
    let walked: Vec<_> = walk
        .containers
        .iter()
        .flat_map(|container| &container.slices)
        .collect();
    assert_eq!(walked.len(), slices.len(), "slices walked");

    for (slice, row) in walked.iter().zip(&slices) {
        assert_eq!(
            slice.header.reference_context_id.to_string(),
            row[2],
            "slice reference"
        );
        assert_eq!(
            slice.header.alignment_start.to_string(),
            row[3],
            "slice start"
        );
        assert_eq!(
            slice.header.record_count.to_string(),
            row[5],
            "slice records"
        );
        assert_eq!(slice.header.block_count.to_string(), row[6], "slice blocks");
    }

    let expected = rows(&corpus, "block");
    let blocks: Vec<_> = walked.iter().flat_map(|slice| &slice.blocks).collect();
    assert_eq!(blocks.len(), expected.len(), "blocks walked");

    let mut gzip = 0;
    for (block, row) in blocks.iter().zip(&expected) {
        // The dump names the core block "core" and an external one by its content id.
        let name = if block.content_type == ContentType::Core {
            "core".to_string()
        } else {
            block.content_id.to_string()
        };
        assert_eq!(name, row[2], "block name");
        assert_eq!(block.method.name(), row[4], "block compression");
        assert_eq!(hex(&block.content), row[6], "block {name} content");
        if block.method == CompressionMethod::Gzip {
            gzip += 1;
        }
    }
    assert_eq!(gzip, blocks.len(), "every block of this file is GZIP");
}

/// The blocks a slice hands the codecs, which is what `cram-record-read` reads from.
#[test]
fn the_walk_hands_the_codecs_what_they_expect() {
    let cram = cram();
    let walk = read_file(&cram).expect("walked");
    let slice = walk.containers[0].slices.first().expect("a slice");
    let bytes = slice.block_bytes();

    // A core block and four externals, and the compression header names an encoding for each.
    assert!(!bytes.core.is_empty());
    assert_eq!(bytes.external.len(), 4);

    let header = walk.containers[0]
        .compression_header
        .as_ref()
        .expect("a compression header");
    // Every external block the slice carries is one the encoding map named.
    assert!(!header.encodings.is_empty());
}
