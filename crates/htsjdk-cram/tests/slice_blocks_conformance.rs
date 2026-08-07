//! Conformance for a slice's blocks, against `htsjdk.samtools.cram.structure.SliceBlocks`.
//!
//! Goldens from `tools/cram-conformance/CramSliceBlocksDump.java` in the pinned oracle.
//!
//! The rows that justify the suite:
//!
//! ```text
//! write  3.0  3,2,1        1,2,3        0005000202aa55...
//! write  3.0  300,2,128    2,128,300    0005000202aa55...
//! read   3.0  core-last    4            1,2
//! err    duplicate  id=4  CRAMException  Attempt to add a duplicate block (id 4 of type EXTERNAL) to compression header encoding map. Existing block is of type EXTERNAL.
//! ```
//!
//! The written order is by content id and not by insertion. The reader takes a count rather than
//! an order, so a stream whose core block comes last reads the same as one whose core block comes
//! first. And a repeated content id is refused rather than overwritten.

use std::io::Read;

use htsjdk_cram::external_codecs::SliceBlockBytes;
use htsjdk_cram::slice_blocks::{from_blocks, read_blocks, write_blocks, SliceBlocksError};

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/cram_slice_blocks.txt.gz");
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

fn major(version: &str) -> u8 {
    version
        .split('.')
        .next()
        .expect("major")
        .parse()
        .expect("major")
}

/// The dump's external block for a content id: its own id and the next byte.
fn external(content_id: i32) -> (i32, Vec<u8>) {
    (content_id, vec![content_id as u8, (content_id + 1) as u8])
}

/// Every stream the reference wrote, with the externals in ascending content id whatever order
/// they were added in.
#[test]
fn the_blocks_are_written_in_the_reference_order() {
    let corpus = corpus();
    let mut compared = 0;
    let mut reordered = 0;
    for row in rows(&corpus, "write") {
        let (version, added, written, expected) = (row[0], row[1], row[2], row[3]);
        let added_ids: Vec<i32> = added
            .split(',')
            .map(|value| value.parse().expect("content id"))
            .collect();

        let blocks = from_blocks(
            vec![0xAA, 0x55],
            added_ids.iter().map(|id| external(*id)).collect(),
        )
        .expect("built");

        // The map's own order is what the reference wrote, and it is not the order given.
        let order: Vec<String> = blocks.external.keys().map(|id| id.to_string()).collect();
        assert_eq!(order.join(","), written, "write {version} {added}");
        if added != written {
            reordered += 1;
        }

        assert_eq!(
            hex(&write_blocks(&blocks, major(version))),
            expected,
            "write {version} {added}"
        );
        compared += 1;
    }
    assert_eq!(compared, 7, "streams written");
    assert!(reordered >= 2, "and at least two were written out of order");
}

/// A stream read back, whether its core block came first or last.
#[test]
fn a_stream_reads_back_whatever_order_it_is_in() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "read") {
        let (version, arrangement, core_size, ids) = (
            row[0],
            row[1],
            row[2].parse::<usize>().expect("core size"),
            row[3],
        );

        // The dump's stream: a four-byte core block and the externals 2 and 1, in that order.
        let core = vec![1u8, 2, 3, 4];
        let externals = [external(2), external(1)];
        let mut stream = Vec::new();
        if arrangement == "core-first" {
            stream.extend_from_slice(&write_one_core(&core, major(version)));
        }
        for (content_id, content) in &externals {
            stream.extend_from_slice(&write_one_external(*content_id, content, major(version)));
        }
        if arrangement == "core-last" {
            stream.extend_from_slice(&write_one_core(&core, major(version)));
        }

        let blocks = read_blocks(&stream, 3, major(version)).expect("read");
        assert_eq!(blocks.core.len(), core_size, "read {version} {arrangement}");
        let order: Vec<String> = blocks.external.keys().map(|id| id.to_string()).collect();
        assert_eq!(order.join(","), ids, "read {version} {arrangement}");
        assert_eq!(blocks.core, core);
        compared += 1;
    }
    assert_eq!(compared, 3, "streams read");
}

/// One block on its own, through the same writer the whole stream uses.
fn write_one_core(content: &[u8], major: u8) -> Vec<u8> {
    write_blocks(
        &SliceBlockBytes {
            core: content.to_vec(),
            external: Default::default(),
        },
        major,
    )
}

fn write_one_external(content_id: i32, content: &[u8], major: u8) -> Vec<u8> {
    let whole = write_blocks(
        &from_blocks(Vec::new(), vec![(content_id, content.to_vec())]).expect("built"),
        major,
    );
    // Drop the empty core block the writer puts first: its header is five bytes, plus a checksum
    // from version 3 on.
    let core_length = 5 + if major >= 3 { 4 } else { 0 };
    whole[core_length..].to_vec()
}

/// What building or reading a slice's blocks refuses.
#[test]
fn the_failures_are_the_reference_failures() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "err") {
        let (what, class, message) = (row[0], row[2], row[3]);
        let error = match what {
            "duplicate" => {
                from_blocks(vec![1], vec![external(4), external(4)]).expect_err("refused")
            }
            "no-core" => {
                let mut stream = write_one_external(1, &[1, 2], 3);
                stream.extend_from_slice(&write_one_external(2, &[2, 3], 3));
                read_blocks(&stream, 2, 3).expect_err("refused")
            }
            "wrong-type" => {
                // A raw compression header block: method 0, type 1, no content id, sizes of 1.
                let mut block = vec![0u8, 1, 0, 1, 1, 1];
                let checksum = htsjdk_cram::compression_header::crc32(&block);
                block.extend_from_slice(&checksum.to_le_bytes());
                read_blocks(&block, 1, 3).expect_err("refused")
            }
            other => panic!("{other}"),
        };
        assert_eq!(error.java_exception(), class, "err {what}");
        assert_eq!(error.message(), message, "err {what}");
        compared += 1;
    }
    assert_eq!(compared, 3, "refusals compared");

    // The core block is looked for after every block has been read, not while reading them, so a
    // stream of externals alone is read in full before it is refused.
    assert!(matches!(
        read_blocks(&write_one_external(1, &[1], 3), 1, 3),
        Err(SliceBlocksError::NoCoreBlock)
    ));
}
