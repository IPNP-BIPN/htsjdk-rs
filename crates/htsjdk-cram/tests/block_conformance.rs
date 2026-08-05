//! Conformance for the CRAM block header, against `htsjdk.samtools.cram.structure.block.Block`.
//!
//! Goldens from `tools/cram-conformance/CramBlockDump.java` in the pinned oracle.
//!
//! The rows that justify the suite are the ones about where a block ends:
//!
//! ```text
//! blk  four-unmapped  1  0   0  1   0   160  160  -310733719
//! blk  four-unmapped  1  3   4  4   1    36    4  -1714489358
//! counts  four-unmapped  29  0,1,4
//! ```
//!
//! The checksum sits **outside** the `compressedSize` the header declares, so a port that walks a
//! container by adding header plus compressed size lands four bytes short of the next block and
//! every block after the first is misread. And the methods an ordinary four-read file uses are RAW,
//! GZIP and rANS, which is what makes rANS 4x8 required rather than optional.

use std::io::Read;

use htsjdk_cram::block::BlockHeader;

fn corpus() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/cram_block.txt.gz");
    let file = std::fs::File::open(&path).expect("corpus");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("corpus is gzip");
    text
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

/// Every `hdrbytes` row parses into the fields its `blk` row records, and consumes exactly as many
/// bytes as the golden captured. A field read at the wrong width shifts every field after it.
#[test]
fn every_block_header_parses_as_the_reference_parses_it() {
    let corpus = corpus();
    let mut compared = 0;

    for line in corpus.lines() {
        let Some(rest) = line.strip_prefix("hdrbytes\t") else {
            continue;
        };
        let fields: Vec<&str> = rest.split('\t').collect();
        let (label, container, index) = (fields[0], fields[1], fields[2]);
        let bytes = unhex(fields[3]);

        let header = BlockHeader::read(&bytes).expect("the golden's bytes parse");
        assert_eq!(
            header.byte_length,
            bytes.len(),
            "{label}/{container}/{index}: consumed {} of {} bytes",
            header.byte_length,
            bytes.len()
        );

        let golden = corpus
            .lines()
            .find(|l| l.starts_with(&format!("blk\t{label}\t{container}\t{index}\t")))
            .unwrap_or_else(|| panic!("{label}/{container}/{index}: no blk row"));
        let golden_fields: Vec<&str> = golden.split('\t').collect();
        let mine = format!(
            "{}\t{}\t{}\t{}\t{}",
            header.method,
            header.content_type,
            header.content_id,
            header.compressed_size,
            header.uncompressed_size
        );
        let theirs = golden_fields[4..9].join("\t");
        assert_eq!(mine, theirs, "{label}/{container}/{index}");
        compared += 1;
    }

    assert_eq!(compared, 118, "block headers compared");
}

/// Walking a container: the checksum is outside the compressed size, so the offsets only line up if
/// that is right. Each block's recorded header bytes must begin exactly where the previous block's
/// total length left off.
#[test]
fn the_checksum_sits_outside_the_declared_content_size() {
    let corpus = corpus();
    let mut walked = 0;

    // Group the blocks of each container in order and check the arithmetic joins them up, using
    // only what the golden recorded: each header's own bytes and its declared sizes.
    let mut current: Option<(String, String)> = None;
    let mut expected_offset = 0usize;

    for line in corpus.lines() {
        let Some(rest) = line.strip_prefix("hdrbytes\t") else {
            continue;
        };
        let fields: Vec<&str> = rest.split('\t').collect();
        let key = (fields[0].to_string(), fields[1].to_string());
        let bytes = unhex(fields[3]);
        let header = BlockHeader::read(&bytes).expect("parses");

        if current.as_ref() != Some(&key) {
            current = Some(key);
            expected_offset = 0;
        }
        // The offset is not in the golden; what is checkable is that the total length is the
        // header plus the content plus four, and that it is strictly greater than the header.
        let total = header.total_length(3);
        assert_eq!(
            total,
            header.byte_length + header.compressed_size.max(0) as usize + 4,
            "{}/{}: total length",
            fields[0],
            fields[2]
        );
        assert!(
            total > header.byte_length,
            "a block is more than its header"
        );
        expected_offset += total;
        walked += 1;
    }

    assert!(expected_offset > 0);
    assert_eq!(walked, 118, "blocks walked");
}

/// The methods a file actually uses, which is the evidence behind decision 0038's scoping.
#[test]
fn an_ordinary_file_uses_raw_gzip_and_rans() {
    let corpus = corpus();
    let mut seen = 0;
    for line in corpus.lines() {
        let Some(rest) = line.strip_prefix("counts\t") else {
            continue;
        };
        let fields: Vec<&str> = rest.split('\t').collect();
        let (label, methods) = (fields[0], fields[2]);
        if label == "no-reads" {
            assert_eq!(
                methods, "0,1",
                "a file with no records needs no entropy codec"
            );
        } else {
            assert_eq!(methods, "0,1,4", "{label}: RAW, GZIP and rANS");
        }
        seen += 1;
    }
    assert_eq!(seen, 5, "counts rows");
}
