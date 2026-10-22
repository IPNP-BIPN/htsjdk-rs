//! Conformance for the CRAM slice header, against `htsjdk.samtools.cram.structure.Slice`.
//!
//! Goldens from `tools/cram-conformance/CramSliceHeaderDump.java` in the pinned oracle.
//!
//! The rows that justify the suite:
//!
//! ```text
//! blockcount  four-unmapped  0  27  27
//! ```
//!
//! The declared block count equals the blocks that follow the header, so a reader that counts the
//! header among them stops one short. htsjdk 5.0.0 writes a slice header with no tags at all,
//! where 4.2.0 wrote six digests.

use std::io::Read;

use htsjdk_cram::slice_header::{
    SliceHeader, EMBEDDED_REFERENCE_ABSENT, MD5_BYTE_SIZE, NO_ALIGNMENT, UNMAPPED_UNPLACED_ID,
};

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/cram_slice_header.txt.gz");
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

/// The constants, measured rather than assumed.
#[test]
fn the_constants_are_the_reference_constants() {
    let corpus = corpus();
    let sizes = rows(&corpus, "sizes");
    assert_eq!(sizes.len(), 1);
    assert_eq!(sizes[0][0], MD5_BYTE_SIZE.to_string());
    assert_eq!(sizes[0][1], EMBEDDED_REFERENCE_ABSENT.to_string());
    assert_eq!(sizes[0][2], NO_ALIGNMENT.to_string());
    assert_eq!(sizes[0][3], NO_ALIGNMENT.to_string());
}

/// Every header the reference wrote parses into the fields it recorded, and writing it back gives
/// the same bytes.
#[test]
fn every_header_round_trips_to_the_bytes_the_reference_wrote() {
    let corpus = corpus();
    let fields: std::collections::HashMap<(String, String), Vec<&str>> = rows(&corpus, "slice")
        .into_iter()
        .map(|row| ((row[0].to_string(), row[1].to_string()), row))
        .collect();

    let mut compared = 0;
    for row in rows(&corpus, "hdrbytes") {
        let (label, index) = (row[0], row[1]);
        let content = unhex(row[2]);
        let header = SliceHeader::read(&content).expect("the reference's own header parses");
        assert_eq!(
            hex(&header.write(3)),
            row[2],
            "{label}/{index}: written back"
        );

        let golden = &fields[&(label.to_string(), index.to_string())];
        let mine = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            header.reference_context_id,
            header.alignment_start,
            header.alignment_span,
            header.record_count,
            header.global_record_counter,
            header.block_count,
            if header.content_ids.is_empty() {
                "-".to_string()
            } else {
                header
                    .content_ids
                    .iter()
                    .map(|id| id.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            },
            header.embedded_reference_content_id,
            hex(&header.reference_md5),
            hex(&header.tags)
        );
        assert_eq!(mine, golden[2..12].join("\t"), "{label}/{index}: fields");
        compared += 1;
    }
    assert_eq!(compared, 5, "headers compared");
}

/// The declared block count is the blocks that follow the header, not including it.
#[test]
fn the_block_count_is_what_follows_the_header() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "blockcount") {
        let (label, index, declared, following) = (row[0], row[1], row[2], row[3]);
        assert_eq!(
            declared, following,
            "{label}/{index}: the count is exactly what follows"
        );
        compared += 1;
    }
    assert_eq!(compared, 5, "block counts compared");
}

/// The alignment context of an unmapped slice, and the two fields it forces to zero.
#[test]
fn an_unmapped_slice_has_no_start_and_no_span() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "hdrbytes") {
        let header = SliceHeader::read(&unhex(row[2])).expect("parses");
        assert!(header.is_unmapped_unplaced(), "{}: unmapped", row[0]);
        assert_eq!(header.reference_context_id, UNMAPPED_UNPLACED_ID);
        assert_eq!(header.alignment_start, NO_ALIGNMENT);
        assert_eq!(header.alignment_span, NO_ALIGNMENT);
        assert_eq!(
            header.embedded_reference_content_id,
            EMBEDDED_REFERENCE_ABSENT
        );
        assert_eq!(header.reference_md5, [0u8; MD5_BYTE_SIZE]);
        compared += 1;
    }
    assert_eq!(compared, 5, "alignment contexts compared");
}

/// htsjdk 5.0.0 writes no slice tags. 4.2.0 wrote six byte arrays per slice, `BD` and `SD` over
/// the bases and the scores and four digests of the empty string, 168 bytes of constant among them.
#[test]
fn a_slice_carries_no_tags() {
    let corpus = corpus();
    assert!(rows(&corpus, "slicetag").is_empty(), "no slice tag rows");
    let slices = rows(&corpus, "slice");
    assert_eq!(slices.len(), 5);
    for row in slices {
        assert_eq!(row.last(), Some(&"-"), "{}: tags", row[0]);
    }
}

/// One slice per file in this corpus, which is what makes the block counts comparable.
#[test]
fn every_file_here_holds_one_slice() {
    let corpus = corpus();
    let counts = rows(&corpus, "counts");
    assert_eq!(counts.len(), 5);
    for row in counts {
        assert_eq!(row[1], "1", "{}: slices", row[0]);
    }
}
