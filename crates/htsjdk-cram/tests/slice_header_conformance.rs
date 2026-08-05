//! Conformance for the CRAM slice header, against `htsjdk.samtools.cram.structure.Slice`.
//!
//! Goldens from `tools/cram-conformance/CramSliceHeaderDump.java` in the pinned oracle.
//!
//! The rows that justify the suite:
//!
//! ```text
//! blockcount  four-unmapped  0  27  27
//! slicetag    four-unmapped  0  B1  byte[]  da39a3ee5e6b4b0d3255bfef95601890afd80709
//! slicetag    four-unmapped  0  BD  byte[]  9efa8928
//! slicetag    tagged         0  BD  byte[]  9efa8928
//! ```
//!
//! The declared block count equals the blocks that follow the header, so a reader that counts the
//! header among them stops one short. Four of the six slice tags are digests of the empty string
//! and are identical in every file; only `BD` and `SD` move with the reads, and they do not move
//! when only the record tags change.

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

/// Four of the six slice tags are digests of the empty string, identical in every file. Only `BD`
/// and `SD` move with the reads, and they do not move when only the record tags change.
#[test]
fn four_of_the_six_slice_tags_digest_nothing() {
    let corpus = corpus();
    let tags = rows(&corpus, "slicetag");
    assert_eq!(tags.len(), 30, "six tags over five files");

    let mut constant: std::collections::HashMap<&str, &str> = std::collections::HashMap::new();
    let mut varying: std::collections::HashMap<(&str, &str), &str> =
        std::collections::HashMap::new();

    for row in &tags {
        let (label, name, value) = (row[0], row[2], row[4]);
        assert_eq!(row[3], "byte[]", "{name} is a byte array");
        match name {
            "B1" | "S1" | "B5" | "S5" => {
                if let Some(seen) = constant.get(name) {
                    assert_eq!(*seen, value, "{name} is the same in every file");
                } else {
                    constant.insert(name, value);
                }
            }
            "BD" | "SD" => {
                varying.insert((label, name), value);
            }
            other => panic!("{other}: an unexpected slice tag"),
        }
    }

    // The SHA-1 and the SHA-512 of the empty string.
    assert_eq!(constant["B1"], "da39a3ee5e6b4b0d3255bfef95601890afd80709");
    assert_eq!(constant["S1"], constant["B1"]);
    assert_eq!(constant["B1"].len(), 40, "twenty bytes");
    assert_eq!(constant["B5"].len(), 128, "sixty-four bytes");
    assert_eq!(constant["S5"], constant["B5"]);

    // Only the reads move BD and SD, and the record tags do not.
    assert_ne!(
        varying[&("four-unmapped", "BD")],
        varying[&("one-unmapped", "BD")],
        "different reads, different digest"
    );
    assert_eq!(
        varying[&("four-unmapped", "BD")],
        varying[&("tagged", "BD")],
        "the same reads with tags added, the same digest"
    );
    assert_eq!(
        varying[&("four-unmapped", "SD")],
        varying[&("tagged", "SD")]
    );

    // 168 bytes of constant per slice.
    let constant_bytes =
        (constant["B1"].len() + constant["S1"].len() + constant["B5"].len() + constant["S5"].len())
            / 2;
    assert_eq!(constant_bytes, 168);
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
