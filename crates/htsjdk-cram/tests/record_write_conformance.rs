//! Conformance for writing a CRAM record, against
//! `htsjdk.samtools.cram.encoding.writer.CramRecordWriter`.
//!
//! Goldens from `tools/cram-conformance/CramRecordWriteDump.java` in the pinned oracle, which
//! writes records through the reference's own writer and hands back both the records, field by
//! field and feature by feature, and every block they landed in.
//!
//! This is the other half of `cram-record-read`: the port rebuilds each record from its row,
//! writes it, and the blocks have to come out byte for byte what the reference produced. A
//! difference in the order, in a delta, or in a single branch shows here as bytes that do not
//! match rather than as an error.
//!
//! The blocks are carried uncompressed, as they are on the read side: what compressor each one got
//! is `cram-block`'s business.

use std::collections::BTreeMap;
use std::io::Read;

use htsjdk_cram::compression_header::CompressionHeader;
use htsjdk_cram::external_codecs::SliceWriteStreams;
use htsjdk_cram::read_features::ReadFeature;
use htsjdk_cram::record_flags::Flags;
use htsjdk_cram::record_read::{CramRecord, SliceContext};
use htsjdk_cram::record_write::{write_record, RecordWriters};

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/cram_record_write.txt.gz");
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

fn reference_context(text: &str) -> i32 {
    match text {
        "UNMAPPED_UNPLACED" => -1,
        "MULTIPLE_REFERENCE" => SliceContext::MULTIPLE_REFERENCE,
        other => other
            .strip_prefix("SINGLE_REFERENCE: ")
            .expect("a reference context")
            .parse()
            .expect("a reference index"),
    }
}

/// A record rebuilt from its row: every field the writer reads is on it.
fn record(row: &[&str]) -> CramRecord {
    CramRecord {
        flags: Flags {
            bam: row[2].parse().expect("bam flags"),
            cram: row[3].parse().expect("cram flags"),
            mate: row[4].parse().expect("mate flags"),
        },
        read_name: Some(row[5].to_string()),
        read_length: row[6].parse().expect("read length"),
        reference_index: row[7].parse().expect("reference index"),
        alignment_start: row[8].parse().expect("alignment start"),
        mate_reference_index: row[9].parse().expect("mate reference"),
        mate_alignment_start: row[10].parse().expect("mate start"),
        template_size: row[11].parse().expect("template size"),
        mapping_quality: row[12].parse().expect("mapping quality"),
        read_group_id: row[13].parse().expect("read group"),
        read_features: features(row[14]),
        read_bases: unhex(row[15]),
        quality_scores: unhex(row[16]),
        records_to_next_fragment: -1,
        tags: Vec::new(),
    }
}

/// The features as the dump prints them: operator, position and payload.
fn features(column: &str) -> Vec<ReadFeature> {
    if column == "-" {
        return Vec::new();
    }
    column
        .split(';')
        .map(|feature| {
            let (head, payload) = feature.rsplit_once(':').expect("a payload");
            let (operator, position) = head.split_once('@').expect("a position");
            let position: i32 = position.parse().expect("a position");
            let operator = operator.as_bytes()[0];
            match operator {
                b'i' => ReadFeature::InsertBase {
                    position,
                    base: unhex(payload)[0],
                },
                b'I' => ReadFeature::Insertion {
                    position,
                    sequence: unhex(payload),
                },
                b'S' => ReadFeature::SoftClip {
                    position,
                    sequence: unhex(payload),
                },
                b'D' => ReadFeature::Deletion {
                    position,
                    length: payload.parse().expect("a length"),
                },
                other => panic!("{}", other as char),
            }
        })
        .collect()
}

/// Every slice the reference wrote, written again from the records its rows carry.
#[test]
fn every_slice_writes_the_reference_blocks() {
    let corpus = corpus();
    let mut labels: Vec<&str> = Vec::new();
    for row in rows(&corpus, "record") {
        if !labels.contains(&row[0]) {
            labels.push(row[0]);
        }
    }
    assert_eq!(labels.len(), 8, "slices measured");

    let mut compared = 0;
    for label in &labels {
        let header_hex = rows(&corpus, "header")
            .into_iter()
            .find(|row| row[0] == *label)
            .expect("a compression header")[1]
            .to_string();
        let header = CompressionHeader::read_block(&unhex(&header_hex), 3)
            .unwrap_or_else(|error| panic!("{label}: {}", error.message()));

        let start_row = rows(&corpus, "start")
            .into_iter()
            .find(|row| row[0] == *label)
            .expect("a slice context");
        let slice = SliceContext {
            reference_context: reference_context(start_row[2]),
            alignment_start: start_row[1].parse().expect("alignment start"),
        };

        let writers = RecordWriters::new(&header.encodings)
            .unwrap_or_else(|error| panic!("{label}: {}", error.message()));
        let mut streams = SliceWriteStreams::new();
        let mut previous_start = slice.alignment_start;

        for row in rows(&corpus, "record")
            .into_iter()
            .filter(|row| row[0] == *label)
        {
            let record = record(&row);
            let tag_ids_index: i32 = row[17].parse().expect("a tag list index");
            write_record(
                &writers,
                &header,
                &slice,
                &mut streams,
                &record,
                tag_ids_index,
                previous_start,
            )
            .unwrap_or_else(|error| panic!("{label}: {}", error.message()));
            previous_start = record.alignment_start;
        }

        // The blocks the reference produced for this slice, each already uncompressed.
        let mut expected: BTreeMap<String, String> = BTreeMap::new();
        for row in rows(&corpus, "block")
            .into_iter()
            .filter(|row| row[0] == *label)
        {
            expected.insert(row[1].to_string(), row[3].to_string());
        }

        let written = streams.finish();
        assert_eq!(hex(&written.core), expected["core"], "{label} core block");
        for (content_id, content) in &written.external {
            // A block the port wrote to must be the block the reference wrote, byte for byte.
            let expected = expected
                .get(&content_id.to_string())
                .unwrap_or_else(|| panic!("{label}: the reference has no block {content_id}"));
            assert_eq!(&hex(content), expected, "{label} block {content_id}");
            compared += 1;
        }

        // And every block the reference wrote to must be one the port wrote.
        for (content_id, content) in &expected {
            if content == "-" || content_id == "core" {
                continue;
            }
            let id: i32 = content_id.parse().expect("a content id");
            assert!(
                written.external.contains_key(&id),
                "{label}: the port left block {content_id} empty and the reference did not"
            );
        }
    }
    assert!(compared >= 8, "blocks compared");
}

/// The order is prescribed and shared with the reader, so a record written and read again is the
/// record that went in. That is the property the two suites hold together.
#[test]
fn what_the_writer_produces_the_reader_consumes() {
    use htsjdk_cram::external_codecs::SliceReadStreams;
    use htsjdk_cram::record_read::{read_record, RecordReaders};

    let corpus = corpus();
    let header_hex = rows(&corpus, "header")
        .into_iter()
        .find(|row| row[0] == "features")
        .expect("the features slice")[1]
        .to_string();
    let header = CompressionHeader::read_block(&unhex(&header_hex), 3).expect("a header");
    let slice = SliceContext {
        reference_context: 0,
        alignment_start: 100,
    };

    let row: Vec<&str> = rows(&corpus, "record")
        .into_iter()
        .find(|row| row[0] == "features")
        .expect("its record");
    let original = record(&row);

    let writers = RecordWriters::new(&header.encodings).expect("writers");
    let mut streams = SliceWriteStreams::new();
    write_record(
        &writers,
        &header,
        &slice,
        &mut streams,
        &original,
        0,
        slice.alignment_start,
    )
    .expect("written");
    let blocks = streams.finish();

    let readers = RecordReaders::new(&header.encodings).expect("readers");
    let mut input = SliceReadStreams::new(&blocks);
    let back = read_record(&readers, &header, &slice, &mut input, slice.alignment_start)
        .expect("read back");

    assert_eq!(back.flags, original.flags);
    assert_eq!(back.read_name, original.read_name);
    assert_eq!(back.alignment_start, original.alignment_start);
    assert_eq!(back.read_features, original.read_features);
    assert_eq!(back.mapping_quality, original.mapping_quality);
}
