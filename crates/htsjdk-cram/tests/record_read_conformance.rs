//! Conformance for reading a CRAM record, against
//! `htsjdk.samtools.cram.encoding.reader.CramRecordReader`.
//!
//! Goldens from `tools/cram-conformance/CramRecordReadDump.java` in the pinned oracle, which
//! writes records through the reference's own writer, hands back the compression header and every
//! block of the slice, and records what its own reader made of them.
//!
//! This is the first suite where the port reads bytes the reference wrote from records, rather
//! than bytes a dump built by hand. The blocks are carried uncompressed: what compressor each one
//! got is another suite's business, and carrying the compressed stream here would make this suite
//! depend on it.
//!
//! The rows that justify the suite:
//!
//! ```text
//! start   delta-negative  100  SINGLE_REFERENCE: 0
//! record  delta-negative  0  0  2  0  r0  10  0  200  -1  0  0  40  -1  -  -
//! record  delta-negative  1  0  2  0  r1  10  0  100  -1  0  0  40  -1  -  -
//! record  features  0  0  2  0  r0  10  0  100  -1  0  0  40  -1  i@4;I@5;D@7;S@8  -
//! ```
//!
//! The alignment start is a delta and the delta may be negative. And a record's read features come
//! last, after the mate block and the tag list, in a series of their own.

use std::collections::BTreeMap;
use std::io::Read;

use htsjdk_cram::compression_header::CompressionHeader;
use htsjdk_cram::external_codecs::{SliceBlockBytes, SliceReadStreams};
use htsjdk_cram::read_features::ReadFeature;
use htsjdk_cram::record_read::{read_record, RecordReaders, SliceContext};

fn corpus() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/cram_record_read.txt.gz");
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

fn unhex(text: &str) -> Vec<u8> {
    if text == "-" {
        return Vec::new();
    }
    (0..text.len())
        .step_by(2)
        .map(|at| u8::from_str_radix(&text[at..at + 2], 16).expect("hex"))
        .collect()
}

/// The reference context the dump printed, as the id the reader compares against.
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

fn features(features: &[ReadFeature]) -> String {
    if features.is_empty() {
        return "-".to_string();
    }
    features
        .iter()
        .map(|feature| {
            let (operator, position) = match feature {
                ReadFeature::ReadBase { position, .. } => ('B', position),
                ReadFeature::Substitution { position, .. } => ('X', position),
                ReadFeature::Insertion { position, .. } => ('I', position),
                ReadFeature::SoftClip { position, .. } => ('S', position),
                ReadFeature::HardClip { position, .. } => ('H', position),
                ReadFeature::Padding { position, .. } => ('P', position),
                ReadFeature::Deletion { position, .. } => ('D', position),
                ReadFeature::RefSkip { position, .. } => ('N', position),
                ReadFeature::InsertBase { position, .. } => ('i', position),
                ReadFeature::BaseQualityScore { position, .. } => ('Q', position),
                ReadFeature::Bases { position, .. } => ('b', position),
                ReadFeature::Scores { position, .. } => ('q', position),
            };
            format!("{operator}@{position}")
        })
        .collect::<Vec<_>>()
        .join(";")
}

/// Every record the reference wrote and read back, read again from the same bytes.
#[test]
fn every_record_reads_back_the_way_the_reference_read_it() {
    let corpus = corpus();
    let mut labels: Vec<&str> = Vec::new();
    for row in rows(&corpus, "record") {
        if !labels.contains(&row[0]) {
            labels.push(row[0]);
        }
    }
    assert_eq!(labels.len(), 6, "slices measured");

    let mut compared = 0;
    for label in &labels {
        let header_hex = rows(&corpus, "header")
            .into_iter()
            .find(|row| row[0] == *label)
            .expect("a compression header")[1]
            .to_string();
        let header = CompressionHeader::read_block(&unhex(&header_hex), 3)
            .unwrap_or_else(|error| panic!("{label}: {}", error.message()));

        // The blocks, each already uncompressed by the dump.
        let mut core = Vec::new();
        let mut external: BTreeMap<i32, Vec<u8>> = BTreeMap::new();
        for row in rows(&corpus, "block")
            .into_iter()
            .filter(|row| row[0] == *label)
        {
            let content = unhex(row[3]);
            if row[1] == "core" {
                core = content;
            } else {
                external.insert(row[1].parse().expect("content id"), content);
            }
        }
        let blocks = SliceBlockBytes { core, external };

        let start_row = rows(&corpus, "start")
            .into_iter()
            .find(|row| row[0] == *label)
            .expect("a slice context");
        let slice = SliceContext {
            reference_context: reference_context(start_row[2]),
            alignment_start: start_row[1].parse().expect("alignment start"),
        };

        let readers = RecordReaders::new(&header.encodings)
            .unwrap_or_else(|error| panic!("{label}: {}", error.message()));
        let mut streams = SliceReadStreams::new(&blocks);
        let mut previous_start = slice.alignment_start;

        for row in rows(&corpus, "record")
            .into_iter()
            .filter(|row| row[0] == *label)
        {
            let index = row[1];
            let record = read_record(&readers, &header, &slice, &mut streams, previous_start)
                .unwrap_or_else(|error| panic!("{label}[{index}]: {}", error.message()));
            previous_start = record.alignment_start;

            assert_eq!(
                record.flags.bam.to_string(),
                row[2],
                "{label}[{index}] bam flags"
            );
            assert_eq!(
                record.flags.cram_flags().to_string(),
                row[3],
                "{label}[{index}] cram flags"
            );
            assert_eq!(
                record.flags.mate_flags().to_string(),
                row[4],
                "{label}[{index}] mate flags"
            );
            assert_eq!(
                record.read_name.clone().unwrap_or_default(),
                row[5],
                "{label}[{index}] read name"
            );
            assert_eq!(
                record.read_length.to_string(),
                row[6],
                "{label}[{index}] length"
            );
            assert_eq!(
                record.reference_index.to_string(),
                row[7],
                "{label}[{index}] reference"
            );
            assert_eq!(
                record.alignment_start.to_string(),
                row[8],
                "{label}[{index}] alignment start"
            );
            assert_eq!(
                record.mate_reference_index.to_string(),
                row[9],
                "{label}[{index}] mate reference"
            );
            assert_eq!(
                record.mate_alignment_start.to_string(),
                row[10],
                "{label}[{index}] mate start"
            );
            assert_eq!(
                record.template_size.to_string(),
                row[11],
                "{label}[{index}] template size"
            );
            assert_eq!(
                record.mapping_quality.to_string(),
                row[12],
                "{label}[{index}] mapping quality"
            );
            assert_eq!(
                record.read_group_id.to_string(),
                row[13],
                "{label}[{index}] read group"
            );
            assert_eq!(
                features(&record.read_features),
                row[14],
                "{label}[{index}] features"
            );
            compared += 1;
        }
    }
    assert_eq!(compared, 10, "records read");
}

/// The alignment start is a delta, and the delta may be negative. Both directions are in the
/// corpus, and reading the second record of each pair proves which one was applied.
#[test]
fn the_alignment_start_is_a_delta_in_both_directions() {
    let corpus = corpus();
    let forwards: Vec<i32> = rows(&corpus, "record")
        .into_iter()
        .filter(|row| row[0] == "delta")
        .map(|row| row[8].parse().expect("start"))
        .collect();
    assert_eq!(forwards, [100, 140]);

    let backwards: Vec<i32> = rows(&corpus, "record")
        .into_iter()
        .filter(|row| row[0] == "delta-negative")
        .map(|row| row[8].parse().expect("start"))
        .collect();
    assert_eq!(backwards, [200, 100], "a negative delta is legal");
}

/// An unmapped record reads no read features at all, and a multi-reference slice reads a reference
/// index per record where a single-reference slice takes the slice's own.
#[test]
fn what_a_record_reads_depends_on_the_slice_and_on_its_own_flags() {
    let corpus = corpus();

    let unmapped: Vec<&str> = rows(&corpus, "record")
        .into_iter()
        .filter(|row| row[0] == "unmapped")
        .map(|row| row[14])
        .collect();
    assert_eq!(unmapped, ["-"], "an unmapped record has no read features");

    // The mixed slice is multi-reference, and its records carry their own reference index: one
    // mapped to 0, one unmapped at -1, one mapped to 0 again.
    let mixed: Vec<&str> = rows(&corpus, "record")
        .into_iter()
        .filter(|row| row[0] == "mixed")
        .map(|row| row[7])
        .collect();
    assert_eq!(mixed, ["0", "-1", "0"]);

    let context = rows(&corpus, "start")
        .into_iter()
        .find(|row| row[0] == "mixed")
        .expect("a context")[2]
        .to_string();
    assert_eq!(context, "MULTIPLE_REFERENCE");
}
