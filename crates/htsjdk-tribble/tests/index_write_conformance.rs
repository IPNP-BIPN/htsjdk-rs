//! Conformance for *building* a `.idx`, against `LinearIndexCreator` and `DynamicIndexCreator`.
//!
//! Goldens from `tools/tribble-conformance/TribbleIndexWriteDump.java` in the pinned oracle.
//!
//! The rows that justify the suite are the ones where the bin width in the file is not the one the
//! creator was given, and differs between contigs of the same file:
//!
//! ```text
//! chr  sparse-linear        chr1  16000    2  600  4  0,45,64
//! chr  sparse-linear        chr2   8000    1   10  1  64,77
//! chr  sparse-dynamic-seek  chr2   2000    1   10  1  64,77
//! chr  mixed-linear         chr2 512000    2   10  2  67293,67311,67335
//! ```
//!
//! And the ones where the same data produces two different index *types*, each written in full:
//!
//! ```text
//! idx  sparse-dynamic-seek  LinearIndex
//! idx  sparse-dynamic-size  IntervalTreeIndex
//! idx  dense-dynamic-seek   IntervalTreeIndex
//! idx  dense-dynamic-size   LinearIndex
//! ```
//!
//! # What is compared, and what is taken as given
//!
//! The header carries the indexed file's URI, size, modification time and MD5, which are facts
//! about the run rather than about the creator. They are read out of the golden and handed back to
//! the writer, and everything else, the properties and every contig record, is produced here. So a
//! byte comparison of the whole file is still a comparison of this port's decisions: the only
//! bytes it is handed are the ones it could not know.

use std::io::Read;

use htsjdk_tribble::index::TribbleIndex;
use htsjdk_tribble::index::{INTERVAL_TREE, LINEAR};
use htsjdk_tribble::index_write::{
    check_ordering, BalanceApproach, BuiltIndex, DynamicIndexCreator, Feature,
    IntervalIndexCreator, LinearIndexCreator, DEFAULT_BIN_WIDTH, DEFAULT_FEATURE_COUNT,
};

/// How each case was built, mirroring the dump's third argument.
#[derive(Clone, Copy)]
enum How {
    Linear,
    IntervalTree,
    Dynamic(BalanceApproach),
}

/// A BED fixture: the lines, exactly as the dump wrote them.
fn fixture(name: &str) -> Vec<String> {
    match name {
        "sparse" => vec![
            "chr1\t100\t110\ta".into(),
            "chr1\t200\t210\tb".into(),
            "chr1\t300\t900\tc".into(),
            "chr1\t20000\t20010\td".into(),
            "chr2\t50\t60\te".into(),
        ],
        "dense" => (0..4000)
            .map(|i| format!("chr1\t{}\t{}\tf{i}", 100 + i * 2, 110 + i * 2))
            .collect(),
        "mixed" => {
            let mut lines: Vec<String> = (0..3000)
                .map(|i| format!("chr1\t{}\t{}\tc1_{i}", 100 + i * 3, 110 + i * 3))
                .collect();
            lines.push("chr2\t100\t110\tc2_a".into());
            lines.push("chr2\t900000\t900010\tc2_b".into());
            lines
        }
        "single" => vec!["chr1\t100\t110\tonly".into()],
        "far" => vec![
            "chr1\t100\t110\tnear".into(),
            "chr1\t5000000\t5000010\tfar".into(),
        ],
        "empty" => Vec::new(),
        "unsorted" => vec![
            "chr1\t500\t510\tsecond".into(),
            "chr1\t100\t110\tfirst".into(),
        ],
        "revisited" => vec![
            "chr1\t100\t110\ta".into(),
            "chr2\t100\t110\tb".into(),
            "chr1\t200\t210\tc".into(),
        ],
        other => panic!("no fixture {other}"),
    }
}

/// The features and their file positions.
///
/// `BEDCodec`'s default `StartOffset.ONE` shifts the start by one on the way in, so the feature the
/// creator sees is not the interval the file states. The position is the byte offset of the line,
/// which is what the creator records as a block start.
fn features(name: &str) -> (Vec<(Feature, i64)>, i64) {
    let mut out = Vec::new();
    let mut offset = 0i64;
    for line in fixture(name) {
        let columns: Vec<&str> = line.split('\t').collect();
        out.push((
            Feature {
                contig: columns[0].to_string(),
                start: columns[1].parse::<i32>().expect("a start") + 1,
                end: columns[2].parse().expect("an end"),
            },
            offset,
        ));
        offset += line.len() as i64 + 1;
    }
    (out, offset)
}

/// Label, fixture, how it was built. The dump's order.
const CASES: &[(&str, &str, How)] = &[
    ("sparse-linear", "sparse", How::Linear),
    ("sparse-interval", "sparse", How::IntervalTree),
    (
        "sparse-dynamic-seek",
        "sparse",
        How::Dynamic(BalanceApproach::ForSeekTime),
    ),
    (
        "sparse-dynamic-size",
        "sparse",
        How::Dynamic(BalanceApproach::ForSize),
    ),
    ("dense-linear", "dense", How::Linear),
    (
        "dense-dynamic-seek",
        "dense",
        How::Dynamic(BalanceApproach::ForSeekTime),
    ),
    (
        "dense-dynamic-size",
        "dense",
        How::Dynamic(BalanceApproach::ForSize),
    ),
    ("mixed-linear", "mixed", How::Linear),
    ("single-linear", "single", How::Linear),
    ("far-linear", "far", How::Linear),
    ("empty-linear", "empty", How::Linear),
    ("unsorted-linear", "unsorted", How::Linear),
    ("revisited-linear", "revisited", How::Linear),
];

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/tribble_index_write.txt.gz");
    let file = std::fs::File::open(&path).expect("corpus");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("corpus is gzip");
    text
}

fn row(text: &str, prefix: &str) -> Option<String> {
    text.lines()
        .find_map(|line| line.strip_prefix(prefix))
        .map(str::to_string)
}

/// What a built case yields: the index itself and the header properties.
type Built = (BuiltIndex, Vec<(String, String)>);

/// Build a case, returning either that or the refusal a row would carry.
fn build(name: &str, how: How) -> Result<Built, String> {
    let (features, final_position) = features(name);
    let plain: Vec<Feature> = features.iter().map(|(f, _)| f.clone()).collect();
    // The path in the message is the one the dump ran with, so the label is the source here and
    // the message is compared by its stable half.
    check_ordering(&plain, "SOURCE").map_err(|error| error.message())?;

    match how {
        How::IntervalTree => {
            let mut creator = IntervalIndexCreator::new(DEFAULT_FEATURE_COUNT);
            for (feature, position) in &features {
                creator.add_feature(feature, *position);
            }
            Ok((
                BuiltIndex::IntervalTree(creator.finalize(final_position)),
                Vec::new(),
            ))
        }
        How::Linear => {
            let mut creator = LinearIndexCreator::new(DEFAULT_BIN_WIDTH);
            for (feature, position) in &features {
                creator.add_feature(feature, *position);
            }
            creator
                .finalize(final_position, Vec::new())
                .map(|contigs| (BuiltIndex::Linear(contigs), Vec::new()))
                .map_err(|error| error.message())
        }
        How::Dynamic(approach) => {
            let mut creator = DynamicIndexCreator::new(approach);
            for (feature, position) in &features {
                creator.add_feature(feature, *position);
            }
            let properties = creator.properties();
            creator
                .finalize(final_position)
                .map(|built| (built, properties))
                .map_err(|error| error.message())
        }
    }
}

#[test]
fn every_index_is_built_as_the_reference_builds_it() {
    let text = corpus();
    let mut compared = 0;
    let mut trees = 0;

    for (label, name, how) in CASES {
        let golden_idx = row(&text, &format!("idx\t{label}\t"));
        let golden_err = row(&text, &format!("err\t{label}\t"));

        match build(name, *how) {
            Err(message) => {
                let golden = golden_err
                    .unwrap_or_else(|| panic!("{label}: the port refused and the golden did not"));
                let (_class, golden_message) =
                    golden.split_once('\t').expect("a class and a message");
                // The upstream message ends with the dump's own working directory, which is not a
                // property of the creator; the half before it is.
                let stable = golden_message
                    .split(", for input source:")
                    .next()
                    .expect("a message");
                // The dump escapes newlines, and one of these two messages carries one.
                let escaped = message.replace('\\', "\\\\").replace('\n', "\\n");
                assert!(
                    escaped.starts_with(stable),
                    "{label}:\n  port:   {escaped}\n  golden: {stable}"
                );
                compared += 1;
            }
            Ok((built, properties)) => {
                let idx = golden_idx.unwrap_or_else(|| {
                    panic!("{label}: the port built a file and the golden did not")
                });
                let fields: Vec<&str> = idx.split('\t').collect();
                let expected_class = match built {
                    BuiltIndex::Linear(_) => "LinearIndex",
                    BuiltIndex::IntervalTree(_) => {
                        trees += 1;
                        "IntervalTreeIndex"
                    }
                };
                assert_eq!(fields[0], expected_class, "{label}: index type");
                let timestamp_offset: usize = fields[1].parse().expect("an offset");
                let golden_bytes = base64_decode(fields[2]);

                // The header facts about the run, taken from the golden and handed back, so what
                // is compared is what this port decides.
                let reference = TribbleIndex::read(&golden_bytes).expect("the golden parses");
                let mine = match built {
                    BuiltIndex::Linear(contigs) => TribbleIndex {
                        index_type: LINEAR,
                        contigs,
                        interval_contigs: Vec::new(),
                        properties,
                        ..reference.clone()
                    },
                    BuiltIndex::IntervalTree(interval_contigs) => TribbleIndex {
                        index_type: INTERVAL_TREE,
                        contigs: Vec::new(),
                        interval_contigs,
                        properties,
                        ..reference.clone()
                    },
                };
                assert_eq!(
                    mine.index_type, reference.index_type,
                    "{label}: the port chose a different layout from the reference"
                );
                let mut bytes = mine.write().expect("an index writes");
                assert!(
                    timestamp_offset + 8 <= bytes.len(),
                    "{label}: the timestamp offset is past the end"
                );
                bytes[timestamp_offset..timestamp_offset + 8].fill(0);

                assert_eq!(
                    bytes.len(),
                    golden_bytes.len(),
                    "{label}: the port wrote {} bytes and the reference wrote {}",
                    bytes.len(),
                    golden_bytes.len()
                );
                let first = bytes
                    .iter()
                    .zip(golden_bytes.iter())
                    .position(|(a, b)| a != b);
                assert_eq!(first, None, "{label}: first differing byte");
                compared += 1;
            }
        }
    }

    // Named rather than implied: a run that compared nothing would otherwise pass.
    assert_eq!(compared, 13, "cases compared byte for byte or by refusal");
    assert_eq!(trees, 3, "cases the reference built as an interval tree");
}

/// The per-contig rows, checked on their own so a divergence names the contig rather than a byte
/// offset. Redundant with the byte comparison and worth the redundancy: a wrong bin width and a
/// wrong block list are one failure in bytes and two here.
#[test]
fn every_contig_carries_the_width_the_optimizer_chose() {
    let text = corpus();
    let mut compared = 0;

    for (label, name, how) in CASES {
        let Ok((BuiltIndex::Linear(contigs), _)) = build(name, *how) else {
            continue;
        };
        for contig in &contigs {
            let golden = row(&text, &format!("chr\t{label}\t{}\t", contig.name))
                .unwrap_or_else(|| panic!("{label}/{}: no golden row", contig.name));
            let positions: Vec<String> = contig
                .blocks
                .iter()
                .map(|block| block.start.to_string())
                .chain(
                    contig
                        .blocks
                        .last()
                        .map(|last| (last.start + last.size).to_string()),
                )
                .collect();
            let mine = format!(
                "{}\t{}\t{}\t{}\t{}",
                contig.bin_width,
                contig.blocks.len(),
                contig.longest_feature,
                contig.n_features,
                positions.join(",")
            );
            assert_eq!(mine, golden, "{label}/{}", contig.name);
            compared += 1;
        }
    }

    assert_eq!(compared, 10, "contig rows compared");
}

/// `Base64.getEncoder()`, decoded without a dependency: the alphabet is fixed and the input is a
/// golden this repository produced.
fn base64_decode(text: &str) -> Vec<u8> {
    const ALPHABET: &[u8] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = Vec::with_capacity(text.len() * 3 / 4);
    let mut accumulator: u32 = 0;
    let mut bits = 0;
    for byte in text.bytes() {
        if byte == b'=' {
            break;
        }
        let value = ALPHABET
            .iter()
            .position(|c| *c == byte)
            .unwrap_or_else(|| panic!("not base64: {byte:?}")) as u32;
        accumulator = (accumulator << 6) | value;
        bits += 6;
        if bits >= 8 {
            bits -= 8;
            out.push((accumulator >> bits) as u8);
        }
    }
    out
}
