//! Conformance for the Tribble `.idx` against htsjdk 4.2.0.
//!
//! Golden from `tools/tribble-conformance/TribbleIndexDump.java`.
//!
//! The dump was written before the port, and this file is what that bought: every number below
//! came out of the reference rather than out of a reading of the format. Three of them would have
//! been wrong the other way round, and two would not have failed loudly — the type identifiers,
//! which the Java defines through a circular reference and cannot be read from source at all; the
//! bin width, which is **per contig** rather than the creator's default; and the header's
//! timestamp, which makes two indexes over the same input differ byte for byte.
//!
//! # The bytes in the golden have their timestamp masked
//!
//! Deliberately, and the offset travels with them. `indexedFileTS` is the source file's
//! modification time, measured to change on every run, so a golden carrying it would fail
//! intermittently — or be "repaired" by regeneration, which is how a suite quietly stops meaning
//! anything. Everything else in the layout is still compared byte for byte.
//!
//! # Only the linear index is parsed
//!
//! The interval-tree chromosome record is a different shape, and this reader **refuses** it rather
//! than guessing. The golden holds its rows too, so the suite records what the reference does for
//! both and says plainly which half the port covers.

use std::io::Read;

use htsjdk_tribble::index::{IndexError, TribbleIndex, INTERVAL_TREE, LINEAR, MAGIC_NUMBER};

fn golden() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/tribble_index.txt.gz");
    let file = std::fs::File::open(&path).expect("golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("golden is gzip");
    text
}

/// The golden carries the `.idx` base64-encoded, so the test decodes rather than pulling in a
/// dependency for four lines.
fn decode_base64(text: &str) -> Vec<u8> {
    const ALPHABET: &[u8] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = Vec::with_capacity(text.len() * 3 / 4);
    let mut accumulator = 0u32;
    let mut bits = 0u32;
    for byte in text.bytes() {
        if byte == b'=' {
            break;
        }
        let Some(value) = ALPHABET.iter().position(|c| *c == byte) else {
            continue;
        };
        accumulator = (accumulator << 6) | value as u32;
        bits += 6;
        if bits >= 8 {
            bits -= 8;
            out.push((accumulator >> bits) as u8);
        }
    }
    out
}

/// One `idx` row: the label, the masked bytes, and where the timestamp was.
fn indexes(text: &str) -> Vec<(String, usize, Vec<u8>)> {
    text.lines()
        .filter_map(|line| line.strip_prefix("idx\t"))
        .map(|rest| {
            let fields: Vec<&str> = rest.split('\t').collect();
            (
                fields[0].to_string(),
                fields[1].parse().expect("an offset"),
                decode_base64(fields[2]),
            )
        })
        .collect()
}

fn rows<'a>(text: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    text.lines()
        .filter_map(|line| line.strip_prefix(kind))
        .filter_map(|rest| rest.strip_prefix('\t'))
        .map(|rest| rest.split('\t').collect())
        .collect()
}

#[test]
fn every_linear_index_parses_to_the_references_numbers() {
    let text = golden();
    let chr_rows = rows(&text, "chr");
    let header_rows = rows(&text, "header");
    let mut parsed = 0usize;
    let mut refused = 0usize;

    for (label, timestamp_offset, bytes) in indexes(&text) {
        let header = header_rows
            .iter()
            .find(|row| row[0] == label)
            .expect("a header row per index");
        assert_eq!(header[1], MAGIC_NUMBER.to_string(), "{label}: magic");

        let index_type: i32 = header[2].parse().expect("a type");
        match TribbleIndex::read(&bytes) {
            Err(IndexError::UnsupportedType { found }) => {
                // The interval-tree index is refused rather than mis-parsed, and the golden says
                // the reference wrote one.
                assert_eq!(found, INTERVAL_TREE, "{label}");
                assert_eq!(index_type, INTERVAL_TREE, "{label}");
                refused += 1;
                continue;
            }
            Err(other) => panic!("{label}: {other:?}"),
            Ok(index) => {
                assert_eq!(index_type, LINEAR, "{label}");
                assert_eq!(index.version.to_string(), header[3], "{label}: version");
                assert_eq!(index.flags.to_string(), header[4], "{label}: flags");
                assert_eq!(
                    index.properties.len().to_string(),
                    header[5],
                    "{label}: properties"
                );
                // The masked field reads back as zero, which is the point of masking it.
                assert_eq!(
                    index.indexed_file_timestamp, 0,
                    "{label}: the golden's timestamp is masked at byte {timestamp_offset}"
                );

                let expected: Vec<&Vec<&str>> =
                    chr_rows.iter().filter(|row| row[0] == label).collect();
                assert_eq!(index.contigs.len(), expected.len(), "{label}: contig count");
                for (contig, row) in index.contigs.iter().zip(&expected) {
                    assert_eq!(contig.name, row[1], "{label}: contig name");
                    assert_eq!(contig.bin_width.to_string(), row[2], "{label}: bin width");
                    assert_eq!(
                        contig.blocks.len().to_string(),
                        row[3],
                        "{label}/{}: block count",
                        contig.name
                    );
                    assert_eq!(
                        contig.longest_feature.to_string(),
                        row[4],
                        "{label}/{}: longest feature",
                        contig.name
                    );
                    assert_eq!(contig.unused.to_string(), row[5], "{label}: unused slot");
                    assert_eq!(
                        contig.n_features.to_string(),
                        row[6],
                        "{label}/{}: feature count",
                        contig.name
                    );
                    // N blocks from N+1 positions: the reconstruction is the thing most easily
                    // got wrong, so the positions are compared rather than the sizes.
                    let mut positions: Vec<String> =
                        contig.blocks.iter().map(|b| b.start.to_string()).collect();
                    if let Some(last) = contig.blocks.last() {
                        positions.push((last.start + last.size).to_string());
                    }
                    assert_eq!(
                        positions.join(","),
                        row[7],
                        "{label}/{}: block positions",
                        contig.name
                    );
                }
                parsed += 1;
            }
        }
    }

    assert_eq!(parsed, 2, "the golden changed size");
    assert_eq!(
        refused, 1,
        "the interval-tree index must still be in the corpus"
    );
    println!("tribble index: {parsed} linear indexes parsed, {refused} refused as interval-tree");
}

#[test]
fn every_query_resolves_to_the_references_blocks() {
    let text = golden();
    let mut compared = 0usize;
    let mut refusals = 0usize;
    let mut empty = 0usize;

    let parsed: Vec<(String, TribbleIndex)> = indexes(&text)
        .into_iter()
        .filter_map(|(label, _, bytes)| TribbleIndex::read(&bytes).ok().map(|i| (label, i)))
        .collect();

    for row in rows(&text, "query") {
        let (label, interval, expected) = (row[0], row[1], row[2]);
        let Some((_, index)) = parsed.iter().find(|(l, _)| l == label) else {
            // An interval-tree index, which this reader does not parse; its query rows are in the
            // golden as a record of the reference rather than as a claim about the port.
            continue;
        };
        let (contig, range) = interval.split_once(':').expect("contig:start-end");
        let (start, end) = range.split_once('-').expect("start-end");
        let start: i32 = start.parse().expect("a start");
        let end: i32 = end.parse().expect("an end");

        match index.blocks(contig, start, end) {
            Err(error) => {
                assert!(
                    expected.starts_with("E:"),
                    "{label} {interval}: the reference answered {expected}"
                );
                assert_eq!(
                    format!("E:{}", error.class()),
                    expected,
                    "{label} {interval}: exception class"
                );
                refusals += 1;
            }
            Ok(blocks) if blocks.is_empty() => {
                assert_eq!(expected, "none", "{label} {interval}");
                empty += 1;
            }
            Ok(blocks) => {
                let rendered: Vec<String> = blocks
                    .iter()
                    .map(|b| format!("{},{}", b.start, b.size))
                    .collect();
                assert_eq!(rendered.join("|"), expected, "{label} {interval}");
                // Never a list: linear blocks are adjacent, so a query merges into one run.
                assert_eq!(blocks.len(), 1, "{label} {interval}: one block or none");
            }
        }
        compared += 1;
    }

    assert_eq!(compared, 12, "the golden changed size");
    assert!(
        refusals > 0 && empty > 0,
        "the corpus must keep both ways of getting nothing: a refused contig ({refusals}) and an \
         empty answer ({empty})"
    );
    println!(
        "tribble index: {compared} queries, {refusals} refused for an unknown contig, {empty} empty"
    );
}
