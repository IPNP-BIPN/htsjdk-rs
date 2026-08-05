//! Conformance for rANS 4x8 order 1, against
//! `htsjdk.samtools.cram.compression.rans.rans4x8.RANS4x8Encode.compressOrder1Way4`.
//!
//! Goldens from `tools/cram-conformance/CramRansOrder1Dump.java` in the pinned oracle.
//!
//! The rows that justify the suite are the ones where order 1 stops being order 0 with more
//! tables:
//!
//! ```text
//! quarter  four-bytes  1  1  2  3  67  71  84
//! norm     four-bytes  0  65  1  1024
//! norm     four-bytes  0  67  1  1024
//! zerofreq 0  65  4096
//! ```
//!
//! The four bytes `ACGT` produce a context 0 holding all four symbols at 1024 apiece, because the
//! table counts the byte at each quarter boundary as if it followed nothing. And a frequency byte
//! of zero, which no writer emits, reads as the whole 4096.

use std::io::Read;

use htsjdk_cram::rans::{self, Order, NUMBER_OF_SYMBOLS, PREFIX_LENGTH};
use htsjdk_cram::rans_order1::{
    calc_frequencies_order1, compress_order1, write_frequencies_order1, Table,
};
use sha2::{Digest, Sha256};

/// The dump prints the whole stream only when it is no longer than this.
const MAX_INLINE_BYTES: usize = 4096;

fn corpus() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/cram_rans_order1.txt.gz");
    let file = std::fs::File::open(&path).expect("corpus");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("corpus is gzip");
    text
}

fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

fn sha256(bytes: &[u8]) -> String {
    hex(&Sha256::digest(bytes))
}

fn repeat(source: &[u8], length: usize) -> Vec<u8> {
    (0..length).map(|i| source[i % source.len()]).collect()
}

fn spread(from: usize, to: usize, count: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity((to - from) * count);
    for _ in 0..count {
        for symbol in from..to {
            out.push(symbol as u8);
        }
    }
    out
}

/// The same inputs the dump built, by the same construction. The `in` rows' digests are what say
/// so.
fn input(label: &str) -> Vec<u8> {
    match label {
        "four-bytes" => b"ACGT".to_vec(),
        "five-bytes" => b"ACGTA".to_vec(),
        "six-bytes" => b"ACGTAC".to_vec(),
        "seven-bytes" => b"ACGTACG".to_vec(),
        "eight-bytes" => b"ACGTACGT".to_vec(),
        "acgt-1000" => repeat(b"ACGT", 1000),
        "acgt-1001" => repeat(b"ACGT", 1001),
        "acgt-1002" => repeat(b"ACGT", 1002),
        "acgt-1003" => repeat(b"ACGT", 1003),
        "uniform-1000" => repeat(b"A", 1000),
        "two-symbols" => {
            let mut out = vec![b'A'; 1000];
            out[900..].fill(b'B');
            out
        }
        "contiguous-run" => spread(5, 13, 100),
        "motif-1000" => repeat(b"ACGTACGTAA", 1000),
        "quality-band" => {
            let mut seed = 0x5DEECE66Du64;
            (0..5000)
                .map(|_| {
                    seed = seed
                        .wrapping_mul(6364136223846793005)
                        .wrapping_add(1442695040888963407);
                    30 + ((seed >> 33) % 8) as u8
                })
                .collect()
        }
        "noise-10000" => {
            let mut seed = 0x1234_5678u64;
            (0..10_000)
                .map(|_| {
                    seed = seed
                        .wrapping_mul(6364136223846793005)
                        .wrapping_add(1442695040888963407);
                    (seed >> 33) as u8
                })
                .collect()
        }
        "all-256-once" => (0..=255u8).collect(),
        other => panic!("{other}: no such input"),
    }
}

fn rows<'a>(corpus: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    let prefix = format!("{kind}\t");
    corpus
        .lines()
        .filter_map(|line| line.strip_prefix(&prefix))
        .map(|rest| rest.split('\t').collect())
        .collect()
}

/// The dump's `digestOf`: every non-zero entry as `context:symbol=frequency;`, then SHA-256.
fn table_digest(table: &Table) -> String {
    let mut text = String::new();
    for (context, row) in table.iter().enumerate() {
        for (symbol, frequency) in row.iter().enumerate() {
            if *frequency != 0 {
                text.push_str(&format!("{context}:{symbol}={frequency};"));
            }
        }
    }
    sha256(text.as_bytes())
}

/// Every input is rebuilt exactly, or nothing below means anything.
#[test]
fn the_inputs_are_the_inputs_the_oracle_compressed() {
    let corpus = corpus();
    let mut checked = 0;
    for row in rows(&corpus, "in") {
        let (label, length, digest) = (row[0], row[1].parse::<usize>().expect("length"), row[2]);
        let input = input(label);
        assert_eq!(input.len(), length, "{label}: length");
        assert_eq!(sha256(&input), digest, "{label}: digest");
        checked += 1;
    }
    assert_eq!(checked, 16, "inputs rebuilt");
}

/// The whole compressed stream, byte for byte where the golden carries it and by digest where it
/// is too long to.
#[test]
fn every_input_compresses_to_the_bytes_the_reference_produced() {
    let corpus = corpus();
    let inline: std::collections::HashMap<&str, &str> = rows(&corpus, "bytes")
        .into_iter()
        .map(|row| (row[0], row[1]))
        .collect();

    let mut compared = 0;
    for row in rows(&corpus, "enc") {
        let (label, length, digest) = (row[0], row[1].parse::<usize>().expect("length"), row[2]);
        let ours = compress_order1(&input(label));
        assert_eq!(ours.len(), length, "{label}: compressed length");
        assert_eq!(sha256(&ours), digest, "{label}: compressed digest");
        if let Some(golden) = inline.get(label) {
            assert_eq!(&hex(&ours), golden, "{label}: compressed bytes");
        }
        compared += 1;
    }
    assert_eq!(compared, 16, "streams compared");

    // `compress` reaches order 1 through the same length rule order 0 declared.
    for label in ["four-bytes", "acgt-1000"] {
        let bytes = input(label);
        assert_eq!(
            rans::compress(&bytes, Order::One).map_err(|e| e.message()),
            Ok(compress_order1(&bytes)),
            "{label}: dispatched"
        );
    }
}

/// The prefix carries order 1, not order 0, on every input here: they are all at least four bytes.
#[test]
fn the_prefix_records_order_one() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "prefix") {
        let label = row[0];
        let ours = compress_order1(&input(label));
        assert_eq!(row[1], "1", "{label}: the reference wrote order 1");
        assert_eq!(
            i32::from(ours[0]).to_string(),
            row[1],
            "{label}: order byte"
        );
        assert_eq!(
            i32::from_le_bytes(ours[1..5].try_into().unwrap()).to_string(),
            row[2],
            "{label}: compressed size"
        );
        assert_eq!(
            i32::from_le_bytes(ours[5..9].try_into().unwrap()).to_string(),
            row[3],
            "{label}: raw size"
        );
        compared += 1;
    }
    assert_eq!(compared, 16, "prefixes compared");
}

/// The three bytes that are counted as if they followed nothing, and where they sit.
#[test]
fn three_quarter_boundaries_are_counted_as_bigrams_that_do_not_exist() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "quarter") {
        let label = row[0];
        let bytes = input(label);
        let quarter = bytes.len() >> 2;
        let mine = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            quarter,
            quarter,
            2 * quarter,
            3 * quarter,
            bytes[quarter],
            bytes[2 * quarter],
            bytes[3 * quarter]
        );
        assert_eq!(mine, row[1..8].join("\t"), "{label}");
        compared += 1;
    }
    assert_eq!(compared, 16, "quarter rows compared");
}

/// The raw per-context totals, which include the three counts the input does not contain.
#[test]
fn every_context_total_includes_the_three_extra_counts() {
    let corpus = corpus();
    let mut current: Option<(String, [i32; NUMBER_OF_SYMBOLS])> = None;
    let mut compared = 0;

    for row in rows(&corpus, "ctxtotal") {
        let label = row[0];
        let context: usize = row[1].parse().expect("context");
        if current.as_ref().map(|(l, _)| l.as_str()) != Some(label) {
            let bytes = input(label);
            let mut totals = [0i32; NUMBER_OF_SYMBOLS];
            let mut last = 0usize;
            for byte in &bytes {
                totals[last] += 1;
                last = *byte as usize;
            }
            totals[0] += 3;
            current = Some((label.to_string(), totals));
        }
        let (_, totals) = current.as_ref().expect("set above");
        assert_eq!(totals[context].to_string(), row[2], "{label}/{context}");
        compared += 1;
    }
    assert_eq!(compared, 583, "context totals compared");
}

/// The floating-point normalisation, context by context and symbol by symbol, for the inputs whose
/// alphabet is narrow enough for the golden to carry the whole table.
#[test]
fn the_normalised_table_is_the_reference_arithmetic() {
    let corpus = corpus();
    let mut current: Option<(String, Table, Table)> = None;
    let mut compared = 0;

    for row in rows(&corpus, "norm") {
        let label = row[0];
        let context: usize = row[1].parse().expect("context");
        let symbol: usize = row[2].parse().expect("symbol");
        if current.as_ref().map(|(l, _, _)| l.as_str()) != Some(label) {
            let bytes = input(label);
            let mut raw: Table = vec![[0i32; NUMBER_OF_SYMBOLS]; NUMBER_OF_SYMBOLS];
            let mut last = 0usize;
            for byte in &bytes {
                raw[last][*byte as usize] += 1;
                last = *byte as usize;
            }
            let quarter = bytes.len() >> 2;
            for multiple in 1..=3 {
                raw[0][bytes[quarter * multiple] as usize] += 1;
            }
            current = Some((label.to_string(), raw, calc_frequencies_order1(&bytes)));
        }
        let (_, raw, normalised) = current.as_ref().expect("set above");
        assert_eq!(
            raw[context][symbol].to_string(),
            row[3],
            "{label}/{context}/{symbol}: raw count"
        );
        assert_eq!(
            normalised[context][symbol].to_string(),
            row[4],
            "{label}/{context}/{symbol}: normalised"
        );
        compared += 1;
    }
    assert_eq!(compared, 149, "table entries compared");
}

/// The whole 256 by 256 table for every input, including the ones too wide to print. It fails on
/// one wrong frequency anywhere in it.
#[test]
fn the_whole_table_matches_even_where_it_is_too_wide_to_print() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "normdigest") {
        let label = row[0];
        let bytes = input(label);
        let table = calc_frequencies_order1(&bytes);
        let contexts = table
            .iter()
            .filter(|row| row.iter().any(|f| *f != 0))
            .count();
        assert_eq!(contexts.to_string(), row[1], "{label}: contexts used");
        assert_eq!(table_digest(&table), row[2], "{label}: table digest");
        compared += 1;
    }
    assert_eq!(compared, 16, "tables compared");
}

/// The frequency table, with both levels of run length.
#[test]
fn the_frequency_table_is_written_byte_for_byte() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "freqtab") {
        let label = row[0];
        let table = calc_frequencies_order1(&input(label));
        let mut ours = Vec::new();
        let size = write_frequencies_order1(&table, &mut ours);
        assert_eq!(size.to_string(), row[1], "{label}: table size");
        let expected = if size <= MAX_INLINE_BYTES {
            hex(&ours)
        } else {
            sha256(&ours)
        };
        assert_eq!(expected, row[2], "{label}: table bytes");
        compared += 1;
    }
    assert_eq!(compared, 16, "frequency tables compared");
}

/// The four states at the head of the blob, after the same two cancelling reversals order 0 has.
#[test]
fn the_four_states_arrive_reversed_twice() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "states") {
        let label = row[0];
        let bytes = input(label);
        let ours = compress_order1(&bytes);
        let table = calc_frequencies_order1(&bytes);
        let mut written = Vec::new();
        let table_size = write_frequencies_order1(&table, &mut written);

        let at = PREFIX_LENGTH + table_size;
        let states: Vec<String> = (0..4)
            .map(|lane| {
                let start = at + lane * 4;
                u32::from_le_bytes(ours[start..start + 4].try_into().expect("four bytes"))
                    .to_string()
            })
            .collect();
        assert_eq!(states.join("\t"), row[1..5].join("\t"), "{label}");
        compared += 1;
    }
    assert_eq!(compared, 16, "state quadruples compared");
}

/// The reference decoded every one of its own streams, and so must this one.
#[test]
fn every_stream_decodes_back_to_its_input() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "roundtrip") {
        let (label, verdict) = (row[0], row[1]);
        assert_eq!(verdict, "ok", "{label}: the reference did not round trip");
        let bytes = input(label);
        assert_eq!(
            rans::uncompress(&compress_order1(&bytes)).map_err(|e| e.message()),
            Ok(bytes),
            "{label}"
        );
        compared += 1;
    }
    assert_eq!(compared, 16, "round trips compared");
}

/// A frequency byte of zero, which no writer emits, reads as the whole table. The golden's row
/// comes from a table built by hand in the dump, because no input can produce one.
#[test]
fn a_zero_frequency_byte_reads_as_four_thousand_and_ninety_six() {
    let corpus = corpus();
    let rows = rows(&corpus, "zerofreq");
    assert_eq!(rows.len(), 1, "zerofreq rows");
    let row = &rows[0];
    assert_eq!((row[0], row[1], row[2]), ("0", "65", "4096"));

    let mut at = 0usize;
    let stats =
        htsjdk_cram::rans_order1::read_stats_order1(&[0x00, 0x41, 0x00, 0x00, 0x00], &mut at)
            .expect("the hand-built table parses");
    assert_eq!(stats.frequency(0, 0x41).to_string(), row[2]);
}
