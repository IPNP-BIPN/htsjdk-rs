//! Conformance for rANS 4x8 order 0, against
//! `htsjdk.samtools.cram.compression.rans.rans4x8.RANS4x8Encode`.
//!
//! Goldens from `tools/cram-conformance/CramRansDump.java` in the pinned oracle.
//!
//! The rows that justify the suite are the arithmetic, not the layout:
//!
//! ```text
//! orderused  3  ONE  0
//! orderused  4  ONE  1
//! norm  skewed-256  254  254  31
//! norm  skewed-256  255  255  152
//! freqtab  uniform-1000  4  41900000
//! states  uniform-1000  8388608  8388608  8388608  8388608
//! ```
//!
//! An order-1 request on three bytes writes a stream that says order 0. One symbol absorbs the
//! whole normalisation residue, and it is a factor of five rather than a nudge. And a uniform
//! input's four states never move, so its blob is the four initial lower bounds.
//!
//! The inputs are rebuilt here rather than carried in the golden, and each one's `in` row records
//! a sha256 that fails if the reconstruction drifts from the Java that produced the corpus.

use std::io::Read;

use htsjdk_cram::rans::{
    self, calc_frequencies_order0, compress_order0, order_used, write_frequencies_order0,
    EncodingSymbol, Order, NUMBER_OF_SYMBOLS, PREFIX_LENGTH, TOTAL_FREQ_SHIFT,
};
use sha2::{Digest, Sha256};

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/cram_rans.txt.gz");
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

/// `source` repeated until exactly `length` bytes have been produced.
fn repeat(source: &[u8], length: usize) -> Vec<u8> {
    (0..length).map(|i| source[i % source.len()]).collect()
}

/// `count` occurrences of each symbol in `[from, to)`, interleaved.
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
        "empty" => Vec::new(),
        "one-byte" => b"A".to_vec(),
        "two-bytes" => b"AB".to_vec(),
        "three-bytes" => b"ACG".to_vec(),
        "four-bytes" => b"ACGT".to_vec(),
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
        "all-256-once" => (0..=255u8).collect(),
        "skewed-256" => {
            let mut out = Vec::new();
            for symbol in 0..NUMBER_OF_SYMBOLS {
                out.extend(std::iter::repeat_n(symbol as u8, symbol));
            }
            out
        }
        "contiguous-run" => spread(5, 13, 100),
        "high-run" => spread(250, 256, 100),
        "zero-heavy" => {
            let mut out = vec![0u8; 1000];
            for i in 0..40usize {
                out[i * 25] = 1 + (i % 3) as u8;
            }
            out
        }
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
    assert_eq!(checked, 18, "inputs rebuilt");
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
        let ours = compress_order0(&input(label));
        assert_eq!(ours.len(), length, "{label}: compressed length");
        assert_eq!(sha256(&ours), digest, "{label}: compressed digest");
        if let Some(golden) = inline.get(label) {
            assert_eq!(&hex(&ours), golden, "{label}: compressed bytes");
        }
        compared += 1;
    }
    assert_eq!(compared, 18, "streams compared");
}

/// The nine-byte prefix: an order byte, then the compressed and raw lengths, little-endian. The
/// order it carries is the order that was used, which is not always the order that was asked for.
#[test]
fn the_prefix_records_the_order_that_was_used() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "prefix") {
        let label = row[0];
        let ours = compress_order0(&input(label));
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
    assert_eq!(compared, 17, "prefixes compared");
}

/// The fixed-point normalisation, symbol by symbol, beside the raw counts it came from. This is
/// where a port diverges first and where a byte-level golden alone would not say why.
#[test]
fn the_normalised_frequencies_are_the_reference_arithmetic() {
    let corpus = corpus();
    let mut current: Option<(String, [i32; NUMBER_OF_SYMBOLS], [i32; NUMBER_OF_SYMBOLS])> = None;
    let mut compared = 0;

    for row in rows(&corpus, "norm") {
        let label = row[0];
        let symbol: usize = row[1].parse().expect("symbol");
        if current.as_ref().map(|(l, _, _)| l.as_str()) != Some(label) {
            let bytes = input(label);
            let mut counts = [0i32; NUMBER_OF_SYMBOLS];
            for byte in &bytes {
                counts[*byte as usize] += 1;
            }
            current = Some((label.to_string(), counts, calc_frequencies_order0(&bytes)));
        }
        let (_, counts, normalised) = current.as_ref().expect("set above");
        assert_eq!(
            counts[symbol].to_string(),
            row[2],
            "{label}/{symbol}: count"
        );
        assert_eq!(
            normalised[symbol].to_string(),
            row[3],
            "{label}/{symbol}: normalised frequency"
        );
        compared += 1;
    }
    assert_eq!(compared, 822, "normalised frequencies compared");
}

/// The frequency table, including its run lengths, which start at the second consecutive symbol
/// and so encode a run of two as a run byte of zero.
#[test]
fn the_frequency_table_is_written_byte_for_byte() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "freqtab") {
        let label = row[0];
        let frequencies = calc_frequencies_order0(&input(label));
        let mut ours = Vec::new();
        let size = write_frequencies_order0(&frequencies, &mut ours);
        assert_eq!(size.to_string(), row[1], "{label}: table size");
        assert_eq!(hex(&ours), row[2], "{label}: table bytes");
        compared += 1;
    }
    assert_eq!(compared, 17, "frequency tables compared");
}

/// The encoding symbols, field by field: the reciprocal, its shift and the bias are the whole of
/// the divide-by-multiplying, and they are the part a port gets subtly wrong.
#[test]
fn every_encoding_symbol_matches_field_for_field() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "sym") {
        let label = row[0];
        let start: i32 = row[2].parse().expect("start");
        let frequency: i32 = row[3].parse().expect("freq");
        let ours = EncodingSymbol::set(start, frequency, TOTAL_FREQ_SHIFT);
        let mine = format!(
            "{}\t{}\t{}\t{}\t{}",
            ours.x_max, ours.rcp_freq, ours.bias, ours.cmpl_freq, ours.rcp_shift
        );
        assert_eq!(mine, row[4..9].join("\t"), "{label}/{}", row[1]);
        compared += 1;
    }
    assert_eq!(compared, 822, "encoding symbols compared");
}

/// The four states at the head of the blob, little-endian, in the order `rans0, rans1, rans2,
/// rans3`. They were written big-endian in the opposite order and the whole blob was then
/// reversed: two reversals that cancel, and a port that performs only one of them puts sixteen
/// plausible bytes in the wrong place.
#[test]
fn the_four_states_arrive_reversed_twice() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "states") {
        let label = row[0];
        let bytes = input(label);
        let ours = compress_order0(&bytes);
        let frequencies = calc_frequencies_order0(&bytes);
        let mut table = Vec::new();
        let table_size = write_frequencies_order0(&frequencies, &mut table);

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
    assert_eq!(compared, 17, "state quadruples compared");
}

/// The reference decoded every one of its own streams. So must this one, and it must decode the
/// reference's bytes rather than only its own.
#[test]
fn every_stream_decodes_back_to_its_input() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "roundtrip") {
        let (label, verdict) = (row[0], row[1]);
        assert_eq!(verdict, "ok", "{label}: the reference did not round trip");
        let bytes = input(label);
        let compressed = compress_order0(&bytes);
        assert_eq!(
            rans::uncompress(&compressed).map_err(|e| e.message()),
            Ok(bytes),
            "{label}"
        );
        compared += 1;
    }
    assert_eq!(compared, 18, "round trips compared");
}

/// The order the writer used against the order it was asked for. Below four bytes the answer is
/// not the question, and an empty input produces no stream at all to carry an answer in.
#[test]
fn an_order_one_request_below_four_bytes_writes_order_zero() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "orderused") {
        let length: usize = row[0].parse().expect("length");
        let requested = match row[1] {
            "ZERO" => Order::Zero,
            "ONE" => Order::One,
            other => panic!("{other}: no such order"),
        };
        let written = if length == 0 {
            "-".to_string()
        } else {
            (order_used(requested, length) as i32).to_string()
        };
        assert_eq!(written, row[2], "{length} bytes, requested {:?}", requested);
        compared += 1;
    }
    assert_eq!(compared, 18, "order decisions compared");
}
