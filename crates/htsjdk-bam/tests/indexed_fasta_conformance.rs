//! Conformance for reading a FASTA through its `.fai`, against the oracle.
//!
//! Golden from `tools/bam-conformance/IndexedFastaDump.java`.
//!
//! The rows that a reader gets wrong by scanning for newlines instead of trusting the index:
//!
//! ```text
//! index  crlf.fasta  chr1  12  7  6  8      the terminator is TWO bytes, and only this says so
//! query  crlf.fasta  chr1:6-7  CG           a query across that boundary jumps two bytes
//! query  lf.fasta    chr1:5-4  (empty)      start == stop + 1 is legal and answers nothing
//! error  lf.fasta    chr1:6-4  Malformed query; start point 6 lies after end point 4
//! ```

use std::io::Read;

use htsjdk_bam::fasta_index::{FastaIndex, FastaIndexError, IndexedFasta};

fn golden() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/indexed_fasta.txt.gz");
    let file = std::fs::File::open(&path).expect("golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("golden is gzip");
    text
}

fn rows(text: &str, kind: &str) -> Vec<Vec<String>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').map(str::to_string).collect::<Vec<_>>())
        .filter(|parts| parts.first().map(String::as_str) == Some(kind))
        .collect()
}

/// The three files the dump wrote, byte for byte.
fn contents(file: &str) -> &'static [u8] {
    match file {
        "lf.fasta" => b">chr1\nACGTAC\nGTACGT\nAC\n>chr2\nTTTTTT\nGG\n",
        "crlf.fasta" => b">chr1\r\nACGTAC\r\nGTACGT\r\n",
        "single.fasta" => b">chrS\nACGTACGTAC\n",
        other => panic!("unknown fixture {other}"),
    }
}

/// The index as the golden reports it, which is the one htsjdk's creator produced.
fn index_of(text: &str, file: &str) -> FastaIndex {
    let mut entries = Vec::new();
    for row in rows(text, "index") {
        if row[1] != file {
            continue;
        }
        entries.push(htsjdk_bam::fasta_index::FastaIndexEntry {
            name: row[2].clone(),
            size: row[3].parse().expect("size"),
            location: row[4].parse().expect("location"),
            bases_per_line: row[5].parse().expect("bases per line"),
            bytes_per_line: row[6].parse().expect("bytes per line"),
        });
    }
    assert!(!entries.is_empty(), "no index rows for {file}");
    FastaIndex { entries }
}

fn parse_label(label: &str) -> (String, i64, i64) {
    let (contig, range) = label.rsplit_once(':').expect("contig:range");
    let (start, stop) = range.split_once('-').expect("start-stop");
    (
        contig.to_string(),
        start.parse().expect("start"),
        stop.parse().expect("stop"),
    )
}

fn escape(bases: &[u8]) -> String {
    String::from_utf8_lossy(bases)
        .replace('\\', "\\\\")
        .replace('\t', "\\t")
        .replace('\r', "\\r")
        .replace('\n', "\\n")
}

#[test]
fn every_query_is_the_reference() {
    let text = golden();
    let queries = rows(&text, "query");
    assert!(queries.len() >= 15, "the golden carries every case");

    for row in &queries {
        let (file, label) = (row[1].as_str(), row[2].as_str());
        let (contig, start, stop) = parse_label(label);
        let mut fasta =
            IndexedFasta::new(std::io::Cursor::new(contents(file)), index_of(&text, file));
        let bases = fasta
            .query(&contig, start, stop)
            .unwrap_or_else(|error| panic!("{file} {label}: {}", error.message()));
        let expected = row.get(3).cloned().unwrap_or_default();
        assert_eq!(escape(&bases), expected, "{file} {label}");
    }
}

#[test]
fn every_refusal_is_the_reference() {
    let text = golden();
    let errors = rows(&text, "error");
    assert_eq!(errors.len(), 3, "two past-the-end and one reversed query");

    for row in &errors {
        let (file, label, message) = (row[1].as_str(), row[2].as_str(), row[4].as_str());
        let (contig, start, stop) = parse_label(label);
        let mut fasta =
            IndexedFasta::new(std::io::Cursor::new(contents(file)), index_of(&text, file));
        let error = fasta
            .query(&contig, start, stop)
            .expect_err("this query is refused");
        assert_eq!(error.message(), message, "{file} {label}");
    }
}

/// The index the creator wrote is the one this port parses back.
#[test]
fn the_index_round_trips_through_the_parser() {
    let text = golden();
    let parsed = FastaIndex::parse("chr1\t14\t6\t6\t7\nchr2\t8\t29\t6\t7\n").expect("parses");
    assert_eq!(parsed, index_of(&text, "lf.fasta"));

    assert!(matches!(
        FastaIndex::parse("chr1\t14\t6\t6\n"),
        Err(FastaIndexError::Malformed(_))
    ));
}
