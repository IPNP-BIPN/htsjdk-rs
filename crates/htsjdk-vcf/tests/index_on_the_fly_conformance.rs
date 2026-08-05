//! Conformance for writing a VCF and its `.idx` in one pass, against
//! `IndexingVariantContextWriter`.
//!
//! Goldens from `tools/vcf-conformance/VcfIndexOnTheFlyDump.java` in the pinned oracle.
//!
//! The rows that justify the suite are the ones where the file's bytes and the index's numbers
//! have to agree, and where the layout was never the caller's choice:
//!
//! ```text
//! vcf   one-record         267    242  1     242
//! idx   one-record         LinearIndex
//! idx   many-records       IntervalTreeIndex
//! idx   header-only        IntervalTreeIndex
//! prop  one-record         flags=0;DICT:chr1=100000;DICT:chr2=200000;FEATURE_LENGTH_MEAN=1.0;...
//! ```
//!
//! 242 is the header's length, so the first record's position is not zero; the dictionary is four
//! properties and not a flag; and the same writer produced two different layouts from two files
//! that differ only in how many records they hold.

use std::io::Read;

use htsjdk_tribble::index::{TribbleIndex, INTERVAL_TREE, LINEAR};
use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType, VcfHeader};
use htsjdk_vcf::indexed_writer::{write_vcf_indexed, SequenceEntry};
use htsjdk_vcf::variant::{Value, VariantContext};

/// `VcfIndexOnTheFlyDump.header()`.
fn header() -> VcfHeader {
    let mut h = VcfHeader::new();
    h.lines.push(HeaderLine::info(
        "DP",
        Cardinality::Fixed(1),
        LineType::Integer,
        "Depth",
    ));
    h.lines.push(HeaderLine::info(
        "NOTE",
        Cardinality::Fixed(1),
        LineType::String,
        "A note",
    ));
    h.lines.push(HeaderLine::contig("chr1", 100000, 0));
    h.lines.push(HeaderLine::contig("chr2", 200000, 1));
    h
}

fn dictionary() -> Vec<SequenceEntry> {
    vec![
        SequenceEntry {
            name: "chr1".into(),
            length: 100000,
        },
        SequenceEntry {
            name: "chr2".into(),
            length: 200000,
        },
    ]
}

fn vc(contig: &str, start: i64, reference: &str, alt: &str) -> VariantContext {
    let mut record = VariantContext::new(
        contig,
        start,
        vec![
            Allele::from_str(reference, true).unwrap(),
            Allele::from_str(alt, false).unwrap(),
        ],
    );
    record.stop = start + reference.len() as i64 - 1;
    record.attributes = vec![("DP".to_string(), Value::Str("10".to_string()))];
    record
}

fn long_info(contig: &str, start: i64) -> VariantContext {
    let mut record = VariantContext::new(
        contig,
        start,
        vec![
            Allele::from_str("A", true).unwrap(),
            Allele::from_str("T", false).unwrap(),
        ],
    );
    record.attributes = vec![("NOTE".to_string(), Value::Str("x".repeat(400)))];
    record
}

/// Label and records, in the dump's order. The two refusals are handled separately because they
/// never reach a file.
fn cases() -> Vec<(&'static str, Vec<VariantContext>)> {
    vec![
        (
            "two-contigs",
            vec![
                vc("chr1", 100, "A", "T"),
                vc("chr1", 20000, "C", "G"),
                vc("chr2", 50, "GG", "G"),
            ],
        ),
        ("one-record", vec![vc("chr1", 100, "A", "T")]),
        (
            "many-records",
            (0..2000)
                .map(|i| vc("chr1", 100 + i * 5, "A", "T"))
                .collect(),
        ),
        ("header-only", Vec::new()),
        ("undeclared-contig", vec![vc("chrX", 9, "A", "T")]),
        (
            "uneven-lines",
            vec![
                vc("chr1", 100, "A", "T"),
                long_info("chr1", 200),
                vc("chr1", 300, "A", "T"),
            ],
        ),
    ]
}

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/vcf_index_on_the_fly.txt.gz");
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

fn decode_base64(text: &str) -> Vec<u8> {
    const ALPHABET: &[u8] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = Vec::with_capacity(text.len() * 3 / 4);
    let mut accumulator = 0u32;
    let mut bits = 0u32;
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

#[test]
fn the_file_and_its_index_agree_with_the_reference() {
    let corpus = corpus();
    let mut compared = 0;
    let mut trees = 0;

    for (label, records) in cases() {
        let golden_vcf = row(&corpus, &format!("vcf\t{label}\t"))
            .unwrap_or_else(|| panic!("{label}: no vcf row"));
        let golden_idx = row(&corpus, &format!("idx\t{label}\t"))
            .unwrap_or_else(|| panic!("{label}: no idx row"));
        let golden_prop = row(&corpus, &format!("prop\t{label}\t"))
            .unwrap_or_else(|| panic!("{label}: no prop row"));

        let idx_fields: Vec<&str> = golden_idx.split('\t').collect();
        let golden_bytes = decode_base64(idx_fields[2]);
        let reference = TribbleIndex::read(&golden_bytes).expect("the golden index parses");

        let written = write_vcf_indexed(
            &header(),
            &records,
            Some(&dictionary()),
            &reference.indexed_path,
        )
        .expect("the fixture writes");

        // The file: its length, the header's length, and every record's offset.
        let header_length = header().write().len();
        let positions = if written.record_positions.is_empty() {
            "-".to_string()
        } else {
            written
                .record_positions
                .iter()
                .map(i64::to_string)
                .collect::<Vec<_>>()
                .join(",")
        };
        assert_eq!(
            format!(
                "{}\t{}\t{}\t{positions}",
                written.text.len(),
                header_length,
                records.len()
            ),
            golden_vcf,
            "{label}: the file's shape"
        );

        // The layout, which the writer chose and the caller did not.
        let expected_class = match written.index.index_type {
            LINEAR => "LinearIndex",
            INTERVAL_TREE => {
                trees += 1;
                "IntervalTreeIndex"
            }
            other => panic!("{label}: unexpected index type {other}"),
        };
        assert_eq!(idx_fields[0], expected_class, "{label}: index type");

        // The properties, with the flags beside them because the dictionary used to live there.
        let rendered = written
            .index
            .properties
            .iter()
            .map(|(key, value)| format!("{key}={value}"))
            .collect::<Vec<_>>()
            .join(";");
        assert_eq!(
            format!("flags={};{rendered}", written.index.flags),
            golden_prop,
            "{label}: the header properties"
        );

        // And the bytes. The path, size, timestamp and MD5 are facts about the run rather than
        // decisions of the writer, so they come from the golden; everything else is produced here.
        let timestamp_offset: usize = idx_fields[1].parse().expect("an offset");
        let mine = TribbleIndex {
            index_type: written.index.index_type,
            properties: written.index.properties.clone(),
            contigs: written.index.contigs.clone(),
            interval_contigs: written.index.interval_contigs.clone(),
            ..reference.clone()
        };
        let mut bytes = mine.write().expect("an index writes");
        bytes[timestamp_offset..timestamp_offset + 8].fill(0);
        assert_eq!(
            bytes.len(),
            golden_bytes.len(),
            "{label}: the port wrote {} index bytes and the reference wrote {}",
            bytes.len(),
            golden_bytes.len()
        );
        assert_eq!(
            bytes
                .iter()
                .zip(golden_bytes.iter())
                .position(|(a, b)| a != b),
            None,
            "{label}: first differing index byte"
        );
        compared += 1;
    }

    assert_eq!(compared, 6, "files compared");
    assert_eq!(trees, 3, "files the writer indexed with an interval tree");
}

/// Indexing is on by default, so both refusals are reached by a caller who asked for nothing
/// unusual.
#[test]
fn the_two_refusals_are_the_reference_s() {
    let corpus = corpus();

    let no_dictionary = write_vcf_indexed(&header(), &[], None, "file:///x.vcf")
        .expect_err("indexing needs a dictionary");
    let golden = row(&corpus, "err\tno-dictionary\t").expect("a row");
    assert_eq!(
        format!("{}\t{}", no_dictionary.class(), no_dictionary.message()),
        golden
    );

    let to_a_stream = write_vcf_indexed(&header(), &[], Some(&dictionary()), "")
        .expect_err("indexing needs a path");
    let golden = row(&corpus, "err\tto-a-stream\t").expect("a row");
    assert_eq!(
        format!("{}\t{}", to_a_stream.class(), to_a_stream.message()),
        golden,
        "documented as an IllegalArgumentException and measured as an NPE"
    );
}

/// Turning indexing off writes the same file and no index, which is what says the index costs the
/// VCF nothing.
#[test]
fn the_index_does_not_change_the_vcf() {
    let corpus = corpus();
    let records = [vc("chr1", 100, "A", "T")];
    let indexed = write_vcf_indexed(
        &header(),
        &records,
        Some(&dictionary()),
        "file:///indexing-off.vcf",
    )
    .expect("writes");
    let plain = htsjdk_vcf::vcf_file::write_vcf(&header(), &records).expect("writes");
    assert_eq!(indexed.text, plain);

    let golden = row(&corpus, "vcf\tindexing-off\t").expect("a row");
    assert_eq!(
        format!(
            "{}\t{}\t1\t{}",
            plain.len(),
            header().write().len(),
            header().write().len()
        ),
        golden
    );
    assert_eq!(
        row(&corpus, "idx\tindexing-off\t").expect("a row"),
        "none\t-\t-",
        "no index is written when the option is unset"
    );
}
