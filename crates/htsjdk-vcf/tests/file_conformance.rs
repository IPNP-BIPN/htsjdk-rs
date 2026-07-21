//! Conformance for whole VCF files against htsjdk's `VariantContextWriter`.
//!
//! Goldens from `tools/vcf-conformance/VcfFileDump.java` in the pinned oracle.
//!
//! The header and the data lines have their own suites. What this covers is the join: the
//! version line the writer substitutes for whatever the header carries, the absence of a
//! separator between the header and the first record, and the trailing newline.

use std::io::Read;

use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType, VcfHeader};
use htsjdk_vcf::variant::{Genotype, VariantContext};
use htsjdk_vcf::vcf_file::write_vcf;

fn corpus() -> Vec<(String, String)> {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/vcf_file.txt.gz");
    let f = std::fs::File::open(&p).expect("corpus");
    let mut s = String::new();
    flate2::read::GzDecoder::new(f)
        .read_to_string(&mut s)
        .expect("corpus is gzip");
    s.lines()
        .filter(|l| !l.is_empty())
        .map(|l| {
            let (name, value) = l.split_once('\t').expect("name TAB value");
            (name.to_string(), unescape(value))
        })
        .collect()
}

fn unescape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut chars = s.chars();
    while let Some(c) = chars.next() {
        if c != '\\' {
            out.push(c);
            continue;
        }
        match chars.next() {
            Some('n') => out.push('\n'),
            Some('t') => out.push('\t'),
            Some('\\') => out.push('\\'),
            Some(other) => {
                out.push('\\');
                out.push(other);
            }
            None => out.push('\\'),
        }
    }
    out
}

fn allele(s: &str, is_ref: bool) -> Allele {
    Allele::from_str(s, is_ref).unwrap()
}

fn header(with_samples: bool) -> VcfHeader {
    let mut h = VcfHeader::new();
    h.lines.push(HeaderLine::info(
        "DP",
        Cardinality::Fixed(1),
        LineType::Integer,
        "Total Depth",
    ));
    h.lines.push(HeaderLine::format(
        "GT",
        Cardinality::Fixed(1),
        LineType::String,
        "Genotype",
    ));
    h.lines.push(HeaderLine::format(
        "GQ",
        Cardinality::Fixed(1),
        LineType::Integer,
        "Genotype Quality",
    ));
    h.lines.push(HeaderLine::filter("q10", "Quality below 10"));
    for i in 0..3 {
        h.lines
            .push(HeaderLine::contig(&format!("chr{}", i + 1), 1000, i));
    }
    if with_samples {
        h.samples = vec!["s1".to_string(), "s2".to_string()];
    }
    h
}

fn vc(contig: &str, start: i64, r: &str, a: &str) -> VariantContext {
    VariantContext::new(contig, start, vec![allele(r, true), allele(a, false)])
}

fn cases() -> Vec<(String, String)> {
    let mut out = Vec::new();
    let mut push = |name: &str, h: &VcfHeader, records: &[VariantContext]| {
        out.push((name.to_string(), write_vcf(h, records).expect(name)));
    };

    push("header_only", &header(false), &[]);
    push("header_only_with_samples", &header(true), &[]);
    push("one_record", &header(false), &[vc("chr1", 100, "A", "T")]);
    push(
        "many_records",
        &header(false),
        &[
            vc("chr1", 100, "A", "T"),
            vc("chr1", 200, "C", "G"),
            vc("chr2", 50, "GG", "G"),
            vc("chr3", 1, "T", "TTTT"),
        ],
    );

    let mut genotyped = vc("chr1", 100, "A", "T");
    let mut g1 = Genotype::new("s1", vec![allele("A", true), allele("T", false)]);
    g1.gq = Some(30);
    let mut g2 = Genotype::new("s2", vec![allele("A", true), allele("A", true)]);
    g2.gq = Some(40);
    genotyped.genotypes = vec![g1, g2];
    push("genotyped", &header(true), &[genotyped]);

    let mut with_format = VcfHeader::new();
    with_format.lines.push(HeaderLine::info(
        "DP",
        Cardinality::Fixed(1),
        LineType::Integer,
        "Total Depth",
    ));
    with_format.lines.push(HeaderLine::Unstructured {
        key: "fileformat".to_string(),
        value: "VCFv4.3".to_string(),
    });
    push(
        "header_declares_its_own_fileformat",
        &with_format,
        &[vc("chr1", 100, "A", "T")],
    );

    push(
        "unsorted_records",
        &header(false),
        &[
            vc("chr3", 500, "A", "T"),
            vc("chr1", 100, "A", "T"),
            vc("chr2", 300, "A", "T"),
        ],
    );

    out
}

#[test]
fn every_file_matches_htsjdks() {
    let golden = corpus();
    let ours = cases();
    assert_eq!(golden.len(), ours.len(), "case lists differ in length");

    let mut mismatches = Vec::new();
    for (i, (name, text)) in ours.iter().enumerate() {
        let (gname, gtext) = &golden[i];
        assert_eq!(gname, name, "case {i}: the lists have drifted");
        if gtext != text {
            let at = text
                .lines()
                .zip(gtext.lines())
                .position(|(a, b)| a != b)
                .unwrap_or(0);
            mismatches.push(format!(
                "{name}: first difference at line {at}\n  htsjdk: {:?}\n  ours  : {:?}",
                gtext.lines().nth(at).unwrap_or("<end>"),
                text.lines().nth(at).unwrap_or("<end>")
            ));
        }
    }
    assert!(
        mismatches.is_empty(),
        "{} of {} files diverge:\n{}",
        mismatches.len(),
        ours.len(),
        mismatches.join("\n")
    );
}

/// The writer does not sort. A file whose records are out of coordinate order comes back out in
/// the order it was given, which is what makes sorting the caller's problem and not the layer's.
#[test]
fn the_writer_preserves_record_order() {
    let golden = corpus();
    let (_, text) = golden
        .iter()
        .find(|(n, _)| n == "unsorted_records")
        .expect("unsorted_records");
    let contigs: Vec<&str> = text
        .lines()
        .filter(|l| !l.starts_with('#'))
        .filter_map(|l| l.split('\t').next())
        .collect();
    assert_eq!(contigs, ["chr3", "chr1", "chr2"]);
}
