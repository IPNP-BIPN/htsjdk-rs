//! Conformance for the VCF header against htsjdk's `VCFWriter.writeHeader`.
//!
//! Goldens from `tools/vcf-conformance/VcfHeaderDump.java` in the pinned oracle.
//!
//! The property this exists for is the line ordering, which is the **opposite** of the SAM
//! header's (decision 0009): a `TreeSet` over the rendered string of each whole line, so plain
//! ASCII order over complete lines. The `shuffled` case is inserted in an order that is neither
//! sorted nor reverse-sorted, so a port preserving insertion order fails it.

use std::io::Read;

use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType, VcfHeader};

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/vcf_header.txt.gz");
    let f = std::fs::File::open(&p).expect("corpus");
    let mut s = String::new();
    flate2::read::GzDecoder::new(f)
        .read_to_string(&mut s)
        .expect("corpus is gzip");
    s
}

/// The harness escapes `\\`, `\t` and `\n`. Undoing only two of the three was a bug that made a
/// backslash-escaped quote look like a divergence.
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

fn goldens() -> Vec<(String, String)> {
    corpus()
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .map(|l| {
            let (name, body) = l.split_once('\t').expect("name TAB body");
            (name.to_string(), unescape(body))
        })
        .collect()
}

fn info(id: &str, n: Cardinality, t: LineType, d: &str) -> HeaderLine {
    HeaderLine::info(id, n, t, d)
}

/// Rebuilds the harness's cases, in order.
fn cases() -> Vec<(String, VcfHeader)> {
    let mut out = Vec::new();

    out.push(("minimal".to_string(), VcfHeader::new()));

    let mut one = VcfHeader::new();
    one.lines = vec![
        info(
            "DP",
            Cardinality::Fixed(1),
            LineType::Integer,
            "Total Depth",
        ),
        HeaderLine::format("GT", Cardinality::Fixed(1), LineType::String, "Genotype"),
        HeaderLine::filter("q10", "Quality below 10"),
        HeaderLine::contig("chr1", 249_250_621, 0),
    ];
    out.push(("one_of_each".to_string(), one));

    let mut shuffled = VcfHeader::new();
    shuffled.lines = vec![
        info("ZZ", Cardinality::Fixed(1), LineType::Integer, "last by id"),
        HeaderLine::filter("zFilter", "z"),
        info("AA", Cardinality::Fixed(1), LineType::String, "first by id"),
        HeaderLine::format("ZQ", Cardinality::Fixed(1), LineType::Float, "z format"),
        HeaderLine::filter("aFilter", "a"),
        HeaderLine::format("AQ", Cardinality::Fixed(1), LineType::Float, "a format"),
    ];
    out.push(("shuffled".to_string(), shuffled));

    let mut quoting = VcfHeader::new();
    quoting.lines = vec![
        info(
            "NOSPACE",
            Cardinality::Fixed(1),
            LineType::String,
            "nospace",
        ),
        info(
            "WITHSPACE",
            Cardinality::Fixed(1),
            LineType::String,
            "with space",
        ),
        info(
            "WITHCOMMA",
            Cardinality::Fixed(1),
            LineType::String,
            "with,comma",
        ),
        info(
            "WITHQUOTE",
            Cardinality::Fixed(1),
            LineType::String,
            "with\"quote",
        ),
        HeaderLine::Unstructured {
            key: "unstructured".into(),
            value: "plain value".into(),
        },
        HeaderLine::Unstructured {
            key: "alsoUnstructured".into(),
            value: "value,with,commas".into(),
        },
    ];
    out.push(("quoting".to_string(), quoting));

    let mut counts = VcfHeader::new();
    counts.lines = vec![
        info("FIXED", Cardinality::Fixed(3), LineType::Integer, "three"),
        info("PERALT", Cardinality::A, LineType::Float, "per alt"),
        info("PERGT", Cardinality::G, LineType::Float, "per genotype"),
        info("PERALLELE", Cardinality::R, LineType::Float, "per allele"),
        info("UNBOUNDED", Cardinality::Unbounded, LineType::String, "any"),
        info("FLAG", Cardinality::Fixed(0), LineType::Flag, "a flag"),
    ];
    out.push(("cardinalities".to_string(), counts));

    let mut samples = VcfHeader::new();
    samples.lines = vec![HeaderLine::format(
        "GT",
        Cardinality::Fixed(1),
        LineType::String,
        "Genotype",
    )];
    samples.samples = vec!["NA12878".into(), "NA12891".into(), "NA12892".into()];
    out.push(("samples".to_string(), samples));

    let mut contigs = VcfHeader::new();
    contigs.lines = (1..=12)
        .map(|i| {
            HeaderLine::contig(
                &format!("chr{i}"),
                250_000_000 - i * 1_000_000,
                (i - 1) as i32,
            )
        })
        .collect();
    out.push(("contigs".to_string(), contigs));

    out
}

/// The gate: every header renders to exactly htsjdk's bytes.
#[test]
fn every_header_matches_htsjdks() {
    let cases = cases();
    let goldens = goldens();
    assert_eq!(cases.len(), goldens.len(), "case lists have drifted");

    let mut failures = Vec::new();
    for ((name, header), (golden_name, expected)) in cases.iter().zip(&goldens) {
        assert_eq!(name, golden_name, "case lists out of step");
        let ours = header.write();
        if &ours != expected {
            let at = ours
                .lines()
                .zip(expected.lines())
                .position(|(a, b)| a != b)
                .unwrap_or(0);
            failures.push(format!(
                "{name}: line {}\n  ours:   {:?}\n  htsjdk: {:?}",
                at + 1,
                ours.lines().nth(at).unwrap_or("<eof>"),
                expected.lines().nth(at).unwrap_or("<eof>")
            ));
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {} headers diverge:\n{}",
        failures.len(),
        cases.len(),
        failures.join("\n\n")
    );
}

/// The `shuffled` case must actually be out of order in the corpus, or it proves nothing about
/// the sort.
#[test]
fn the_shuffled_case_is_genuinely_unsorted_on_input() {
    let (_, golden) = goldens()
        .into_iter()
        .find(|(n, _)| n == "shuffled")
        .expect("the shuffled case");
    let meta: Vec<&str> = golden
        .lines()
        .filter(|l| l.starts_with("##") && !l.contains("fileformat"))
        .collect();
    let mut sorted = meta.clone();
    sorted.sort();
    assert_eq!(meta, sorted, "htsjdk's output must itself be sorted");
    // And the input order in `cases()` is not that order.
    let input_first = "##INFO=<ID=ZZ";
    assert!(
        !meta[0].starts_with(input_first),
        "the first line inserted must not be the first line written"
    );
}
