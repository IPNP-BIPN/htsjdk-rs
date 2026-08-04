//! Conformance for `VariantContextComparator` against htsjdk 4.2.0.
//!
//! Golden from `tools/vcf-conformance/VariantComparatorDump.java`.
//!
//! The class is thirty lines of Java and the suite is here because three of its behaviours are
//! decisions rather than consequences, and a port that tidied any of them would still compile and
//! still sort:
//!
//!  * **the two constructors refuse different things, and word the empty case differently.** This
//!    port returned the contig-list message for both until this golden said otherwise. Nothing in
//!    the class explains the second sentence; someone wrote it, and reproducing it is the job;
//!  * **an unknown contig throws.** htsjdk's own comment says it does so "happily", so the
//!    exception is the behaviour and not an oversight;
//!  * **`compare` returns a subtraction**, so the magnitude is observable. A comparator built from
//!    header lines with sparse indexes returns the difference of those *indexes*, not of their
//!    positions in the collection.

use std::io::Read;

use htsjdk_vcf::comparator::{ComparatorError, VariantContextComparator};
use htsjdk_vcf::header::HeaderLine;
use htsjdk_vcf::variant::VariantContext;

fn golden() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/variant_comparator.txt.gz");
    let file = std::fs::File::open(&path).expect("golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("golden is gzip");
    text
}

fn contig_line(id: &str, index: i32) -> HeaderLine {
    HeaderLine::Contig {
        index,
        fields: vec![
            ("ID".to_string(), id.to_string()),
            ("length".to_string(), "1000".to_string()),
        ],
    }
}

/// The inputs each labelled constructor case was given. A label is a configuration and the row
/// carries nothing to derive it from, so it is written beside the dump that produced it.
fn constructor(label: &str) -> Result<VariantContextComparator, ComparatorError> {
    let names =
        |names: &[&str]| -> Vec<String> { names.iter().map(|n| (*n).to_string()).collect() };
    match label {
        "one" => VariantContextComparator::from_contigs(&names(&["chr1"])),
        "three" => VariantContextComparator::from_contigs(&names(&["chr1", "chr2", "chr3"])),
        "empty" => VariantContextComparator::from_contigs(&[]),
        "duplicate-name" => VariantContextComparator::from_contigs(&names(&["chr1", "chr1"])),
        "lines-two" => VariantContextComparator::from_header_lines(&[
            contig_line("chr1", 0),
            contig_line("chr2", 1),
        ]),
        "lines-sparse" => VariantContextComparator::from_header_lines(&[
            contig_line("chr1", 5),
            contig_line("chr2", 9),
        ]),
        "lines-reversed" => VariantContextComparator::from_header_lines(&[
            contig_line("chr1", 1),
            contig_line("chr2", 0),
        ]),
        "lines-empty" => VariantContextComparator::from_header_lines(&[]),
        "lines-duplicate-name" => VariantContextComparator::from_header_lines(&[
            contig_line("chr1", 0),
            contig_line("chr1", 1),
        ]),
        "lines-shared-index" => VariantContextComparator::from_header_lines(&[
            contig_line("chr1", 0),
            contig_line("chr2", 0),
        ]),
        // The `cmp` rows use their own labels: `list` is the four-contig comparator the dump
        // built once and reused for every pair.
        "list" => VariantContextComparator::from_contigs(&names(&["chr1", "chr2", "chr3", "chr4"])),
        other => panic!("{other} is in the golden but not configured here"),
    }
}

/// `chr1:100` and the rest, as the dump names them.
fn variant(key: &str) -> VariantContext {
    let (contig, start) = key.split_once(':').expect("contig:start");
    VariantContext::new(contig, start.parse().expect("a start"), Vec::new())
}

#[test]
fn every_constructor_accepts_or_refuses_as_the_reference_does() {
    let text = golden();
    let mut compared = 0usize;
    let mut refusals = 0usize;

    for line in text.lines() {
        let Some(rest) = line.strip_prefix("ctor\t") else {
            continue;
        };
        let (label, expected) = rest.split_once('\t').expect("label and outcome");
        match (constructor(label), expected) {
            (Ok(_), "ok") => {}
            (Err(error), expected) => {
                let body = expected
                    .strip_prefix("E:")
                    .unwrap_or_else(|| panic!("{label}: the reference accepted, the port refused"));
                let (class, message) = body.split_once(':').expect("class:message");
                assert_eq!(error.class(), class, "{label}: exception class");
                assert_eq!(error.message(), message, "{label}: exception message");
                refusals += 1;
            }
            (Ok(_), expected) => {
                panic!("{label}: the reference said {expected}, the port accepted")
            }
        }
        compared += 1;
    }

    assert_eq!(compared, 10, "the golden changed size");
    // Five refusals across four distinct messages: the two empty cases word it differently, and
    // that difference is the finding this suite exists for.
    assert_eq!(
        refusals, 5,
        "the golden lost a refusal; the two constructors between them refuse in five places"
    );
    println!("constructors: {compared} cases, {refusals} refusals, each with its own message");
}

#[test]
fn compare_returns_the_references_number_and_not_merely_its_sign() {
    let text = golden();
    let mut compared = 0usize;
    let mut thrown = 0usize;

    for line in text.lines() {
        let Some(rest) = line.strip_prefix("cmp\t") else {
            continue;
        };
        let fields: Vec<&str> = rest.split('\t').collect();
        let (label, first, second, expected) = (fields[0], fields[1], fields[2], fields[3]);
        let comparator = constructor(label).expect("the comparator this label built");

        match comparator.compare(&variant(first), &variant(second)) {
            Ok(value) => {
                assert_eq!(
                    value.to_string(),
                    expected,
                    "{label}: compare({first}, {second}) is a subtraction, so the magnitude counts"
                );
            }
            Err(unknown) => {
                assert_eq!(
                    expected, "E:java.lang.NullPointerException",
                    "{label}: {} is unknown here but the reference ordered it",
                    unknown.contig
                );
                thrown += 1;
            }
        }
        compared += 1;
    }

    assert_eq!(compared, 53, "the golden changed size");
    assert!(
        thrown > 0,
        "the golden lost the unknown-contig rows, which are the case htsjdk throws on deliberately"
    );
    println!("compare: {compared} pairs, {thrown} refused for an unknown contig");
}

#[test]
fn compatibility_needs_the_index_and_not_only_the_name() {
    let text = golden();
    let comparator = constructor("lines-sparse").expect("the sparse comparator");
    let mut compared = 0usize;

    for line in text.lines() {
        let Some(rest) = line.strip_prefix("compat\t") else {
            continue;
        };
        let (label, expected) = rest.split_once('\t').expect("label and outcome");
        let lines = match label {
            "same" => vec![contig_line("chr1", 5)],
            "other-index" => vec![contig_line("chr1", 0)],
            "unknown-name" => vec![contig_line("chrX", 5)],
            "empty" => vec![],
            other => panic!("{other} is in the golden but not configured here"),
        };
        assert_eq!(
            comparator.is_compatible(&lines).to_string(),
            expected,
            "{label}"
        );
        compared += 1;
    }

    assert_eq!(compared, 4, "the golden changed size");
    println!("isCompatible: {compared} cases");
}
