//! Conformance for `VCFUtils.smartMergeHeaders` against htsjdk 4.2.0.
//!
//! Golden from `tools/vcf-conformance/SmartMergeDump.java`.
//!
//! The suite exists because three of the method's behaviours are not what reading it suggests, and
//! this port had all three wrong before the dump ran:
//!
//!  * **every output carries a `fileformat` line no source wrote.** `getMetaDataInSortedOrder`
//!    prepends one, so two headers holding one `INFO` line apiece merge to **two** lines;
//!  * **the version is a field, not a line.** A header assembled in memory has none whatever its
//!    `##fileformat` line says, so it prepends `VCFv4.2` and never reaches the version policy;
//!  * **both promotion arms are no-ops.** An Integer seen first stays Integer and a Float seen
//!    first stays Float; the Java's `put` writes back what the map already holds.
//!
//! The rendered lines are compared as strings, because that is what the merge is for: the result
//! is written into a header, and a difference that does not reach the rendering is not a
//! difference in the file.

use std::io::Read;

use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType, VcfHeader};
use htsjdk_vcf::merge::{smart_merge_headers, MergeError, Source, VCF4_2, VCF4_3};

fn golden() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/smart_merge_headers.txt.gz");
    let file = std::fs::File::open(&path).expect("golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("golden is gzip");
    text
}

fn info(id: &str, number: Cardinality, line_type: LineType, description: &str) -> HeaderLine {
    HeaderLine::Compound {
        key: "INFO".to_string(),
        id: id.to_string(),
        number,
        line_type,
        description: description.to_string(),
        extra: Vec::new(),
    }
}

fn format(id: &str, number: Cardinality, line_type: LineType, description: &str) -> HeaderLine {
    HeaderLine::Compound {
        key: "FORMAT".to_string(),
        id: id.to_string(),
        number,
        line_type,
        description: description.to_string(),
        extra: Vec::new(),
    }
}

fn filter(id: &str, description: &str) -> HeaderLine {
    HeaderLine::Filter {
        id: id.to_string(),
        description: description.to_string(),
    }
}

fn contig(id: &str, index: i32) -> HeaderLine {
    HeaderLine::Contig {
        index,
        fields: vec![
            ("ID".to_string(), id.to_string()),
            ("length".to_string(), "1000".to_string()),
        ],
    }
}

fn unstructured(key: &str, value: &str) -> HeaderLine {
    HeaderLine::Unstructured {
        key: key.to_string(),
        value: value.to_string(),
    }
}

fn header(lines: Vec<HeaderLine>) -> VcfHeader {
    VcfHeader {
        lines,
        samples: Vec::new(),
    }
}

/// The two sources each labelled case was given, and the version each carried.
///
/// A label is a configuration and the row carries nothing to derive it from, so it is written here
/// beside the dump that produced it. The `version` is `None` for every case built from a set of
/// lines, which is the finding: htsjdk leaves the field null there.
fn sources(label: &str) -> (Vec<VcfHeader>, Vec<Option<&'static str>>) {
    let dp_int = info("DP", Cardinality::Fixed(1), LineType::Integer, "depth");
    let dp_other = info("DP", Cardinality::Fixed(1), LineType::Integer, "read depth");
    let af_one = info("AF", Cardinality::Fixed(1), LineType::Float, "af");
    let af_two = info("AF", Cardinality::Fixed(2), LineType::Float, "af");
    let x_int = info("X", Cardinality::Fixed(1), LineType::Integer, "x");
    let x_float = info("X", Cardinality::Fixed(1), LineType::Float, "x");
    let x_string = info("X", Cardinality::Fixed(1), LineType::String, "x");
    let gq = format("GQ", Cardinality::Fixed(1), LineType::Integer, "gq");

    let two =
        |a: Vec<HeaderLine>, b: Vec<HeaderLine>| (vec![header(a), header(b)], vec![None, None]);

    match label {
        "identical" => two(vec![dp_int.clone()], vec![dp_int]),
        "disjoint" => two(vec![dp_int], vec![gq]),
        "number-differs" => two(vec![af_one], vec![af_two]),
        "number-differs-reversed" => two(vec![af_two], vec![af_one]),
        "int-then-float" => two(vec![x_int], vec![x_float]),
        "float-then-int" => two(vec![x_float], vec![x_int]),
        "int-then-string" => two(vec![x_int], vec![x_string]),
        "description-differs" => two(vec![dp_int], vec![dp_other]),
        "filter-description-differs" => two(
            vec![filter("LowQual", "low quality")],
            vec![filter("LowQual", "poor")],
        ),
        "unstructured-conflict" => two(
            vec![unstructured("source", "one")],
            vec![unstructured("source", "two")],
        ),
        "order-sorted-within-source" => two(
            vec![
                info("B", Cardinality::Fixed(1), LineType::Integer, "b"),
                info("A", Cardinality::Fixed(1), LineType::Integer, "a"),
            ],
            vec![info("C", Cardinality::Fixed(1), LineType::Integer, "c")],
        ),
        "contigs" => two(
            vec![contig("chr2", 1), contig("chr1", 0)],
            vec![contig("chr3", 2)],
        ),
        // The version cases: the dump built these through `new VCFHeader(set)`, which leaves the
        // version field null however the set is labelled. That is the whole point of the case.
        "version-43-and-42" | "version-43-and-43" | "version-42-and-41" => {
            two(vec![dp_int.clone()], vec![dp_int])
        }
        "version-43-and-empty" => two(vec![dp_int], vec![]),
        other => panic!("{other} is in the golden but not configured here"),
    }
}

fn run(label: &str) -> Result<String, MergeError> {
    let (headers, versions) = sources(label);
    let sources: Vec<Source> = headers
        .iter()
        .zip(&versions)
        .map(|(header, version)| Source {
            header,
            version: *version,
        })
        .collect();
    smart_merge_headers(&sources, true).map(|(merged, _)| {
        merged
            .iter()
            .map(HeaderLine::render)
            .collect::<Vec<_>>()
            .join("|")
    })
}

#[test]
fn every_merge_renders_as_the_reference_does() {
    let text = golden();
    let mut compared = 0usize;
    let mut refused = 0usize;

    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("merged\t") {
            let fields: Vec<&str> = rest.split('\t').collect();
            let (label, count, rendered) = (fields[0], fields[1], fields[2]);
            let ours = run(label).unwrap_or_else(|error| {
                panic!(
                    "{label}: the reference merged, the port refused: {}",
                    error.message()
                )
            });
            assert_eq!(ours, rendered, "{label}");
            assert_eq!(
                ours.split('|').count().to_string(),
                count,
                "{label}: line count"
            );
            compared += 1;
        } else if let Some(rest) = line.strip_prefix("err\t") {
            let (label, expected) = rest.split_once('\t').expect("label and error");
            let error = run(label).expect_err("the reference refused");
            let (class, message) = expected.split_once(':').expect("class:message");
            assert_eq!(error.class(), class, "{label}: exception class");
            assert_eq!(error.message(), message, "{label}: exception message");
            refused += 1;
            compared += 1;
        }
    }

    assert_eq!(compared, 16, "the golden changed size");
    assert!(refused > 0, "the golden lost the collision that throws");
    println!("smartMergeHeaders: {compared} merges, {refused} refused");
}

/// The prepended line is the finding, so it gets its own assertion rather than riding along inside
/// a rendered string.
#[test]
fn the_prepended_version_line_is_in_every_output_and_belongs_to_no_source() {
    let empty = header(vec![]);

    // No version: 4.2, and it is the only line even though the source had none.
    let (merged, _) = smart_merge_headers(
        &[Source {
            header: &empty,
            version: None,
        }],
        true,
    )
    .expect("merges");
    assert_eq!(merged.len(), 1);
    assert_eq!(merged[0].render(), format!("fileformat={VCF4_2}"));

    // A 4.3 header prepends 4.3 instead.
    let (merged, _) = smart_merge_headers(
        &[Source {
            header: &empty,
            version: Some(VCF4_3),
        }],
        true,
    )
    .expect("merges");
    assert_eq!(merged[0].render(), format!("fileformat={VCF4_3}"));
}

/// The version policy fires only when a version was actually set, which is why every `version-*`
/// case in the golden merges rather than throwing.
#[test]
fn the_version_policy_needs_a_version_that_was_set() {
    let a = header(vec![]);
    let b = header(vec![]);

    let error = smart_merge_headers(
        &[
            Source {
                header: &a,
                version: Some(VCF4_3),
            },
            Source {
                header: &b,
                version: Some(VCF4_2),
            },
        ],
        true,
    )
    .expect_err("refuses");
    assert_eq!(error.class(), "java.lang.IllegalArgumentException");

    // The same pair with no version set is what the golden's `version-43-and-42` actually is.
    assert!(run("version-43-and-42").is_ok());
}
