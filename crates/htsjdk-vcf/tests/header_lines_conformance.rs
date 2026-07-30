//! Conformance for typing one `##` line, against `VCFInfoHeaderLine`, `VCFFormatHeaderLine`,
//! `VCFFilterHeaderLine` and `VCFContigHeaderLine`.
//!
//! Goldens from `tools/vcf-conformance/VcfHeaderLineDump.java` in the pinned oracle.
//!
//! Four exception classes come out of this one layer, and the golden is what says which:
//!
//! ```text
//! Number=x      java.lang.NumberFormatException          nothing wraps Integer.parseInt
//! Number=-1     TribbleException$InvalidHeader           "Count < 0 for fixed size ..."
//! Number=0      java.lang.IllegalArgumentException       from validate(), non-Flag only
//! Type=integer  htsjdk.tribble.TribbleException          no "malformed header" prefix
//! ```
//!
//! The successful rows are compared as the **rendered line**, because that is what a writer emits
//! and it carries the quoting rule along with the values.

use std::io::Read;

use htsjdk_vcf::header::HeaderLine;
use htsjdk_vcf::header_lines::parse_meta_line;
use htsjdk_vcf::header_parse::VcfVersion;

/// The dump's cases, in its order: label, the whole `##` line, the version, and the contig index
/// the codec's counter would be at.
const CASES: &[(&str, &str, VcfVersion, i32)] = &[
    (
        "info-int",
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">",
        VcfVersion::Vcf4_2,
        0,
    ),
    (
        "info-a",
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">",
        VcfVersion::Vcf4_2,
        0,
    ),
    (
        "info-r",
        "##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Depths\">",
        VcfVersion::Vcf4_2,
        0,
    ),
    (
        "info-g",
        "##INFO=<ID=PL,Number=G,Type=Integer,Description=\"Likelihoods\">",
        VcfVersion::Vcf4_2,
        0,
    ),
    (
        "info-unbounded",
        "##INFO=<ID=X,Number=.,Type=String,Description=\"Any\">",
        VcfVersion::Vcf4_2,
        0,
    ),
    (
        "info-flag",
        "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"In dbSNP\">",
        VcfVersion::Vcf4_2,
        0,
    ),
    (
        "info-flag-count-2",
        "##INFO=<ID=DB,Number=2,Type=Flag,Description=\"In dbSNP\">",
        VcfVersion::Vcf4_2,
        0,
    ),
    (
        "info-source-42",
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\",Source=\"s\",Version=\"3\">",
        VcfVersion::Vcf4_2,
        0,
    ),
    (
        "info-source-41",
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\",Source=\"s\">",
        VcfVersion::Vcf4_1,
        0,
    ),
    (
        "info-no-description",
        "##INFO=<ID=DP,Number=1,Type=Integer>",
        VcfVersion::Vcf4_2,
        0,
    ),
    (
        "info-number-not-a-number",
        "##INFO=<ID=DP,Number=x,Type=Integer,Description=\"d\">",
        VcfVersion::Vcf4_2,
        0,
    ),
    (
        "info-number-negative",
        "##INFO=<ID=DP,Number=-1,Type=Integer,Description=\"d\">",
        VcfVersion::Vcf4_2,
        0,
    ),
    (
        "info-number-zero",
        "##INFO=<ID=DP,Number=0,Type=Integer,Description=\"d\">",
        VcfVersion::Vcf4_2,
        0,
    ),
    (
        "info-type-lowercase",
        "##INFO=<ID=DP,Number=1,Type=integer,Description=\"d\">",
        VcfVersion::Vcf4_2,
        0,
    ),
    ("info-no-number", "##INFO=<ID=DP>", VcfVersion::Vcf4_2, 0),
    (
        "info-id-angle",
        "##INFO=<ID=\"a<b\",Number=1,Type=Integer,Description=\"d\">",
        VcfVersion::Vcf4_2,
        0,
    ),
    (
        "info-id-equals",
        "##INFO=<ID=\"a=b\",Number=1,Type=Integer,Description=\"d\">",
        VcfVersion::Vcf4_2,
        0,
    ),
    (
        "format-int",
        "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Quality\">",
        VcfVersion::Vcf4_2,
        0,
    ),
    (
        "format-flag",
        "##FORMAT=<ID=DB,Number=0,Type=Flag,Description=\"In dbSNP\">",
        VcfVersion::Vcf4_2,
        0,
    ),
    (
        "filter-plain",
        "##FILTER=<ID=LowQual,Description=\"Low quality\">",
        VcfVersion::Vcf4_2,
        0,
    ),
    (
        "filter-no-description",
        "##FILTER=<ID=LowQual>",
        VcfVersion::Vcf4_2,
        0,
    ),
    (
        "filter-wrong-order",
        "##FILTER=<Description=\"d\",ID=LowQual>",
        VcfVersion::Vcf4_2,
        0,
    ),
    (
        "filter-no-id",
        "##FILTER=<Description=\"d\">",
        VcfVersion::Vcf4_2,
        0,
    ),
    (
        "filter-extra-tag",
        "##FILTER=<ID=LowQual,Description=\"d\",Extra=1>",
        VcfVersion::Vcf4_2,
        0,
    ),
    (
        "contig-plain",
        "##contig=<ID=chr1,length=1000>",
        VcfVersion::Vcf4_2,
        0,
    ),
    (
        "contig-reordered",
        "##contig=<length=1000,ID=chr1>",
        VcfVersion::Vcf4_2,
        3,
    ),
    (
        "contig-extra",
        "##contig=<ID=chr1,length=1000,assembly=b37,md5=abc>",
        VcfVersion::Vcf4_2,
        7,
    ),
    (
        "contig-no-id",
        "##contig=<length=1000>",
        VcfVersion::Vcf4_2,
        0,
    ),
    (
        "contig-negative-index",
        "##contig=<ID=chr1,length=1000>",
        VcfVersion::Vcf4_2,
        -1,
    ),
];

fn corpus() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/vcf_header_lines.txt.gz");
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

/// The Java class the reference builds for each key, which the dump reports alongside the line.
fn java_class(line: &HeaderLine) -> &'static str {
    match line {
        HeaderLine::Compound { key, .. } if key == "INFO" => "htsjdk.variant.vcf.VCFInfoHeaderLine",
        HeaderLine::Compound { key, .. } if key == "FORMAT" => {
            "htsjdk.variant.vcf.VCFFormatHeaderLine"
        }
        HeaderLine::Structured { key, .. } if key == "FILTER" => {
            "htsjdk.variant.vcf.VCFFilterHeaderLine"
        }
        HeaderLine::Contig { .. } => "htsjdk.variant.vcf.VCFContigHeaderLine",
        _ => "htsjdk.variant.vcf.VCFHeaderLine",
    }
}

#[test]
fn every_header_line_types_as_the_reference_types_it() {
    let text = corpus();

    for (label, line, version, contig_index) in CASES {
        match parse_meta_line(line, *version, *contig_index) {
            Ok(parsed) => {
                let ours = format!("{}\t{}", java_class(&parsed), parsed.render());
                let expected = row(&text, &format!("hline\t{label}\t")).unwrap_or_else(|| {
                    panic!("{label}: the port built a line, the golden has no `hline` row")
                });
                assert_eq!(ours, expected, "{label}");
            }
            Err(error) => {
                let expected = row(&text, &format!("hlineerror\t{label}\t")).unwrap_or_else(|| {
                    panic!("{label}: the port refused, the golden has no `hlineerror` row")
                });
                assert_eq!(
                    format!("{}\t{}", error.class(), error.message()),
                    expected,
                    "{label}"
                );
            }
        }
    }

    let in_golden = text
        .lines()
        .filter(|line| line.starts_with("hline\t") || line.starts_with("hlineerror\t"))
        .count();
    assert_eq!(
        in_golden,
        CASES.len(),
        "the golden and the test disagree on the case list"
    );
    println!("{} header lines identical", CASES.len());
}

/// The three rows that are values rather than errors, named so a regression in any of them is a
/// failing test with a reason attached rather than one line in a list of thirty.
#[test]
fn the_silent_rewrites_are_the_reference_rewrites() {
    let text = corpus();

    // A Flag with a count of 2 is stored as 0, not refused.
    let flag = row(&text, "hline\tinfo-flag-count-2\t").expect("the golden carries it");
    assert!(
        flag.ends_with("Number=0,Type=Flag,Description=\"In dbSNP\">"),
        "{flag}"
    );

    // Under 4.1 the Source tag is parsed and then dropped, with no error at all.
    let v41 = row(&text, "hline\tinfo-source-41\t").expect("the golden carries it");
    assert!(!v41.contains("Source"), "{v41}");

    // A contig written with its fields the other way round renders with ID first.
    let contig = row(&text, "hline\tcontig-reordered\t").expect("the golden carries it");
    assert!(contig.ends_with("contig=<ID=chr1,length=1000>"), "{contig}");
}
