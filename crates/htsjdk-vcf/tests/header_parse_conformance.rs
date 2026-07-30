//! Conformance for *reading* a VCF header, against `VCFHeaderLineTranslator` and `VCFCodec`.
//!
//! Goldens from `tools/vcf-conformance/VcfHeaderParseDump.java` in the pinned oracle.
//!
//! Everything else in this crate writes VCF. This is the first half of a reader, and the rows that
//! justify it are the ones a port written from the VCF specification gets wrong:
//!
//! ```text
//! line  angle-open-inside     ID=A|Foo=xy      the unquoted '<' is dropped, not kept
//! line  no-trailing-angle     ID=A             the last field is lost when the line has no '>'
//! line  quote-mid-value       ID=A|Desc=ab,cd  a mid-value quote makes the comma stop separating
//! line  repeated-key          ID=B|Number=1    last value, first position
//! ```
//!
//! The refusal messages are compared too, because they are what distinguish two failures of the
//! same file: a missing version line and a missing `#CHROM` line are both "malformed header" and a
//! caller can act on which.

use std::io::Read;

use htsjdk_vcf::header_parse::{parse_structured_value, read_header_frame};

/// The tags an INFO line must carry, in order, and the two that may follow them.
const INFO_TAGS: [&str; 4] = ["ID", "Number", "Type", "Description"];
const INFO_RECOMMENDED: [&str; 2] = ["Source", "Version"];

/// The same cases as the dump's, in the same order: label, value, and whether the INFO tag order is
/// enforced. A label in one and not the other is a failure, so neither side can quietly drop a case.
const LINES: &[(&str, &str, bool)] = &[
    (
        "plain",
        "<ID=DP,Number=1,Type=Integer,Description=\"Approximate depth\">",
        false,
    ),
    ("comma-in-quotes", "<ID=A,Description=\"one, two\">", false),
    (
        "escaped-quote",
        "<ID=A,Description=\"say \\\"hi\\\"\">",
        false,
    ),
    ("backslash-n", "<ID=A,Description=\"path\\nnext\">", false),
    ("double-backslash", "<ID=A,Description=\"a\\\\b\">", false),
    ("unclosed-quote", "<ID=A,Description=\"unterminated>", false),
    ("angle-open-inside", "<ID=A,Foo=x<y>", false),
    ("angle-close-inside", "<ID=A,Foo=x>y>", false),
    ("no-trailing-angle", "<ID=A,Number=1", false),
    ("spaces", "< ID = A , Number = 1 >", false),
    ("empty-value", "<ID=,Number=1>", false),
    ("quote-mid-value", "<ID=A,Desc=a\"b,c\"d>", false),
    ("repeated-key", "<ID=A,Number=1,ID=B>", false),
    ("empty-brackets", "<>", false),
    ("no-brackets", "ID=A,Number=1", false),
    (
        "info-ok",
        "<ID=DP,Number=1,Type=Integer,Description=\"d\">",
        true,
    ),
    (
        "info-wrong-order",
        "<Number=1,ID=DP,Type=Integer,Description=\"d\">",
        true,
    ),
    (
        "info-unexpected-tag",
        "<ID=DP,Foo=1,Type=Integer,Description=\"d\">",
        true,
    ),
    (
        "info-recommended-early",
        "<ID=DP,Source=\"s\",Number=1,Type=Integer,Description=\"d\">",
        true,
    ),
    (
        "info-recommended-late",
        "<ID=DP,Number=1,Type=Integer,Description=\"d\",Source=\"s\">",
        true,
    ),
    ("info-no-tags", "<>", true),
    (
        "info-extra-trailing",
        "<ID=DP,Number=1,Type=Integer,Description=\"d\",Whatever=1>",
        true,
    ),
];

/// The header texts, again in the dump's order.
const FRAMES: &[(&str, &str)] = &[
    (
        "minimal",
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    ),
    (
        "with-samples",
        "##fileformat=VCFv4.2\n##source=probe\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA1\tNA2\n",
    ),
    (
        "v40",
        "##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    ),
    (
        "v43",
        "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    ),
    (
        "v33",
        "##fileformat=VCFv3.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    ),
    (
        "v44",
        "##fileformat=VCFv4.4\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    ),
    (
        "unknown-version",
        "##fileformat=VCFv9.9\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    ),
    (
        "version-line-extra-equals",
        "##fileformat=VCFv4.2=x\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    ),
    (
        "no-version",
        "##source=probe\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    ),
    ("no-column-line", "##fileformat=VCFv4.2\n"),
    (
        "data-before-column-line",
        "##fileformat=VCFv4.2\nchr1\t1\t.\tA\tT\t.\t.\t.\n",
    ),
    (
        "too-few-columns",
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\n",
    ),
    (
        "misspelled-column",
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTR\tINFO\n",
    ),
    (
        "swapped-columns",
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tINFO\tFILTER\n",
    ),
    (
        "format-without-samples",
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n",
    ),
    (
        "ninth-column-not-format",
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tNA1\n",
    ),
    (
        "repeated-sample",
        "##fileformat=VCFv4.2\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA1\tNA1\tNA2\n",
    ),
];

/// `TribbleException.InvalidHeader`, which is the only exception either layer raises.
const EXCEPTION: &str = "htsjdk.tribble.TribbleException$InvalidHeader";

fn corpus() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/vcf_header_parse.txt.gz");
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

/// The dump flattens tabs and newlines out of a message so a row stays one row.
fn one_line(message: &str) -> String {
    message.replace(['\n', '\t'], " ")
}

#[test]
fn every_structured_line_parses_as_the_reference_parses_it() {
    let text = corpus();

    for (label, value, validated) in LINES {
        let expected_tags: Option<&[&str]> = if *validated { Some(&INFO_TAGS) } else { None };
        let outcome = parse_structured_value(value, expected_tags, &INFO_RECOMMENDED);

        match outcome {
            Ok(pairs) => {
                let ours = pairs
                    .iter()
                    .map(|(key, value)| format!("{key}={value}"))
                    .collect::<Vec<_>>()
                    .join("|");
                let expected = row(&text, &format!("line\t{label}\t")).unwrap_or_else(|| {
                    panic!("{label}: the port parsed the line, the golden has no `line` row")
                });
                assert_eq!(ours, expected, "{label}");
            }
            Err(error) => {
                let expected = row(&text, &format!("lineerror\t{label}\t")).unwrap_or_else(|| {
                    panic!("{label}: the port refused, the golden has no `lineerror` row")
                });
                assert_eq!(
                    format!("{EXCEPTION}\t{}", one_line(&error.message())),
                    expected,
                    "{label}"
                );
            }
        }
    }

    let in_golden = text
        .lines()
        .filter(|line| line.starts_with("line\t") || line.starts_with("lineerror\t"))
        .count();
    assert_eq!(
        in_golden,
        LINES.len(),
        "the golden and the test disagree on the case list"
    );
    println!("{} structured lines identical", LINES.len());
}

#[test]
fn every_header_frame_reads_as_the_reference_reads_it() {
    let text = corpus();

    for (label, header) in FRAMES {
        match read_header_frame(header) {
            Ok(frame) => {
                let ours = format!(
                    "{}\t{}\t{}",
                    frame.version.version_string(),
                    frame.meta_keys().join(","),
                    frame.samples.join(",")
                );
                let expected = row(&text, &format!("frame\t{label}\t")).unwrap_or_else(|| {
                    panic!("{label}: the port read the header, the golden has no `frame` row")
                });
                assert_eq!(ours, expected, "{label}");
            }
            Err(error) => {
                let expected = row(&text, &format!("frameerror\t{label}\t")).unwrap_or_else(|| {
                    panic!("{label}: the port refused, the golden has no `frameerror` row")
                });
                assert_eq!(
                    format!("{EXCEPTION}\t{}", one_line(&error.message())),
                    expected,
                    "{label}"
                );
            }
        }
    }

    let in_golden = text
        .lines()
        .filter(|line| line.starts_with("frame\t") || line.starts_with("frameerror\t"))
        .count();
    assert_eq!(
        in_golden,
        FRAMES.len(),
        "the golden and the test disagree on the case list"
    );
    println!("{} header frames identical", FRAMES.len());
}
