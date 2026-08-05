//! Conformance for *reading* one VCF data line, against `AbstractVCFCodec.decodeLine`.
//!
//! Goldens from `tools/vcf-conformance/VcfRecordParseDump.java` in the pinned oracle.
//!
//! The rows that justify the suite are the ones where the header, and not the line, decides what
//! the record contains:
//!
//! ```text
//! rec  flag-equals-zero   ...  (no attributes)   a declared Flag written =0 is dropped
//! rec  undeclared-zero    ...  XX=0              the same text, undeclared, is kept as a string
//! rec  bare-non-flag      ...  DP=.              a bare key typed Integer becomes "."
//! rec  undeclared-bare    ...  XX=true           a bare key with no declaration is a flag
//! ```
//!
//! And the two that would look like the same failure to a port with one error type:
//!
//! ```text
//! recerror  pos-not-a-number   htsjdk.tribble.TribbleException   ... at approximately line number 9
//! recerror  qual-not-a-number  java.lang.NumberFormatException   For input string: "x"
//! ```

use std::io::Read;

use htsjdk_vcf::header_lines::parse_meta_line;
use htsjdk_vcf::header_parse::{read_header_frame, VcfVersion};
use htsjdk_vcf::record_parse::decode_line;
use htsjdk_vcf::variant::{Value, VariantContext};
use htsjdk_vcf::VcfHeader;

/// The dump's header, with two samples, declaring the types the INFO parsing consults.
const HEADER: &str = "##fileformat=VCFv4.2\n\
    ##INFO=<ID=DB,Number=0,Type=Flag,Description=\"In dbSNP\">\n\
    ##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
    ##INFO=<ID=AF,Number=A,Type=Float,Description=\"Frequency\">\n\
    ##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n\
    ##FILTER=<ID=LowQual,Description=\"Low quality\">\n\
    ##contig=<ID=chr1,length=100000>\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA1\tNA2\n";

/// The same header with no samples, which changes the column count and the declared types.
const SITES_ONLY_HEADER: &str = "##fileformat=VCFv4.2\n\
    ##INFO=<ID=DB,Number=0,Type=Flag,Description=\"In dbSNP\">\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// Label, data line, and whether the sites-only header is used. Same order as the dump.
const CASES: &[(&str, &str, bool)] = &[
    (
        "snp",
        "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\tGT\t0/1\t1/1",
        false,
    ),
    (
        "multiallelic",
        "chr1\t100\trs1\tA\tT,C\t50\tPASS\tDP=10;AF=0.5,0.25\tGT\t0/1\t1/2",
        false,
    ),
    (
        "deletion",
        "chr1\t100\t.\tACGT\tA\t50\t.\tDP=10\tGT\t0/1\t0/0",
        false,
    ),
    (
        "no-qual",
        "chr1\t100\t.\tA\tT\t.\tPASS\tDP=10\tGT\t0/1\t1/1",
        false,
    ),
    (
        "qual-minus-one",
        "chr1\t100\t.\tA\tT\t-1\tPASS\tDP=10\tGT\t0/1\t1/1",
        false,
    ),
    (
        "qual-minus-one-point-zero",
        "chr1\t100\t.\tA\tT\t-1.0\tPASS\tDP=10\tGT\t0/1\t1/1",
        false,
    ),
    (
        "unfiltered",
        "chr1\t100\t.\tA\tT\t50\t.\tDP=10\tGT\t0/1\t1/1",
        false,
    ),
    (
        "two-filters",
        "chr1\t100\t.\tA\tT\t50\tLowQual;q10\tDP=10\tGT\t0/1\t1/1",
        false,
    ),
    (
        "end-key",
        "chr1\t100\t.\tA\t<DEL>\t50\tPASS\tEND=200\tGT\t0/1\t1/1",
        false,
    ),
    (
        "end-before-start",
        "chr1\t100\t.\tA\t<DEL>\t50\tPASS\tEND=50\tGT\t0/1\t1/1",
        false,
    ),
    (
        "flag-bare",
        "chr1\t100\t.\tA\tT\t50\tPASS\tDB\tGT\t0/1\t1/1",
        false,
    ),
    (
        "flag-equals-zero",
        "chr1\t100\t.\tA\tT\t50\tPASS\tDB=0\tGT\t0/1\t1/1",
        false,
    ),
    (
        "flag-equals-one",
        "chr1\t100\t.\tA\tT\t50\tPASS\tDB=1\tGT\t0/1\t1/1",
        false,
    ),
    (
        "bare-non-flag",
        "chr1\t100\t.\tA\tT\t50\tPASS\tDP\tGT\t0/1\t1/1",
        false,
    ),
    (
        "undeclared-bare",
        "chr1\t100\t.\tA\tT\t50\tPASS\tXX\tGT\t0/1\t1/1",
        false,
    ),
    (
        "undeclared-zero",
        "chr1\t100\t.\tA\tT\t50\tPASS\tXX=0\tGT\t0/1\t1/1",
        false,
    ),
    (
        "empty-value",
        "chr1\t100\t.\tA\tT\t50\tPASS\tDP=\tGT\t0/1\t1/1",
        false,
    ),
    (
        "empty-info",
        "chr1\t100\t.\tA\tT\t50\tPASS\t.\tGT\t0/1\t1/1",
        false,
    ),
    (
        "repeated-key",
        "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10;DP=20\tGT\t0/1\t1/1",
        false,
    ),
    (
        "lowercase-ref",
        "chr1\t100\t.\ta\tT\t50\tPASS\tDP=10\tGT\t0/1\t1/1",
        false,
    ),
    (
        "alt-missing",
        "chr1\t100\t.\tA\t.\t50\tPASS\tDP=10\tGT\t0/1\t1/1",
        false,
    ),
    (
        "alt-star",
        "chr1\t100\t.\tA\tT,*\t50\tPASS\tDP=10\tGT\t0/1\t1/1",
        false,
    ),
    (
        "symbolic-alt",
        "chr1\t100\t.\tA\t<NON_REF>\t50\tPASS\tDP=10\tGT\t0/1\t1/1",
        false,
    ),
    (
        "leading-tab",
        "\tchr1\t100\t.\tA\tT\t50\tPASS\tDP=10\tGT\t0/1\t1/1",
        false,
    ),
    (
        "header-line",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        false,
    ),
    (
        "pos-not-a-number",
        "chr1\tx\t.\tA\tT\t50\tPASS\tDP=10\tGT\t0/1\t1/1",
        false,
    ),
    (
        "empty-id",
        "chr1\t100\t\tA\tT\t50\tPASS\tDP=10\tGT\t0/1\t1/1",
        false,
    ),
    (
        "qual-not-a-number",
        "chr1\t100\t.\tA\tT\tx\tPASS\tDP=10\tGT\t0/1\t1/1",
        false,
    ),
    (
        "filter-zero",
        "chr1\t100\t.\tA\tT\t50\t0\tDP=10\tGT\t0/1\t1/1",
        false,
    ),
    (
        "info-with-space",
        "chr1\t100\t.\tA\tT\t50\tPASS\tDP=1 0\tGT\t0/1\t1/1",
        false,
    ),
    (
        "info-empty-string",
        "chr1\t100\t.\tA\tT\t50\tPASS\t\tGT\t0/1\t1/1",
        false,
    ),
    (
        "end-not-a-number",
        "chr1\t100\t.\tA\t<DEL>\t50\tPASS\tEND=x\tGT\t0/1\t1/1",
        false,
    ),
    (
        "ref-missing",
        "chr1\t100\t.\t.\tT\t50\tPASS\tDP=10\tGT\t0/1\t1/1",
        false,
    ),
    (
        "ref-symbolic",
        "chr1\t100\t.\t<DEL>\tT\t50\tPASS\tDP=10\tGT\t0/1\t1/1",
        false,
    ),
    (
        "ref-bad-base",
        "chr1\t100\t.\tQ\tT\t50\tPASS\tDP=10\tGT\t0/1\t1/1",
        false,
    ),
    (
        "alt-breakend",
        "chr1\t100\t.\tA\tA[chr2:200[\t50\tPASS\tDP=10\tGT\t0/1\t1/1",
        false,
    ),
    (
        "vcf3-deletion",
        "chr1\t100\t.\tA\tD2\t50\tPASS\tDP=10\tGT\t0/1\t1/1",
        false,
    ),
    ("too-few-columns", "chr1\t100\t.\tA\tT\t50\tPASS", false),
    (
        "sites-only-eight",
        "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10",
        true,
    ),
    (
        "sites-only-nine",
        "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\tGT",
        true,
    ),
    (
        "sites-only-undeclared-dp",
        "chr1\t100\t.\tA\tT\t50\tPASS\tDP",
        true,
    ),
];

fn corpus() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/vcf_record_parse.txt.gz");
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

/// Build the header the codec would have, through the two slices that came before this one.
fn header(text: &str) -> (VcfHeader, usize, VcfVersion) {
    let frame = read_header_frame(text).expect("the fixture header parses");
    let mut header = VcfHeader::new();
    header.samples = frame.samples.clone();
    let mut contigs = 0;
    for line in &frame.meta_lines {
        if let Ok(parsed) = parse_meta_line(line, frame.version, contigs) {
            if matches!(parsed, htsjdk_vcf::HeaderLine::Contig { .. }) {
                contigs += 1;
            }
            header.lines.push(parsed);
        }
    }
    // The codec's line counter after the header: every `##` line plus the `#CHROM` line. The
    // version travels with it, because it chooses the transformer every INFO value goes through.
    (header, frame.meta_lines.len() + 1, frame.version)
}

/// The dump's rendering of a decoded record, field for field.
fn render(variant: &VariantContext, samples: &[String]) -> String {
    let alleles = variant
        .alleles
        .iter()
        .map(|allele| {
            format!(
                "{}{}",
                allele.display_string(),
                if allele.is_reference() { "*" } else { "" }
            )
        })
        .collect::<Vec<_>>()
        .join(",");

    let qual = if variant.has_log10_p_error() {
        java_double(variant.log10_p_error)
    } else {
        "none".to_string()
    };

    let filters = match &variant.filters {
        None => "unfiltered".to_string(),
        Some(list) if list.is_empty() => "PASS".to_string(),
        Some(list) => {
            let mut sorted = list.clone();
            sorted.sort();
            sorted.join(",")
        }
    };

    // Sorted: upstream they live in a HashMap, so the order is the JDK's and not the file's.
    let mut attributes: Vec<String> = variant
        .attributes
        .iter()
        .map(|(key, value)| format!("{key}={}", java_value(value)))
        .collect();
    attributes.sort();

    format!(
        "{}\t{}\t{}\t{}\t{alleles}\t{qual}\t{filters}\t{}\t{}",
        variant.contig,
        variant.start,
        variant.stop,
        variant.id,
        attributes.join(";"),
        samples.join(",")
    )
}

/// `Double.toString`, for the two shapes QUAL takes here.
fn java_double(value: f64) -> String {
    if value == value.trunc() && value.abs() < 1e7 {
        format!("{value:.1}")
    } else {
        value.to_string()
    }
}

/// `String.valueOf(Object)` over the attribute values the decoder stores.
fn java_value(value: &Value) -> String {
    match value {
        Value::Str(text) => text.clone(),
        Value::Bool(flag) => flag.to_string(),
        Value::List(items) => format!(
            "[{}]",
            items.iter().map(java_value).collect::<Vec<_>>().join(",")
        ),
        other => other.format().unwrap_or_default(),
    }
}

#[test]
fn every_record_decodes_as_the_reference_decodes_it() {
    let text = corpus();
    let (with_samples, samples_line_no, samples_version) = header(HEADER);
    let (sites_only, sites_line_no, sites_version) = header(SITES_ONLY_HEADER);

    for (label, line, sites) in CASES {
        let (header, line_no, version) = if *sites {
            (&sites_only, sites_line_no, sites_version)
        } else {
            (&with_samples, samples_line_no, samples_version)
        };

        match decode_line(line, header, line_no, version) {
            Ok(None) => {
                assert!(
                    text.lines().any(|l| l == format!("recnull\t{label}")),
                    "{label}: the port skipped the line, the golden did not"
                );
            }
            Ok(Some(decoded)) => {
                let ours = render(&decoded.variant, &header.samples);
                let expected = row(&text, &format!("rec\t{label}\t")).unwrap_or_else(|| {
                    panic!("{label}: the port decoded the line, the golden has no `rec` row")
                });
                assert_eq!(ours, expected, "{label}");
            }
            Err(error) => {
                let expected = row(&text, &format!("recerror\t{label}\t")).unwrap_or_else(|| {
                    panic!("{label}: the port refused, the golden has no `recerror` row")
                });
                let message = error.message().replace('\t', " ");
                assert_eq!(format!("{}\t{message}", error.class()), expected, "{label}");
            }
        }
    }

    let in_golden = text
        .lines()
        .filter(|line| {
            line.starts_with("rec\t")
                || line.starts_with("recerror\t")
                || line.starts_with("recnull\t")
        })
        .count();
    assert_eq!(
        in_golden,
        CASES.len(),
        "the golden and the test disagree on the case list"
    );
    println!("{} records identical", CASES.len());
}

/// The four rows where the header, not the line, decides what the record contains. Named so a
/// regression carries its reason rather than being one line in a list of forty.
#[test]
fn the_header_decides_what_the_info_field_means() {
    let text = corpus();

    // A declared Flag written `=0` leaves no attribute at all.
    let dropped = row(&text, "rec\tflag-equals-zero\t").expect("the golden carries it");
    assert!(!dropped.contains("DB"), "{dropped}");

    // The same text under a key the header does not declare survives as a string.
    let kept = row(&text, "rec\tundeclared-zero\t").expect("the golden carries it");
    assert!(kept.contains("XX=0"), "{kept}");

    // A bare key the header types as Integer becomes the missing value, not a flag.
    let bare = row(&text, "rec\tbare-non-flag\t").expect("the golden carries it");
    assert!(bare.contains("DP=."), "{bare}");

    // A bare key with no declaration is a flag.
    let flag = row(&text, "rec\tundeclared-bare\t").expect("the golden carries it");
    assert!(flag.contains("XX=true"), "{flag}");
}
