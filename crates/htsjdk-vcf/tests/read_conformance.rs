//! Conformance for reading a whole VCF file, against `VCFCodec` driven the way a feature reader
//! drives it.
//!
//! Goldens from `tools/vcf-conformance/VcfReadDump.java` in the pinned oracle.
//!
//! The rows that justify the suite are the ones about state carried between lines, which no
//! single-line suite can see:
//!
//! ```text
//! err  short-line-first             ... Line 12: there aren't enough columns ...
//! pct  % 1                          ... at approximately line number 13: ...
//! ```
//!
//! Two files with identical twelve-line headers and one data line each, and the same position
//! reported as 12 by one check and 13 by the other.
//!
//! And the ones about the header not being the header:
//!
//! ```text
//! hdr  fileformat-4-0    VCFv4.0  none  fileformat=VCFv4.2 | source=mine
//! hdr  repair-wrong-type VCFv4.2  none  fileformat=VCFv4.2 | INFO=<ID=DP,...,Description="Approximate read depth; ...">
//! ```
//!
//! The first file declared v4.0 and the second declared its own description; neither survives.

use std::io::Read;

use htsjdk_vcf::header::HeaderLine;
use htsjdk_vcf::header_parse::VcfVersion;
use htsjdk_vcf::reader::{read_vcf, VcfFile};
use htsjdk_vcf::text_transformer::TextTransformer;
use htsjdk_vcf::variant::{Value, VariantContext};
use htsjdk_vcf::vcf_file::{reject_vcf_v43_headers, write_vcf};

/// `VcfReadDump.META`, the metadata every whole-file case shares.
const META: &str = concat!(
    "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n",
    "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Frequency\">\n",
    "##INFO=<ID=NOTE,Number=1,Type=String,Description=\"A note\">\n",
    "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"In dbSNP\">\n",
    "##FILTER=<ID=LowQual,Description=\"Low quality\">\n",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
    "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n",
    "##FORMAT=<ID=SB,Number=1,Type=String,Description=\"A string\">\n",
    "##contig=<ID=chr1,length=100000>\n",
    "##contig=<ID=chr2,length=100000>\n",
);

const COLUMNS_SITES: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
const COLUMNS_SAMPLES: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA1\tNA2\n";

fn sites(version: &str, body: &str) -> String {
    format!("##fileformat={version}\n{META}{COLUMNS_SITES}{body}")
}

fn samples(version: &str, body: &str) -> String {
    format!("##fileformat={version}\n{META}{COLUMNS_SAMPLES}{body}")
}

/// Every `read` case of the dump, in its order.
fn read_cases() -> Vec<(&'static str, String)> {
    let encoded_body = concat!(
        "chr1\t100\t.\tA\tT\t50\tPASS\tNOTE=a%3Ab%3Bc\n",
        "chr1\t200\t.\tC\tG\t50\tPASS\tNOTE=100%25\n"
    );
    let encoded_genotypes = "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\tGT:SB\t0/1:a%3Ab\t1/1:c%3Ad\n";
    let encoded_flag = "chr1\t100\t.\tA\tT\t50\tPASS\tDB=%30\n";

    vec![
        (
            "sites-only",
            sites(
                "VCFv4.2",
                concat!(
                    "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\n",
                    "chr1\t200\trs1\tC\tG\t.\t.\tDP=20;AF=0.5\n",
                    "chr2\t1\t.\tGG\tG\t30\tLowQual\tDB\n"
                ),
            ),
        ),
        (
            "genotyped",
            samples(
                "VCFv4.2",
                concat!(
                    "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\tGT:GQ\t0/1:30\t1|1:40\n",
                    "chr1\t200\t.\tC\tG\t50\tPASS\tDP=20\tGT\t./.\t0/0\n"
                ),
            ),
        ),
        ("header-only", sites("VCFv4.2", "")),
        (
            "comment-in-body",
            sites(
                "VCFv4.2",
                concat!(
                    "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\n",
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
                    "chr1\t200\t.\tC\tG\t50\tPASS\tDP=20\n"
                ),
            ),
        ),
        (
            "hash-comment-in-body",
            sites(
                "VCFv4.2",
                concat!(
                    "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\n",
                    "# a comment\n",
                    "chr1\t200\t.\tC\tG\t50\tPASS\tDP=20\n"
                ),
            ),
        ),
        (
            "blank-line-in-body",
            sites(
                "VCFv4.2",
                concat!(
                    "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\n",
                    "\n",
                    "chr1\t200\t.\tC\tG\t50\tPASS\tDP=20\n"
                ),
            ),
        ),
        (
            "short-line-reports-its-number",
            sites(
                "VCFv4.2",
                concat!(
                    "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\n",
                    "chr1\t200\t.\tC\tG\t50\tPASS\tDP=20\n",
                    "chr1\t300\t.\tC\tG\t50\tPASS\tDP=30\n",
                    "chr1\t400\t.\tC\tG\t50\tPASS\n"
                ),
            ),
        ),
        (
            "bad-qual-reports-its-number",
            sites(
                "VCFv4.2",
                concat!(
                    "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\n",
                    "chr1\t200\t.\tC\tG\t50\tPASS\tDP=20\n",
                    "chr1\t300\t.\tC\tG\t50\tPASS\tDP=30\n",
                    "chr1\t400\t.\tC\tG\tx\tPASS\tDP=40\n"
                ),
            ),
        ),
        (
            "short-line-first",
            sites("VCFv4.2", "chr1\t100\t.\tA\tT\t50\tPASS\n"),
        ),
        (
            "crlf",
            sites(
                "VCFv4.2",
                concat!(
                    "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\r\n",
                    "chr1\t200\t.\tC\tG\t50\tPASS\tNOTE=x\r\n"
                ),
            ),
        ),
        (
            "no-trailing-newline",
            sites("VCFv4.2", "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10"),
        ),
        (
            "no-fileformat",
            format!("{META}{COLUMNS_SITES}chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\n"),
        ),
        (
            "no-chrom-line",
            format!("##fileformat=VCFv4.2\n{META}chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\n"),
        ),
        ("empty-file", String::new()),
        (
            "fileformat-with-two-equals",
            format!("##fileformat=VCFv4.2=x\n{META}{COLUMNS_SITES}"),
        ),
        (
            "unsupported-version",
            format!("##fileformat=VCFv3.3\n{META}{COLUMNS_SITES}"),
        ),
        ("percent-under-4-2", sites("VCFv4.2", encoded_body)),
        ("percent-under-4-3", sites("VCFv4.3", encoded_body)),
        (
            "percent-genotype-4-2",
            samples("VCFv4.2", encoded_genotypes),
        ),
        (
            "percent-genotype-4-3",
            samples("VCFv4.3", encoded_genotypes),
        ),
        ("percent-flag-zero-4-2", sites("VCFv4.2", encoded_flag)),
        ("percent-flag-zero-4-3", sites("VCFv4.3", encoded_flag)),
    ]
}

/// The `read` cases the dump emits *after* its `hdr` block, kept in that order because the
/// comparison is line by line.
fn read_cases_after_headers() -> Vec<(&'static str, String)> {
    vec![
        (
            "duplicate-samples",
            format!(
                "##fileformat=VCFv4.2\n{META}\
                 #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA1\tNA1\n\
                 chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\tGT\t0/1\t1/1\n"
            ),
        ),
        (
            "unsorted-and-undeclared-contig",
            sites(
                "VCFv4.2",
                concat!(
                    "chr2\t500\t.\tA\tT\t50\tPASS\tDP=10\n",
                    "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\n",
                    "chrX\t9\t.\tA\tT\t50\tPASS\tDP=10\n"
                ),
            ),
        ),
    ]
}

/// Every `hdr` case of the dump, in its order.
fn header_cases() -> Vec<(&'static str, String)> {
    let one =
        |version: &str, line: &str| format!("##fileformat={version}\n{line}\n{COLUMNS_SITES}");
    vec![
        ("repair-none", one("VCFv4.2", "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth; some reads may have been filtered\">")),
        ("repair-wrong-type", one("VCFv4.2", "##INFO=<ID=DP,Number=1,Type=Float,Description=\"Depth\">")),
        ("repair-wrong-count", one("VCFv4.2", "##INFO=<ID=DP,Number=2,Type=Integer,Description=\"Depth\">")),
        ("repair-wrong-count-type", one("VCFv4.2", "##INFO=<ID=DP,Number=A,Type=Integer,Description=\"Depth\">")),
        ("repair-wrong-description-only", one("VCFv4.2", "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"my own depth\">")),
        ("repair-format-gq", one("VCFv4.2", "##FORMAT=<ID=GQ,Number=1,Type=String,Description=\"Genotype Quality\">")),
        ("repair-not-a-standard-id", one("VCFv4.2", "##INFO=<ID=XX,Number=2,Type=Float,Description=\"Mine\">")),
        ("repair-under-4-3", one("VCFv4.3", "##INFO=<ID=DP,Number=1,Type=Float,Description=\"Depth\">")),
        ("fileformat-4-0", one("VCFv4.0", "##source=mine")),
        ("fileformat-4-1", one("VCFv4.1", "##source=mine")),
        ("fileformat-4-2", one("VCFv4.2", "##source=mine")),
        ("fileformat-4-3", one("VCFv4.3", "##source=mine")),
        ("fileformat-twice", format!("##fileformat=VCFv4.2\n##source=mine\n##fileformat=VCFv4.1\n{COLUMNS_SITES}")),
    ]
}

/// The raw texts of the `pct` rows, in the dump's order.
const PERCENT_CASES: &[&str] = &[
    "plain", "%41", "%3D%41", "a%3Ab", "%", "x%", "%4", "x%4", "%4G", "%G4", "%+1", "%-1", "%%41",
    "%09", "%00", "%7e", "%7E", "100%25", "%zz", "% 1",
];

/// The `trip` cases, in the dump's order.
fn trip_cases() -> Vec<(&'static str, String)> {
    vec![
        (
            "trip-sites",
            sites(
                "VCFv4.2",
                concat!(
                    "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\n",
                    "chr1\t200\trs1\tC\tG\t.\t.\tDP=20;AF=0.5\n"
                ),
            ),
        ),
        (
            "trip-genotyped",
            samples(
                "VCFv4.2",
                "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\tGT:GQ\t0/1:30\t1|1:40\n",
            ),
        ),
        (
            "trip-4-3",
            sites("VCFv4.3", "chr1\t100\t.\tA\tT\t50\tPASS\tNOTE=a%3Ab\n"),
        ),
    ]
}

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/vcf_read.txt.gz");
    let file = std::fs::File::open(&path).expect("corpus");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("corpus is gzip");
    text
}

/// Every golden line, so a row the port never produces is caught as well as one it gets wrong.
fn golden_rows(text: &str) -> Vec<&str> {
    text.lines().filter(|line| !line.starts_with('#')).collect()
}

fn version_text(version: Option<VcfVersion>) -> String {
    match version {
        None => "none".to_string(),
        Some(version) => version.version_string().to_string(),
    }
}

/// The dump's escaping: backslash, tab, newline, carriage return, then `\uXXXX` for anything
/// outside printable ASCII.
fn escape(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for c in text.chars() {
        match c {
            '\\' => out.push_str("\\\\"),
            '\t' => out.push_str("\\t"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            c if (c as u32) < 0x20 || (c as u32) > 0x7e => {
                out.push_str(&format!("\\u{:04x}", c as u32));
            }
            c => out.push(c),
        }
    }
    out
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

fn sorted_attributes(attributes: &[(String, Value)]) -> String {
    let mut rendered: Vec<String> = attributes
        .iter()
        .map(|(key, value)| format!("{key}={}", java_value(value)))
        .collect();
    rendered.sort();
    if rendered.is_empty() {
        "-".to_string()
    } else {
        rendered.join(";")
    }
}

fn render_record(label: &str, index: usize, variant: &VariantContext) -> String {
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

    let genotypes: Vec<String> = variant
        .genotypes
        .iter()
        .map(|g| {
            let separator = if g.phased { "|" } else { "/" };
            let called = g
                .alleles
                .iter()
                .map(|a| a.display_string())
                .collect::<Vec<_>>()
                .join(separator);
            format!(
                "{}:{called}:GQ={}:{}",
                g.sample_name,
                // `Genotype.getGQ()` returns -1 when there is none.
                g.gq.unwrap_or(-1),
                sorted_attributes(&g.extended)
            )
        })
        .collect();

    format!(
        "rec\t{label}\t{index}\t{}\t{}\t{}\t{}\t{alleles}\t{qual}\t{filters}\t{}\t{}",
        variant.contig,
        variant.start,
        variant.stop,
        variant.id,
        sorted_attributes(&variant.attributes),
        if genotypes.is_empty() {
            "-".to_string()
        } else {
            genotypes.join(" ")
        }
    )
}

fn render_read(label: &str, text: &str) -> Vec<String> {
    let mut rows = Vec::new();
    match read_vcf(text) {
        Ok(VcfFile {
            codec_version,
            header_version,
            ref records,
            ref skipped,
            ref header,
            ..
        }) => {
            for index in skipped {
                rows.push(format!("null\t{label}\t{index}"));
            }
            rows.push(format!(
                "file\t{label}\t{}\t{}\t{}\t{}",
                codec_version.version_string(),
                version_text(header_version),
                if header.samples.is_empty() {
                    "-".to_string()
                } else {
                    header.samples.join(",")
                },
                records.len()
            ));
            for (index, record) in records.iter().enumerate() {
                rows.push(render_record(label, index, record));
            }
        }
        Err(failure) => {
            rows.push(format!(
                "file\t{label}\taborted\taborted\t-\t{}",
                failure.records.len()
            ));
            for (index, record) in failure.records.iter().enumerate() {
                rows.push(render_record(label, index, record));
            }
            rows.push(format!(
                "err\t{label}\t{}\t{}",
                failure.error.class(),
                escape(&failure.error.message())
            ));
        }
    }
    rows
}

/// One INFO value read under both versions, which is how the dump measures the transformer: the
/// codec chooses it, and choosing it is half the behaviour.
fn note_value(version: &str, raw: &str) -> String {
    let text = sites(
        version,
        &format!("chr1\t100\t.\tA\tT\t50\tPASS\tNOTE={raw}\n"),
    );
    match read_vcf(&text) {
        Ok(file) => java_value(&file.records[0].attributes[0].1),
        Err(failure) => format!(
            "{}: {}",
            // The dump prints `getSimpleName()`, so the package and the outer class go.
            failure.error.class().rsplit('.').next().unwrap(),
            escape(&failure.error.message())
        ),
    }
}

fn render_trip(label: &str, text: &str) -> String {
    let file = match read_vcf(text) {
        Ok(file) => file,
        Err(failure) => {
            return format!(
                "err\t{label}\t{}\t{}",
                failure.error.class(),
                escape(&failure.error.message())
            )
        }
    };
    // The writer refuses a 4.3 header, and it can only do that because the header still knows it
    // is 4.3 while a 4.2 one no longer knows what it is.
    if let Err(message) = reject_vcf_v43_headers(file.header_version) {
        return format!(
            "err\t{label}\tjava.lang.IllegalArgumentException\t{}",
            escape(&message)
        );
    }
    let rewritten = write_vcf(&file.header, &file.records).expect("the fixture writes");
    let mut at: i64 = -1;
    for (index, (a, b)) in text.bytes().zip(rewritten.bytes()).enumerate() {
        if a != b {
            at = index as i64;
            break;
        }
    }
    if at == -1 && text.len() != rewritten.len() {
        at = text.len().min(rewritten.len()) as i64;
    }
    format!(
        "trip\t{label}\t{}\t{}\t{}",
        if at == -1 { "same" } else { "differs" },
        if at == -1 {
            "-".to_string()
        } else {
            at.to_string()
        },
        escape(&rewritten)
    )
}

#[test]
fn every_file_reads_as_the_reference_reads_it() {
    let text = corpus();
    let golden = golden_rows(&text);

    let mut produced: Vec<String> = Vec::new();
    for (label, file) in read_cases() {
        produced.extend(render_read(label, &file));
    }
    for (label, file) in header_cases() {
        produced.push(render_header(label, &file));
    }
    for (label, file) in read_cases_after_headers() {
        produced.extend(render_read(label, &file));
    }
    for raw in PERCENT_CASES {
        produced.push(format!(
            "pct\t{}\t{}\t{}",
            escape(raw),
            escape(&note_value("VCFv4.3", raw)),
            escape(&note_value("VCFv4.2", raw))
        ));
    }
    for (label, file) in trip_cases() {
        produced.push(render_trip(label, &file));
    }

    assert_eq!(
        produced.len(),
        golden.len(),
        "the port produced {} rows and the golden has {}",
        produced.len(),
        golden.len()
    );
    for (index, (mine, theirs)) in produced.iter().zip(golden.iter()).enumerate() {
        assert_eq!(mine, theirs, "row {index}");
    }
}

fn render_header(label: &str, text: &str) -> String {
    match read_vcf(text) {
        Ok(file) => {
            let lines: Vec<String> = file
                .meta_data_in_input_order()
                .iter()
                .map(HeaderLine::render)
                .collect();
            format!(
                "hdr\t{label}\t{}\t{}\t{}",
                file.codec_version.version_string(),
                version_text(file.header_version),
                lines.join(" | ")
            )
        }
        Err(failure) => format!(
            "err\t{label}\t{}\t{}",
            failure.error.class(),
            escape(&failure.error.message())
        ),
    }
}

/// The transformer applies to genotype values through a different call site than INFO values, and
/// it runs before any key is looked at, so it reaches the GT string itself.
#[test]
fn the_transformer_reaches_genotype_values_and_the_gt_string() {
    assert_eq!(
        TextTransformer::for_version(VcfVersion::Vcf4_3).decode("0%2F1"),
        "0/1"
    );
    assert_eq!(
        TextTransformer::for_version(VcfVersion::Vcf4_2).decode("0%2F1"),
        "0%2F1"
    );
}
