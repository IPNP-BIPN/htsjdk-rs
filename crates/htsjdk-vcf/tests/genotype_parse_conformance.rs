//! Conformance for the genotype columns, against `AbstractVCFCodec.createGenotypeMap`.
//!
//! Goldens from `tools/vcf-conformance/VcfGenotypeParseDump.java` in the pinned oracle.
//!
//! Four rows carry the weight:
//!
//! ```text
//! gt  backslash-separator  A,T          a backslash separates alleles like a slash
//! gt  ad-not-a-number      AD none      a malformed AD disappears
//! gterror  dp-not-a-number  NumberFormatException   a malformed DP aborts the record
//! gt  gq-negative-half     GQ -2 / none  Math.round is floor(x + 0.5), and -1 means "absent"
//! ```

use std::io::Read;

use htsjdk_vcf::genotype_parse::{parse_genotypes, GenotypeContext};
use htsjdk_vcf::header_lines::parse_meta_line;
use htsjdk_vcf::header_parse::{read_header_frame, VcfVersion};
use htsjdk_vcf::record_parse::{decode_line, split_condensed};
use htsjdk_vcf::variant::Genotype;
use htsjdk_vcf::{HeaderLine, VcfHeader};

const COMMON_META: &str = "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
    ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
    ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Quality\">\n\
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Depths\">\n\
    ##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Likelihoods\">\n\
    ##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Genotype filter\">\n\
    ##FORMAT=<ID=XX,Number=1,Type=String,Description=\"Anything\">\n\
    ##FILTER=<ID=LowQual,Description=\"Low quality\">\n";

const COLUMNS_TWO: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA1\tNA2\n";
const COLUMNS_ONE: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA1\n";

/// Which header a case uses, matching the dump's three.
#[derive(Clone, Copy, PartialEq)]
enum Which {
    Two,
    V40,
    V41,
}

/// Label, the genotype block, and which header. Same order as the dump.
const CASES: &[(&str, &str, Which)] = &[
    ("plain", "GT\t0/1\t1/1", Which::Two),
    ("phased", "GT\t0|1\t1|1", Which::Two),
    ("haploid", "GT\t0\t1", Which::Two),
    ("no-call", "GT\t./.\t0/1", Which::Two),
    (
        "full-format",
        "GT:GQ:DP:AD:PL\t0/1:99:30:10,20:100,0,200\t1/1:50:10:0,10:255,30,0",
        Which::Two,
    ),
    ("backslash-separator", "GT\t0\\1\t1\\1", Which::Two),
    ("doubled-separator", "GT\t0//1\t1/1", Which::Two),
    ("leading-separator", "GT\t/0/1\t1/1", Which::Two),
    ("trailing-separator", "GT\t0/1/\t1/1", Which::Two),
    ("missing-gq", "GT:GQ\t0/1:.\t1/1:50", Which::Two),
    (
        "short-value-list",
        "GT:GQ:DP\t0/1:99\t1/1:50:10",
        Which::Two,
    ),
    ("gq-minus-one", "GT:GQ\t0/1:-1\t1/1:50", Which::Two),
    ("gq-half", "GT:GQ\t0/1:2.5\t1/1:3.5", Which::Two),
    ("gq-negative-half", "GT:GQ\t0/1:-2.5\t1/1:-1.5", Which::Two),
    ("ad-not-a-number", "GT:AD\t0/1:1,x\t1/1:0,10", Which::Two),
    (
        "pl-not-a-number",
        "GT:PL\t0/1:1,x,3\t1/1:0,10,20",
        Which::Two,
    ),
    ("dp-not-a-number", "GT:DP\t0/1:x\t1/1:10", Which::Two),
    ("ft-pass", "GT:FT\t0/1:PASS\t1/1:LowQual", Which::Two),
    ("ft-missing", "GT:FT\t0/1:.\t1/1:PASS", Which::Two),
    ("ft-two", "GT:FT\t0/1:LowQual;q10\t1/1:PASS", Which::Two),
    ("extended-key", "GT:XX\t0/1:hello\t1/1:world", Which::Two),
    ("no-gt", "GQ\t99\t50", Which::Two),
    ("gt-not-first", "GQ:GT\t99:0/1\t50:1/1", Which::Two),
    ("too-many-values", "GT\t0/1:99\t1/1", Which::Two),
    ("too-few-columns", "GT\t0/1", Which::Two),
    ("allele-index-out-of-range", "GT\t0/2\t1/1", Which::Two),
    ("allele-index-not-a-number", "GT\tx/1\t1/1", Which::Two),
    ("no-gt-v40", "GQ\t99", Which::V40),
    ("no-gt-v41", "GQ\t99", Which::V41),
];

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/vcf_genotype_parse.txt.gz");
    let file = std::fs::File::open(&path).expect("corpus");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("corpus is gzip");
    text
}

fn header_text(which: Which) -> String {
    match which {
        Which::Two => format!("##fileformat=VCFv4.2\n{COMMON_META}{COLUMNS_TWO}"),
        Which::V40 => format!("##fileformat=VCFv4.0\n{COMMON_META}{COLUMNS_ONE}"),
        Which::V41 => format!("##fileformat=VCFv4.1\n{COMMON_META}{COLUMNS_ONE}"),
    }
}

/// The header, the version, and the codec's line counter after reading it.
fn header(which: Which) -> (VcfHeader, VcfVersion, usize) {
    let text = header_text(which);
    let frame = read_header_frame(&text).expect("the fixture header parses");
    let mut header = VcfHeader::new();
    header.samples = frame.samples.clone();
    let mut contigs = 0;
    for line in &frame.meta_lines {
        if let Ok(parsed) = parse_meta_line(line, frame.version, contigs) {
            if matches!(parsed, HeaderLine::Contig { .. }) {
                contigs += 1;
            }
            header.lines.push(parsed);
        }
    }
    (header, frame.version, frame.meta_lines.len() + 1)
}

/// The dump's rendering of one genotype.
fn render(label: &str, genotype: &Genotype) -> String {
    let alleles = genotype
        .alleles
        .iter()
        .map(|allele| allele.display_string())
        .collect::<Vec<_>>()
        .join(",");

    // `-1` is `Genotype`'s own sentinel for "no GQ", so a GQ that rounds to -1 reads as absent.
    let gq = match genotype.gq {
        Some(value) if value != -1 => value.to_string(),
        _ => "none".to_string(),
    };
    let dp = genotype
        .dp
        .map_or_else(|| "none".to_string(), |value| value.to_string());
    let ad = ints(genotype.ad.as_deref());
    let pl = ints(genotype.pl.as_deref());
    // An empty filter list is `PASS`, which the accessor reports as no filters at all.
    let filters = match genotype.filters.as_deref() {
        None | Some("") => "unfiltered".to_string(),
        Some(text) => text.to_string(),
    };
    let mut extended: Vec<String> = genotype
        .extended
        .iter()
        .map(|(key, value)| format!("{key}={}", value.format().unwrap_or_default()))
        .collect();
    extended.sort();

    format!(
        "gt\t{label}\t{}\t{alleles}\t{}\t{gq}\t{dp}\t{ad}\t{pl}\t{filters}\t{}",
        genotype.sample_name,
        genotype.phased,
        extended.join(";")
    )
}

fn ints(values: Option<&[i32]>) -> String {
    match values {
        None => "none".to_string(),
        Some(list) => list
            .iter()
            .map(i32::to_string)
            .collect::<Vec<_>>()
            .join(","),
    }
}

#[test]
fn every_genotype_decodes_as_the_reference_decodes_it() {
    let text = corpus();
    let mut rows = 0;

    for (label, block, which) in CASES {
        let (header, version, line_no) = header(*which);
        let line = format!("chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\t{block}");

        let outcome = decode_line(&line, &header, line_no, version).and_then(|decoded| {
            let record = decoded.expect("a data line");
            let site_parts = split_condensed(&line, '\t', 9, true);
            parse_genotypes(
                record.genotype_block.as_deref().unwrap_or(""),
                &record.variant.alleles,
                &GenotypeContext {
                    site_parts: &site_parts,
                    header: &header,
                    version,
                    contig: &record.variant.contig,
                    pos: record.variant.start,
                    line_number: line_no + 1,
                },
            )
        });

        match outcome {
            Ok(genotypes) => {
                for genotype in &genotypes {
                    let ours = render(label, genotype);
                    let expected = text
                        .lines()
                        .filter(|line| line.starts_with(&format!("gt\t{label}\t")))
                        .find(|line| line.split('\t').nth(2) == Some(genotype.sample_name.as_str()))
                        .unwrap_or_else(|| {
                            panic!(
                                "{label}/{}: no `gt` row in the golden",
                                genotype.sample_name
                            )
                        });
                    assert_eq!(ours, expected, "{label}/{}", genotype.sample_name);
                    rows += 1;
                }
            }
            Err(error) => {
                let expected = text
                    .lines()
                    .find_map(|line| line.strip_prefix(&format!("gterror\t{label}\t")))
                    .unwrap_or_else(|| {
                        panic!("{label}: the port refused, the golden has no `gterror` row")
                    });
                let message = error.message().replace('\t', " ");
                assert_eq!(format!("{}\t{message}", error.class()), expected, "{label}");
                rows += 1;
            }
        }
    }

    let in_golden = text
        .lines()
        .filter(|line| line.starts_with("gt\t") || line.starts_with("gterror\t"))
        .count();
    assert_eq!(
        in_golden, rows,
        "the golden and the test disagree on how many rows there are"
    );
    println!("{rows} genotype rows identical");
}

/// The four rows that separate this suite from a plausible reimplementation.
#[test]
fn the_separators_and_the_silent_drops_are_the_reference() {
    let text = corpus();
    let row = |prefix: &str| {
        text.lines()
            .find(|line| line.starts_with(prefix))
            .unwrap_or_else(|| panic!("the golden carries {prefix}"))
            .to_string()
    };

    // A backslash separates alleles, so this is two alleles and not one token.
    assert!(row("gt\tbackslash-separator\tNA1\t").contains("\tA,T\t"));
    // A malformed AD leaves no AD at all rather than refusing the record.
    let ad = row("gt\tad-not-a-number\tNA1\t");
    assert_eq!(ad.split('\t').nth(7), Some("none"), "{ad}");
    // A malformed DP under the very same rules aborts the whole record.
    assert!(row("gterror\tdp-not-a-number\t").contains("java.lang.NumberFormatException"));
    // Math.round is floor(x + 0.5): -2.5 becomes -2, and -1.5 becomes -1, which reads as absent.
    let gq = row("gt\tgq-negative-half\tNA1\t");
    assert_eq!(gq.split('\t').nth(5), Some("-2"), "{gq}");
    let absent = row("gt\tgq-negative-half\tNA2\t");
    assert_eq!(absent.split('\t').nth(5), Some("none"), "{absent}");
}
