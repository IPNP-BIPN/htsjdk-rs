//! Conformance for what the genotype columns are written from, against the oracle.
//!
//! Golden from `tools/vcf-conformance/VcfLazyGenotypesDump.java`.
//!
//! The rule the golden pins, and the reason a port gets it wrong by writing the obvious Rust:
//!
//! ```text
//! lazy  unsorted-format  copy   true    <- read and written untouched
//! line  unsorted-format  copy   ... GT:GQ:DP  0/1:60:10 ...   the FILE's order
//! lazy  unsorted-format  touch  false   <- one genotype was READ
//! line  unsorted-format  touch  ... GT:DP:GQ  0/1:10:60 ...   the computed order
//! ```
//!
//! A reader that decodes eagerly and throws the text away can only produce the second, and every
//! tool that copies records through then rewrites their FORMAT column.

use std::io::Read;

use htsjdk_vcf::encoder::VcfEncoder;
use htsjdk_vcf::reader::read_vcf;
use htsjdk_vcf::variant::Genotype;

fn golden() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/vcf_lazy_genotypes.txt.gz");
    let file = std::fs::File::open(&path).expect("golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("golden is gzip");
    text
}

fn rows(text: &str, kind: &str) -> Vec<Vec<String>> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').map(str::to_string).collect::<Vec<_>>())
        .filter(|parts| parts.first().map(String::as_str) == Some(kind))
        .collect()
}

fn escape(text: &str) -> String {
    text.replace('\\', "\\\\")
        .replace('\t', "\\t")
        .replace('\n', "\\n")
}

/// The dump's own header, which every case shares.
fn header_text() -> String {
    concat!(
        "##fileformat=VCFv4.2\n",
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
        "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Quality\">\n",
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n",
        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Depths\">\n",
        "##contig=<ID=chr1,length=100000>\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\n",
    )
    .to_string()
}

/// The four data lines the dump reads, by case name.
fn line_of(case: &str) -> &'static str {
    match case {
        "unsorted-format" => "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\tGT:GQ:DP\t0/1:60:10\t1/1:70:20",
        "sorted-format" => "chr1\t200\t.\tC\tG\t50\tPASS\tDP=10\tGT:DP:GQ\t0/1:10:60\t1/1:20:70",
        "trailing-missing" => "chr1\t300\t.\tG\tA\t50\tPASS\tDP=10\tGT:GQ:AD\t0/1:60:5,6\t1/1:.:.",
        "short-format" => "chr1\t400\t.\tT\tC\t50\tPASS\tDP=10\tGT\t0/1\t1/1",
        other => panic!("unknown case {other}"),
    }
}

/// One case in one mode: read the file, do what the mode says, and encode.
fn run(case: &str, mode: &str) -> (bool, String) {
    let text = header_text() + line_of(case) + "\n";
    let file = read_vcf(&text).expect("the dump's inputs all parse");
    let mut record = file.records.into_iter().next().expect("one record");

    match mode {
        "copy" => {}
        // A READ, and nothing else. `getGenotypes().get(0)` decodes in the reference, and
        // `decode()` is what drops the text.
        "touch" => {
            let _ = record.genotypes.first();
        }
        // The genotypes replaced by an equal list, which is what a builder does.
        "rebuild" => {
            let rebuilt: Vec<Genotype> = record.genotypes.to_vec();
            record.genotypes = rebuilt.into();
        }
        other => panic!("unknown mode {other}"),
    }

    let lazy = record.genotypes.unparsed().is_some();
    let header = file.header.clone();
    let encoder = VcfEncoder::new(&header);
    let encoded = encoder.encode(&record).expect("the record encodes");
    (lazy, encoded)
}

#[test]
fn every_state_is_the_reference() {
    let text = golden();
    let lazy_rows = rows(&text, "lazy");
    let line_rows = rows(&text, "line");
    assert_eq!(lazy_rows.len(), 12, "four cases in three modes");
    assert_eq!(line_rows.len(), lazy_rows.len());

    for (state, line) in lazy_rows.iter().zip(&line_rows) {
        let (case, mode) = (state[1].as_str(), state[2].as_str());
        let (lazy, encoded) = run(case, mode);
        assert_eq!(
            lazy.to_string(),
            state[3],
            "whether the context still answers with its text: {case}/{mode}"
        );
        assert_eq!(escape(&encoded), line[3], "the encoded line: {case}/{mode}");
    }
}

/// The one that costs a consumer its output: a copied record keeps the file's key order.
#[test]
fn a_copied_record_keeps_the_files_own_format_order() {
    let (lazy, encoded) = run("unsorted-format", "copy");
    assert!(lazy);
    assert!(encoded.contains("GT:GQ:DP"), "{encoded}");

    // And looking at one genotype is enough to lose it.
    let (lazy, encoded) = run("unsorted-format", "touch");
    assert!(!lazy);
    assert!(encoded.contains("GT:DP:GQ"), "{encoded}");
}

/// `size()` and `isEmpty()` are answered without decoding, which is htsjdk's own optimisation.
#[test]
fn counting_the_genotypes_does_not_decode_them() {
    let text = header_text() + line_of("unsorted-format") + "\n";
    let file = read_vcf(&text).expect("it parses");
    let record = &file.records[0];
    assert_eq!(record.genotypes.len(), 2);
    assert!(!record.genotypes.is_empty());
    assert!(
        record.genotypes.unparsed().is_some(),
        "counting must not be what rewrites the FORMAT column"
    );
}
