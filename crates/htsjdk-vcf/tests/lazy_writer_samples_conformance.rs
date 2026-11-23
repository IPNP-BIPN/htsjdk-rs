//! Conformance for a copied record's genotype columns under a writer whose samples are not the
//! file's, against the oracle.
//!
//! Golden from `tools/vcf-conformance/VcfLazyWriterSamplesDump.java`.
//!
//! `VCFEncoder` copies a lazy record's text without asking the writer's header anything, so the
//! columns come out as the file had them whatever samples the header names:
//!
//! ```text
//! header    #CHROM ... FORMAT  s1         s0
//! copy      ...        GT:GQ:DP  0/1:60:10  1/1:70:20   <- s0's genotype under s1, verbatim
//! touch     ...        GT:DP:GQ  1/1:20:70  0/1:10:60   <- decoded, written by sample name
//! ```
//!
//! A file whose samples are not sorted is decoded at the read, so it is written by name either way.

use std::io::Read;

use htsjdk_vcf::reader::read_vcf;
use htsjdk_vcf::vcf_file::write_vcf;

fn golden() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/vcf_lazy_writer_samples.txt.gz");
    let file = std::fs::File::open(&path).expect("golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("golden is gzip");
    text
}

fn escape(text: &str) -> String {
    text.replace('\\', "\\\\")
        .replace('\t', "\\t")
        .replace('\n', "\\n")
}

/// The dump's input file, its samples named in `samples` order.
fn input(samples: &str) -> String {
    format!(
        concat!(
            "##fileformat=VCFv4.2\n",
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
            "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Quality\">\n",
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n",
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Depths\">\n",
            "##contig=<ID=chr1,length=100000>\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n",
            "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\tGT:GQ:DP\t0/1:60:10\t1/1:70:20\n",
            "chr1\t200\t.\tG\tA\t50\tPASS\tDP=10\tGT:GQ:AD\t0/1:60:5,6\t1/1:.:.\n",
        ),
        samples
    )
}

/// One row: read the input, touch every record or not, and write under `samples`.
fn run(input_name: &str, samples: &str, mode: &str) -> String {
    let text = match input_name {
        "sorted" => input("s0\ts1"),
        "unsorted" => input("s1\ts0"),
        other => panic!("unknown input {other}"),
    };
    let file = read_vcf(&text).expect("the dump's inputs all parse");
    let mut header = file.header.clone();
    header.samples = samples
        .split(',')
        .filter(|name| !name.is_empty())
        .map(str::to_string)
        .collect();
    match mode {
        "copy" => {}
        // A READ, which is `getGenotypes().get(0)` in the reference and decodes.
        "touch" => {
            for record in &file.records {
                let _ = record.genotypes.first();
            }
        }
        other => panic!("unknown mode {other}"),
    }
    write_vcf(&header, &file.records).expect("every case writes")
}

#[test]
fn every_writer_header_is_the_reference() {
    let text = golden();
    let rows: Vec<Vec<&str>> = text
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').collect())
        .collect();
    assert_eq!(
        rows.len(),
        16,
        "six headers and two for the unsorted input, two modes each"
    );
    for row in &rows {
        let (input_name, samples, mode) = (row[1], row[2], row[3]);
        assert_eq!(row[4], "ok", "{input_name}/{samples}/{mode}");
        assert_eq!(
            escape(&run(input_name, samples, mode)),
            row[5],
            "the written file: {input_name}/{samples}/{mode}"
        );
    }
}

/// The one this suite exists for: reordering the header's samples does not reorder the columns.
#[test]
fn a_copied_record_keeps_its_columns_under_a_reordered_header() {
    let written = run("sorted", "s1,s0", "copy");
    assert!(
        written.contains(
            "FORMAT\ts1\ts0\nchr1\t100\t.\tA\tT\t50\tPASS\tDP=10\tGT:GQ:DP\t0/1:60:10\t1/1:70:20\n"
        ),
        "{written}"
    );
}
