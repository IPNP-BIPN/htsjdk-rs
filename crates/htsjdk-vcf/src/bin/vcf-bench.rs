//! The VCF half of the I/O floor benchmark: reading a file and writing it back.
//!
//! Issue #78 names "VCF encode and decode" as the third path under everything. The two modes here
//! are not two implementations of one thing, they are the two the library actually has:
//!
//!   * **copy** encodes records nothing has looked at, which takes `VCFEncoder`'s verbatim branch
//!     and writes the genotype columns straight from the file's own text;
//!   * **touched** reads one genotype per record first, which is `decode()` in the reference and
//!     makes the encoder rebuild the FORMAT column and every sample column under it.
//!
//! Every tool that copies records through pays the first and every tool that inspects them pays the
//! second, so the gap between the two rows is what laziness is worth.
//!
//! It writes the corpus to a file and reads it back, so `tools/benchmark/CodecBench.java` can time
//! htsjdk's `VCFFileReader` and `VCFWriter` on the very same bytes, and the runner compares the
//! files the two sides write before it prints either duration.
//!
//! Usage: vcf-bench [records] [reps] [samples] [dir]

use std::path::{Path, PathBuf};

use htsjdk_vcf::reader::{read_vcf, VcfFile};
use htsjdk_vcf::vcf_file::write_vcf;

fn corpus_text(records: usize, samples: usize) -> String {
    let mut text = String::from("##fileformat=VCFv4.2\n");
    text.push_str("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n");
    text.push_str("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Frequency\">\n");
    text.push_str("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    text.push_str("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Quality\">\n");
    text.push_str("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n");
    text.push_str("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Depths\">\n");
    text.push_str("##contig=<ID=chr1,length=250000000>\n");
    text.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
    if samples > 0 {
        text.push_str("\tFORMAT");
        for sample in 0..samples {
            text.push_str(&format!("\ts{sample}"));
        }
    }
    text.push('\n');
    let mut seed = 12345u64;
    let mut next = move || {
        seed = seed
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        seed >> 33
    };
    for index in 0..records {
        let position = 1 + index as u64 * 37 % 249_000_000;
        let bases = *b"ACGT";
        // ALT is offset from REF so the two never coincide: htsjdk rejects a duplicate allele.
        let base = (next() % 4) as usize;
        let reference = bases[base] as char;
        let alternate = bases[(base + 1 + (next() % 3) as usize) % 4] as char;
        // The FORMAT is deliberately NOT in the order the encoder computes, which is what makes
        // the two modes measurably different files as well as different speeds.
        text.push_str(&format!(
            "chr1\t{position}\trs{index}\t{reference}\t{alternate}\t{}\tPASS\tDP={};AF={:.3}{}",
            next() % 100,
            next() % 200,
            (next() % 1000) as f64 / 1000.0,
            if samples > 0 { "\tGT:GQ:DP:AD" } else { "" },
        ));
        for sample in 0..samples {
            text.push_str(&format!(
                "\t{}/{}:{}:{}:{},{}",
                sample % 2,
                (sample + 1) % 2,
                next() % 99,
                next() % 100,
                next() % 50,
                next() % 50
            ));
        }
        text.push('\n');
    }
    text
}

fn main() {
    let mut args = std::env::args().skip(1);
    let records: usize = args
        .next()
        .map(|a| a.parse().expect("records"))
        .unwrap_or(100_000);
    let reps: usize = args.next().map(|a| a.parse().expect("reps")).unwrap_or(3);
    // The sample count is the axis the genotype columns live on: reading a sites-only file and
    // reading the same sites with four samples separates the site parser's cost from theirs.
    let samples: usize = args
        .next()
        .map(|a| a.parse().expect("samples"))
        .unwrap_or(4);
    let dir = PathBuf::from(args.next().unwrap_or_else(|| ".".to_string()));

    let corpus = dir.join("corpus.vcf");
    let text = corpus_text(records, samples);
    std::fs::write(&corpus, &text).expect("writes the corpus");
    println!(
        "vcf_records={records} samples={samples} vcf_text_bytes={}",
        text.len()
    );
    assert_eq!(read(&corpus, false).records.len(), records);

    // `lazy` is what htsjdk's reader leaves behind, and the port's too; `decoded` reaches one
    // genotype per record, which is `LazyGenotypesContext.decode()` in the reference.
    for (mode, decode) in [("lazy", false), ("decoded", true)] {
        for run in 0..reps {
            let start = std::time::Instant::now();
            std::hint::black_box(read(&corpus, decode));
            println!(
                "rust_vcf_read_{mode}_run{run}_secs={:.6}",
                start.elapsed().as_secs_f64()
            );
        }
    }

    for (mode, decode) in [("copy", false), ("touched", true)] {
        if samples == 0 && decode {
            // Nothing to touch, and nothing the mode would measure.
            break;
        }
        let out = dir.join(format!("port_{mode}.vcf"));
        for run in 0..reps {
            // A fresh read each run, untimed: decoding is what makes the encoder leave the verbatim
            // branch, so a write must not see records an earlier one touched.
            let file = read(&corpus, decode);
            let start = std::time::Instant::now();
            let written = write_vcf(&file.header, &file.records).expect("encodes");
            std::fs::write(&out, written).expect("writes");
            println!(
                "rust_vcf_write_{mode}_run{run}_secs={:.6}",
                start.elapsed().as_secs_f64()
            );
        }
    }
}

fn read(path: &Path, decode: bool) -> VcfFile {
    let file = read_vcf(&std::fs::read_to_string(path).expect("reads")).expect("parses");
    if decode {
        for record in &file.records {
            let _ = record.genotypes.first();
        }
    }
    file
}
