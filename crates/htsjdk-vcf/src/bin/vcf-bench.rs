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
//! Usage: vcf-bench [records] [reps]

use htsjdk_vcf::encoder::VcfEncoder;
use htsjdk_vcf::reader::read_vcf;

fn corpus(records: usize, samples: usize) -> String {
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
        let reference = bases[(next() % 4) as usize] as char;
        let alternate = bases[((next() % 3) as usize + 1) % 4] as char;
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

    let text = corpus(records, samples);
    let megabytes = text.len() as f64 / (1024.0 * 1024.0);
    println!(
        "vcf_records={records} samples={samples} text_bytes={}",
        text.len()
    );

    let file = read_vcf(&text).expect("the corpus parses");
    assert_eq!(file.records.len(), records);
    for run in 0..reps {
        let start = std::time::Instant::now();
        let parsed = read_vcf(&text).expect("parses");
        let seconds = start.elapsed().as_secs_f64();
        println!(
            "rust_vcf_read_run{run}_mbps={:.2} recs_per_sec={:.0}",
            megabytes / seconds,
            records as f64 / seconds
        );
        std::hint::black_box(parsed);
    }

    for mode in ["copy", "touched"] {
        for run in 0..reps {
            // A fresh read each run: encoding is destructive of the laziness, so a second pass
            // over the same records would measure the other mode.
            let file = read_vcf(&text).expect("parses");
            let encoder = VcfEncoder::new(&file.header);
            let mut records_out = file.records;
            if mode == "touched" {
                for record in &mut records_out {
                    let _ = record.genotypes.first();
                }
            }
            if samples == 0 && mode == "touched" {
                // Nothing to touch, and nothing the mode would measure.
                break;
            }
            let mut out = String::with_capacity(text.len());
            let start = std::time::Instant::now();
            for record in &records_out {
                encoder.encode_into(record, &mut out).expect("encodes");
                out.push('\n');
            }
            let seconds = start.elapsed().as_secs_f64();
            println!(
                "rust_vcf_encode_{mode}_run{run}_mbps={:.2} recs_per_sec={:.0} bytes={}",
                out.len() as f64 / (1024.0 * 1024.0) / seconds,
                records as f64 / seconds,
                out.len()
            );
            std::hint::black_box(out);
        }
    }
}
