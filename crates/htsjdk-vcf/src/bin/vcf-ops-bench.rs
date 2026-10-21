//! The `htsjdk-vcf` and `htsjdk-tribble` operations `vcf-bench` does not time, on the inputs
//! `tools/benchmark/VcfOpsBench.java` times htsjdk on: BCF write and read, a VCF with its `.idx`
//! written on the fly, a `.vcf.gz` with its `.tbi`, reading a `.vcf.gz`, and BED decoding.
//!
//! The output protocol is `bam-ops-bench`'s: an `op` line, `rust_<op>_run<r>_secs` lines, and
//! `port_<op><suffix>` files the runner compares with htsjdk's `java_<op><suffix>`. `check=idx`
//! compares a Tribble index with the indexed file's path and modification time taken out, the two
//! fields that name the file rather than describe it.
//!
//! Usage: vcf-ops-bench dir reps   (dir holds vcf-bench's corpus.vcf)

/// The allocator picard-rs links, so the rows measure the library as it ships there.
#[cfg(feature = "mimalloc")]
#[global_allocator]
static ALLOCATOR: mimalloc::MiMalloc = mimalloc::MiMalloc;

use std::fmt::Write as _;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use htsjdk_bgzf::{decompress_all, BgzfWriter};
use htsjdk_tribble::bed::{decode as decode_bed, StartOffset};
use htsjdk_tribble::tabix::{FeatureRef, TabixFormat, TabixIndexCreator};
use htsjdk_vcf::bcf2::codec::read_bcf;
use htsjdk_vcf::bcf2::writer::{Bcf2Writer, Origin};
use htsjdk_vcf::encoder::VcfEncoder;
use htsjdk_vcf::indexed_writer::{write_vcf_indexed, SequenceEntry};
use htsjdk_vcf::reader::{read_vcf, VcfFile};
use htsjdk_vcf::vcf_file::write_vcf;

struct Bench {
    dir: PathBuf,
    reps: usize,
}

impl Bench {
    fn op(&self, name: &str, bytes: usize, records: usize, check: &str, suffixes: &str) {
        println!("op name={name} bytes={bytes} records={records} check={check} files={suffixes}");
    }

    /// One untimed run, whose answer is the check, then `reps` timed ones, each on a fresh
    /// `setup()` made outside the timing.
    fn time<S, T>(&self, name: &str, setup: impl Fn() -> S, f: impl Fn(S) -> T) -> T {
        let out = f(setup());
        for run in 0..self.reps {
            let input = setup();
            let start = Instant::now();
            let answer = std::hint::black_box(f(input));
            let seconds = start.elapsed().as_secs_f64();
            // The answer is freed after the clock stops, as the JVM's collector frees htsjdk's later.
            drop(answer);
            println!("rust_{name}_run{run}_secs={:.6}", seconds);
        }
        out
    }

    fn path(&self, op: &str, suffix: &str) -> PathBuf {
        self.dir.join(format!("port_{op}{suffix}"))
    }
}

fn read(path: &Path) -> VcfFile {
    read_vcf(&std::fs::read_to_string(path).expect("reads")).expect("parses")
}

/// The dictionary of `vcf-bench`'s corpus, whose one `##contig` line this is.
fn dictionary() -> Vec<SequenceEntry> {
    vec![SequenceEntry {
        name: "chr1".to_string(),
        length: 250_000_000,
    }]
}

fn main() {
    let mut args = std::env::args().skip(1);
    let dir = PathBuf::from(args.next().expect("dir"));
    let reps: usize = args.next().map(|a| a.parse().expect("reps")).unwrap_or(3);
    let b = Bench { dir, reps };
    let corpus = b.dir.join("corpus.vcf");
    let size = std::fs::metadata(&corpus).expect("exists").len() as usize;
    let n = read(&corpus).records.len();

    // BCF: htsjdk writes it uncompressed, as BCF2Writer does to a file.
    b.op("bcf_write", size, n, "bytes", ".bcf");
    let bcf = b.path("bcf_write", ".bcf");
    b.time(
        "bcf_write",
        || read(&corpus),
        |file| {
            let mut writer = Bcf2Writer::new(&file.header, false).expect("a writer");
            let mut out = BufWriter::new(std::fs::File::create(&bcf).expect("creates"));
            out.write_all(&writer.header_bytes()).expect("writes");
            for vc in &file.records {
                out.write_all(&writer.encode(vc, Origin::VcfCodec).expect("encodes"))
                    .expect("writes");
            }
            out.flush().expect("flushes");
        },
    );

    // Both read the BCF the port wrote; the check is the records written back out as VCF text.
    b.op("bcf_read", size, n, "decoded", ".vcf");
    let file = b.time(
        "bcf_read",
        || (),
        |()| read_bcf(&std::fs::read(&bcf).expect("reads")).expect("decodes"),
    );
    let encoder = VcfEncoder::new(&file.header);
    let mut text = file.header.write();
    for vc in &file.records {
        encoder.encode_into(vc, &mut text).expect("encodes");
        text.push('\n');
    }
    std::fs::write(b.path("bcf_read", ".vcf"), text).expect("writes");

    // The VCF and its `.idx` in one pass, as every GATK tool writing a `.vcf` does.
    let dict = dictionary();
    b.op("vcf_idx_write", size, n, "idx", ".vcf,.vcf.idx");
    let out = b.path("vcf_idx_write", ".vcf");
    let uri = format!("file://{}", out.display());
    b.time(
        "vcf_idx_write",
        || read(&corpus),
        |file| {
            let written =
                write_vcf_indexed(&file.header, &file.records, Some(&dict), &uri).expect("writes");
            std::fs::write(&out, &written.text).expect("writes");
            // The two fields the writer leaves to its caller, filled as picard-rs's `VcfWriter`
            // fills them: from the file it just wrote.
            let mut index = written.index;
            index.indexed_file_size = written.text.len() as i64;
            index.indexed_file_timestamp = std::fs::metadata(&out)
                .and_then(|m| m.modified())
                .expect("a mtime")
                .duration_since(std::time::UNIX_EPOCH)
                .expect("after 1970")
                .as_millis() as i64;
            let idx = index.write().expect("serialises");
            std::fs::write(b.path("vcf_idx_write", ".vcf.idx"), idx).expect("writes");
        },
    );

    b.op("vcf_gz_write", size, n, "bytes", ".vcf.gz,.vcf.gz.tbi");
    let gz = b.path("vcf_gz_write", ".vcf.gz");
    b.time(
        "vcf_gz_write",
        || read(&corpus),
        |file| write_vcf_gz(&file, &dict, &gz),
    );

    // Reading it back: the inflate, then the same parse as a `.vcf`.
    b.op("vcf_gz_read", size, n, "decoded", ".vcf");
    let file = b.time(
        "vcf_gz_read",
        || (),
        |()| {
            let text = decompress_all(&std::fs::read(&gz).expect("reads")).expect("inflates");
            read_vcf(std::str::from_utf8(&text).expect("text")).expect("parses")
        },
    );
    std::fs::write(
        b.path("vcf_gz_read", ".vcf"),
        write_vcf(&file.header, &file.records).expect("encodes"),
    )
    .expect("writes");

    bed(&b);
}

/// `VariantContextWriterBuilder` on a `.vcf.gz`: BGZF through the default deflater, a
/// `TabixIndexCreator(dict, TabixFormat.VCF)` fed each record's virtual offset before its line, and
/// the `.tbi` written at close through the same deflater. picard-rs's `VcfWriter` is this.
fn write_vcf_gz(file: &VcfFile, dict: &[SequenceEntry], path: &Path) {
    let mut creator = TabixIndexCreator::new(TabixFormat::VCF);
    let mut bgzf = BgzfWriter::new(BufWriter::new(
        std::fs::File::create(path).expect("creates"),
    ));
    bgzf.write_all(file.header.write().as_bytes())
        .expect("writes");
    let encoder = VcfEncoder::new(&file.header);
    let mut line = String::new();
    for vc in &file.records {
        let at = bgzf.file_pointer().expect("a position");
        let length = dict
            .iter()
            .find(|s| s.name == vc.contig)
            .map_or(0, |s| s.length as i32);
        creator
            .add_feature(
                FeatureRef {
                    contig: &vc.contig,
                    start: vc.start as i32,
                    end: vc.stop as i32,
                    description: &vc.contig,
                    sequence_length: length,
                },
                at as i64,
            )
            .expect("in order");
        line.clear();
        encoder.encode_into(vc, &mut line).expect("encodes");
        line.push('\n');
        bgzf.write_all(line.as_bytes()).expect("writes");
    }
    // `VCFWriter.close()` flushes before the indexer asks for the final position.
    bgzf.flush().expect("flushes");
    let end = bgzf.file_pointer().expect("a position") as i64;
    bgzf.into_inner()
        .expect("finishes")
        .flush()
        .expect("flushes");
    let index = creator.finish(end).expect("an index");
    let mut tbi = path.as_os_str().to_os_string();
    tbi.push(".tbi");
    std::fs::write(tbi, index.write()).expect("writes");
}

/// `BEDCodec.decode` on every line of a six-column BED.
fn bed(b: &Bench) {
    let mut seed = 99u64;
    let mut next = move |bound: u64| {
        seed = seed
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (seed >> 33) % bound
    };
    let mut text = String::new();
    for i in 0..500_000 {
        let start = next(100_000_000);
        writeln!(
            text,
            "chr{}\t{start}\t{}\tfeature_{i}\t{}\t{}",
            1 + next(22),
            start + 1 + next(5000),
            next(1000),
            ["+", "-", "."][next(3) as usize]
        )
        .unwrap();
    }
    let input = b.dir.join("corpus.bed");
    std::fs::write(&input, &text).expect("writes");
    b.op("bed_decode", text.len(), 500_000, "decoded", ".txt");
    let features = b.time(
        "bed_decode",
        || (),
        |()| {
            std::fs::read_to_string(&input)
                .expect("reads")
                .lines()
                .map(|l| decode_bed(l, StartOffset::One).expect("decodes"))
                .collect::<Vec<_>>()
        },
    );
    let mut out = String::new();
    for f in features.into_iter().flatten() {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}",
            f.contig,
            f.start,
            f.end,
            f.name,
            f.score.to_bits() as i32,
            f.strand.name()
        )
        .unwrap();
    }
    std::fs::write(b.path("bed_decode", ".txt"), out).expect("writes");
}
