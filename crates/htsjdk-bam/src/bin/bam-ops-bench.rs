//! Every `htsjdk-bam` operation that has an htsjdk counterpart, timed on the inputs
//! `tools/benchmark/BamOpsBench.java` times htsjdk on.
//!
//! `record-bench` covers BAM read and write; this covers the rest of the crate: SAM text, sorting,
//! the BAI and SBI indexes, block-copy reheader and gather, FASTA and its `.fai`/`.dict`, indexed
//! FASTA queries, FASTQ, interval lists, liftover, MD/NM and the overlap detector.
//!
//! Each operation prints one `op` line, then `rust_<op>_run<r>_secs` per timed run. The untimed
//! first run is both the warm-up and the output the runner checks: `port_<op><suffix>` against the
//! `java_<op><suffix>` htsjdk wrote, byte for byte. `check=bytes` means those files are the
//! operation's own output; `check=decoded` means the operation's output is in memory (records, a
//! sorted list, query answers) and the files are a rendering of it, written outside the timing.
//!
//! Usage: bam-ops-bench dir reps megabytes   (dir holds record-bench's corpus.bam)

use std::fmt::Write as _;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use htsjdk_bam::fasta::read_fasta_file;
use htsjdk_bam::fasta_index::IndexedFasta;
use htsjdk_bam::fasta_writer::FastaReferenceWriter;
use htsjdk_bam::fastq::{as_fastq_record, parse_fastq, write_record, FastqRecord};
use htsjdk_bam::header::SequenceRecord;
use htsjdk_bam::interval::{Interval, IntervalList};
use htsjdk_bam::liftover::LiftOver;
use htsjdk_bam::md_nm::calculate_md_and_nm;
use htsjdk_bam::overlap::OverlapDetector;
use htsjdk_bam::reader::parse_header_text;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::{read_sam_with, write_sam};
use htsjdk_bam::sbi::SbiIndexWriter;
use htsjdk_bam::text_parse::ValidationStringency;
use htsjdk_bam::{
    build_bam_index, coordinate, gather_bam_files, query_name, reheader_bam, BamReader, BamWriter,
    SamHeader,
};
use htsjdk_bgzf::decompress_all;
use md5::{Digest, Md5};

struct Bench {
    dir: PathBuf,
    reps: usize,
}

impl Bench {
    /// The `op` line the runner keys every row on.
    fn op(&self, name: &str, bytes: usize, records: usize, check: &str, suffixes: &str) {
        println!("op name={name} bytes={bytes} records={records} check={check} files={suffixes}");
    }

    /// One untimed run, whose answer is returned for the check, then `reps` timed ones.
    fn time<T>(&self, name: &str, mut f: impl FnMut() -> T) -> T {
        let out = f();
        for run in 0..self.reps {
            let start = Instant::now();
            std::hint::black_box(f());
            println!(
                "rust_{name}_run{run}_secs={:.6}",
                start.elapsed().as_secs_f64()
            );
        }
        out
    }

    fn path(&self, op: &str, suffix: &str) -> PathBuf {
        self.dir.join(format!("port_{op}{suffix}"))
    }

    fn input(&self, name: &str) -> PathBuf {
        self.dir.join(name)
    }
}

/// The same 64-bit LCG the other benches use.
struct Lcg(u64);

impl Lcg {
    fn below(&mut self, bound: u64) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 33) % bound
    }
}

fn md5_hex(bytes: &[u8]) -> String {
    Md5::digest(bytes)
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}

fn read_bam(path: &Path) -> (SamHeader, Vec<BamRecord>) {
    let data = decompress_all(&std::fs::read(path).expect("reads")).expect("inflates");
    let reader = BamReader::new(&data).expect("a BAM");
    let header = reader.header.text.clone();
    (header, reader.map(|r| r.expect("decodes")).collect())
}

fn write_bam<'a>(
    path: &Path,
    header: &SamHeader,
    records: impl IntoIterator<Item = &'a BamRecord>,
) {
    let file = std::fs::File::create(path).expect("creates");
    let mut writer = BamWriter::new(BufWriter::new(file), header).expect("writes");
    for record in records {
        writer.write(record).expect("encodes");
    }
    writer.finish().expect("finishes").flush().expect("flushes");
}

fn main() {
    let mut args = std::env::args().skip(1);
    let dir = PathBuf::from(args.next().expect("dir"));
    let reps: usize = args.next().map(|a| a.parse().expect("reps")).unwrap_or(3);
    let megabytes: usize = args.next().map(|a| a.parse().expect("MB")).unwrap_or(64);
    let b = Bench { dir, reps };

    let (header, records) = read_bam(&b.input("corpus.bam"));
    sam(&b, &header, &records);
    sorting(&b, &header, &records);
    indexes(&b);
    block_copy(&b);
    fastq(&b, &records);
    let genome = fasta(&b, megabytes);
    md_nm(&b, &records, &genome[0].1);
    intervals(&b, &genome);
}

fn sam(b: &Bench, header: &SamHeader, records: &[BamRecord]) {
    let text = write_sam(header, records).expect("renders");
    b.op("sam_write", text.len(), records.len(), "bytes", ".sam");
    let out = b.path("sam_write", ".sam");
    b.time("sam_write", || {
        let text = write_sam(header, records).expect("renders");
        std::fs::write(&out, text).expect("writes");
    });

    // Both sides read the file the port just wrote, and the check is what they read, as SAM again.
    b.op("sam_read", text.len(), records.len(), "decoded", ".sam");
    let (h, read) = b.time("sam_read", || {
        let text = std::fs::read_to_string(&out).expect("reads");
        read_sam_with(&text, ValidationStringency::Silent).expect("parses")
    });
    std::fs::write(
        b.path("sam_read", ".sam"),
        write_sam(&h, &read).expect("renders"),
    )
    .expect("writes");
}

/// Queryname order from the corpus's own order, then coordinate order from queryname order, which
/// is what SortSam does to an aligner's name-grouped output. The corpus is nearly in coordinate
/// order already, and sorting it again would time a merge of four runs.
fn sorting(b: &Bench, header: &SamHeader, records: &[BamRecord]) {
    let bytes = std::fs::metadata(b.input("corpus.bam"))
        .expect("exists")
        .len() as usize;
    // A vector of references, which is what `new ArrayList<>(records)` copies there.
    let mut input: Vec<&BamRecord> = records.iter().collect();
    for name in ["queryname", "coordinate"] {
        let op = format!("sort_{name}");
        b.op(&op, bytes, records.len(), "decoded", ".bam");
        let sorted = b.time(&op, || {
            let mut sorted = input.clone();
            if name == "coordinate" {
                sorted.sort_by(|a, b| coordinate::compare(a, b));
            } else {
                sorted.sort_by(|a, b| query_name::compare(a, b));
            }
            sorted
        });
        let mut h = header.clone();
        h.attributes.set("SO", name);
        write_bam(&b.path(&op, ".bam"), &h, sorted.iter().copied());
        input = sorted;
    }
}

fn indexes(b: &Bench) {
    let sorted = b.path("sort_coordinate", ".bam");
    let bytes = std::fs::metadata(&sorted).expect("exists").len() as usize;
    let (_, records) = read_bam(&sorted);

    b.op("bai_build", bytes, records.len(), "bytes", ".bai");
    let out = b.path("bai_build", ".bai");
    b.time("bai_build", || {
        let bai = build_bam_index(&std::fs::read(&sorted).expect("reads")).expect("indexes");
        std::fs::write(&out, bai).expect("writes");
    });

    // The writer alone, over a synthetic offset sequence `BamOpsBench` reproduces: the port has no
    // BAM walker that hands out virtual offsets, so the file-level `SBIIndexer` is not comparable.
    let count = 10_000_000u64;
    b.op(
        "sbi_write",
        (count * 8) as usize,
        count as usize,
        "bytes",
        ".sbi",
    );
    let out = b.path("sbi_write", ".sbi");
    b.time("sbi_write", || {
        let mut writer = SbiIndexWriter::new(4096);
        for i in 0..count {
            writer.process_record(sbi_offset(i)).expect("increasing");
        }
        let sbi = writer
            .finish(sbi_offset(count), 1 << 40, None, None)
            .expect("finishes");
        std::fs::write(&out, sbi).expect("writes");
    });
}

/// A virtual offset per record, increasing: forty records to a block of 9000 compressed bytes.
fn sbi_offset(i: u64) -> u64 {
    ((i / 40) * 9000) << 16 | ((i % 40) * 1500)
}

fn block_copy(b: &Bench) {
    let input = b.path("sort_coordinate", ".bam");
    let (mut header, records) = read_bam(&input);
    header.comments.push("@CO\tbenchmark".to_string());
    let bytes = std::fs::metadata(&input).expect("exists").len() as usize;
    b.op("bam_reheader", bytes, records.len(), "bytes", ".bam");
    let out = b.path("bam_reheader", ".bam");
    b.time("bam_reheader", || {
        let bam = reheader_bam(&header, &std::fs::read(&input).expect("reads")).expect("copies");
        std::fs::write(&out, bam).expect("writes");
    });

    let inputs = [
        b.input("corpus.bam"),
        b.path("sort_coordinate", ".bam"),
        b.path("sort_queryname", ".bam"),
    ];
    let bytes: usize = inputs
        .iter()
        .map(|p| std::fs::metadata(p).expect("exists").len() as usize)
        .sum();
    b.op("bam_gather", bytes, records.len() * 3, "bytes", ".bam");
    let out = b.path("bam_gather", ".bam");
    b.time("bam_gather", || {
        let data: Vec<Vec<u8>> = inputs
            .iter()
            .map(|p| std::fs::read(p).expect("reads"))
            .collect();
        let slices: Vec<&[u8]> = data.iter().map(Vec::as_slice).collect();
        std::fs::write(&out, gather_bam_files(&slices).expect("gathers")).expect("writes");
    });
}

fn fastq(b: &Bench, records: &[BamRecord]) {
    let fastq: Vec<FastqRecord> = records.iter().map(as_fastq_record).collect();
    let size: usize = fastq.iter().map(|r| write_record(r).len()).sum();
    b.op("fastq_write", size, fastq.len(), "bytes", ".fastq");
    let out = b.path("fastq_write", ".fastq");
    b.time("fastq_write", || {
        let mut w = BufWriter::new(std::fs::File::create(&out).expect("creates"));
        for record in &fastq {
            w.write_all(write_record(record).as_bytes())
                .expect("writes");
        }
        w.flush().expect("flushes");
    });

    b.op("fastq_read", size, fastq.len(), "decoded", ".fastq");
    let read = b.time("fastq_read", || {
        parse_fastq(&std::fs::read_to_string(&out).expect("reads"), false).expect("parses")
    });
    let text: String = read.iter().map(write_record).collect();
    std::fs::write(b.path("fastq_read", ".fastq"), text).expect("writes");
}

/// A genome of eight contigs: mixed-case ACGT with runs of N, the shape of a real assembly.
fn genome(megabytes: usize) -> Vec<(String, Vec<u8>)> {
    let mut rng = Lcg(424242);
    let per = megabytes * 1024 * 1024 / 8;
    (1..=8)
        .map(|c| {
            let mut bases = Vec::with_capacity(per);
            while bases.len() < per {
                if rng.below(1000) == 0 {
                    let run = 1 + rng.below(500) as usize;
                    bases.extend(std::iter::repeat_n(b'N', run.min(per - bases.len())));
                } else {
                    bases.push(b"ACGTacgt"[rng.below(8) as usize]);
                }
            }
            (format!("chr{c}"), bases)
        })
        .collect()
}

fn fasta(b: &Bench, megabytes: usize) -> Vec<(String, Vec<u8>)> {
    let genome = genome(megabytes);
    let size: usize = genome.iter().map(|(_, s)| s.len()).sum();
    b.op(
        "fasta_write",
        size,
        genome.len(),
        "bytes",
        ".fasta,.fasta.fai,.dict",
    );
    let fasta = b.path("fasta_write", ".fasta");
    b.time("fasta_write", || {
        let mut writer = FastaReferenceWriter::new(60, true).expect("a writer");
        for (name, bases) in &genome {
            writer.start_sequence(name).expect("starts");
            writer.append_bases(bases).expect("appends");
        }
        let out = writer.close().expect("closes");
        std::fs::write(&fasta, out.fasta).expect("writes");
        std::fs::write(b.path("fasta_write", ".fasta.fai"), out.index).expect("writes");
        std::fs::write(b.path("fasta_write", ".dict"), out.dictionary).expect("writes");
    });

    b.op("fasta_read", size, genome.len(), "decoded", ".txt");
    let read = b.time("fasta_read", || read_fasta_file(&fasta).expect("reads"));
    let mut text = String::new();
    for s in &read {
        writeln!(text, "{}\t{}\t{}", s.name, s.bases.len(), md5_hex(&s.bases)).unwrap();
    }
    std::fs::write(b.path("fasta_read", ".txt"), text).expect("writes");

    // Read-length queries at random places, through the `.fai`.
    let queries = queries(&genome, 200_000, 151);
    b.op(
        "fasta_query",
        queries.len() * 151,
        queries.len(),
        "decoded",
        ".txt",
    );
    write_queries(&b.input("fasta_queries.txt"), &queries);
    let digest = b.time("fasta_query", || {
        let mut fasta = IndexedFasta::open(&fasta).expect("indexed");
        let mut md5 = Md5::new();
        for (contig, start, end) in &queries {
            md5.update(
                fasta
                    .query(contig, *start as i64, *end as i64)
                    .expect("in range"),
            );
        }
        md5.finalize()
    });
    let hex: String = digest.iter().map(|b| format!("{b:02x}")).collect();
    std::fs::write(b.path("fasta_query", ".txt"), hex + "\n").expect("writes");
    genome
}

fn queries(genome: &[(String, Vec<u8>)], count: usize, length: i32) -> Vec<(String, i32, i32)> {
    let mut rng = Lcg(777);
    (0..count)
        .map(|_| {
            let (name, bases) = &genome[rng.below(genome.len() as u64) as usize];
            let start = 1 + rng.below((bases.len() as i32 - length) as u64) as i32;
            (name.clone(), start, start + length - 1)
        })
        .collect()
}

fn write_queries(path: &Path, queries: &[(String, i32, i32)]) {
    let mut text = String::new();
    for (c, s, e) in queries {
        writeln!(text, "{c}\t{s}\t{e}").unwrap();
    }
    std::fs::write(path, text).expect("writes");
}

/// `calculateMdAndNmTags(record, chr1, true, true)` for every corpus record.
fn md_nm(b: &Bench, records: &[BamRecord], chr1: &[u8]) {
    let size: usize = records.iter().map(|r| r.read_bases.len()).sum();
    b.op("md_nm", size, records.len(), "decoded", ".txt");
    let tags = b.time("md_nm", || {
        records
            .iter()
            .map(|r| calculate_md_and_nm(r.alignment_start, &r.cigar, &r.read_bases, chr1))
            .collect::<Vec<_>>()
    });
    let mut text = String::new();
    for (md, nm) in tags {
        writeln!(text, "{md}\t{nm}").unwrap();
    }
    std::fs::write(b.path("md_nm", ".txt"), text).expect("writes");
}

fn intervals(b: &Bench, genome: &[(String, Vec<u8>)]) {
    let mut header = SamHeader::default();
    for (name, bases) in genome {
        header
            .sequences
            .push(SequenceRecord::new(name, bases.len() as i32));
    }
    let dictionary: Vec<String> = genome.iter().map(|(n, _)| n.clone()).collect();
    let mut rng = Lcg(31337);
    let mut list = IntervalList::new(dictionary.clone());
    for i in 0..500_000 {
        let (name, bases) = &genome[rng.below(genome.len() as u64) as usize];
        let start = 1 + rng.below(bases.len() as u64 - 2000) as i32;
        let end = start + rng.below(1000) as i32;
        let negative = rng.below(2) == 1;
        let label = format!("target_{i}");
        list.intervals.push(Interval::with_strand_and_name(
            name,
            start,
            end,
            negative,
            Some(&label),
        ));
    }
    let n = list.intervals.len();
    let render = |h: &SamHeader, l: &IntervalList| h.encode() + &l.write_body();

    let text = render(&header, &list);
    b.op(
        "interval_list_write",
        text.len(),
        n,
        "bytes",
        ".interval_list",
    );
    let out = b.path("interval_list_write", ".interval_list");
    b.time("interval_list_write", || {
        std::fs::write(&out, render(&header, &list)).expect("writes")
    });

    b.op(
        "interval_list_read",
        text.len(),
        n,
        "decoded",
        ".interval_list",
    );
    let (h, read) = b.time("interval_list_read", || {
        let text = std::fs::read_to_string(&out).expect("reads");
        let header_end: usize = text
            .lines()
            .take_while(|l| l.starts_with('@'))
            .map(|l| l.len() + 1)
            .sum();
        let h = parse_header_text(&text[..header_end]);
        let names = h.sequences.iter().map(|s| s.name.clone()).collect();
        let list = IntervalList::parse_body(names, &text[header_end..]).expect("parses");
        (h, list)
    });
    std::fs::write(
        b.path("interval_list_read", ".interval_list"),
        render(&h, &read),
    )
    .expect("writes");

    // `sorted()` sets SO:coordinate on htsjdk's header; `uniqued()` keeps the input's.
    let mut sorted_header = header.clone();
    sorted_header.attributes.set("SO", "coordinate");
    b.op(
        "interval_list_sort",
        text.len(),
        n,
        "decoded",
        ".interval_list",
    );
    let sorted = b.time("interval_list_sort", || list.sorted());
    std::fs::write(
        b.path("interval_list_sort", ".interval_list"),
        render(&sorted_header, &sorted),
    )
    .expect("writes");
    b.op(
        "interval_list_uniq",
        text.len(),
        n,
        "decoded",
        ".interval_list",
    );
    let uniqued = b.time("interval_list_uniq", || list.uniqued(true));
    std::fs::write(
        b.path("interval_list_uniq", ".interval_list"),
        render(&header, &uniqued),
    )
    .expect("writes");

    overlaps(b, &list);
    liftover(b, genome);
}

/// `OverlapDetector`: 50,000 targets loaded, 100,000 intervals queried against them. Sized down:
/// the port scans a contig's targets linearly (decision 0020) where htsjdk walks an interval tree.
fn overlaps(b: &Bench, list: &IntervalList) {
    let targets = &list.intervals[..50_000];
    let queries = &list.intervals[..100_000];
    b.op(
        "overlap_detector",
        0,
        targets.len() + queries.len(),
        "decoded",
        ".txt",
    );
    let counts = b.time("overlap_detector", || {
        let mut detector = OverlapDetector::new(0, 0);
        for (index, i) in targets.iter().enumerate() {
            detector.add(&i.contig, i.start, i.end, index);
        }
        queries
            .iter()
            .map(|i| detector.get_overlaps(&i.contig, i.start, i.end).len())
            .collect::<Vec<_>>()
    });
    let text: String = counts.iter().map(|c| format!("{c}\n")).collect();
    std::fs::write(b.path("overlap_detector", ".txt"), text).expect("writes");
}

/// A chain per contig onto a renamed build, blocks of 1-5 kb separated by small gaps on both sides,
/// and 200,000 intervals lifted through it.
fn liftover(b: &Bench, genome: &[(String, Vec<u8>)]) {
    let mut rng = Lcg(2718);
    let mut chain = String::new();
    for (id, (name, bases)) in genome.iter().enumerate() {
        let mut blocks = Vec::new();
        let (mut t, mut q) = (0u64, 0u64);
        let limit = bases.len() as u64 - 10_000;
        loop {
            let size = 1000 + rng.below(4000);
            if t + size > limit {
                blocks.push((size, None));
                t += size;
                q += size;
                break;
            }
            let (dt, dq) = (rng.below(50), rng.below(50));
            blocks.push((size, Some((dt, dq))));
            t += size + dt;
            q += size + dq;
        }
        writeln!(
            chain,
            "chain 1000 {name} {} + 0 {t} {name}_new {} + 0 {q} {}",
            bases.len(),
            q + 10_000,
            id + 1
        )
        .unwrap();
        for (size, gap) in blocks {
            match gap {
                Some((dt, dq)) => writeln!(chain, "{size}\t{dt}\t{dq}").unwrap(),
                None => writeln!(chain, "{size}\n").unwrap(),
            }
        }
    }
    std::fs::write(b.input("liftover.chain"), &chain).expect("writes");
    let queries = queries(genome, 200_000, 500);
    write_queries(&b.input("liftover_intervals.txt"), &queries);
    b.op("liftover", 0, queries.len(), "decoded", ".txt");
    let lifted = b.time("liftover", || {
        let lift = LiftOver::load(&chain).expect("a chain");
        queries
            .iter()
            .map(|(c, s, e)| lift.lift_over(&Interval::new(c, *s, *e)))
            .collect::<Vec<_>>()
    });
    let mut text = String::new();
    for l in lifted {
        match l {
            Some(i) => writeln!(text, "{}\t{}\t{}\t{}", i.contig, i.start, i.end, i.strand()),
            None => writeln!(text, "null"),
        }
        .unwrap();
    }
    std::fs::write(b.path("liftover", ".txt"), text).expect("writes");
}
