//! The record half of the I/O floor benchmark: BAM record decode and encode.
//!
//! Issue #78 names "BAM record decode, which every read filter and every walker pays per record"
//! as one of the four paths under everything, and it had no number. This is that number, in the
//! shape `tools/benchmark/run.sh` already uses: `name=value` lines, and a digest of what was
//! encoded so the runner can compare BYTES before it compares seconds.
//!
//! The corpus is synthetic and deterministic -- the same LCG the conformance harnesses use -- so a
//! run on one machine is comparable with a run on another, and the digests say whether an
//! optimisation moved a byte.
//!
//! Given a directory, it also writes the corpus there as a BAM file and times the whole-file paths
//! htsjdk has an answer for: reading every record of it (inflate and decode), and writing those
//! records to a new BAM (encode and deflate). `tools/benchmark/CodecBench.java` times htsjdk on the
//! same file, and the runner compares the two written files byte for byte before either duration.
//!
//! Usage: record-bench [records] [reps] [dir]

/// The allocator picard-rs links, so the rows measure the library as it ships there.
#[cfg(feature = "mimalloc")]
#[global_allocator]
static ALLOCATOR: mimalloc::MiMalloc = mimalloc::MiMalloc;

use std::io::{BufWriter, Write};
use std::path::Path;

use htsjdk_bam::cigar::{Cigar, CigarElement, Op};
use htsjdk_bam::header::{ReadGroup, SequenceRecord};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue, Tags};
use htsjdk_bam::writer::BamWriter;
use htsjdk_bam::{BamReader, SamHeader};
use htsjdk_bgzf::decompress_all;
use md5::{Digest, Md5};

fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

fn md5(bytes: &[u8]) -> String {
    let mut hasher = Md5::new();
    hasher.update(bytes);
    hex(&hasher.finalize())
}

/// The same 64-bit LCG the conformance harnesses use.
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0
    }

    fn below(&mut self, bound: u64) -> u64 {
        (self.next() >> 33) % bound
    }
}

/// One record of the shape a walker actually sees: a name, a cigar of a few elements, 151 bases
/// with qualities, and the tags an aligner leaves behind.
fn record(index: usize, rng: &mut Lcg) -> BamRecord {
    let read_len = 151usize;
    let clip = rng.below(8) as u32;
    let cigar = Cigar::new(vec![
        CigarElement {
            length: clip.max(1),
            op: Op::S,
        },
        CigarElement {
            length: read_len as u32 - clip.max(1),
            op: Op::M,
        },
    ]);
    let bases: Vec<u8> = (0..read_len)
        .map(|_| b"ACGTN"[rng.below(5) as usize])
        .collect();
    let quals: Vec<u8> = (0..read_len).map(|_| rng.below(41) as u8).collect();

    let mut tags = Tags::default();
    tags.insert(Tag::new(b"RG"), TagValue::Str("rg1".to_string()));
    tags.insert(Tag::new(b"NM"), TagValue::Int(rng.below(5) as i64));
    tags.insert(Tag::new(b"AS"), TagValue::Int(rng.below(151) as i64));
    tags.insert(
        Tag::new(b"MD"),
        TagValue::Str(format!("{}A{}", rng.below(100), rng.below(50))),
    );

    BamRecord {
        read_name: format!("HWI:1:FCX:1:{}:{}:{}", index % 8, index % 2048, index),
        flags: if index.is_multiple_of(2) { 99 } else { 147 },
        reference_index: 0,
        alignment_start: 1 + (index as i32 % 100_000),
        mapping_quality: 60,
        cigar,
        mate_reference_index: 0,
        mate_alignment_start: 1 + (index as i32 % 100_000),
        inferred_insert_size: 300,
        read_bases: bases,
        base_qualities: quals,
        tags,
    }
}

fn main() {
    let mut args = std::env::args().skip(1);
    let count: usize = args
        .next()
        .map(|a| a.parse().expect("records"))
        .unwrap_or(200_000);
    let reps: usize = args.next().map(|a| a.parse().expect("reps")).unwrap_or(3);
    let dir = args.next();

    let mut rng = Lcg(12345);
    let records: Vec<BamRecord> = (0..count).map(|i| record(i, &mut rng)).collect();

    // The encoded corpus, and its digest: an optimisation that changes a byte is a defect, and
    // this is where the run says so rather than where a suite says so an hour later.
    let mut encoded: Vec<u8> = Vec::new();
    for record in &records {
        encoded.extend_from_slice(&record.encode().expect("the corpus encodes"));
    }
    println!(
        "records={count} encoded_bytes={} md5={}",
        encoded.len(),
        md5(&encoded)
    );

    let megabytes = encoded.len() as f64 / (1024.0 * 1024.0);

    // Decode.
    let decoded = decode_all(&encoded);
    assert_eq!(decoded, count, "the corpus round trips");
    for run in 0..reps {
        let start = std::time::Instant::now();
        let n = decode_all(&encoded);
        let seconds = start.elapsed().as_secs_f64();
        println!(
            "rust_record_decode_run{run}_mbps={:.2} recs_per_sec={:.0}",
            megabytes / seconds,
            n as f64 / seconds
        );
    }

    // Encode.
    for run in 0..reps {
        let start = std::time::Instant::now();
        let mut out: Vec<u8> = Vec::with_capacity(encoded.len());
        for record in &records {
            record.encode_into(&mut out).expect("encodes");
        }
        let seconds = start.elapsed().as_secs_f64();
        assert_eq!(out.len(), encoded.len());
        println!(
            "rust_record_encode_run{run}_mbps={:.2} recs_per_sec={:.0}",
            megabytes / seconds,
            count as f64 / seconds
        );
        std::hint::black_box(out);
    }

    if let Some(dir) = dir {
        files(Path::new(&dir), &records, reps);
    }
}

/// Reads every record of a BAM file: the inflate and the decode, as a caller pays them.
fn read_bam(path: &Path) -> (SamHeader, Vec<BamRecord>, usize) {
    let data = decompress_all(&std::fs::read(path).expect("reads")).expect("inflates");
    let reader = BamReader::new(&data).expect("a BAM");
    let header = reader.header.text.clone();
    let records = reader.map(|r| r.expect("decodes")).collect();
    (header, records, data.len())
}

fn write_bam(path: &Path, header: &SamHeader, records: &[BamRecord]) {
    let file = std::fs::File::create(path).expect("creates");
    let mut writer = BamWriter::new(BufWriter::new(file), header).expect("writes");
    for record in records {
        writer.write(record).expect("encodes");
    }
    writer.finish().expect("finishes").flush().expect("flushes");
}

/// The file rows. Durations only: the runner divides the counts printed here, so the htsjdk column
/// is over the same bytes and records.
fn files(dir: &Path, records: &[BamRecord], reps: usize) {
    let mut header = SamHeader::default();
    header
        .sequences
        .push(SequenceRecord::new("chr1", 250_000_000));
    let mut group = ReadGroup::new("rg1");
    group.attributes.set("SM", "s1");
    header.read_groups.push(group);
    let corpus = dir.join("corpus.bam");
    write_bam(&corpus, &header, records);

    let (header, decoded, inflated) = read_bam(&corpus);
    println!(
        "bam_records={} bam_inflated_bytes={inflated}",
        decoded.len()
    );
    for run in 0..reps {
        let start = std::time::Instant::now();
        let answer = std::hint::black_box(read_bam(&corpus));
        let seconds = start.elapsed().as_secs_f64();
        // The answer is freed after the clock stops, as the JVM's collector frees htsjdk's later.
        drop(answer);
        println!("rust_bam_decode_run{run}_secs={:.6}", seconds);
    }

    let out = dir.join("port_out.bam");
    for run in 0..reps {
        let start = std::time::Instant::now();
        write_bam(&out, &header, &decoded);
        println!(
            "rust_bam_encode_run{run}_secs={:.6}",
            start.elapsed().as_secs_f64()
        );
    }
}

fn decode_all(bytes: &[u8]) -> usize {
    let mut offset = 0usize;
    let mut count = 0usize;
    while let Some((record, used)) = BamRecord::decode(&bytes[offset..]).expect("decodes") {
        std::hint::black_box(&record);
        offset += used;
        count += 1;
        if offset >= bytes.len() {
            break;
        }
    }
    count
}
