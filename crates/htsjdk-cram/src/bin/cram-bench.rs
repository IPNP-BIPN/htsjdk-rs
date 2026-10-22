//! CRAM and its codecs, on the inputs `tools/benchmark/CramBench.java` times htsjdk on.
//!
//! The corpus is aligned reads drawn from a synthetic reference with one substitution in a hundred,
//! which is what makes reference compression do work. It is written as a BAM for htsjdk to read
//! and as a FASTA with its `.fai` for both sides to compress against.
//!
//! CRAM write is checked byte for byte: both sides write CRAM 3.1 with htsjdk 5.0.0's NORMAL
//! profile, deflating GZIP blocks with libdeflate, 5.0.0's default (the container turns it off for
//! goldens, and `run.sh` turns it back on for `CramBench`). Both name the file `cram_write.cram`,
//! since the name is the file definition's 20-byte id. CRAM read is both sides decoding the port's
//! file, checked as SAM.
//!
//! Each codec is timed both ways over the corpus's own streams cut into 1 MiB blocks, the size a
//! CRAM slice's blocks run to. Compression is checked byte for byte; decompression decodes the
//! port's streams and is checked against the input. fqzcomp also takes each block's record lengths
//! and flags, which are written beside the blocks for `CramBench`.
//!
//! The output protocol is `bam-ops-bench`'s.
//!
//! Usage: cram-bench dir reps [records]

/// The allocator picard-rs links, so the rows measure the library as it ships there.
#[cfg(feature = "mimalloc")]
#[global_allocator]
static ALLOCATOR: mimalloc::MiMalloc = mimalloc::MiMalloc;

use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use htsjdk_bam::cigar::{Cigar, CigarElement, Op};
use htsjdk_bam::fasta_writer::FastaReferenceWriter;
use htsjdk_bam::header::{ReadGroup, SequenceRecord};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::sam_file::write_sam;
use htsjdk_bam::tag::{Tag, TagValue, Tags};
use htsjdk_bam::{BamWriter, SamHeader};
use htsjdk_cram::reader::{CramReader, FastaReference};
use htsjdk_cram::writer::write_cram;
use htsjdk_cram::{fqzcomp, name_tokeniser, range, rans, rans_nx16};

struct Bench {
    dir: PathBuf,
    reps: usize,
}

impl Bench {
    fn op(&self, name: &str, bytes: usize, records: usize, check: &str, suffixes: &str) {
        println!("op name={name} bytes={bytes} records={records} check={check} files={suffixes}");
    }

    fn time<T>(&self, name: &str, mut f: impl FnMut() -> T) -> T {
        let out = f();
        for run in 0..self.reps {
            let start = Instant::now();
            let answer = std::hint::black_box(f());
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

const CONTIG_LENGTH: usize = 2_000_000;
const READ_LENGTH: usize = 151;

/// Coordinate-sorted pairs off one contig: a few soft clips, one substitution in a hundred, and a
/// quality string that drifts the way an instrument's does.
fn corpus(count: usize, reference: &[u8]) -> Vec<BamRecord> {
    let mut rng = Lcg(5150);
    let mut start = 1usize;
    (0..count)
        .map(|index| {
            start += rng.below(20) as usize;
            if start + READ_LENGTH >= CONTIG_LENGTH {
                start = 1;
            }
            let clip = if rng.below(10) == 0 {
                1 + rng.below(10) as u32
            } else {
                0
            };
            let mut elements = Vec::new();
            if clip > 0 {
                elements.push(CigarElement {
                    length: clip,
                    op: Op::S,
                });
            }
            elements.push(CigarElement {
                length: READ_LENGTH as u32 - clip,
                op: Op::M,
            });
            let mut bases = reference[start - 1..start - 1 + READ_LENGTH].to_vec();
            let mut mismatches = 0;
            for base in bases.iter_mut() {
                if rng.below(100) == 0 {
                    *base = b"ACGT"[rng.below(4) as usize];
                    mismatches += 1;
                }
            }
            let mut q = 35i64;
            let quals: Vec<u8> = (0..READ_LENGTH)
                .map(|_| {
                    q = (q + rng.below(5) as i64 - 2).clamp(2, 41);
                    q as u8
                })
                .collect();
            let mut tags = Tags::default();
            tags.insert(Tag::new(b"RG"), TagValue::Str("rg1".to_string()));
            tags.insert(Tag::new(b"NM"), TagValue::Int(mismatches));
            let first = index.is_multiple_of(2);
            BamRecord {
                read_name: format!("HWI:1:FCX:1:{}:{}:{}", index % 8, index % 2048, index / 2),
                flags: if first { 99 } else { 147 },
                reference_index: 0,
                alignment_start: start as i32,
                mapping_quality: 60,
                cigar: Cigar::new(elements),
                mate_reference_index: 0,
                mate_alignment_start: start as i32,
                inferred_insert_size: if first { 151 } else { -151 },
                read_bases: bases,
                base_qualities: quals,
                tags,
            }
        })
        .collect()
}

fn main() {
    let mut args = std::env::args().skip(1);
    let dir = PathBuf::from(args.next().expect("dir"));
    let reps: usize = args.next().map(|a| a.parse().expect("reps")).unwrap_or(3);
    let count: usize = args
        .next()
        .map(|a| a.parse().expect("records"))
        .unwrap_or(200_000);
    let b = Bench { dir, reps };

    let mut rng = Lcg(8080);
    let reference: Vec<u8> = (0..CONTIG_LENGTH)
        .map(|_| b"ACGT"[rng.below(4) as usize])
        .collect();
    let fasta = b.dir.join("cram_ref.fasta");
    let mut writer = FastaReferenceWriter::new(60, true).expect("a writer");
    writer.start_sequence("chr1").expect("starts");
    writer.append_bases(&reference).expect("appends");
    let written = writer.close().expect("closes");
    std::fs::write(&fasta, written.fasta).expect("writes");
    std::fs::write(b.dir.join("cram_ref.fasta.fai"), written.index).expect("writes");

    let mut header = SamHeader::default();
    header.attributes.set("SO", "coordinate");
    header
        .sequences
        .push(SequenceRecord::new("chr1", CONTIG_LENGTH as i32));
    let mut group = ReadGroup::new("rg1");
    group.attributes.set("SM", "s1");
    header.read_groups.push(group);
    let records = corpus(count, &reference);
    write_bam(&b.dir.join("cram_input.bam"), &header, &records);

    cram(&b, &header, &records, &fasta);
    codecs(&b, &records);
}

fn write_bam(path: &Path, header: &SamHeader, records: &[BamRecord]) {
    let file = std::fs::File::create(path).expect("creates");
    let mut writer = BamWriter::new(BufWriter::new(file), header).expect("writes");
    for record in records {
        writer.write(record).expect("encodes");
    }
    writer.finish().expect("finishes").flush().expect("flushes");
}

fn read_cram(path: &Path, fasta: &Path) -> (SamHeader, Vec<BamRecord>) {
    let bytes = std::fs::read(path).expect("reads");
    let mut reader = CramReader::new(&bytes).expect("a CRAM");
    let records = reader
        .read_all(&mut FastaReference::open(fasta))
        .expect("decodes");
    (reader.header.clone(), records)
}

fn cram(b: &Bench, header: &SamHeader, records: &[BamRecord], fasta: &Path) {
    let size = write_sam(header, records).expect("renders").len();
    b.op("cram_write", size, records.len(), "bytes", ".cram");
    let out = b.path("cram_write", ".cram");
    b.time("cram_write", || {
        let mut reference = FastaReference::open(fasta);
        let cram =
            write_cram(header, records, "cram_write.cram", &mut reference).expect("compresses");
        std::fs::write(&out, cram).expect("writes");
    });

    b.op("cram_read", size, records.len(), "decoded", ".sam");
    let (h, decoded) = b.time("cram_read", || read_cram(&out, fasta));
    std::fs::write(
        b.path("cram_read", ".sam"),
        write_sam(&h, &decoded).expect("renders"),
    )
    .expect("writes");
}

const BLOCK: usize = 1 << 20;

/// A codec direction: the block's index and the block in, one block out.
type Transform<'a> = &'a dyn Fn(usize, &[u8]) -> Vec<u8>;

fn frame(parts: &[Vec<u8>]) -> Vec<u8> {
    let mut out = Vec::new();
    for p in parts {
        out.extend_from_slice(&(p.len() as u32).to_le_bytes());
        out.extend_from_slice(p);
    }
    out
}

/// One codec over a stream, both ways. `blocks` are what the codec is handed; the port's
/// compressed blocks are framed as `u32` length + bytes in `port_<codec>_compress.bin`, which is
/// what both sides decompress.
fn codec(b: &Bench, name: &str, blocks: &[Vec<u8>], compress: Transform, uncompress: Transform) {
    let size: usize = blocks.iter().map(Vec::len).sum();
    let op = format!("{name}_compress");
    b.op(&op, size, blocks.len(), "bytes", ".bin");
    let packed = b.time(&op, || {
        blocks
            .iter()
            .enumerate()
            .map(|(i, x)| compress(i, x))
            .collect::<Vec<_>>()
    });
    std::fs::write(b.path(&op, ".bin"), frame(&packed)).expect("writes");

    let op = format!("{name}_uncompress");
    b.op(&op, size, blocks.len(), "decoded", ".bin");
    let unpacked = b.time(&op, || {
        packed
            .iter()
            .enumerate()
            .map(|(i, x)| uncompress(i, x))
            .collect::<Vec<_>>()
    });
    assert_eq!(unpacked, blocks, "{name} round trips");
    std::fs::write(b.path(&op, ".bin"), frame(&unpacked)).expect("writes");
}

fn codecs(b: &Bench, records: &[BamRecord]) {
    let quals: Vec<u8> = records
        .iter()
        .flat_map(|r| r.base_qualities.iter().copied())
        .collect();
    let quality_blocks: Vec<Vec<u8>> = quals.chunks(BLOCK).map(<[u8]>::to_vec).collect();

    let rans0 = |_, x: &[u8]| rans::compress(x, rans::Order::Zero).expect("compresses");
    let rans1 = |_, x: &[u8]| rans::compress(x, rans::Order::One).expect("compresses");
    let unrans = |_, x: &[u8]| rans::uncompress(x).expect("decompresses");
    codec(b, "rans4x8_o0", &quality_blocks, &rans0, &unrans);
    codec(b, "rans4x8_o1", &quality_blocks, &rans1, &unrans);

    let nx0 = |_, x: &[u8]| rans_nx16::compress(x, 0);
    let nx1 = |_, x: &[u8]| rans_nx16::compress(x, 1);
    let unnx = |_, x: &[u8]| rans_nx16::uncompress(x).expect("decompresses");
    codec(b, "ransnx16_o0", &quality_blocks, &nx0, &unnx);
    codec(b, "ransnx16_o1", &quality_blocks, &nx1, &unnx);

    let range0 = |_, x: &[u8]| range::compress(x, 0).expect("compresses");
    let range1 = |_, x: &[u8]| range::compress(x, 1).expect("compresses");
    let unrange = |_, x: &[u8]| range::uncompress(x).expect("decompresses");
    codec(b, "range_o0", &quality_blocks, &range0, &unrange);
    codec(b, "range_o1", &quality_blocks, &range1, &unrange);

    // Read names, NUL-terminated as a CRAM block holds them, ten thousand to a block.
    let name_blocks: Vec<Vec<u8>> = records
        .chunks(10_000)
        .map(|chunk| {
            let mut block = Vec::new();
            for r in chunk {
                block.extend_from_slice(r.read_name.as_bytes());
                block.push(0);
            }
            block
        })
        .collect();
    let names = |_, x: &[u8]| name_tokeniser::compress(x, 0).expect("compresses");
    let unnames = |_, x: &[u8]| name_tokeniser::uncompress(x, 0).expect("decompresses");
    codec(b, "name_tok", &name_blocks, &names, &unnames);

    // fqzcomp, per block of whole records: their lengths and flags are part of its input, and go
    // to `CramBench` as little-endian `i32` arrays framed like the blocks.
    let chunks: Vec<&[BamRecord]> = records.chunks(BLOCK / READ_LENGTH).collect();
    let fqz_blocks: Vec<Vec<u8>> = chunks
        .iter()
        .map(|c| {
            c.iter()
                .flat_map(|r| r.base_qualities.iter().copied())
                .collect()
        })
        .collect();
    let ints = |f: &dyn Fn(&BamRecord) -> i32| -> Vec<Vec<i32>> {
        chunks.iter().map(|c| c.iter().map(f).collect()).collect()
    };
    let lengths = ints(&|r| r.base_qualities.len() as i32);
    let flags = ints(&|r| r.flags as i32);
    for (name, values) in [("lengths", &lengths), ("flags", &flags)] {
        let parts: Vec<Vec<u8>> = values
            .iter()
            .map(|v| v.iter().flat_map(|x| x.to_le_bytes()).collect())
            .collect();
        std::fs::write(b.dir.join(format!("fqzcomp_{name}.bin")), frame(&parts)).expect("writes");
    }
    let fqz = |i: usize, x: &[u8]| fqzcomp::compress(x, &lengths[i], Some(&flags[i]));
    let unfqz = |_, x: &[u8]| fqzcomp::uncompress(x).expect("decompresses");
    codec(b, "fqzcomp", &fqz_blocks, &fqz, &unfqz);
}
