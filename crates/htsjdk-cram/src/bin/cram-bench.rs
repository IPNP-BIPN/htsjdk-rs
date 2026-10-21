//! CRAM and its codecs, on the inputs `tools/benchmark/CramBench.java` times htsjdk on.
//!
//! The corpus is aligned reads drawn from a synthetic reference with one substitution in a hundred,
//! which is what makes reference compression do work. It is written as a BAM for htsjdk to read
//! and as a FASTA with its `.fai` for both sides to compress against.
//!
//! CRAM write is checked on decoded records, not bytes: this port writes htsjdk 5.0.0's CRAM 3.1
//! and the container's htsjdk 4.2.0 writes CRAM 3.0, so each side decodes its own file and the two
//! renderings as SAM must be identical. htsjdk 4.2.0 cannot read CRAM 3.1 ("CRAM version 3.1 is not
//! supported"), so CRAM read is both sides decoding htsjdk's CRAM 3.0 file, in a second invocation
//! with `read` that runs after `CramBench.java` has written it. In both checks the port's records
//! lose the `MD` tag its reader generated: htsjdk 5.0.0's reader restores `MD` and `NM` where a
//! record lacks them, as htslib does, and 4.2.0's does not. The corpus carries its own `NM`.
//!
//! Each codec is timed both ways over the corpus's own streams cut into 1 MiB blocks, the size a
//! CRAM slice's blocks run to. Compression is checked byte for byte; decompression decodes the
//! port's streams and is checked against the input. fqzcomp is decompression only: htsjdk 4.2.0's
//! `FQZCompEncode.compress` throws "FQZComp compression is not implemented".
//!
//! The output protocol is `bam-ops-bench`'s.
//!
//! Usage: cram-bench dir reps [records]   then   cram-bench dir reps records read

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
    if args.next().as_deref() == Some("read") {
        let fasta = b.dir.join("cram_ref.fasta");
        let size = std::fs::metadata(b.dir.join("port_cram_write.sam"))
            .expect("exists")
            .len();
        b.op("cram_read", size as usize, count, "decoded", ".sam");
        let input = b.dir.join("java_cram_write.cram");
        let (h, decoded) = b.time("cram_read", || read_cram(&input, &fasta));
        write_check(&b.path("cram_read", ".sam"), &h, decoded);
        return;
    }

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
    b.op("cram_write", size, records.len(), "decoded", ".sam");
    let out = b.path("cram_write", ".cram");
    b.time("cram_write", || {
        let mut reference = FastaReference::open(fasta);
        let cram = write_cram(header, records, "port_cram_write.cram", &mut reference)
            .expect("compresses");
        std::fs::write(&out, cram).expect("writes");
    });
    let (h, decoded) = read_cram(&out, fasta);
    write_check(&b.path("cram_write", ".sam"), &h, decoded);
}

/// The decoded records as SAM, less the `MD` the port's reader generated (see the module doc).
fn write_check(path: &Path, header: &SamHeader, mut records: Vec<BamRecord>) {
    for r in &mut records {
        r.tags.remove(Tag::new(b"MD"));
    }
    std::fs::write(path, write_sam(header, &records).expect("renders")).expect("writes");
}

const BLOCK: usize = 1 << 20;

/// A codec direction: one block in, one block out.
type Transform<'a> = &'a dyn Fn(&[u8]) -> Vec<u8>;

/// One codec over a stream, both ways. `blocks` are what the codec is handed; the port's
/// compressed blocks are framed as `u32` length + bytes in `port_<codec>_compress.bin`, which is
/// what both sides decompress. `same_encoder` is false where the port follows a later htsjdk
/// encoder than the container's, whose bytes legitimately differ.
fn codec(
    b: &Bench,
    name: &str,
    same_encoder: bool,
    blocks: &[Vec<u8>],
    compress: Option<Transform>,
    uncompress: Transform,
    compressed: Option<Vec<Vec<u8>>>,
) {
    let size: usize = blocks.iter().map(Vec::len).sum();
    let frame = |parts: &[Vec<u8>]| {
        let mut out = Vec::new();
        for p in parts {
            out.extend_from_slice(&(p.len() as u32).to_le_bytes());
            out.extend_from_slice(p);
        }
        out
    };
    let packed = match compress {
        Some(compress) => {
            let op = format!("{name}_compress");
            if same_encoder {
                b.op(&op, size, blocks.len(), "bytes", ".bin");
            } else {
                // The check is htsjdk decoding these streams back to the input (`.dec`).
                b.op(&op, size, blocks.len(), "decoded", ".dec");
                std::fs::write(b.path(&op, ".dec"), frame(blocks)).expect("writes");
            }
            let packed = b.time(&op, || {
                blocks.iter().map(|x| compress(x)).collect::<Vec<_>>()
            });
            std::fs::write(b.path(&op, ".bin"), frame(&packed)).expect("writes");
            packed
        }
        None => {
            let packed = compressed.expect("a stream to decompress");
            std::fs::write(b.path(&format!("{name}_compress"), ".bin"), frame(&packed))
                .expect("writes");
            packed
        }
    };
    let op = format!("{name}_uncompress");
    b.op(&op, size, blocks.len(), "decoded", ".bin");
    let unpacked = b.time(&op, || {
        packed.iter().map(|x| uncompress(x)).collect::<Vec<_>>()
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

    let rans0 = |x: &[u8]| rans::compress(x, rans::Order::Zero).expect("compresses");
    let rans1 = |x: &[u8]| rans::compress(x, rans::Order::One).expect("compresses");
    let unrans = |x: &[u8]| rans::uncompress(x).expect("decompresses");
    codec(
        b,
        "rans4x8_o0",
        true,
        &quality_blocks,
        Some(&rans0),
        &unrans,
        None,
    );
    codec(
        b,
        "rans4x8_o1",
        true,
        &quality_blocks,
        Some(&rans1),
        &unrans,
        None,
    );

    let nx0 = |x: &[u8]| rans_nx16::compress(x, 0);
    let nx1 = |x: &[u8]| rans_nx16::compress(x, 1);
    let unnx = |x: &[u8]| rans_nx16::uncompress(x).expect("decompresses");
    codec(
        b,
        "ransnx16_o0",
        true,
        &quality_blocks,
        Some(&nx0),
        &unnx,
        None,
    );
    codec(
        b,
        "ransnx16_o1",
        true,
        &quality_blocks,
        Some(&nx1),
        &unnx,
        None,
    );

    let range0 = |x: &[u8]| range::compress(x, 0).expect("compresses");
    let range1 = |x: &[u8]| range::compress(x, 1).expect("compresses");
    let unrange = |x: &[u8]| range::uncompress(x).expect("decompresses");
    codec(
        b,
        "range_o0",
        true,
        &quality_blocks,
        Some(&range0),
        &unrange,
        None,
    );
    codec(
        b,
        "range_o1",
        true,
        &quality_blocks,
        Some(&range1),
        &unrange,
        None,
    );

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
    // htsjdk 5.0.0's encoder, which codes names differently from the container's 4.2.0.
    let names = |x: &[u8]| name_tokeniser::compress(x, 0).expect("compresses");
    let unnames = |x: &[u8]| name_tokeniser::uncompress(x, 0).expect("decompresses");
    codec(
        b,
        "name_tok",
        false,
        &name_blocks,
        Some(&names),
        &unnames,
        None,
    );

    // fqzcomp, per block of whole records: the record lengths are part of its input.
    let per_block = BLOCK / READ_LENGTH;
    let mut fqz_blocks = Vec::new();
    let mut fqz_packed = Vec::new();
    for chunk in records.chunks(per_block) {
        let data: Vec<u8> = chunk
            .iter()
            .flat_map(|r| r.base_qualities.iter().copied())
            .collect();
        let lengths: Vec<i32> = chunk
            .iter()
            .map(|r| r.base_qualities.len() as i32)
            .collect();
        let flags: Vec<i32> = chunk.iter().map(|r| r.flags as i32).collect();
        fqz_packed.push(fqzcomp::compress(&data, &lengths, Some(&flags)));
        fqz_blocks.push(data);
    }
    let unfqz = |x: &[u8]| fqzcomp::uncompress(x).expect("decompresses");
    codec(
        b,
        "fqzcomp",
        true,
        &fqz_blocks,
        None,
        &unfqz,
        Some(fqz_packed),
    );
}
