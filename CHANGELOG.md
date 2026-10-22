# Changelog

All notable changes to this repository, newest first. Versions are the workspace's: every crate
here moves together, because they are one port of one library and a consumer pins the set.

The format is [Keep a Changelog](https://keepachangelog.com/en/1.1.0/); this project does not yet
follow semantic versioning, because 0.x is what an unreleased port is.

## [Unreleased]

The reference is htsjdk 5.0.0 (decision 0047). The oracle moved first and 18 of 87 suites moved
with it; each difference is ported and its golden re-frozen from the x86-64 runner.

### Changed

- **The oracle is htsjdk 5.0.0**, with its runtime dependencies pinned by sha256 and libdeflate
  off image-wide, so every suite runs against 5.0.0 (#271).
- **CRAM data series are on htslib's content ids**, 11 to 46, where 4.2.0 numbered them 1 to 32.
- **A tag's compressor is a trial of GZIP against rANS Nx16 order 0**, the strictly smaller
  winning, and every `Z` tag is tab-stopped. `negotiation::trial_choice` replaces
  `best_compressor`, and `Compressor::Rans` is `Compressor::RansNx16`.
- **A CRAM template length is htslib's**, the span of both ends: a pair at 100 and 200 with
  50-base reads is 150, not 101. `record_flags::template_length` is the shared rule.
- `ByteArrayStopCodec::read` returns a `Result`.

### Fixed

- **An unterminated string tag is refused as 5.0.0 refuses it**:
  `SAMFormatException: Null-terminated string tag value is not null terminated.`, where 4.2.0
  ran off the buffer with a bare `BufferUnderflowException`.
- **Reading past the end of an external block is `IllegalStateException`** naming the bytes
  asked for and left, and a stop codec whose block ends before its stop byte is a `CRAMException`
  naming it, where 4.2.0 returned what it had.

## [0.4.0] - 2026-10-22

Speed, measured against htsjdk itself. A benchmark now times 44 operations against htsjdk 4.2.0,
each only after both sides wrote the same output, and the rows htsjdk was ahead on were made faster
without moving a byte. picard-rs v0.2.0's speed-ups rest on these.

### Added

- `BgzfWriter` deflates full blocks on every core (rayon) and writes them in order, with the
  serial writer's bytes for every deflater and level: a block's payload depends only on its own
  bytes, the level and the deflater, and the block boundaries do not move. `parallel(false)` turns
  it off (#249, #257).
- `reheader::block_copy_to` and `gather_bam_files_to`: the block copy and the gather from seekable
  inputs to a writer, in constant memory, with the bytes of `block_copy_with` and
  `gather_bam_files_with`, which now run through them. `ReheaderError::Write` is a failed write.
  picard-rs's GatherBamFiles goes from 703 MB to 4 MB on the benchmark BAM (#258).
- `reheader::reheader_bam_to`, the same for a reheader; `reheader_bam_with` copies once instead of
  twice, and block copies move 1 MiB at a time (#269).
- **A benchmark of every operation against htsjdk 4.2.0, outputs first** (#266, #268, #269): BGZF,
  BAM, VCF and 44 more rows across every crate, run in the pinned container on x86-64. `run.sh`
  fails on any output difference before it prints a rate. The README's Throughput section is
  rebuilt from one run. The bench binaries allocate through mimalloc, as picard-rs does, through
  an opt-in `mimalloc` feature on each crate; the libraries do not choose an allocator.

### Changed

- `BgzfWriter::file_pointer` takes `&mut self` and returns `io::Result<u64>`: it writes the blocks
  in flight before answering, so the offset stays exact. `BamWriter` asks only when indexing.
- **Faster, same bytes.** Each change carries a test against the code it replaced. Ratios are
  against htsjdk on one x86-64 run (#269), before and after:
  - `decompress_all` walks the block headers, then inflates the blocks on every core: BGZF inflate
    0.96x to 2.94x (text, level 5), BAM read 0.74x to 2.20x, FASTQ read 0.70x to 3.19x.
  - `OverlapDetector` answers from a per-contig index sorted by start, built lazily, in insertion
    order: 0.04x to 3.97x.
  - `IntervalList::sorted` looks each contig up once and `uniqued` moves rather than clones: sort
    0.58x to 3.94x, uniq 0.63x to 4.16x. `IndexedFasta::query` reads the span in one read (0.33x
    to 0.86x), and `calculate_md_and_nm` writes its counts in place (0.64x to 1.84x).
  - `gkl-deflate` hashes with a slicing-by-4 CRC32C and compares matches eight bytes at a time:
    10% fewer cycles for picard-rs's CleanSam, below Picard's own (#265).
  - A base maps to its nibble through a table, which doubles BAM record encoding (#263). A VCF
    genotype block is read without copying it (#264). The phred error table is built once, which
    was most of the fingerprint tools' CPU (#262).
  - VCF header lines are deduplicated with a set: 4.34 s to 35 ms on a header with 3,400 contig
    lines (#254). BCF2 looks header lines up by ID once per header (#253).
- `noodles-bam` 0.97.0, `noodles-sam` 0.92.0 and `noodles-vcf` 0.95.0, used by the test oracles
  only; `libz-sys` 1.1.30, `libdeflate-sys` 1.26.1 and `sha2` 0.11.1.

### Fixed

- **A BAM record nothing was set on is written from the bytes it was read from**, as htsjdk's
  `BAMRecord` does: the tags keep the file's order and integer types. picard-rs's SortSam on a
  samtools BAM differed from Picard in 643 of 2,994 records before this (#256).
- **`COMPRESSION_LEVEL=0` under GKL stores each block**, as Picard does, instead of panicking (#255).
- **A record whose alleles repeat is refused**, as `VariantContext.makeAlleles` does, with
  htsjdk's message (#267).
- `gkl-deflate` no longer panics on overflow in a debug build when the window slides with the last
  match below it; release builds already wrapped, so no byte changes.

## [0.3.0] - 2026-10-17

Java's number printing, reconstructed from behaviour, closes the last formatting gap. CRAM 3.1,
BCF2 and libdeflate bring the formats and the deflater Picard 3.5.0 writes; the rest are VCF bytes
and a speed bound that gatk-rs's tools reached.

### Added

- `jmath::to_string`: `Double.toString` and `Float.toString` as Java 17 prints them, including
  where they are not the shortest decimal. Reconstructed from the oracle's measured behaviour, not
  from the JDK's GPL-2.0 source (decision 0046).
- **CRAM 3.1, read and written as htsjdk 5.0.0 does** (#252): the rANS Nx16, arithmetic,
  fqzcomp and name-tokeniser codecs, `reader::CramReader` with htsjdk's slice normalisation,
  `writer::write_cram` (NORMAL profile, per-tag GZIP/rANS trial) and the `.bai` of a CRAM.
  Writing a bzip2 block is not ported, so an `SA` or `XA` tag that would need one is refused.
- **`htsjdk_vcf::bcf2`**: htsjdk's BCF2 codec and writer, BCF 2.1 uncompressed, with INFO in
  Java `HashMap` order and the oracle's refusals (#251).
- **`Deflater::Libdeflate`**: libdeflate v1.25, which htsjdk 5.0.0 uses whenever jlibdeflate
  loads, so `USE_JDK_DEFLATER=true` in Picard 3.5.0 is libdeflate's bytes (#248).
- `_with` variants of `write_bam_header_block`, `block_copy`, `reheader_bam`, `gather_bam_files`
  and `TabixIndex::write` that take the level and the deflater; the existing functions keep
  their signatures (#248).

### Fixed

- `FormatUtil` formatting now matches on all 41,678 corpus values (was 99.84%, 66 quarantined),
  the SAM text writer on all 67 records (`Float.MIN_VALUE` is `1.4E-45`), and `%f` / `%e` on all
  127,803 swept doubles (was 99.85%). Nothing formatting-related is licence-blocked any more.
- **A written GQ is capped at 99**, as `VCFEncoder` does through `GQAccessor`. The genotype keeps
  its real GQ, which GATK's `GQ_MEAN` reads; only the text is capped.
- **A file whose sample columns are not sorted has every record decoded**, as
  `AbstractVCFCodec.parseVCFLine` does when `samplesWereAlreadySorted()` is false, so it is written
  back re-encoded rather than copied verbatim. The order is `String.compareTo`'s (#242).
- **A Tribble index's double properties are written as `Double.toString` writes them**: a
  feature-length variance below a thousandth is `3.111081185151864E-4`, not `0.000311...` (#250).

### Changed

- `jmath` `log` and `log10` answer from a table, an exact reduction and a short polynomial whenever
  twice the error bound rounds to one double, and fall back to the double-double series otherwise:
  277 ns per call becomes 11.6, with 200 million inputs checked and none differing (decision 0045).
- `noodles-bam` 0.96.0, `noodles-sam` 0.91.0, `noodles-vcf` 0.94.0 and `noodles-core` 0.21.0, used
  by the test oracles only.

## [0.2.0] - 2026-10-01

Fixes found by gatk-rs's covering arrays, which run whole tools against the reference and so reach
corners of the library its own suites did not. Each one is a byte the reference writes and this
port did not.

### The claim

- **86 conformance suites, all oracle-backed**, unchanged in number: the fixes below are each
  pinned by a unit test or a frozen refusal rather than a new suite.

### Added

- `jmath::strict_log::log`: `StrictMath.log`, transcribed from FDLIBM 5.3's `e_log.c` in the C
  function's operation order. It is not `Math.log`, which is correctly rounded: the two differ on
  186 of the corpus's 44,996 points, and the test asserts both that this port matches the
  `StrictMath` column on all of them and that the two columns differ somewhere. gatk-rs needs it
  for `Random.nextGaussian`, which `QualByDepth` uses to jitter every QD above 35 (decision 0044).

### Fixed

- **GKL levels 1 and 2 reach ISA-L as themselves.** Decision 0031 had sent both to ISA-L level 1,
  read off four fixtures that happen to compress identically at both levels. A 15 KB VCF header
  written at GATKConfig's level 2 does not: `IntelDeflater` gives 5303 bytes, which is ISA-L level
  2, where level 1 gives 5321. Decision 0043 records the correction.
- **A sample held twice in one record answers with its later genotype**, as htsjdk's
  `sampleNameToOffset` does. CombineGVCFs produces that shape when a block left open on one contig
  is merged with the same sample's block on the next.
- **A QUAL of zero is written as `0`, not `-0`.** `getPhredScaledQual` adds `0.0` to the product,
  which turns `-0.0` into `0.0`; the port returned the negative zero and the encoder wrote it.
- **A FASTA reference closed with no sequence in it is refused** with htsjdk's
  `IllegalStateException`, after the three files are written, as `FastaReferenceWriter.close()`
  does. The port answered three empty files and success.

### Changed

- `noodles-bam` 0.95.0, `noodles-sam` 0.90.0 and `noodles-vcf` 0.93.0, used by the test oracles
  only.

## [0.1.0] - 2026-09-06

The first version worth naming: every format the consumers reach is ported, and every claim about
those bytes is re-derived by CI in a pinned container on a real x86-64 runner.

### The claim

- **86 conformance suites, all oracle-backed.** No committed golden is unchecked: CI regenerates
  each one on every push and compares it against the reference build, in the digest-pinned
  `linux/amd64` container on JDK 17.
- Formats reproduced byte for byte: BGZF, the BAM record codec and file writer, SAM text, the BAI
  and Tribble and Tabix indexes, interval lists, VCF headers, records and whole files, the CRAM
  container model with its codecs and CRAI, and Picard's metrics file layout and number formatting.
- Where a value cannot be matched exactly it is quarantined and reported with its measured
  divergence rate, and the output is called **bio-identical** rather than bit-identical. The two
  cases are `Math.exp` (decisions 0005, 0014, 0025) and the last 112 of 41,678 metrics doubles,
  which are licence-blocked (decisions 0013 and 0040).

### Added

- `htsjdk_vcf::genotypes_context::GenotypesContext`: htsjdk's `LazyGenotypesContext`, so a record
  nobody looked at is written from the file's own text, FORMAT column included. A *read* is what
  drops that text, which is the part a port gets wrong by assuming only a mutation matters.
- `htsjdk_bam::fasta_index`: the `.fai` parser and the indexed FASTA reader, whose query arithmetic
  is `getSubsequenceAt`'s: a seek computed from the index rather than a scan for newlines.
- `record-bench` and `vcf-bench` beside `bgzf-bench`: the I/O floor's three paths, each printing
  the digest of what it encoded so a speed change that moves a byte says so in the run.

### Changed

- The BGZF reader and writer keep one buffer and one zlib state per stream rather than per block,
  which is what htsjdk does. Inflate went from 933 MB/s to 3527 MB/s locally, and 1.4x to 3.2x
  against htsjdk on a real x86-64 runner.
- The BAM record codec stopped cloning the tags and the cigar of every record for the sake of a
  branch almost nothing takes. Decode 680 → 1185 MB/s, encode 175 → 262 MB/s, same bytes.
- The VCF line splitter borrows instead of building a `Vec<char>` and a `String` per column, and
  the genotype parser no longer takes the site columns it reads only when refusing. Reading went
  from 95k to 330k records/s.
- `gatk-engine` reads `.bai` and FASTA through this repository rather than through `noodles`, which
  leaves `noodles`, `rayon` and `crossbeam` out of its dependency tree entirely (decision 0041).

### Fixed

- A record read from one file and written under another file's header copied its genotype columns
  verbatim, putting one sample's genotype in another sample's column. The lazy text is now only an
  answer under the sample order it was written in.

### Known limitations

- The reference version is htsjdk 4.2.0 and stays there until the tool ports are done: all three
  repositories name one set of pins (decision 0042).
- Levels 1 and 2 of GKL's deflater need a linked ISA-L, and refuse rather than substituting zlib's
  bytes when the build cannot reach it (decision 0034).
