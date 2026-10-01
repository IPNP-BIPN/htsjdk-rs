# Changelog

All notable changes to this repository, newest first. Versions are the workspace's: every crate
here moves together, because they are one port of one library and a consumer pins the set.

The format is [Keep a Changelog](https://keepachangelog.com/en/1.1.0/); this project does not yet
follow semantic versioning, because 0.x is what an unreleased port is.

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
