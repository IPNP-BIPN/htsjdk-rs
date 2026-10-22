# htsjdk-rs

Pure-Rust port of [htsjdk](https://github.com/samtools/htsjdk), targeting **byte-identical**
output against a pinned reference build. Work in progress.

> **This is not the official htsjdk.** It is an independent reimplementation, not affiliated
> with or endorsed by the Broad Institute or the samtools organization.

## Why this exists, and why not `rust-htslib`

[`rust-htslib`](https://github.com/rust-bio/rust-htslib) provides FFI bindings to **HTSlib**,
the **C** library from samtools. htsjdk is the **Java** library, an independent implementation
of the same file formats. They are not two names for one codebase.

That distinction is the whole reason this crate exists. BGZF leaves block size, deflate level,
and block-fill policy to the implementation, so byte equality between HTSlib and htsjdk output
is guaranteed by nothing. GATK and Picard are built on htsjdk, so reproducing their output
byte-for-byte requires porting **htsjdk specifically**.

`rust-htslib` remains excellent for reading, and this project uses it as an independent
cross-check oracle.

## Reference version

Measured against htsjdk `5.0.0`, the version GATK 4.7.0.0 and Picard 3.5.0 pin. Most modules were
transcribed from 4.2.0 and say so; every suite now runs against 5.0.0, and the 18 behaviours that
moved are ported
([decision 0047](docs/decisions/0047-the-reference-moves-to-htsjdk-5.0.0.md), which supersedes the
pin of [0042](docs/decisions/0042-the-reference-version-is-pinned-by-the-consumers-not-by-this-repository.md)).
gatk-rs, still at GATK 4.6.2.0, pins this repository by revision and keeps 4.2.0's behaviour
until it moves that pin.

## Scope

| Module | Status |
|---|---|
| BGZF reader/writer (`BlockCompressedOutputStream` semantics) | **byte-identical**, 35 goldens |
| BAM record codec (byte level, incl. `bin` and tag type promotion) | **byte-identical**, 136 goldens |
| SAM text header (`SAMTextHeaderCodec.encode`) | **byte-identical**, 9 goldens |
| BAM file writer (header + dictionary + framing) | **byte-identical**, 4 whole files |
| BAM file reader (header parse + records) | reads htsjdk-produced files |
| SAM text record writer | **byte-identical**, 67 of 67, decision 0046 |
| SAM text record parser | BAM → SAM → BAM byte-identical |
| Whole SAM text files | **byte-identical**, 5 goldens; BAM → SAM matches htsjdk |
| `Histogram` (20 of 44 metrics tools) | **byte-identical**, 338 statistics |
| Metrics number formatting (`FormatUtil`) | **100%** over 41,678 doubles re-derived by the oracle on every run; Java 17's digits reconstructed from behaviour, decisions 0040 and 0046 |
| `MetricsFile` layout | ported (`htsjdk-metrics::file`); no suite of its own, exercised end to end by every picard-rs metrics golden |
| BAM index (`.bai`) | **byte-identical**, 9 goldens |
| Interval lists (`IntervalList`, `.interval_list`) | **byte-identical**, 11 goldens |
| VCF header | **byte-identical**, 7 goldens, decision 0016 |
| VCF records (`VCFEncoder`, `Allele`, genotypes) | **byte-identical**, 54 records + 29 double cases, decisions 0017, 0018 |
| JVM `%f` / `%e` number formatting | **100%** over 127,803 doubles, decisions 0017 and 0046 |
| Whole VCF files (`VariantContextWriter`) | **byte-identical**, 7 goldens, decision 0019 |
| Tribble index | **byte-identical**, both layouts, read and written: 5 suites (`tribble-index`, `tribble-index-write`, `tribble-bed`, `tribble-interval-list`, `tribble-vcf-candecode`) |
| Tabix index (`.tbi`) | **byte-identical**, written: 19 feature streams as body and as block-compressed file (`tabix-index`) |
| CRAM (container model, codecs, negotiation, CRAI) | **byte-identical**, 27 suites; CRAM 3.1 is the extension surface of decision 0039 |

## Throughput

Measured, not claimed, and only where the bytes are pinned first: `tools/benchmark/run.sh` times
every operation the port has an htsjdk counterpart for against htsjdk's, on the same inputs, in
the same pinned container, and asserts that every stream and file the two write is identical
**before** it prints a duration. A speed number for a path with no golden is not a measurement of this port
([#78](https://github.com/IPNP-BIPN/htsjdk-rs/issues/78)).

```sh
tools/benchmark/run.sh 64 3
```

The figures come from the `benchmark` job, which runs it on a real x86-64 runner; a run on a
developer machine goes through an emulated container and measures the emulator. Both sides are
library calls, not command-line tools: htsjdk 5.0.0 on JDK 17, the port at the same commit,
allocating through mimalloc as picard-rs does (each crate's `mimalloc` feature sets it as the
benches' global allocator; the libraries themselves leave the allocator to the binary).

htsjdk 5.0.0 deflates with libdeflate by default, and the oracle image turns that off so every
golden records the deflater its dump chose. The BGZF, BAM and VCF rows below therefore compare
zlib with zlib, the port's `Deflater::Jdk`; the CRAM rows run htsjdk with libdeflate on, its
default, because the port's CRAM writer deflates its GZIP blocks as 5.0.0 does by default.

Every table below is three runs of the same commit, [39760001813](https://github.com/IPNP-BIPN/htsjdk-rs/actions/runs/39760001813)
(AMD EPYC 9V45), [39760014679](https://github.com/IPNP-BIPN/htsjdk-rs/actions/runs/39760014679) (EPYC 7763)
and [39760028056](https://github.com/IPNP-BIPN/htsjdk-rs/actions/runs/39760028056) (EPYC 9V74), each a
median of three repetitions. GitHub's runners are not one machine, and one run on one of them can
mislead by a factor of two on the rows that allocate most, so each row shows the run whose ratio is
the median of the three, its rates, and the range of the ratio over all three. Each answer is freed
after its clock stops, as the JVM's collector frees htsjdk's outside the timed region. The VCF and
BCF readers run deferred, as htsjdk's do: a record's genotype columns are parsed when something first
looks at them (decision 0048).

### BGZF

64 MB payloads, **after** the run asserted that every framed stream the two sides wrote was
byte-identical on all ten keys:

| operation | payload | level | htsjdk | port | ratio | range |
|---|---|---:|---:|---:|---:|---:|
| deflate | text | 1 | 748 MB/s | 1,624 MB/s | 2.17x | 2.09 to 2.23 |
| deflate | text | 5 | 460 MB/s | 1,238 MB/s | 2.69x | 2.55 to 3.21 |
| deflate | random | 5 | 43.9 MB/s | 131 MB/s | 2.98x | 2.96 to 3.01 |
| inflate | text | 1 | 1,667 MB/s | 15,528 MB/s | 9.32x | 8.75 to 14.70 |
| inflate | text | 5 | 2,206 MB/s | 8,258 MB/s | 3.74x | 3.48 to 3.87 |
| inflate | random | 5 | 390 MB/s | 1,463 MB/s | 3.75x | 3.39 to 3.81 |

Deflate is stock zlib on both sides at every level here, so the bytes match by construction. The
`gkl-deflate` crate, a Rust port of Intel GKL's deflater that GATK installs in htsjdk's place, is
not on this path. Both ratios are parallelism, not faster code: the port's BGZF writer deflates
full blocks on every core, and `decompress_all` walks the block headers once and then inflates the
blocks on every core, each into its own slice of the output. htsjdk does both on the calling
thread. The bytes, and for a damaged stream the error, are the serial ones.

### BAM and VCF files

400,000 BAM records of 151 bases (122 MB inflated) and 200,000 VCF records with four samples
(23 MB), **after** the run asserted that the BAM and both VCFs the two sides wrote were
byte-identical. Rates are over the same bytes and records for both columns.

| path | htsjdk | port | ratio | range |
|---|---:|---:|---:|---:|
| BAM read, every record decoded | 223 MB/s, 730,000 records/s | 444 MB/s, 1,455,000 records/s | 1.99x | 1.99 to 2.00 |
| BAM write | 28 MB/s, 90,000 records/s | 84 MB/s, 276,000 records/s | 3.07x | 3.01 to 3.52 |
| VCF read, genotypes left lazy | 46 MB/s, 390,000 records/s | 112 MB/s, 956,000 records/s | 2.45x | 2.13 to 3.24 |
| VCF read, genotypes decoded | 27 MB/s, 232,000 records/s | 53 MB/s, 452,000 records/s | 1.95x | 1.78 to 2.36 |
| VCF write, records copied through | 74 MB/s, 633,000 records/s | 290 MB/s, 2,487,000 records/s | 3.93x | 3.23 to 5.23 |
| VCF write, genotypes read first | 28 MB/s, 244,000 records/s | 47 MB/s, 405,000 records/s | 1.66x | 1.49 to 2.29 |

htsjdk reads the BAM with `EAGERLY_DECODE` and `SILENT` validation, the closest it has to the
port's reader, which decodes every field and validates nothing. Its default `STRICT` would add a
validation pass on top. BAM read starts with the whole-file inflate above, so it gains the same
parallelism; BAM write, like BGZF deflate, owes most of its ratio to deflating on every core. The
VCF reader here is the deferred one (decision 0048): it keeps each record's genotype columns as text
and parses them on first access, as `LazyGenotypesContext` does, so the lazy row parses none and
the decoded row parses each record's once. `read_vcf` itself still parses at the read unless asked
not to, so that no accessor of a caller that has not opted in can fail. A genotype's alleles share
their record's bases rather than copy them, which is most of what keeps the parse cheap when it runs.

The two VCF write rows are one number and its absence: a record nothing looked at is written from
the file's own text, and a record whose genotypes were read is rebuilt column by column. That is
htsjdk's `LazyGenotypesContext`, and reproducing it made copying records through about three
times faster as a side effect of making it correct.

The BAM record codec alone, in memory with no BGZF, has no htsjdk column because htsjdk's codec
writes into a stream rather than answering with bytes: 871 MB/s and 2,854,341 records/s to
decode, 601 MB/s and 1,969,539 records/s to encode (the median run of each).

### Every other operation

The same run then times the rest of what the crates implement (`bam-ops-bench`, `vcf-ops-bench`,
`cram-bench` and `metrics-bench` against `BamOpsBench`, `VcfOpsBench`, `CramBench` and
`MetricsBench`), and asserts all 45 outputs identical before printing a rate. **bytes**: the
files the two sides wrote are identical. **decoded**: the operation's result lives in memory, so
what is compared is a rendering of it (records as SAM, sequences as name, length and MD5, query
answers). **idx**: identical except the indexed file's path and modification time, which name the
file rather than describe it.

| operation | check | over | htsjdk | port | ratio | range |
|---|---|---|---:|---:|---:|---:|
| SAM write | bytes | 153 MB of SAM, 400,000 records | 126 MB/s | 139 MB/s | 1.11x | 1.04 to 1.23 |
| SAM read | decoded | same | 150 MB/s | 517 MB/s | 3.45x | 3.42 to 4.16 |
| sort, queryname | decoded | 400,000 records | 2,438,000 rec/s | 4,701,000 rec/s | 1.93x | 1.77 to 2.54 |
| sort, coordinate (from queryname order) | decoded | same | 2,362,000 rec/s | 5,481,000 rec/s | 2.32x | 2.15 to 3.07 |
| BAI from a sorted BAM | bytes | 75 MB of BAM | 105 MB/s | 134 MB/s | 1.27x | 1.12 to 1.33 |
| SBI writer, 10 M offsets | bytes | 76 MB of offsets | 3,743 MB/s | 13,364 MB/s | 3.57x | 3.52 to 3.58 |
| BAM reheader, block copy, file to file | bytes | 75 MB of BAM | 200 MB/s | 247 MB/s | 1.24x | 0.99 to 2.26 |
| BAM gather, block copy, file to file | bytes | 226 MB of BAM | 199 MB/s | 190 MB/s | 0.96x | 0.82 to 1.12 |
| FASTQ write | bytes | 128 MB of FASTQ | 162 MB/s | 190 MB/s | 1.17x | 1.14 to 1.22 |
| FASTQ read | decoded | same | 641 MB/s | 1,642 MB/s | 2.56x | 2.48 to 2.90 |
| FASTA write with `.fai` and `.dict` | bytes | 64 MB of bases | 18.7 MB/s | 154 MB/s | 8.24x | 8.20 to 13.89 |
| FASTA read | decoded | same | 592 MB/s | 746 MB/s | 1.26x | 1.19 to 1.67 |
| indexed FASTA, 200,000 queries of 151 bp | decoded | 29 MB of bases | 104 MB/s | 120 MB/s | 1.15x | 1.13 to 1.19 |
| MD and NM | decoded | 400,000 reads | 510,000 rec/s | 926,000 rec/s | 1.81x | 1.48 to 1.96 |
| interval list write | bytes | 500,000 intervals | 211 MB/s | 350 MB/s | 1.66x | 1.25 to 2.16 |
| interval list read | decoded | same | 19.4 MB/s | 294 MB/s | 15.15x | 14.06 to 17.87 |
| interval list `sorted()` | decoded | same | 1,443,000 rec/s | 3,685,000 rec/s | 2.55x | 1.75 to 2.68 |
| interval list `uniqued()` | decoded | same | 1,183,000 rec/s | 3,498,000 rec/s | 2.96x | 2.05 to 3.20 |
| `OverlapDetector`, 50,000 targets, 100,000 queries | decoded | 150,000 intervals | 3,013,000 rec/s | 8,106,000 rec/s | 2.69x | 2.39 to 2.93 |
| liftover, chain load and 200,000 intervals | decoded | 200,000 intervals | 122,000 rec/s | 1,012,000 rec/s | 8.29x | 7.48 to 10.12 |
| CRAM 3.1 write | bytes | 200,000 records | 86,000 rec/s | 91,000 rec/s | 1.06x | 1.02 to 1.08 |
| CRAM 3.1 read | decoded | same | 122,000 rec/s | 115,000 rec/s | 0.95x | 0.93 to 0.98 |
| rANS 4x8 order 0, compress | bytes | 29 MB of qualities | 244 MB/s | 253 MB/s | 1.04x | 1.04 to 1.07 |
| rANS 4x8 order 0, uncompress | decoded | same | 227 MB/s | 339 MB/s | 1.49x | 1.41 to 1.52 |
| rANS 4x8 order 1, compress | bytes | same | 154 MB/s | 209 MB/s | 1.36x | 1.24 to 1.41 |
| rANS 4x8 order 1, uncompress | decoded | same | 171 MB/s | 209 MB/s | 1.22x | 1.20 to 1.24 |
| rANS Nx16 order 0, compress | bytes | same | 216 MB/s | 302 MB/s | 1.40x | 1.26 to 1.41 |
| rANS Nx16 order 0, uncompress | decoded | same | 171 MB/s | 275 MB/s | 1.61x | 1.43 to 1.61 |
| rANS Nx16 order 1, compress | bytes | same | 199 MB/s | 269 MB/s | 1.36x | 1.29 to 1.48 |
| rANS Nx16 order 1, uncompress | decoded | same | 221 MB/s | 246 MB/s | 1.11x | 1.04 to 1.17 |
| range coder order 0, compress | bytes | same | 30.8 MB/s | 42.2 MB/s | 1.37x | 1.18 to 1.39 |
| range coder order 0, uncompress | decoded | same | 33.0 MB/s | 41.9 MB/s | 1.27x | 1.15 to 1.28 |
| range coder order 1, compress | bytes | same | 37.0 MB/s | 56.2 MB/s | 1.52x | 1.29 to 1.57 |
| range coder order 1, uncompress | decoded | same | 33.9 MB/s | 47.0 MB/s | 1.39x | 1.31 to 1.55 |
| name tokeniser, compress | bytes | 4.6 MB of names | 16.1 MB/s | 24.0 MB/s | 1.49x | 1.12 to 1.51 |
| name tokeniser, uncompress | decoded | same | 104 MB/s | 90.6 MB/s | 0.87x | 0.80 to 0.89 |
| fqzcomp, compress | bytes | 29 MB of qualities | 19.1 MB/s | 23.4 MB/s | 1.23x | 1.19 to 1.26 |
| fqzcomp, uncompress | decoded | same | 18.7 MB/s | 21.5 MB/s | 1.15x | 1.15 to 1.21 |
| metrics file write | bytes | 200,000 rows, two histograms | 31.6 MB/s | 61.5 MB/s | 1.94x | 1.64 to 3.50 |
| BCF write | bytes | 23 MB of VCF, 200,000 records | 17.0 MB/s | 19.5 MB/s | 1.15x | 1.04 to 1.37 |
| BCF read | decoded | same | 140 MB/s | 224 MB/s | 1.60x | 1.54 to 2.29 |
| VCF write with `.idx` on the fly | idx | same | 76.9 MB/s | 257 MB/s | 3.35x | 2.27 to 3.68 |
| `.vcf.gz` write with `.tbi` | bytes | same | 18.4 MB/s | 22.4 MB/s | 1.22x | 1.16 to 1.29 |
| `.vcf.gz` read | decoded | same | 81.8 MB/s | 136 MB/s | 1.66x | 1.21 to 2.12 |
| BED decode | decoded | 21 MB, 500,000 lines | 23.8 MB/s | 180 MB/s | 7.57x | 6.12 to 8.72 |

CRAM is compared byte for byte: both sides write CRAM 3.1 with 5.0.0's NORMAL profile, fqzcomp
for the qualities and the name tokeniser for the names, and both name the file the same, since the
name is the file definition's 20-byte id. Both read the port's file back, `MD` restored on both
sides.

The rows where htsjdk is ahead or level:

* **Name tokeniser decompression** (0.87x): htsjdk 5.0.0's decoder is twice as fast as 4.2.0's,
  which the port beat; not profiled yet.
* **CRAM read** (0.95x) and **BAM gather** (0.96x, 0.82 to 1.12 across runners): within a few
  percent, not profiled.

**BCF write** is 1.15x where #272 measured 1.53x because it now pays what htsjdk pays: its input
is read deferred, so the genotype parse that `BCF2Writer.add` triggers lands inside the timed write
on both sides.

## Bit-identity contract

Output is compared byte-for-byte against goldens produced by the pinned reference running in
a digest-pinned `linux/amd64` container on JDK 17, on real x86-64 CI. Fields that are
legitimately allowed to vary (timestamps, version strings, command lines in headers) are
canonicalized under explicitly declared rules, and every comparison records which fields were
compared raw and which canonicalized.

Where a value cannot be matched exactly, it is quarantined and reported with its measured
divergence rate, and the affected output is described as **bio-identical** rather than
**bit-identical**. That vocabulary comes from
[broadinstitute/gatk#9384](https://github.com/broadinstitute/gatk/pull/9384), which
established that Java `Math.log/exp/pow` differ by roughly 1 ULP across CPU architectures
while only `StrictMath` is portable.

## Decisions

Each one records a place where the obvious implementation is *valid* and *wrong*, with
the measurement that settled it. They are the only barrier against these mistakes: the
compiler will never catch any of them.

| # | Finding |
|---:|---|
| [0001](docs/decisions/0001-deflate-backend.md) | Deflate backend must be zlib, never miniz_oxide |
| [0002](docs/decisions/0002-sorting-collection-tie-order.md) | Sort ties must be stable and run-indexed; no heap pin is required |
| [0003](docs/decisions/0003-deflate-fallback-is-a-status-not-a-length.md) | The BGZF no-compression fallback is a status test, not a length test |
| [0004](docs/decisions/0004-oracle-platform.md) | The emulated linux/amd64 container is a valid oracle; the plan's AVX assumption was wrong |
| [0005](docs/decisions/0005-java-math-has-three-implementations.md) | Java has three incompatible math libraries, and the port must track which one each call site used |
| [0006](docs/decisions/0006-correct-rounding-is-the-target-for-log-and-log10.md) | `log` and `log10` are correctly rounded, so round rather than port the intrinsic |
| [0007](docs/decisions/0007-pow-may-not-be-portable-across-x86-cpus.md) | `Math.pow` is deferred: its intrinsic depends on an approximate hardware instruction |
| [0008](docs/decisions/0008-the-bam-record-has-four-places-a-correct-encoder-diverges.md) | The BAM record has four places where a correct encoder still diverges |
| [0009](docs/decisions/0009-header-attribute-order-is-a-property-of-a-java-collection.md) | The SAM header's byte order is decided by a Java collection choice |
| [0010](docs/decisions/0010-the-bai-index-is-where-a-divergence-is-least-likely-to-be-noticed.md) | The BAI index is where a divergence is least likely to be noticed |
| [0011](docs/decisions/0011-metrics-number-formatting-depends-on-the-jvm-locale.md) | Metrics number formatting depends on the JVM's locale, and nothing pins it |
| [0012](docs/decisions/0012-nan-sign-bits-are-chosen-by-the-fpu.md) | NaN sign bits are chosen by the FPU, so the port is not bit-identical to itself across architectures |
| [0013](docs/decisions/0013-the-last-divergences-are-blocked-by-a-licence-not-by-difficulty.md) | The last formatting divergences are blocked by a licence, not by difficulty |
| [0014](docs/decisions/0014-math-exp-was-withdrawn-it-was-a-gpl2-transcription.md) | `Math.exp` is withdrawn: it was a transcription of GPL2-only source |
| [0015](docs/decisions/0015-htsjdks-writer-emits-records-its-own-reader-rejects.md) | htsjdk's SAM writer emits records its own reader rejects |
| [0016](docs/decisions/0016-the-vcf-header-comparator-is-not-a-total-order.md) | The VCF header comparator is not a total order |
| [0017](docs/decisions/0017-the-jvm-formats-a-short-decimal-not-the-value.md) | The JVM's `%f` rounds a short decimal, not the value |
| [0018](docs/decisions/0018-four-ordering-rules-in-one-library.md) | Four ordering rules in one library, none of them written down |
| [0021](docs/decisions/0021-the-coordinate-comparator-is-total-so-byte-identity-needs-a-stable-sort.md) | The coordinate comparator is a total order, so byte-identity needs a stable in-memory sort (resolves R5) |
| [0020](docs/decisions/0020-the-interval-tree-need-not-be-byte-exact.md) | The interval tree need not be byte-exact: the overlap set's order does not escape |
| [0019](docs/decisions/0019-re-sorting-a-sorted-header-is-stable-measured-not-proven.md) | Re-sorting an already-sorted VCF header is stable: measured, not proven |
| [0022](docs/decisions/0022-the-format-corpus-was-never-checked-against-the-oracle.md) | The format corpus was never checked against the oracle |
| [0025](docs/decisions/0025-fdlibm-is-portable-and-is-the-worse-stand-in-for-the-intrinsic.md) | FDLIBM is portable, and it is the *worse* stand-in for `Math.exp` |
| [0040](docs/decisions/0040-the-format-corpus-came-back-green.md) | The format corpus came back green, so the 99.73% is the oracle's number (closes 0022) |
| [0041](docs/decisions/0041-a-patch-sent-upstream-is-not-a-dependency-removed.md) | A patch sent upstream is not a dependency removed: Milestone U closed from this side |
| [0042](docs/decisions/0042-the-reference-version-is-pinned-by-the-consumers-not-by-this-repository.md) | The reference version is pinned by the consumers: 4.2.0 until the tool ports are done (superseded by 0047) |
| [0043](docs/decisions/0043-gkl-passes-levels-one-and-two-through.md) | GKL passes levels 1 and 2 through to ISA-L, correcting 0031 |
| [0044](docs/decisions/0044-strictmath-log-is-fdlibm-and-is-not-math-log.md) | `StrictMath.log` is FDLIBM, and it is not `Math.log`: 186 corpus points tell them apart |
| [0047](docs/decisions/0047-the-reference-moves-to-htsjdk-5.0.0.md) | The reference moves to htsjdk 5.0.0: 18 of 87 suites moved, each ported |

## What "finished" means

Milestone H's entries are ticked and the formats are done, so the open question is what is left.
[`docs/htsjdk-classes-a-ported-call-site-reaches.md`](docs/htsjdk-classes-a-ported-call-site-reaches.md)
answers it by walking the consumers rather than the library, the way decision 0023 answered it for
`jmath`: of the 43 htsjdk classes `picard-rs` and `gatk-rs` name, 14 were already ported here, 15
are named only as a type or a message, one is not an htsjdk class at all, and **12 were htsjdk
classes ported inside the repository that consumes them**. Those twelve, and the one method of
`CigarUtil` that was in the same position, **are here now**, each with a suite measured against the
reference class rather than through the tool that drives it. What is left of each move is the other
half of its acceptance test: the consumer calls the class here, drops its copy, and its own goldens
still pass.

## Conformance suites

Every suite is declared once in [`tools/conformance/manifest.json`](tools/conformance/manifest.json):
the harness that regenerates it, the golden it produces, and how the two sides are compared. The
oracle jobs of `.github/workflows/ci.yml` are generated from it, and the same run happens locally:

```sh
python3 tools/conformance/run_suite.py --list
python3 tools/conformance/run_suite.py --suites bgzf
```

A suite's `status` is part of the claim it supports: **oracle-backed** means CI re-derives the
golden in the pinned container on every run; **unchecked** means it has never been re-derived.
Decision 0022 found one corpus in the second category and decision 0040 closed it: `format` has been
re-derived on every push since and agrees, so **all 86 suites are oracle-backed** and no committed
golden here is unchecked. The zlib vectors decision 0001 rests on joined them in the same pass:
they had been confirmed inside the pinned container once, by hand, and are now regenerated on every
push like everything else.

## Part of a three-repository program

| Repo | Ports | Depends on |
|---|---|---|
| `htsjdk-rs` | htsjdk 5.0.0 | (none) |
| `picard-rs` | Picard 3.5.0 | `htsjdk-rs` |
| `gatk-rs` | GATK 4.6.2.0 | `picard-rs`, `htsjdk-rs` |

The topology mirrors upstream, including the direction of the dependencies.

## Relationship to `fulcrumgenomics/riker`

[riker](https://github.com/fulcrumgenomics/riker) is an independent Rust reimplementation of
Picard's QC tools, MIT-licensed, from the maintainers of Picard and htsjdk. It overlaps this
program's Picard layer, and the distinction is worth stating plainly because it decides how the
two relate.

**riker targets functional equivalence; this program targets byte equivalence.** riker's own
README says it "is not intended to be a drop-in replacement for Picard" — lowercase `snake_case`
headers, no metadata lines, and "bug fixes that yield slightly different outputs in specific edge
cases". So riker is a better tool; this is a byte-for-byte reproduction of the existing one,
including the bugs.

That makes riker a **source of divergence candidates, never a source to port from**. Its
`ERRATA.md` documents exactly where a careful reimplementer chooses to differ from Picard, and
every entry is a place this port must *not* differ. Those entries are pinned as conformance
cases in `picard-rs`, measured against the reference rather than trusted. Where riker's reading
of Picard and this port's agree, that is an independent cross-check of the reading; where riker's
errata is silent on something byte comparison surfaces (its list does not mention Picard's
alignment-block cycle binning, for instance), that is evidence for the method rather than against
riker.

## Commit attribution

Commits are co-authored with the model that wrote them. On 2026-07-21 the history of all three
repositories was rewritten to add that trailer uniformly, at the maintainer's request. The
rewrite changed every commit SHA. `picard-rs` and `gatk-rs` pin their dependencies by SHA, so
their **current** pins were updated to the rewritten commits, but pins recorded in *historical*
commits now name pre-rewrite SHAs that no longer exist. Historical builds before that date are
therefore no longer bit-reproducible; current and future ones are. This was a deliberate trade,
made with the cost understood.

## License

MIT, matching the htsjdk sources this ports. See `LICENSE`.

One caveat, taken from htsjdk's own README: htsjdk is **not uniformly MIT**. Licensing is
per-file, and notably the CRAM code is Apache License 2.0. Since CRAM is in scope for this
port, code derived from it will carry Apache 2.0 rather than MIT, and will say so in the file
header. Check the notice on the reference file before porting anything new.
