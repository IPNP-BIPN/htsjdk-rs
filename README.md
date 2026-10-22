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

Every table below is three runs of the same commit, [39752603417](https://github.com/IPNP-BIPN/htsjdk-rs/actions/runs/39752603417)
(AMD EPYC 9V74), [39753783551](https://github.com/IPNP-BIPN/htsjdk-rs/actions/runs/39753783551) and
[39753793628](https://github.com/IPNP-BIPN/htsjdk-rs/actions/runs/39753793628) (EPYC 7763), each a median
of three repetitions. GitHub's runners are not one machine, and one run on one of them can mislead
by a factor of two on the rows that allocate most, so each row shows the run whose ratio is the
median of the three, its rates, and the range of the ratio over all three. Each answer is freed
after its clock stops, as the JVM's collector frees htsjdk's outside the timed region.

### BGZF

64 MB payloads, **after** the run asserted that every framed stream the two sides wrote was
byte-identical on all ten keys:

| operation | payload | level | htsjdk | port | ratio | range |
|---|---|---:|---:|---:|---:|---:|
| deflate | text | 1 | 742 MB/s | 1,618 MB/s | 2.18x | 2.10 to 2.21 |
| deflate | text | 5 | 184 MB/s | 585 MB/s | 3.18x | 2.57 to 3.26 |
| deflate | random | 5 | 44.0 MB/s | 131 MB/s | 2.97x | 2.93 to 2.98 |
| inflate | text | 1 | 1,660 MB/s | 14,882 MB/s | 8.96x | 8.76 to 9.90 |
| inflate | text | 5 | 1,721 MB/s | 6,394 MB/s | 3.72x | 3.42 to 3.92 |
| inflate | random | 5 | 272 MB/s | 940 MB/s | 3.46x | 3.24 to 3.86 |

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
| BAM read, every record decoded | 182 MB/s, 596,000 records/s | 349 MB/s, 1,145,000 records/s | 1.92x | 1.92 to 1.99 |
| BAM write | 25 MB/s, 82,000 records/s | 89 MB/s, 291,000 records/s | 3.53x | 3.01 to 3.62 |
| VCF read, genotypes left lazy | 72 MB/s, 615,000 records/s | 55 MB/s, 469,000 records/s | 0.76x | 0.65 to 0.84 |
| VCF read, genotypes decoded | 31 MB/s, 264,000 records/s | 54 MB/s, 460,000 records/s | 1.74x | 1.67 to 1.87 |
| VCF write, records copied through | 75 MB/s, 643,000 records/s | 200 MB/s, 1,718,000 records/s | 2.67x | 2.01 to 2.90 |
| VCF write, genotypes read first | 31 MB/s, 270,000 records/s | 58 MB/s, 501,000 records/s | 1.86x | 1.51 to 2.05 |

htsjdk reads the BAM with `EAGERLY_DECODE` and `SILENT` validation, the closest it has to the
port's reader, which decodes every field and validates nothing. Its default `STRICT` would add a
validation pass on top. BAM read starts with the whole-file inflate above, so it gains the same
parallelism; BAM write, like BGZF deflate, owes most of its ratio to deflating on every core. The
port's VCF reader parses every genotype up front: it reproduces `LazyGenotypesContext`'s
observable effect on the output, not its deferral, which is why the two read rows are one number
for the port and two for htsjdk. A genotype's alleles share their record's bases rather than copy
them, which is most of what keeps the parse that cheap.

The two VCF write rows are one number and its absence: a record nothing looked at is written from
the file's own text, and a record whose genotypes were read is rebuilt column by column. That is
htsjdk's `LazyGenotypesContext`, and reproducing it made copying records through about three
times faster as a side effect of making it correct.

The BAM record codec alone, in memory with no BGZF, has no htsjdk column because htsjdk's codec
writes into a stream rather than answering with bytes: 706 MB/s and 2,315,252 records/s to
decode, 475 MB/s and 1,558,163 records/s to encode (the median run of each).

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
| SAM write | bytes | 153 MB of SAM, 400,000 records | 129 MB/s | 144 MB/s | 1.12x | 1.02 to 1.40 |
| SAM read | decoded | same | 146 MB/s | 504 MB/s | 3.45x | 3.36 to 3.48 |
| sort, queryname | decoded | 400,000 records | 2,286,000 rec/s | 4,581,000 rec/s | 2.00x | 1.99 to 2.16 |
| sort, coordinate (from queryname order) | decoded | same | 1,767,000 rec/s | 4,888,000 rec/s | 2.77x | 2.04 to 3.67 |
| BAI from a sorted BAM | bytes | 75 MB of BAM | 91.4 MB/s | 110 MB/s | 1.21x | 1.15 to 1.28 |
| SBI writer, 10 M offsets | bytes | 76 MB of offsets | 3,302 MB/s | 11,765 MB/s | 3.56x | 3.54 to 3.62 |
| BAM reheader, block copy, file to file | bytes | 75 MB of BAM | 1,252 MB/s | 1,477 MB/s | 1.18x | 0.90 to 1.37 |
| BAM gather, block copy, file to file | bytes | 226 MB of BAM | 410 MB/s | 376 MB/s | 0.92x | 0.92 to 0.96 |
| FASTQ write | bytes | 128 MB of FASTQ | 344 MB/s | 421 MB/s | 1.22x | 1.20 to 1.36 |
| FASTQ read | decoded | same | 650 MB/s | 1,359 MB/s | 2.09x | 2.03 to 2.86 |
| FASTA write with `.fai` and `.dict` | bytes | 64 MB of bases | 18.0 MB/s | 254 MB/s | 14.12x | 8.18 to 14.20 |
| FASTA read | decoded | same | 806 MB/s | 996 MB/s | 1.24x | 1.23 to 1.27 |
| indexed FASTA, 200,000 queries of 151 bp | decoded | 29 MB of bases | 95.6 MB/s | 112 MB/s | 1.17x | 1.16 to 1.18 |
| MD and NM | decoded | 400,000 reads | 540,000 rec/s | 948,000 rec/s | 1.75x | 1.56 to 1.81 |
| interval list write | bytes | 500,000 intervals | 214 MB/s | 338 MB/s | 1.58x | 1.55 to 1.60 |
| interval list read | decoded | same | 15.7 MB/s | 222 MB/s | 14.10x | 13.63 to 15.06 |
| interval list `sorted()` | decoded | same | 1,199,000 rec/s | 3,095,000 rec/s | 2.58x | 2.38 to 2.66 |
| interval list `uniqued()` | decoded | same | 1,331,000 rec/s | 3,984,000 rec/s | 2.99x | 2.88 to 3.05 |
| `OverlapDetector`, 50,000 targets, 100,000 queries | decoded | 150,000 intervals | 2,531,000 rec/s | 8,190,000 rec/s | 3.24x | 2.68 to 3.58 |
| liftover, chain load and 200,000 intervals | decoded | 200,000 intervals | 122,000 rec/s | 997,000 rec/s | 8.19x | 7.37 to 9.09 |
| CRAM 3.1 write | bytes | 200,000 records | 83,000 rec/s | 90,000 rec/s | 1.08x | 1.06 to 1.12 |
| CRAM 3.1 read | decoded | same | 124,000 rec/s | 119,000 rec/s | 0.96x | 0.94 to 1.01 |
| rANS 4x8 order 0, compress | bytes | 29 MB of qualities | 215 MB/s | 223 MB/s | 1.04x | 1.03 to 1.09 |
| rANS 4x8 order 0, uncompress | decoded | same | 234 MB/s | 339 MB/s | 1.45x | 1.43 to 1.46 |
| rANS 4x8 order 1, compress | bytes | same | 159 MB/s | 210 MB/s | 1.32x | 1.29 to 1.38 |
| rANS 4x8 order 1, uncompress | decoded | same | 178 MB/s | 217 MB/s | 1.22x | 1.19 to 1.22 |
| rANS Nx16 order 0, compress | bytes | same | 211 MB/s | 306 MB/s | 1.45x | 1.42 to 1.47 |
| rANS Nx16 order 0, uncompress | decoded | same | 170 MB/s | 278 MB/s | 1.64x | 1.61 to 1.67 |
| rANS Nx16 order 1, compress | bytes | same | 104 MB/s | 133 MB/s | 1.28x | 0.89 to 1.49 |
| rANS Nx16 order 1, uncompress | decoded | same | 127 MB/s | 145 MB/s | 1.14x | 1.05 to 1.15 |
| range coder order 0, compress | bytes | same | 30.0 MB/s | 41.9 MB/s | 1.40x | 1.28 to 1.50 |
| range coder order 0, uncompress | decoded | same | 33.1 MB/s | 41.6 MB/s | 1.26x | 1.25 to 1.29 |
| range coder order 1, compress | bytes | same | 41.5 MB/s | 64.9 MB/s | 1.57x | 1.53 to 1.62 |
| range coder order 1, uncompress | decoded | same | 32.6 MB/s | 46.7 MB/s | 1.43x | 1.37 to 1.48 |
| name tokeniser, compress | bytes | 4.6 MB of names | 17.0 MB/s | 25.3 MB/s | 1.49x | 1.28 to 1.55 |
| name tokeniser, uncompress | decoded | same | 54.5 MB/s | 46.1 MB/s | 0.85x | 0.84 to 0.93 |
| fqzcomp, compress | bytes | 29 MB of qualities | 18.6 MB/s | 23.0 MB/s | 1.24x | 1.21 to 1.29 |
| fqzcomp, uncompress | decoded | same | 22.5 MB/s | 26.9 MB/s | 1.20x | 1.16 to 1.20 |
| metrics file write | bytes | 200,000 rows, two histograms | 38.4 MB/s | 69.1 MB/s | 1.80x | 1.76 to 2.12 |
| BCF write | bytes | 23 MB of VCF, 200,000 records | 21.6 MB/s | 32.9 MB/s | 1.53x | 1.51 to 1.94 |
| BCF read | decoded | same | 137 MB/s | 137 MB/s | 1.00x | 0.99 to 1.53 |
| VCF write with `.idx` on the fly | idx | same | 69.5 MB/s | 130 MB/s | 1.87x | 1.86 to 2.55 |
| `.vcf.gz` write with `.tbi` | bytes | same | 22.1 MB/s | 25.4 MB/s | 1.15x | 1.13 to 1.21 |
| `.vcf.gz` read | decoded | same | 104 MB/s | 69.5 MB/s | 0.67x | 0.66 to 0.95 |
| BED decode | decoded | 21 MB, 500,000 lines | 28.0 MB/s | 175 MB/s | 6.27x | 5.83 to 7.23 |

CRAM is compared byte for byte: both sides write CRAM 3.1 with 5.0.0's NORMAL profile, fqzcomp
for the qualities and the name tokeniser for the names, and both name the file the same, since the
name is the file definition's 20-byte id. Both read the port's file back, `MD` restored on both
sides.

The rows where htsjdk is ahead or level:

* **VCF read with genotypes left lazy** (0.76x) and **`.vcf.gz` read** (0.67x), which is the same
  parse after an inflate: the port's reader parses every genotype up front, where htsjdk keeps the
  sample columns as text until something asks. With the genotypes decoded the port is ahead (1.74x).
  Deferring the parse would move where a malformed genotype is refused, from the read to the first
  look, so it is a decision about the API and not a speed change.
* **BCF read** (1.00x, 0.99 to 1.53): the same choice, for BCF's genotype block.
* **Name tokeniser decompression** (0.85x): htsjdk 5.0.0's decoder is twice as fast as 4.2.0's,
  which the port beat; not profiled yet.
* **BAM gather** (0.92x) and **CRAM read** (0.96x): within a few percent, not profiled.

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
