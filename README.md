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

Ported from htsjdk `4.2.0`, the version pinned by GATK 4.6.2.0's `build.gradle`. GATK 4.7.0.0
pins htsjdk 5.0.0, and the target stays here until the tool ports are done: all three
repositories name one set of pins, and moving one alone reproduces a combination nobody ships
([decision 0042](docs/decisions/0042-the-reference-version-is-pinned-by-the-consumers-not-by-this-repository.md)).

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
the port's BGZF, BAM and VCF codecs against htsjdk's on the same inputs, in the same pinned
container, and asserts that every stream and file the two write is byte-identical **before** it
prints a duration. A speed number for a path with no golden is not a measurement of this port
([#78](https://github.com/IPNP-BIPN/htsjdk-rs/issues/78)).

```sh
tools/benchmark/run.sh 64 3
```

The figures come from the `benchmark` job, which runs it on a real x86-64 runner; a run on a
developer machine goes through an emulated container and measures the emulator. Both sides are
library calls, not command-line tools: htsjdk 4.2.0 on JDK 17 with its defaults, the port at the
same commit.

### BGZF

Measured on an AMD EPYC 7763 (4 vCPUs), 64 MB payloads, median of three, **after** the run asserted
that every framed stream the two sides wrote was byte-identical on all ten keys:

| operation | payload | level | htsjdk | port | ratio |
|---|---|---:|---:|---:|---:|
| deflate | text | 1 | 757 MB/s | 1612 MB/s | 2.13x |
| deflate | text | 5 | 182 MB/s | 597 MB/s | 3.27x |
| deflate | random | 5 | 36.2 MB/s | 105 MB/s | 2.90x |
| inflate | text | 1 | 1319 MB/s | 4170 MB/s | 3.16x |
| inflate | text | 5 | 1406 MB/s | 2043 MB/s | 1.45x |
| inflate | random | 5 | 275 MB/s | 286 MB/s | 1.04x |

Deflate is stock zlib on both sides at every level: htsjdk's default deflater factory hands out
`java.util.zip.Deflater`, and the port's `Deflater::Jdk` calls the same zlib, so the bytes match
by construction. The `gkl-deflate` crate, a Rust port of Intel GKL's deflater that GATK installs
in htsjdk's place, is not on this path. The deflate ratio is parallelism, not a faster compressor:
the port's BGZF writer deflates full blocks on every core and htsjdk's on the calling thread, and
the bytes are the same either way. Inflate is single-threaded on both sides, and 3.16x on the
compressible payload is what removing a per-block allocation and a per-block inflater is worth
(#226).

### BAM and VCF files

400,000 BAM records of 151 bases (122 MB inflated) and 200,000 VCF records with four samples
(23 MB), median of three, **after** the run asserted that the BAM and both VCFs the two sides wrote
were byte-identical. Rates are over the same bytes and records for both columns.

| path | htsjdk | port | ratio |
|---|---:|---:|---:|
| BAM read, every record decoded | 171 MB/s, 559,000 records/s | 147 MB/s, 481,000 records/s | 0.86x |
| BAM write | 25 MB/s, 83,000 records/s | 88 MB/s, 289,000 records/s | 3.49x |
| VCF read, genotypes left lazy | 64 MB/s, 546,000 records/s | 43 MB/s, 371,000 records/s | 0.68x |
| VCF read, genotypes decoded | 27 MB/s, 235,000 records/s | 43 MB/s, 370,000 records/s | 1.57x |
| VCF write, records copied through | 68 MB/s, 581,000 records/s | 157 MB/s, 1,346,000 records/s | 2.32x |
| VCF write, genotypes read first | 31 MB/s, 261,000 records/s | 47 MB/s, 406,000 records/s | 1.55x |

htsjdk reads the BAM with `EAGERLY_DECODE` and `SILENT` validation, the closest it has to the
port's reader, which decodes every field and validates nothing. Its default `STRICT` would add a
validation pass on top. The two rows where htsjdk is ahead are real. The port reads a BAM by
inflating the whole file before decoding it, where htsjdk streams. And the port's VCF reader
parses every genotype up front: it reproduces `LazyGenotypesContext`'s observable effect on the
output, not its deferral, so leaving the genotypes alone saves htsjdk the parse and saves the port
nothing. BAM write, like BGZF deflate, owes most of its ratio to deflating on every core.

The two VCF write rows are one number and its absence: a record nothing looked at is written from
the file's own text, and a record whose genotypes were read is rebuilt column by column. That is
htsjdk's `LazyGenotypesContext`, and reproducing it made copying records through 2.7x faster as a
side effect of making it correct.

The BAM record codec alone, in memory with no BGZF, has no htsjdk column because htsjdk's codec
writes into a stream rather than answering with bytes: 658 MB/s and 2,158,000 records/s to
decode, 446 MB/s and 1,462,000 records/s to encode.

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
| [0042](docs/decisions/0042-the-reference-version-is-pinned-by-the-consumers-not-by-this-repository.md) | The reference version is pinned by the consumers: 4.2.0 until the tool ports are done |
| [0043](docs/decisions/0043-gkl-passes-levels-one-and-two-through.md) | GKL passes levels 1 and 2 through to ISA-L, correcting 0031 |
| [0044](docs/decisions/0044-strictmath-log-is-fdlibm-and-is-not-math-log.md) | `StrictMath.log` is FDLIBM, and it is not `Math.log`: 186 corpus points tell them apart |

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
| `htsjdk-rs` | htsjdk 4.2.0 | (none) |
| `picard-rs` | Picard 3.4.0 | `htsjdk-rs` |
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
