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
developer machine goes through an emulated container and measures the emulator. The tables below
were measured before decision 0047, against the 4.2.0 container, and have not been re-run against
5.0.0. Both sides are
library calls, not command-line tools: htsjdk 4.2.0 on JDK 17 with its defaults, the port at the
same commit, allocating through mimalloc as picard-rs does (each crate's `mimalloc` feature sets
it as the benches' global allocator; the libraries themselves leave the allocator to the binary).

Every table below is one run, [39651610628](https://github.com/IPNP-BIPN/htsjdk-rs/actions/runs/39651610628),
on an AMD EPYC 9V45, median of three. GitHub's runners are not one machine: an earlier run on an
EPYC 7763 measured both sides about a third slower, so compare ratios across runs, never rates.
Each answer is freed after its clock stops, as the JVM's collector frees htsjdk's outside the timed
region.

### BGZF

64 MB payloads, **after** the run asserted that every framed stream the two sides wrote was
byte-identical on all ten keys:

| operation | payload | level | htsjdk | port | ratio |
|---|---|---:|---:|---:|---:|
| deflate | text | 1 | 1,308 MB/s | 2,890 MB/s | 2.21x |
| deflate | text | 5 | 459 MB/s | 1,180 MB/s | 2.57x |
| deflate | random | 5 | 51.2 MB/s | 158 MB/s | 3.09x |
| inflate | text | 1 | 1,997 MB/s | 19,295 MB/s | 9.66x |
| inflate | text | 5 | 2,142 MB/s | 6,289 MB/s | 2.94x |
| inflate | random | 5 | 374 MB/s | 1,447 MB/s | 3.87x |

Deflate is stock zlib on both sides at every level: htsjdk's default deflater factory hands out
`java.util.zip.Deflater`, and the port's `Deflater::Jdk` calls the same zlib, so the bytes match
by construction. The `gkl-deflate` crate, a Rust port of Intel GKL's deflater that GATK installs
in htsjdk's place, is not on this path. Both ratios are parallelism, not faster code: the port's
BGZF writer deflates full blocks on every core, and `decompress_all` walks the block headers once
and then inflates the blocks on every core, each into its own slice of the output. htsjdk does
both on the calling thread. The bytes, and for a damaged stream the error, are the serial ones.

### BAM and VCF files

400,000 BAM records of 151 bases (122 MB inflated) and 200,000 VCF records with four samples
(23 MB), **after** the run asserted that the BAM and both VCFs the two sides wrote were
byte-identical. Rates are over the same bytes and records for both columns.

| path | htsjdk | port | ratio |
|---|---:|---:|---:|
| BAM read, every record decoded | 250 MB/s, 818,000 records/s | 549 MB/s, 1,800,000 records/s | 2.20x |
| BAM write | 32 MB/s, 104,000 records/s | 99 MB/s, 325,000 records/s | 3.12x |
| VCF read, genotypes left lazy | 86 MB/s, 733,000 records/s | 94 MB/s, 803,000 records/s | 1.10x |
| VCF read, genotypes decoded | 54 MB/s, 465,000 records/s | 95 MB/s, 811,000 records/s | 1.74x |
| VCF write, records copied through | 67 MB/s, 578,000 records/s | 207 MB/s, 1,777,000 records/s | 3.08x |
| VCF write, genotypes read first | 38 MB/s, 329,000 records/s | 86 MB/s, 734,000 records/s | 2.23x |

htsjdk reads the BAM with `EAGERLY_DECODE` and `SILENT` validation, the closest it has to the
port's reader, which decodes every field and validates nothing. Its default `STRICT` would add a
validation pass on top. BAM read starts with the whole-file inflate above, so it gains the same
parallelism; BAM write, like BGZF deflate, owes most of its ratio to deflating on every core. The
port's VCF reader parses every genotype up front: it reproduces `LazyGenotypesContext`'s
observable effect on the output, not its deferral, which is why the two read rows are one number
for the port and two for htsjdk.

The two VCF write rows are one number and its absence: a record nothing looked at is written from
the file's own text, and a record whose genotypes were read is rebuilt column by column. That is
htsjdk's `LazyGenotypesContext`, and reproducing it made copying records through about three
times faster as a side effect of making it correct.

The BAM record codec alone, in memory with no BGZF, has no htsjdk column because htsjdk's codec
writes into a stream rather than answering with bytes: 1,151 MB/s and 3,773,000 records/s to
decode, 745 MB/s and 2,441,000 records/s to encode.

### Every other operation

The same run then times the rest of what the crates implement (`bam-ops-bench`, `vcf-ops-bench`,
`cram-bench` and `metrics-bench` against `BamOpsBench`, `VcfOpsBench`, `CramBench` and
`MetricsBench`), and asserts all 44 outputs identical before printing a rate. **bytes**: the
files the two sides wrote are identical. **decoded**: the operation's result lives in memory, or
its bytes legitimately differ, so what is compared is a rendering of it (records as SAM, sequences
as name, length and MD5, query answers). **idx**: identical except the indexed file's path and
modification time, which name the file rather than describe it.

| operation | check | over | htsjdk | port | ratio |
|---|---|---|---:|---:|---:|
| SAM write | bytes | 153 MB of SAM, 400,000 records | 142 MB/s | 155 MB/s | 1.09x |
| SAM read | decoded | same | 170 MB/s | 645 MB/s | 3.80x |
| sort, queryname | decoded | 400,000 records | 3,090,000 rec/s | 6,247,000 rec/s | 2.02x |
| sort, coordinate (from queryname order) | decoded | same | 2,640,000 rec/s | 6,766,000 rec/s | 2.56x |
| BAI from a sorted BAM | bytes | 75 MB of BAM | 139 MB/s | 160 MB/s | 1.15x |
| SBI writer, 10 M offsets | bytes | 76 MB of offsets | 4,194 MB/s | 3,927 MB/s | 0.94x |
| BAM reheader, block copy, file to file | bytes | 75 MB of BAM | 235 MB/s | 251 MB/s | 1.07x |
| BAM gather, block copy, file to file | bytes | 226 MB of BAM | 204 MB/s | 205 MB/s | 1.01x |
| FASTQ write | bytes | 128 MB of FASTQ | 125 MB/s | 202 MB/s | 1.62x |
| FASTQ read | decoded | same | 744 MB/s | 2,374 MB/s | 3.19x |
| FASTA write with `.fai` and `.dict` | bytes | 64 MB of bases | 20.9 MB/s | 178 MB/s | 8.51x |
| FASTA read | decoded | same | 837 MB/s | 1,504 MB/s | 1.80x |
| indexed FASTA, 200,000 queries of 151 bp | decoded | 29 MB of bases | 131 MB/s | 114 MB/s | 0.86x |
| MD and NM | decoded | 400,000 reads | 665,000 rec/s | 1,221,000 rec/s | 1.84x |
| interval list write | bytes | 500,000 intervals | 112 MB/s | 196 MB/s | 1.74x |
| interval list read | decoded | same | 23.3 MB/s | 436 MB/s | 18.67x |
| interval list `sorted()` | decoded | same | 1,428,000 rec/s | 5,622,000 rec/s | 3.94x |
| interval list `uniqued()` | decoded | same | 1,209,000 rec/s | 5,024,000 rec/s | 4.16x |
| `OverlapDetector`, 50,000 targets, 100,000 queries | decoded | 150,000 intervals | 3,360,000 rec/s | 13,338,000 rec/s | 3.97x |
| liftover, chain load and 200,000 intervals | decoded | 200,000 intervals | 200,000 rec/s | 1,334,000 rec/s | 6.67x |
| CRAM write | decoded | 200,000 records | 149,000 rec/s | 131,000 rec/s | 0.88x |
| CRAM read | decoded | same | 281,000 rec/s | 343,000 rec/s | 1.22x |
| rANS 4x8 order 0, compress / uncompress | bytes / decoded | 29 MB of qualities | 220 / 291 MB/s | 356 / 472 MB/s | 1.62x / 1.62x |
| rANS 4x8 order 1, compress / uncompress | bytes / decoded | same | 185 / 195 MB/s | 272 / 278 MB/s | 1.47x / 1.43x |
| rANS Nx16 order 0, compress / uncompress | bytes / decoded | same | 179 / 188 MB/s | 405 / 366 MB/s | 2.27x / 1.95x |
| rANS Nx16 order 1, compress / uncompress | bytes / decoded | same | 170 / 170 MB/s | 262 / 228 MB/s | 1.54x / 1.34x |
| range coder order 0, compress / uncompress | bytes / decoded | same | 47 / 35 MB/s | 54 / 53 MB/s | 1.14x / 1.51x |
| range coder order 1, compress / uncompress | bytes / decoded | same | 49 / 39 MB/s | 71 / 60 MB/s | 1.46x / 1.54x |
| name tokeniser, compress / uncompress | decoded / decoded | 4.6 MB of names | 3.7 / 45 MB/s | 50 / 89 MB/s | 13.44x / 1.98x |
| fqzcomp, uncompress | decoded | 29 MB of qualities | 24.7 MB/s | 31.7 MB/s | 1.28x |
| metrics file write | bytes | 200,000 rows, two histograms | 37.4 MB/s | 107 MB/s | 2.86x |
| BCF write | bytes | 23 MB of VCF, 200,000 records | 24.2 MB/s | 35.4 MB/s | 1.46x |
| BCF read | decoded | same | 154 MB/s | 100 MB/s | 0.65x |
| VCF write with `.idx` on the fly | idx | same | 73.0 MB/s | 114 MB/s | 1.56x |
| `.vcf.gz` write with `.tbi` | bytes | same | 25.4 MB/s | 31.2 MB/s | 1.23x |
| `.vcf.gz` read | decoded | same | 92.3 MB/s | 89.7 MB/s | 0.97x |
| BED decode | decoded | 21 MB, 500,000 lines | 41.6 MB/s | 340 MB/s | 8.18x |

Where the bytes differ, the reason is a version, never a defect. The port writes htsjdk 5.0.0's
CRAM 3.1 and the container's htsjdk 4.2.0 writes CRAM 3.0, and 4.2.0 cannot read 3.1, so CRAM
write compares each side's file decoded by its own reader and CRAM read is both sides on htsjdk's
3.0 file. The port's reader restores `MD` where a record lacks it, as 5.0.0's and htslib's do and
4.2.0's does not, so the comparison drops that tag. The port's name tokeniser is 5.0.0's encoder,
so its row checks that htsjdk 4.2.0 decodes the port's streams to the input. fqzcomp has no
compress row because 4.2.0's `FQZCompEncode.compress` throws "not implemented".

The rows where htsjdk is still ahead:

* **BCF read** (0.65x): not found yet. The samples are spread over the decoder with no single
  hot spot once allocation is out of the way.
* **CRAM write** (0.88x) is not like for like: 3.1 with 5.0.0's codec choices, most of its time
  in the range coder, against 3.0, which never runs it.
* **Indexed FASTA** (0.86x): the span is one read now, but every query still seeks an unbuffered
  file and allocates its answer.
* **SBI writer** (0.94x) and **`.vcf.gz` read** (0.97x): within the run-to-run spread.

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
