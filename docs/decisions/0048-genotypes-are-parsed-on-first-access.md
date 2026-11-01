# 0048. Genotypes are parsed on first access

**Status:** accepted, staged: deferral is opt-in and `Deref` stays
**Date:** 2026-10-22
**Would supersede:** the eager parse that `genotypes_context.rs` describes under "What this port
does differently, and why", and the comment in `reader::read_vcf` that keeps the refusal at the read
**Measured on:** htsjdk 5.0.0's source (`AbstractVCFCodec`, `LazyGenotypesContext`,
`BCF2LazyGenotypesDecoder`, `TribbleIndexedFeatureReader`), picard-rs `eee5dcc`, gatk-rs `e8322292`

## Why this is on the table

Two reasons, and the second matters more than the first.

**Speed.** Over three benchmark runs, reading a VCF without touching its genotypes is 0.76x
htsjdk (range 0.65 to 0.84), a `.vcf.gz` 0.67x, a BCF 1.00x. On the benchmark corpus (four
samples) the genotype parse is 57% of the read: the site columns alone, with the genotype block
kept as text, decode at 330 MB/s on a laptop, and 140 MB/s with the parse. The share grows with
the sample count, so on a cohort VCF it is most of the read.

**Fidelity.** htsjdk never parses the genotypes of a record nothing looks at. A tool that copies
records through (`MergeVcfs`, `SortVcf`, `GatherVcfs`, `SplitVcfs`, `VcfFormatConverter`,
`UpdateVcfSequenceDictionary`, `RenameSampleInVcf`) copies a malformed sample column verbatim and
succeeds. The port refuses the file at the read. Today's eager parse is therefore a deviation in
outcome, not only in cost; the comment that justifies it ("which call reports a broken file")
understates it.

## What htsjdk 5.0.0 does

* `parseVCFLine` wraps `parts[8]` in a `LazyGenotypesContext` and decodes it on the spot **only**
  when the header's samples are not already sorted (`!samplesWereAlreadySorted()`). `BCF2Codec`
  follows the same rule, and the port already mirrors it in both readers.
* Every accessor decodes first, except `size()` and `isEmpty()`, which answer from
  `nUnparsedGenotypes`: the header's sample count, not the number of columns the line has.
  `VCFEncoder` writes the unparsed text as it stands; a writer whose header has other samples
  reaches the genotypes by sample name, which decodes, so a write can be where the parse runs.
  (Measured below: it does not; the text is copied whatever the writer's sample order.)
* The parse is `codec.createGenotypeMap`, run with the **codec's state at the time of the
  access**:
  * `lineNo` is the codec's counter then, not the record's line;
  * `parts`, which one message quotes ("There are too many keys for the sample ..., values =
    `parts[genotypeOffset]`"), is the split of the line the codec read last.
* `VCFFileReader` iterates through `TribbleIndexedFeatureReader.WFIterator`, which reads one record
  ahead: when a tool first touches record N, the codec has already split record N+1. A late
  refusal therefore names N+1's line number and quotes N+1's columns. A tool that read the whole
  file into a list first gets the last line's counter and columns.
* Two refusals carry no line at all and are unaffected: `InternalCodecException` for a GT allele
  index that is not a number, and for one not defined in REF/ALT.
* The exception is thrown from the getter, so it surfaces in the tool's own loop, after whatever
  the tool had already written.

## The design

### Three states for `GenotypesContext`

| state | holds | made by |
|---|---|---|
| built | genotypes | callers, `new`, `collect` (unchanged) |
| parsed, text kept | genotypes and the text; drops the text on any look | today's `lazy`, `lazy_bcf` (unchanged) |
| **deferred** | the text block or BCF `LazyData`, a parse context, a `OnceLock<Result<Vec<Genotype>, GenotypeError>>` | the deferring readers |

A deferred VCF context keeps:

* the genotype block (kept today already);
* the record's alleles, contig and start, captured when the line is read, as
  `LazyVCFGenotypesParser` captures them. With `Allele` holding `Arc<[u8]>` the alleles cost a
  reference count each;
* an `Arc` to the file's parse context (header, codec version, sample order);
* an `Arc` to the reader's **cursor**: the codec's current line counter and the line it split
  last. The reader updates it as it reads; the parse consults it when it runs. That is what makes
  the line number and the quoted columns htsjdk's.

A deferred BCF context keeps the block (`LazyData`, kept today) and an `Arc` to the codec state
`decode_genotypes` reads (dictionary, header, sample order).

### The fallible accessor

```rust
impl GenotypesContext {
    /// `getGenotypes()`: parses a deferred block on first use and drops the text, as `decode()`.
    pub fn try_get(&self) -> Result<&[Genotype], GenotypeError>;
    pub fn try_get_mut(&mut self) -> Result<&mut Vec<Genotype>, GenotypeError>;
    pub fn try_into_vec(self) -> Result<Vec<Genotype>, GenotypeError>;
    /// `size()` / `isEmpty()`: never parse; the header's sample count while deferred.
    pub fn len(&self) -> usize;
    pub fn is_empty(&self) -> bool;
    /// Unchanged: the text, for a writer with the same samples. Never parses.
    pub fn unparsed_for(&self, samples: &[String]) -> Option<&str>;
}

/// The refusal as htsjdk throws it from the getter: class and message, the message worded with
/// the cursor's line number and columns at the moment of the parse.
pub struct GenotypeError { pub class: &'static str, pub message: String }
```

`VcfEncoder::encode_into` and the BCF writer gain an `EncodeError::Genotypes(GenotypeError)` for
the write that has to parse (a writer whose samples differ, as in `MergeVcfs`).

### What happens to `Deref`

`GenotypesContext` derefs to `Vec<Genotype>` today, and most callers read genotypes through it.
A `Deref` cannot return an error. Two ways:

1. **Strict.** Remove `Deref`, `DerefMut` and `IntoIterator for &GenotypesContext`. The compiler
   finds every read, and no input can panic anything. Every caller migrates (counts below).
2. **Staged.** Keep `Deref`, as `try_get().expect(...)`: it panics only on a deferred context
   whose block is malformed. Deferral is opt-in per reader (`read_vcf_with`, `read_bcf_with`, and a
   `GenotypesContext::deferred_vcf` constructor for picard-rs's own reader), and `read_vcf` keeps its
   eager default, so a caller that does not opt in sees no change at all. A caller that opts in
   must use `try_*`; nothing checks that it does.

**Recommendation: staged now, strict as the end state.** Staged lands the speed and the fidelity
where they matter first (picard-rs's copy-through tools) without forcing gatk-rs to move with
htsjdk-rs; strict follows once both consumers read through `try_*`. A grep in picard-rs's CI for
genotype reads that do not go through `try_*` would cover the gap between the two.

### Where refusals move

| case | port today | htsjdk 5.0.0 | deferred |
|---|---|---|---|
| sorted samples, record never touched (copy-through) | the read refuses | no refusal; column copied verbatim | **no refusal; copied verbatim** |
| sorted samples, genotypes touched | the read refuses, record's own line | the first access throws; cursor's line and columns | **the first `try_*` returns `Err`, cursor's line and columns** |
| sorted samples, written under other samples | the read refuses | no refusal; copied verbatim | `EncodeError::Genotypes` (the port's encoder decodes here; see the measurements) |
| samples not sorted | the read refuses | the read refuses | the read refuses (unchanged) |
| bad GT allele index (`InternalCodecException`) | the read refuses | the first access throws, no line | the first `try_*`, no line |
| BCF genotype block | the read refuses | the first access throws | the first `try_*` |
| `len()` of a malformed record | (never reached) | header sample count | header sample count |

Where a tool streams, a refusal that moves from the read to the tool's loop leaves behind what the
tool had already written, which is what Picard leaves.

## What picard-rs would change

picard-rs reads through its own `vcf_io` reader, which decodes lines ahead on every core and calls
`parse_genotypes` itself.

* **Reader** (`vcf_io::decode_at`): build `GenotypesContext::deferred_vcf` instead of calling
  `parse_genotypes` when the samples are sorted (unsorted stays eager). Keep a cursor per input
  and move it as records are **handed out**, one record ahead of the one handed out, not as the
  parallel read-ahead decodes them: the read-ahead is not htsjdk's, the iterator's one record is.
  `Record::lazy_genotypes` stays: it is how picard-rs copies a block verbatim under a header whose
  samples are ordered differently, which the port's encoder does not do (see the measurements).
* **Error plumbing**: a `GenotypeError` comes from the tool's own code, not from the reader. It
  ranks as the tool's own exception in `vcf_io`'s whole-file precedence (gzip, then UTF-8, then
  header, then first record, then the tool's exception), since htsjdk throws it from the tool's
  loop after the reads before it succeeded.
* **Call sites**: 39 lines in 20 files read or change genotypes: `LiftoverVcf` 6, `vcf_io` 4,
  `GenotypeConcordance`, `FindMendelianViolations` and `filter_vcf` 3 each, two each in
  `collect_variant_calling_metrics`, `MergePedIntoVcf`, `MakeSitesOnlyVcf`,
  `CombineGenotypingArrayVcfs` and `CollectArraysVariantCallingMetrics`, and one each in
  `CollectVariantCallingMetrics`, `CollectIndependentReplicateMetrics`, `VcfToAdpc`, `GtcToVcf`,
  `FixVcfHeader`, `CompareGtcFiles`, `fingerprint`, `extract_run`, `infinium` and `javascript`.
  Each read becomes a `try_*` with the error routed as above.
* **Behaviour change, by design**: the copy-through tools stop refusing a malformed sample column
  and copy it, as Picard does. New error-row goldens: the same malformed VCF through `MergeVcfs`
  (Picard succeeds, bytes compared) and through `FilterVcf` (message compared, line number
  included).

## What gatk-rs would change

gatk-rs uses htsjdk-rs's `VariantContext` directly and pins htsjdk-rs by revision, so nothing moves
until its pin does.

* **Staged**: nothing. `read_vcf` stays eager by default and `Deref` stays. A tool can opt in when
  its goldens are ready for it.
* **Strict**: 316 lines in 57 files touch `.genotypes` (gatk-tools 201, gatk-cli 66,
  gatk-annotation 28, gatk-engine 21): about 102 reads (`iter`, indexing, `get`, `first`), 52
  mutations, about 100 chains and `&vc.genotypes` arguments that rely on `Deref` coercion to
  `&[Genotype]`, and 52 `len`/`is_empty` calls that need nothing. Each read and mutation becomes a
  `try_*`, reported as GATK reports a `TribbleException` from a getter: unwrapped, exit code 3
  (measured below). The 79 `read_vcf`
  call sites change only if they opt in.

## What htsjdk-rs itself would change

* `genotypes_context.rs`: the deferred state, `try_*`, `len` from the header's sample count.
* `reader.rs` and `bcf2/codec.rs`: the opt-in readers and the cursor; `decode_genotypes` moved
  behind the BCF context.
* About 35 internal reads move to `try_*`: `variant.rs` 8, `bcf2/writer.rs` 7, `reader.rs` 6,
  `encoder.rs` 4, `genotype_type.rs` 3, `bcf2/decode.rs` 3, `chromosome_counts.rs` 2, `vcf_file.rs`
  and `bcf2/codec.rs` 1 each. `EncodeError` and the BCF writer's error gain the genotype refusal.
* Conformance: the existing genotype-parse golden (decode right after the read) and the
  lazy-genotypes golden stay as they are. A new dump pins what this design depends on:
  * a whole-file read, then the genotypes of record k: the line number and the quoted columns;
  * an iteration that touches each record as it is handed out: the one-ahead line number;
  * a copy-through write of a malformed block: succeeds, bytes;
  * the same block written under a writer with other samples: the write's refusal.
* The benchmark's lazy rows read through the deferred reader.

## Measured before building

All three were measured on 2026-10-22, before any code moved.

**The line a late refusal names** (`VcfDeferredGenotypesDump`, htsjdk 5.0.0, a record at file
line 9 with three values for two keys, followed by a good record at line 10, POS 300):

| how the record was read | line number | columns quoted |
|---|---|---|
| `VCFFileReader`, touched as handed out | 10 | line 10's (`keys = GT:GQ:DP ...`, `values = 300`) |
| `VCFFileReader`, collected first, then touched | 13, the last line | the last line's |
| `VCFCodec.decode`, touched before the next line | 9 | its own |

The count check names the record's own `chr1:400` with the cursor's line number; the allele-index
refusals carry no line; a file whose samples are not sorted is refused at the read, with `, for
input source: file://...` appended by the reader, which a late refusal never gets.

One thing the design above assumed is wrong: **a writer whose header orders the samples
differently does not decode.** `VCFEncoder` copies the unparsed text verbatim, columns not
reordered, so the write succeeds. The port's `unparsed_for` decodes in that case, which is an
existing difference from htsjdk, not one this change introduces; picard-rs already works around
it with its own copy of the block (`Record::lazy_genotypes`), which therefore has to stay. Left for
its own change.

**How the tools report it** (the same file through the pinned oracles):

| tool | exit | outcome |
|---|---|---|
| Picard 3.5.0 `SortVcf`, `MergeVcfs` | 0 | the malformed column copied verbatim |
| Picard 3.5.0 `FilterVcf` | 1 | `Exception in thread "main" htsjdk.tribble.TribbleException: ...line number 10...`, the record written before it kept |
| Picard 3.5.0 `MakeSitesOnlyVcf` | 1 | the same exception, nothing written |
| GATK 4.6.2.0 `CountVariants`, `SelectVariants` | 0 | counted; copied verbatim |
| GATK 4.6.2.0 `SelectVariants -sn s0`, `VariantsToTable -GF GT` | 3 | the same `TribbleException`, unwrapped (no "A USER ERROR"), the rows before it written |

**What deferral saves** (a laptop, the reader's decode without file I/O):

| samples | sites, block kept | sites and genotypes | the genotype parse's share |
|---:|---:|---:|---:|
| 4 | 330 MB/s | 140 MB/s | 57% |
| 100 | 3,693 MB/s | 141 MB/s | 96% |

htsjdk's lazy read of the 100-sample file, file I/O included, ran at 78 MB/s on the same laptop.

## What was built (htsjdk-rs)

* `reader::GenotypeParsing::{Eager, Deferred}`, `read_vcf_with` and `bcf2::codec::read_bcf_with`.
  `read_vcf` and `read_bcf` are the eager ones, unchanged.
* `GenotypesContext::deferred_vcf(text, samples, DeferredVcfGenotypes)`, for a reader of its own;
  `genotype_parse::VcfGenotypeSource::new(header, version)` is the file's codec and
  `advance(line_number, line)` moves its cursor; `DeferredVcfGenotypes::new(&source, alleles,
  contig, start)` is one record's parser.
* `try_get`, `try_get_mut`, `try_into_vec`, and `try_peek`, which parses without marking the
  context decoded. `GenotypeError { class, message }`. A refused parse leaves the context lazy.
* `EncodeError::Genotypes`; the BCF writer and `fully_decode` return a refusal as a `Bcf2Error` of
  the same class and message.
* `Deref`, `DerefMut`, `peek`, `into_vec` and `IntoIterator` parse as needed and panic on a
  refused block, with a message that names `try_get`.
* The deferred and eager readers produce the same records and the same bytes, copied or touched,
  for every well-formed input the BCF and VCF suites hold; a malformed block is refused with the
  same exception either way, at the read eagerly and on the first `try_*` deferred. The suite
  `vcf-deferred-genotypes` replays the dump's 38 rows.

## Not chosen: validate at the read, build on access

Scan every block at the read for anything the parse would refuse, without allocating, and build
the genotypes on first access. The refusal stays where it is today and no accessor can fail, so no
caller changes. It saves only the allocations, not the scan, and it keeps the deviation that
matters more: the copy-through tools would still refuse a file Picard copies.
