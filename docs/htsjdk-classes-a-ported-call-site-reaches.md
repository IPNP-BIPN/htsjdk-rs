# The htsjdk classes a ported call site reaches

Milestone H's entries are ticked, and "the htsjdk port is finished" is still not a statement anyone
can check. This is the same problem decision 0023 had for `jmath`, and the answer is the same one:
**walk the call sites rather than the library**. A port of htsjdk is finished when every htsjdk
class a ported call site reaches is reproduced, and every one that is not is named.

The list below is built from the two consumers. Every `htsjdk.<package>.<Class>` named in a Rust
source file of `picard-rs` or `gatk-rs` is collected, and each is answered with where its Rust
equivalent actually lives.

```sh
grep -rhoE 'htsjdk\.[a-z0-9_.]*[A-Z][A-Za-z0-9]*' picard-rs/crates gatk-rs/crates --include='*.rs'
```

43 distinct classes are named. 14 of them are named here as well and are ported here. The other 29
split four ways -- 15, 1, 12 and 1 -- and only the third group is work.

The fourth is one row, and it is the reason every row was checked against the reference tree rather
than taken from the comment that named it: `CircularByteBuffer` is **not an htsjdk class**. It is
`picard.util.CircularByteBuffer`, and `picard-rs`'s `fifo_buffer` says `htsjdk.samtools.util` in its
header. A misattribution in a comment becomes a work item in a list built from comments, so the
list is built from the comments and then answered from `htsjdk/src/main/java`: each of the twelve
below has a file there, and this one does not.

## 1. Named as a type or a message, not as behaviour (15)

`SAMException`, `SAMFormatException`, `RuntimeIOException`, `CodecLineParsingException`,
`SAMTextReader`, `SamReaderFactory`, `BAMRecord`, `Cigar`, `SAMFileHeader`, `SAMSequenceDictionary`,
`SAMSequenceRecord`, `SAMReadGroupRecord`, `Feature`, `Strand` and `StringUtil` appear in a comment
that explains a refusal message, a header field, or which Java type a record came from. The
behaviour behind each is either in `htsjdk-bam`/`htsjdk-tribble` already, or is the JVM's rather
than htsjdk's.

The exception types are the clearest case: a port reproduces the *text* htsjdk throws, because that
text is compared, and there is nothing else to port. `Exception in thread "main"
htsjdk.samtools.SAMException: Mate CIGAR (Tag MC) not found` is a string this programme must
produce; `SAMException` is not a class it must have.

## 2. Ported here under a different name (partial, 1)

`CigarUtil` is here (`htsjdk-bam::cigar`) for `softClipEndOfRead`, `clipEndOfRead` and
`mergeClippingCigarElement`, and **not** here for `softClip3PrimeEndOfRead`, which
`picard-rs`'s `merge_bam_alignment_clip` ports locally. One class, two homes: the adapter clip
belongs beside its siblings.

## 3. htsjdk classes living in a repository that consumes htsjdk (12, in 9 rows)

This is the list that makes "finished" checkable, and every row is the same shape: a class that is
htsjdk's, ported inside `picard-rs` or `gatk-rs` because this crate set did not offer it.

| htsjdk class | ported in | why it belongs here |
|---|---|---|
| `util.QualityUtil` | **moved here** (`htsjdk-bam::quality_util`); `picard-rs` still has its three copies until it bumps the pin | two call sites already, in two unrelated tools, and GATK reaches it as well |
| `ConstantMemoryDownsamplingIterator` | **moved here** (`htsjdk-bam::downsampling`) | it is an htsjdk iterator; `DownsampleSam` only drives it, and GATK's downsamplers reach the same family |
| `DuplicateScoringStrategy` | **moved here** (`htsjdk-bam::duplicate_scoring`) | scoring is htsjdk's, and both `MarkDuplicates` and GATK's `MarkDuplicatesSpark` use it |
| `filter.AlignedFilter` | **moved here** (`htsjdk-bam::filter`) | `htsjdk.samtools.filter` is a package of reusable predicates |
| `filter.ReadNameFilter` | **moved here** (`htsjdk-bam::filter`) | as above |
| `filter.TagFilter` | **moved here** (`htsjdk-bam::filter`) | as above |
| `SamFileValidator` | `picard-rs`: `validate_sam_file` | `ValidateSamFile` is a thin wrapper around it; the validation rules are htsjdk's |
| `QueryInterval`, `Chunk` **moved here** (`htsjdk-bam::query`), with the bins a region reaches and the span a query reads; `BAMIteratorFilter` still to come | `gatk-rs`: `gatk-engine::reads` | the **read** side of the BAI. This crate builds an index (`htsjdk-bam::build_index`) and parses one; it cannot answer a query with one |
| `SBIIndexWriter`, `TextualBAMIndexWriter` | `gatk-rs`: two tools | index formats htsjdk writes, with no home here |

The index query is the one with visible cost already. Nothing here turns a `.bai` plus an interval
list into the records that overlap it, so a consumer that needs it writes its own overlap loop:
`picard-rs`'s `MergeSamFiles` port does exactly that for its `INTERVALS` path, and its first run
was five records out because `queryOverlapping` returns placed-but-unmapped reads and a hand-written
filter did not. That is the failure mode this list exists to predict: a reimplementation is not
wrong because it is duplicated, it is wrong because it is unverified.

## What this changes

Nothing about the formats. BAM, SAM, BGZF, CRAM, VCF, Tribble and the index *writer* are ported and
oracle-backed, and no consumer reimplements any of them.

What it changes is the meaning of the milestone. "Finish htsjdk-rs" is these twelve classes plus one
method of `CigarUtil`, and each row is a move with a suite attached rather than an open-ended
"more of htsjdk". A row is done when the class lives here, its consumer calls it, and the consumer's
own goldens still pass, which is the cheapest possible acceptance test, because it already exists.

The first two moves are `QualityUtil` and the three `filter` predicates, and both repaid the trip.

The filters repaid it in the shape of the port rather than in a number: `filterOut` has a **pair**
form as well as a single one, and it is not the single form applied twice. `AlignedFilter` with
`includeAligned` keeps a pair only when both ends are mapped and without it keeps a pair when
either end is unmapped; `ReadNameFilter` needs both names listed or neither; `TagFilter` keeps a
pair when either end carries a listed value and drops one only when both do. Three classes, three
different asymmetries, and a port can agree on every single-record row while getting all three
wrong. The suite dumps both forms for that reason.

`DuplicateScoringStrategy` repaid it in a number that is not a sign. `compare` returns
`String.compareTo` verbatim when the scores tie, and `String.compareTo` answers the difference of
the first differing code units: comparing `read2` with `read5` is **-3**, not -1. A port that
normalized to a sign sorts identically and disagrees with the reference on every tie-break. Its
`RANDOM` strategy and its `compare` had no port at all before this: `picard-rs` needed the two
scores and stopped there.

`QualityUtil` repaid it immediately: the three copies in `picard-rs`
each closed on `f64::round`, which rounds half away from zero where `Math.round` rounds half up, so
every negative argument, meaning every observed error rate above one, was a different integer. The
version here goes through `jmath` and is compared against the reference's own answers by the
`quality-util` suite.

## 4. Named as htsjdk and belonging to something else (1)

`CircularByteBuffer` is `picard.util.CircularByteBuffer`. `picard-rs` ports it in `fifo_buffer`,
where it belongs, and the module's header names the wrong package. The row is here rather than
deleted because the check that caught it is the point: a list built from comments inherits the
comments' mistakes unless every entry is answered from the reference tree.

## How to regenerate this

The command at the top produces the raw list; the three groups are judgement, and the judgement is
in this file rather than in a script, because "named in a message" and "ported behaviour" cannot be
told apart mechanically. Re-run it when a consumer gains a tool: a new name in the third group is a
new row, and a new name in the first is not.
