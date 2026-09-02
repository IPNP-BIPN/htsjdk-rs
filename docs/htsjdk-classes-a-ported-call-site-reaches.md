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

A fifth group follows the four, and it is not from the grep at all: a class a consumer neither
ported nor named is invisible to a list built from names, whether it reimplemented the class badly
or refused the work outright. Both rows there cost a covering-array row.

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

`CigarUtil` was here (`htsjdk-bam::cigar`) for `softClipEndOfRead`, `clipEndOfRead` and
`mergeClippingCigarElement`, and not for `softClip3PrimeEndOfRead`, which `picard-rs`'s
`merge_bam_alignment_clip` ported locally. **It is here now**, and the half picard-rs's copy left
out is the half that is not about cigars: on a negative strand the alignment start moves by the
reference span the clip removed, a record with nothing aligned left becomes unmapped and loses its
coordinates, and `NM`/`MD`/`UQ` are dropped whenever the reference length changed at all.

The row also produced the sort of finding this list exists for: a cigar that is **all clipping**
comes back **unclipped**, not unmapped. `clip3PrimeEndOfRead` opens with
`if (!isValidCigar(rec, cigar, true)) return;`, and `Cigar.isValid` requires at least one real
operator, so `50S` leaves the method having done nothing at all.

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
| `SamFileValidator` | **moved here** (`htsjdk-bam::validation`) | `ValidateSamFile` is a thin wrapper around it; the validation rules are htsjdk's |
| `QueryInterval`, `Chunk`, `BAMIteratorFilter` | **all three moved here** (`htsjdk-bam::query` and `::iterator_filter`), with the bins a region reaches, the span a query reads and the per-record decision that ends it | the **read** side of the BAI. This crate builds an index (`htsjdk-bam::build_index`) and parses one; it cannot answer a query with one |
| `SBIIndexWriter`, `TextualBAMIndexWriter` | **both moved here** (`htsjdk-bam::sbi` and `::textual_index`) | index formats htsjdk writes, with no home here |

The index query is the one with visible cost already. Nothing here turns a `.bai` plus an interval
list into the records that overlap it, so a consumer that needs it writes its own overlap loop:
`picard-rs`'s `MergeSamFiles` port does exactly that for its `INTERVALS` path, and its first run
was five records out because `queryOverlapping` returns placed-but-unmapped reads and a hand-written
filter did not. That is the failure mode this list exists to predict: a reimplementation is not
wrong because it is duplicated, it is wrong because it is unverified.

## What this changes

Nothing about the formats. BAM, SAM, BGZF, CRAM, VCF, Tribble and the index *writer* are ported and
oracle-backed, and no consumer reimplements any of them.

What it changed is the meaning of the milestone. "Finish htsjdk-rs" was these twelve classes plus
one method of `CigarUtil`, each a move with a suite attached rather than an open-ended "more of
htsjdk". **All thirteen are here now**, each with a suite of its own measured against the reference
class rather than through the tool that drives it.

What is left of the moves is the other half of a row's acceptance test: the consumer calls the class
here and drops its copy, and its own goldens still pass. That is a pin bump in `picard-rs` and in
`gatk-rs`, and it is where a move stops being a duplicate and becomes a move.

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

## 5. Named by nothing, because there was no copy to name (2)

The list is built from names, so it can only see a reimplementation that says what it is
reimplementing. `SamFiles.findIndex` was reached by no comment and no import: `gatk-rs`'s
`CountReads` runner opened the index with `path.with_extension("bam.bai")`, four words that name
neither htsjdk nor the class whose rule they are half of.

Half, because htsjdk writes the index of `reads.bam` as `reads.bai` at least as often, and
`findIndex` looks for that name **first**. A covering-array row over the corpus's BAM therefore
counted two reads in the reference and none in the port: the interval query found no index, and no
index is not an error (IPNP-BIPN/gatk-rs#1020).

It is here now (`htsjdk-bam::sam_files`), and the measurement it arrived with says the rule is not
the one the method's own summary gives. The extension is **replaced** before it is **appended**, so
a `.csi` beside `reads.bam` beats a `reads.bam.bai` beside the same file. A CRAM has no replaced
`.csi` and no replaced `.bai`: `reads.crai`, then `reads.cram.crai`, then the shared fallthrough.
`Files.isRegularFile` refuses a directory, so a directory named `reads.bai` is stepped over rather
than returned. And a search that finds nothing resolves the path's symbolic links and runs again at
the real location, which means an index beside the link wins and an index beside the target is
still found.

The row is what the grep cannot produce, and it is the reason this file is judgement rather than a
script: the entry that costs a consumer a wrong answer is the one that never named the class.

The second row is the same blind spot reached from the other side. `TabixIndexCreator` and
`TabixIndex` were named by nothing because nobody reimplemented them: `gatk-rs`'s
`IndexFeatureFile` writes the Tribble index a plain VCF gets and **refuses** a block-compressed one,
in a message whose own words say the refusal is the port's and not GATK's. A refusal is the honest
form of a gap and it is still a gap, and it was the whole of the last row of that tool's covering
array (IPNP-BIPN/gatk-rs#1030).

It is here now (`htsjdk-tribble::tabix`), and half of it was here already: the creator drives the
same `BinningIndexBuilder` a `.bai` is built with, so what had to be ported is the layer that
decides the layout. Four things there are not what the class's summary suggests. A feature cannot
be indexed until the **next** one arrives, because the next feature's start is what closes the
previous one's chunk, and the last feature waits for `finalizeIndex`'s own position. The bin comes
from `regionToBin(start - 1, end)`, because `getIndexingBin` returns null and the builder computes
it from a zero-based half-open region. A feature with no end is one base for the **bin** and a
*shifted window* for the **linear index**, which are two different rules reached by one feature. And
the name block's declared size counts its null terminators, the names themselves being written a
low byte at a time rather than encoded.

The suite dumps each case twice, as the little-endian body and as the block-compressed file, so the
composition is measured as well as the arithmetic: a `.tbi` is a BGZF stream, and its bytes depend
on the deflate pin like every other block-compressed golden.

## How to regenerate this

The command at the top produces the raw list; the three groups are judgement, and the judgement is
in this file rather than in a script, because "named in a message" and "ported behaviour" cannot be
told apart mechanically. Re-run it when a consumer gains a tool: a new name in the third group is a
new row, and a new name in the first is not.

Group 5 has no command. It grows when a consumer's own measurement disagrees with the reference and
the reason turns out to be a rule that belongs to a class nobody named, which is to say from a
failure rather than from a sweep.
