# 0037. noodles reads htsjdk's CRAM, and 0036 was written too broadly

**Status:** accepted; corrects [0036](0036-noodles-is-a-reader-to-disagree-with.md), which stands on its evidence and overreached in its title
**Date:** 2026-08-05
**Follows:** [0036](0036-noodles-is-a-reader-to-disagree-with.md)

## What 0036 got right and what it did not

0036 measured two things and both hold. `noodles-cram`'s **primitives** are unreachable:
`io::reader::num` is `pub(crate)` and `codecs::rans_4x8::{encode, decode}` are `pub(crate)` inside a
`pub mod`, so calling them is `error[E0603]`. And where the ITF8 reader *is* reachable it disagrees
with htsjdk on truncation, returning `UnexpectedEof` where htsjdk returns `-1`.

From those two facts it concluded that noodles "cannot be more" than a reader to disagree with, and
that the port must hand-write everything. **That conclusion does not follow from that evidence**, and
it is wrong. The primitives being private says nothing about the layer above them, and the layer
above them is fully public: `read_file_definition`, `read_file_header`, `read_container`, `records`,
`seek`. 0036 measured the floor and generalised to the building.

## The measurement it should have made

Written with htsjdk in the pinned container — four unmapped reads, no reference — and read back with
`noodles_cram::io::Reader` 0.98.0:

| | htsjdk wrote | noodles read |
|---|---|---|
| file | 1646 bytes, CRAM 3.0 | version 3.0, 1 reference sequence |
| `read0` | flags 4, `ACGTACGT` | flags 4, `ACGTACGT` |
| `read1` | flags 4, `ACGTACGTACGT` | flags 4, `ACGTACGTACGT` |
| `read2` | flags 4, `ACGTACGTACGTACGT` | flags 4, `ACGTACGTACGTACGT` |
| `read3` | flags 4, `ACGTACGTACGTACGTACGT` | flags 4, `ACGTACGTACGTACGTACGT` |

Four of four. **noodles reads the CRAM htsjdk writes.**

## The rule, restated on the axis that actually matters

The dividing line is not "third-party crate or not". It is **reading against writing**, and it is the
same line decision 0001 draws for deflate.

**Writing is not reusable, and cannot become so.** The goal is byte-identity with htsjdk. Two
implementations of a compressing format produce different bytes for the same input, so any crate
that *emits* CRAM, BAM or VCF is disqualified by construction rather than by quality. That is
0034's lesson and it is unchanged.

**Reading may be reusable, per format, once measured.** A reader's output is records, not bytes, and
records are what the format defines. Where an independent reader returns the same records as htsjdk,
using it is a saving with no loss; where it does not, the difference is a finding. Which of the two
holds is an experiment, not a preference — and the experiment above is the first one, on four reads.

**The gate is which side of the tool the format sits on.** A GATK tool that reads CRAM and writes
BAM needs htsjdk's *records* on the way in and htsjdk's *bytes* on the way out. If the records agree,
the read side is free.

## What this does not license yet

Four unmapped reads with no reference is the smallest CRAM there is. It exercises no reference-based
compression, no mapped records, no substitution or insertion codes, no multi-container file, no tag
of any type but one integer, and no slice boundary. The claim earned here is "the container model
and the record loop agree on the simplest case", which is worth having and is not "the read side is
free".

The experiment that would license that is a corpus: a CRAM with mapped reads against a real
reference, several containers, every tag type and the quality-score paths, read by both and compared
record for record. Until it is run, the port keeps building its own reader, and noodles keeps its
0036 role of second opinion — which it has already earned twice, agreeing with the five-byte ITF8
nibble and disagreeing on truncation.

## The record of how this happened

0036 was written and merged in the same pass that measured the private primitives, and its
conclusion was reached by reasoning from them rather than by testing it. The test took one probe and
five minutes, and it was prompted by being asked why the crate was not reusable. The question was
the right one and the answer in 0036 was not.
