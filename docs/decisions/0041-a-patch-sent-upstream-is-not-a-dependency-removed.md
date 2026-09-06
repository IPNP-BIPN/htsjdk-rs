# 0041: A patch sent upstream is not a dependency removed

**Status.** Accepted. Closes Milestone U from this side.

## What Milestone U was

Six things this program wanted from `noodles`, each sent upstream as a pull request rather than
worked around:

| # | Wanted | Upstream |
|---|---|---|
| U.1 (#110) | CRAM's ITF-8, LTF-8 and uint7 codecs public | zaeleus/noodles#411 |
| U.2 (#111) | rANS 4x8 and Nx16 `encode`/`decode` callable | zaeleus/noodles#412 |
| U.3 (#112) | A VCF header round trip that is idempotent | zaeleus/noodles#413, #414 |
| U.4 (#113) | The raw blocks of a CRAM slice reachable | zaeleus/noodles#415 |
| U.5 (#115) | BGZF without a mandatory `rayon` and `crossbeam` | zaeleus/noodles#416 |
| U.6 (#116) | A `.bai` parser that does not cost `noodles-sam` | (blocked on U.5's shape) |

All five open pull requests are still open. That is not a complaint: they are an unpaid
maintainer's queue, and the patches were offered without a deadline attached.

## The mistake the milestone made

It tracked **the patches** rather than **what this program needed**, and those are not the same
thing. Every one of the six is either already answered inside this repository or is a dependency
this repository chose and can un-choose:

- **U.1, U.2 and U.4 were already answered here.** `htsjdk-cram` does not depend on `noodles` at
  all. Its ITF-8 and LTF-8 codecs are `htsjdk_cram::varint`, its rANS is `htsjdk_cram::rans` and
  `rans_order1`, and its block headers are `htsjdk_cram::block` -- all public, all measured by the
  27 CRAM suites. The upstream patches were offered for *other* users of `noodles`, and holding
  our own issues open for them confused a contribution with a gap.

- **U.3 is a defect in `noodles`' writer, and this port does not have it.** A VCF header here is
  one list of lines in the order they were read, sorted only where htsjdk sorts (decision 0016),
  so reading and writing is the identity -- `crates/htsjdk-vcf/tests/header_round_trip.rs` asserts
  it on a header whose FILTER, INFO and FORMAT lines are deliberately interleaved.

- **U.5 and U.6 were a dependency, and the dependency is gone.** `gatk-engine` parsed `.bai`
  through `noodles-bam` and read FASTA through `noodles-fasta`. It now uses
  `htsjdk_bam::index::read_bai` (nine goldens) and `htsjdk_bam::fasta_index` (the `indexed-fasta`
  suite), and `noodles`, `rayon` and `crossbeam` have left its dependency tree entirely.

## Why the FASTA reader had to be ported rather than borrowed

`getSubsequenceAt` answers a query by seeking: the first base's byte offset is computed from the
`.fai`'s bases-per-line and bytes-per-line, and every line boundary crossed is a jump over a
terminator whose length is the difference between those two columns. Nothing scans for a newline.
Two readers of the format therefore agree until a file's index disagrees with its own lines, and
the first version of the conformance dump proved it: a hand-written `.fai` with the CRLF offset one
byte wrong made the reference return a *terminator byte* among the bases, reported as an answer
rather than as an error.

## What `noodles` is here now

What decision 0036 said it should be: an independent implementation to check against, never a
source of bytes. It remains a **dev-dependency** of `htsjdk-vcf` and `htsjdk-bam`, which is the
role that cannot mislead -- a second opinion in a test, not a component of the answer.

## The rule this leaves

A patch offered upstream is a good thing to do and a bad thing to wait on. When a dependency's
behaviour is in the way, the question is not "when will they merge it" but "is this ours to
own" -- and for a byte-identical port of htsjdk, a reader of an htsjdk format always is.
