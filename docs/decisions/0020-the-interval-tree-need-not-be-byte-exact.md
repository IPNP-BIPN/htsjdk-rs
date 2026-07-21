# 0020. The interval tree need not be byte-exact, because the overlap set's order does not escape

**Status:** accepted; measured in the oracle before a line was ported
**Date:** 2026-07-21

## The question that decides 1273 lines

`CollectRnaSeqMetrics` maps each read to the genes it overlaps through
`OverlapDetector.getOverlaps`, which is backed by `htsjdk.samtools.util.IntervalTree`, a 1273-line
augmented red-black tree. The reflex is that reproducing Picard byte-for-byte requires reproducing
that tree exactly, including its rotation and balancing, because a different tree shape yields a
different traversal order.

That reflex is expensive and, here, wrong. The question is not whether the tree's *shape* matches
but whether its *order reaches the output*. If the consumer folds the overlap set with a
commutative operation, the order is invisible and any correct overlap-set structure suffices.

## Reading the consumer

`getOverlaps` returns a `Set<Gene>` — a `HashSet`. `RnaSeqMetricsCollector` iterates it and, per
overlapping transcript, does two things:

1. `Gene.Transcript.assignLocusFunctionForRange`, which writes a per-base `LocusFunction` only
   when the new value's ordinal is **higher** (`Gene.java:161`). That is a max reduction:
   commutative and associative.
2. `addCoverageCounts`, which accumulates into **that transcript's own** array. Per-key
   summation: also commutative.

Both foldings are order-independent, so the `HashSet`'s iteration order cannot reach the output.

## Measured, not merely read

The reading is confirmed in the pinned oracle by `picard-rs`'s `tools/rnaseq-conformance/OverlapOrderProbe.java`. The probe exercises `CollectRnaSeqMetrics`, a *Picard* class, so it lives and runs in picard-rs CI where the oracle contains Picard, not here.
`Gene` does not override `hashCode`, so it inherits identity hash, whose value is not stable across
constructions; a `HashSet<Gene>` therefore iterates in a different order on each run. The probe
runs `CollectRnaSeqMetrics` twice, over an input with **fourteen genes overlapping the same
region** with differing locus functions — the only situation where a last-writer-wins fold would
diverge from a max fold — and compares the bytes:

```
TWO_RUNS_IDENTICAL=true
```

If the order reached the output, identity-hash variation between the two runs would have changed
it. It did not. The order does not escape.

## Decision

`OverlapDetector` is ported against a **correct** overlap structure, not a byte-exact port of
`IntervalTree`'s red-black balancing. The requirement is that `get_overlaps` returns the right
*set*; the order within it is unobservable downstream and is not part of the claim.

This is bounded, and the bound is stated so it is not quietly overrun. It holds for consumers that
fold the overlap set commutatively, which is what `CollectRnaSeqMetrics` does. A future consumer
that iterates the overlaps into an order-dependent computation — writing them to a file in
traversal order, say — would reintroduce the dependency, and this decision would have to be
revisited with a probe against *that* consumer. Each consumer's use is checked; the data structure
is not granted a blanket exemption.

## Why this is worth a decision and not just a code comment

It reverses the default effort. Without this, the RnaSeq port begins with a multi-week transcription
of a red-black tree whose every rotation must match. With it, the tree is an implementation detail
and the effort goes where the output actually comes from. A wrong answer to "must this be
byte-exact?" is expensive in exactly one direction, and the measurement is cheap, so it is taken
first.
