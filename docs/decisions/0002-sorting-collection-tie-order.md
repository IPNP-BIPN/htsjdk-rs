# 0002. Sort ties must be stable and run-indexed; no heap pin is required

**Status:** accepted
**Date:** 2026-07-21
**Resolves:** program plan open question 2, and downgrades risk R5
**Source:** htsjdk 4.2.0, commit `4cc010022ac038fb30f26e6f9717fabff3e808c1`

## Question

Is `SAMRecordCoordinateComparator` a total order? If not, does coordinate-sorted output depend
on `MAX_RECORDS_IN_RAM` and heap size, which would force those settings into the bit-identity
claim as pinned parameters?

The plan assumed the answer was "not total, therefore pin the heap". Both halves turned out to
need correcting.

## Finding 1: the comparator is not a total order

`SAMRecordCoordinateComparator.compare` chains: reference index, alignment start, negative
strand flag, read name, flags, mapping quality, mate reference index, mate alignment start,
inferred insert size. It then returns.

There is **no tie-break on index-in-file**. This is a deliberate contrast with
`picard.sam.markduplicates.MarkDuplicates.ReadEndsMDComparator`, which ends with
`read1IndexInFile` / `read2IndexInFile` and carries the in-source comment that this is
"arbitrary and is only included for completeness". Two records agreeing on all nine fields
compare equal, which is reachable in practice (exact duplicate records, and unmapped reads
where `fileOrderCompare` returns 0 immediately).

## Finding 2: the output is deterministic anyway, and heap-independent

The determinism comes from `SortingCollection`, not from the comparator. The chain:

1. Records are added in input order, so each RAM batch is a **contiguous run** of the input.
2. `spillToDisk()` sorts with `Arrays.parallelSort(ramRecords, 0, numRecordsInRam, comparator)`
   (l.247). The object overload of `Arrays.parallelSort` is **specified stable**, so ties keep
   input order within a batch. Parallelism does not affect the result.
3. Spill files are numbered `n` in creation order, so file index increases with input position.
4. `doneAdding()` spills the RAM remainder as the last numbered file (or, if nothing ever
   spilled, `iterator()` returns `InMemoryIterator` over a single stable sort). Either way all
   records are treated uniformly.
5. The merge is a **single** k-way merge, not multi-pass. `PeekFileRecordIteratorComparator`
   (l.654-656) delegates to the record comparator and, on tie, returns `lhs.n - rhs.n`: the
   lower-numbered spill file wins.

Therefore, for any two records comparing equal, the one earlier in the **input** is emitted
first, **regardless of where the batch boundaries fell**.

## Decision

**No heap pin is required.** `MAX_RECORDS_IN_RAM`, `-Xmx`, and the number of spill files are
not part of the bit-identity claim for coordinate-sorted output. They affect performance and
temp-file count only. This is a strictly better outcome than pinning them, because a pinned
heap would have leaked an environment detail into every downstream claim.

**The port must replicate stability, not just ordering.** Two specific requirements:

- The in-memory sort must be **stable**. Rust's `sort_unstable_by` is faster and is the wrong
  choice here: it reorders equal elements and silently diverges on any input containing
  tie-comparing records. Use `sort_by`.
- The external merge must break ties by **run index**, ascending, exactly as `lhs.n - rhs.n`
  does. A generic k-way merge on a binary heap without this tie-break is not equivalent.

## Why this matters beyond SortSam

This is the same failure shape as decision 0001. The wrong choice (`sort_unstable_by` there,
`miniz_oxide` here) produces output that is *valid and plausible*: correctly sorted, accepted
by samtools, round-trips fine. Only byte comparison against the reference catches it, and only
on inputs that happen to contain ties. Both traps are silent, both are the ergonomic default in
Rust, and both must be caught by convention rather than by the type system.

## Follow-up

`Arrays.parallelSort` stability is taken here from its specification. Confirm empirically in
the differential corpus with inputs deliberately containing tie-comparing records, including
exact duplicate records and multiple unmapped reads, at several `MAX_RECORDS_IN_RAM` values so
that batch boundaries land differently across runs.
