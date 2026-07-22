# 0021. The coordinate comparator is a total order, so byte-identity needs a stable in-memory sort

**Status:** accepted; resolves plan risk R5/R11
**Date:** 2026-07-22

## The question

The program's risk register flagged `SAMRecordCoordinateComparator` as *suspected non-total*: if it
were not a total order, the output of a coordinate-sorting tool (SortSam, MarkDuplicates) would
depend on heap size and spill boundaries, and no byte-identity claim could stand without pinning
`-Xmx` and `MAX_RECORDS_IN_RAM`. This is the fifth distinct ordering rule the port has had to pin
down (after the SAM header's `LinkedHashMap`, the VCF header `TreeSet`, VCF INFO/FORMAT, and the
queryname comparator), so it was worth settling rather than assuming.

## What reading it settles

It is a **total order**. `compare` is a fixed chain of tie-breaks, evaluated in order:

1. `fileOrderCompare`: reference index, then alignment start. An unmapped read (reference index
   `-1`) sorts **after** every mapped read.
2. strand: at the same position, forward before reverse.
3. within the same strand: read name, then flags, then mapping quality, then mate reference index,
   then mate alignment start, then inferred insert size.

Every step compares deterministic integers or a string, and the chain is consistent and transitive
(the strand step partitions one position into a forward group entirely before a reverse group, each
internally ordered by the field chain). Two records that agree on **all** of the above compare
**equal** — the comparator returns 0.

So the comparator is not the source of any nondeterminism. What is left is the treatment of records
it calls equal: their relative order is decided by the **sort**, not the comparator.

## The consequence for the port

Byte-identity for coordinate-sorted output requires a **stable, in-memory sort**: stable so equal
records keep their input order, in-memory so no spill-and-merge can reorder them. Rust's
`slice::sort_by` is stable, which satisfies the first half directly. The second half is a property
of the *oracle*, not the port: htsjdk's `SortingCollection` spills to disk past
`MAX_RECORDS_IN_RAM`, and its spill-merge is not guaranteed to preserve the input order of equal
records. The oracle contract therefore pins `MAX_RECORDS_IN_RAM` high enough that the conformance
corpora sort entirely in memory, and any tool whose real inputs exceed that is reported with the
pin stated as part of its claim, exactly as the risk register's mitigation anticipated.

This is the coordinate analogue of decision 0002 (sorting-collection tie order): the tie-break is
made total where it can be, and where genuinely-equal records remain, the sort's stability and the
in-memory bound carry the guarantee.
