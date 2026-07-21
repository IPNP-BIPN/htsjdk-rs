# 0019. Re-sorting an already-sorted VCF header is stable: measured, not proven

**Status:** accepted; hypothesis tested and refuted, guard kept
**Date:** 2026-07-21
**Follows:** [0016](0016-the-vcf-header-comparator-is-not-a-total-order.md)

## The hypothesis

Decision 0016 established that `VCFHeaderLine.compareTo` is not consistent: `VCFContigHeaderLine`
compares to other contigs by index and to everything else by rendered string, which admits a
cycle, so a `VCFHeader`'s output depends on the order its lines were inserted.

Porting the file layer surfaced a place where htsjdk feeds a header's own sorted order straight
back into a new sort. `VCFWriter.setHeader`, when `DO_NOT_WRITE_GENOTYPES` is set:

```java
this.mHeader = doNotWriteGenotypes ? new VCFHeader(header.getMetaDataInSortedOrder()) : header;
```

`getMetaDataInSortedOrder()` returns a `TreeSet`'s iteration order, and the `VCFHeader`
constructor pours it into a **new** `TreeSet` with the same comparator. Under a consistent
comparator that round trip is the identity. Under an inconsistent one there is no such
guarantee: a red-black tree's shape depends on the order of insertion and on comparisons whose
results need not be coherent.

So the hypothesis was that the same header, written the two ways, would produce two different
files, and that decision 0016's divergence would therefore reach a second code path.

## What the oracle says

It does not. `tools/vcf-conformance/ResortProbe.java` writes the same header both ways and
compares the metadata lines, over three inputs:

| input | stable? |
|---|---|
| the three contigs from decision 0016, whose comparison cycles | **yes** |
| contigs whose index order and string order agree (control) | yes |
| a realistic dictionary, `chr1 chr2 chr3 chr10 chr11 chrX` | **yes** |

`CYCLIC_RESORT_STABLE=true`. The hypothesis is refuted.

## Why the wording of this record matters

The honest statement is **measured stable, not proven stable**. Re-inserting an already-ordered
sequence into a `TreeSet` reproduces that order here, over these inputs, on this JDK. It is not
a theorem: with an inconsistent comparator, whether the round trip is the identity depends on
which comparisons the tree happens to make along each insertion path, and that depends on the
tree's shape, which depends on the input. A different contig set could behave differently.

Writing "re-sorting is a no-op" in the port would be a claim the evidence does not support.
Writing "re-sorting is stable for every header measured, and CI re-measures it" is what the
evidence supports, so that is what the port says and what CI checks. If a future htsjdk or a
future JDK changes `TreeMap`'s balancing, the probe fails and this record is reopened rather
than the divergence appearing silently in someone's file.

## Consequence for the port

None, which is the useful part. `write_vcf` writes the header through the same total order it
always uses, and no second ordering path is needed. A negative result that removes a suspected
divergence is worth as much as one that finds it, and costs a probe either way.

The record exists because the *absence* of a divergence here is not obvious, and the next
person to read `setHeader` will have the same suspicion. Without this they would re-derive it;
with it they can read the probe and the CI step.
