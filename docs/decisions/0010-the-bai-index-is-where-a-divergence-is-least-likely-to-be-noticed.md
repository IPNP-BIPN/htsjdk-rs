# 0010. The BAI index is where a divergence is least likely to be noticed

**Status:** accepted, measured against htsjdk
**Date:** 2026-07-21
**Follows:** [0009](0009-header-attribute-order-is-a-property-of-a-java-collection.md)

## Finding

The `.bai` is a pure side file. Nothing in it changes what a BAM contains, and every consumer
treats it as a hint: a wrong index degrades to a slower query, or to a query that re-checks
more records than it needed to, but not to an error. Delete it and tools rebuild it. That
combination makes it the part of the format where a divergence is least likely to be noticed
and most likely to persist once introduced.

Three of its properties are htsjdk's choices rather than the format's.

### 1. The linear index is back-filled, and htsjdk says so in a comment

A 16 kb window with no read overlapping it is not written as zero, and not as absent. It is
written as **the last non-empty window's offset**:

```java
// C (samtools index) also fills in intermediate 0's with values.  This seems unnecessary, but safe
long lastNonZeroOffset = 0;
```

htsjdk's own comment calls it unnecessary and does it anyway, for samtools compatibility. A
port that left empty windows at zero produces an index of exactly the same size that no reader
complains about.

The sentinel for "not yet set" is `-1`, not `0`, and that is load-bearing: `0` is a legitimate
virtual file pointer, so a zero-initialised array cannot distinguish "no read here" from "a
read at the very start of the file".

### 2. "Adjacent blocks" means the block address plus one

`Bin.addChunk` coalesces a new chunk onto the previous one when
`areInSameOrAdjacentBlocks(lastChunk.end, newChunk.start)`, which is:

```java
return (block1 == block2 || block1 + 1 == block2);
```

Block addresses are *byte offsets* of the block in the compressed stream, so `block1 + 1` would
require a one-byte BGZF block, which cannot exist. The second clause never fires on real files.
The rule is effectively "same block", and the port reproduces the literal form because being
right for the right reason costs nothing here.

### 3. Pseudo-bin 37450 holds four statistics disguised as two chunks

`writeChunkMetaData` writes bin number `MAX_BINS` (37450) with `n_chunk = 2`, then four
64-bit values that are not chunk boundaries at all: first offset, last offset, aligned record
count, unaligned record count.

The counting is asymmetric in two ways that are easy to get backwards. `firstOffset` starts at
`-1` while `lastOffset` starts at `0`, and both are written as they stand. And the pseudo-bin
is **counted in `n_bin`** while being **excluded from the loop** that writes the real bins:

```java
codec.writeInt(size + ((metaData != null) ? 1 : 0));
for (final Bin bin : bins) {
    if (bin.getBinNumber() == GenomicIndexUtil.MAX_BINS) continue;
    writeBin(bin);
}
```

A record with no coordinates is counted in a separate running total and is otherwise skipped;
a record that is unmapped but *placed* is both indexed and counted as unaligned. Those two
counters answer different questions and conflating them is the natural mistake.

## Verification

`tools/bam-conformance/BaiDump.java` emits 9 cases in the pinned oracle container, each as a
BAM **and** its BAI: an empty reference, a single read, two reads sharing a window and a block,
sparse reads leaving empty windows, reads at all six bin levels, three references with a gap in
the middle, mapped / placed-unmapped / unplaced-unmapped, 20,000 reads crossing block
boundaries, and reads sitting on a window boundary.

All 9 indices are byte-identical. The BAM of each case is checked first and separately, because
an index that matches for a file that does not is a coincidence rather than a result: the index
is made of virtual file pointers into that exact byte stream.

Sabotage:

| sabotage | indices that diverged | first differing byte |
|---|---:|---|
| back-fill replaced by zero | 2 of 9 | 136 and 280, inside the linear index |
| file pointer taken after the write instead of around it | 8 of 9 | 20, the first chunk |

Two things are worth noting in that table. The sizes were **identical in every diverging case**,
so nothing short of byte comparison would have caught either. And the back-fill is detectable
by only 2 of 9 cases: like the tag-ordering trap in decision 0008, it is invisible to most of
the corpus and would have been missed by a corpus assembled without it in mind.

## Consequence

`BgzfWriter` gained `file_pointer()`, ported from
`BlockCompressedOutputStream.getFilePointer()`. Indexing must be enabled before the first
record, because the chunk of a record already written cannot be recovered; `BamWriter::with_index`
is therefore a constructor-stage choice rather than a flag that can be flipped later.
