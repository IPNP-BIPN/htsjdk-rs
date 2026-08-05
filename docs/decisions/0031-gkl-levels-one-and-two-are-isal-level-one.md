# 0031. GKL's levels 1 and 2 are ISA-L level 1, and the buffer size is part of the answer

**Status:** accepted
**Date:** 2026-08-05
**Follows:** [0029](0029-only-levels-one-and-two-are-igzip.md), [0030](0030-gkls-default-bam-bytes-depend-on-sse42.md)

## The gap

Decision 0029 read the level branch out of `libgkl_compression.so` and found `isal_deflate_stateless`
on levels 1 and 2. That names a function. igzip's output also depends on the ISA-L level, on the
size of the token buffer it is given, and on whether the stream is closed in one call, and none of
those are visible in the branch.

`crates/gkl-deflate` covers levels 3 to 9 exactly and refuses 1 and 2. This closes the gap on what
they are; the port of them is separate work.

## The answer

ISA-L 2.30.0, `isal_deflate_stateless`, with:

```c
s.level          = 1;
s.level_buf_size = ISAL_DEF_LVL1_DEFAULT;   /* 282624 */
s.end_of_stream  = 1;
```

All four fixtures, both Java levels, byte for byte against the column GKL produced in the pinned
container:

| fixture | bytes | |
|---|---|---|
| acgt | 19044 | identical |
| runs | 327 | identical |
| random | 60005 | identical |
| acgt-2blocks | 63311 | identical |

It also explains a result 0028 recorded without accounting for: **Java levels 1 and 2 produce
identical bytes**. GKL does not pass the level through to ISA-L, so both land on ISA-L level 1.
That was the measurement 0029 used to argue those levels are not zlib, and it now has a mechanism
rather than only a signature.

## The near miss worth recording

The obvious reading of the disassembly is wrong, and it is wrong in the way that survives testing.
The level-1 branch contains

```asm
55b8: movl $0x141d0, %esi     ; 82384
55bd: movl $0x1, %edi
55c2: callq calloc@plt
```

which looks exactly like a level buffer being allocated. Feeding ISA-L 82384 bytes reproduces GKL
on **three of the four fixtures** and misses the fourth by 62 bytes out of 63373.

Three of four is what a smaller corpus would have called a match. The fixture that caught it is the
200 KB one, which exists only because a corpus of single-BGZF-block inputs would never fill a token
buffer and so could never see the buffer size at all. A sweep of every buffer size from
`ISAL_DEF_LVL1_MIN` to 300000 in steps of 256 found no size that reproduced GKL, which is what sent
the search back to the assumption instead of onward to the ISA-L version.

## What is still open

**The port.** `deflate_gkl` refuses levels 1 and 2 and will keep refusing them until igzip's level 1
is reproduced in Rust. That is a larger piece than `deflate_medium` was: igzip does not emit deflate
directly but into an intermediate compressed format, then builds Huffman tables from the token
histogram and encodes from those.

What this decision removes from that work is every unknown except the code. The target is fixed, the
configuration is known, and the reference implementation is BSD-3-Clause.

**The CPU question, again.** ISA-L dispatches on CPU features exactly as Intel's zlib does, and
decision 0030 showed that dispatch changes bytes on the zlib side at the default level. Whether it
does on the igzip side is not measured here. It is measured the same way 0030 measured it: force the
feature detection and diff.
