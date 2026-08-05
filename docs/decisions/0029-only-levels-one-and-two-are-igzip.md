# 0029. Only levels 1 and 2 are igzip; 3 to 9 are Intel's zlib fork

**Status:** accepted; corrects [0028](0028-gkl-is-igzip-below-seven-and-zlib-above.md)
**Date:** 2026-08-05

## What 0028 got wrong

Decision 0028 measured GKL against the JDK at every level and read the result as:

> levels 1 to 6 are igzip, levels 7 to 9 are zlib

Every number in that table is still correct. The **explanation attached to them was inferred, not
measured**, and it is wrong. The inference was that a level differing from the JDK must be igzip,
which is only sound if the two candidate backends are igzip and *the JDK's own zlib*. They are not.

## What the code does

`Java_com_intel_gkl_compression_IntelDeflater_resetNative` branches on the level before it
branches on anything else:

```asm
5415: movl %eax, %r13d      ; the level, read through FID_level
5418: subl $0x1, %eax
5422: cmpl $0x1, %eax
542e: jbe  0x5560           ; (level - 1) <= 1, so level is 1 or 2
```

The taken branch calls `isal_deflate_stateless_init`. The fall-through calls `deflateInit2_` with
`method=8`, `memLevel=8`, `strategy=0` and `windowBits = nowrap ? -15 : 15`. So:

| Java level | backend |
|---|---|
| 1, 2 | ISA-L igzip |
| 3 to 9 | zlib, inside the same shared object |

That zlib is **not** the JDK's. Two things separate them:

- **Version.** The GKL library says `deflate 1.2.13 Copyright 1995-2022`; the JDK in the pinned
  image says `deflate 1.3.2 Copyright 1995-2026`.
- **Patches.** The GKL copy exports `deflate_medium`, `slide_hash_sse` and a global
  `longest_match`, none of which exist in stock zlib. `deflate_medium` is Intel's addition and
  covers levels 4 to 6, which is exactly where 0028's table shows a difference that levels 7 to 9
  do not have.

So levels 3 to 6 differ from the JDK because they are a *different, patched zlib*, and levels 7 to
9 agree because on the `deflate_slow` path the patches change nothing and 1.2.13 and 1.3.2 emit the
same bytes.

## The measurement that already said so

0028's own data corroborates the disassembly and was misread at the time. **Levels 1 and 2 produce
byte-identical output on all four fixtures**, while 3, 4, 5 and 6 differ from each other. zlib
cannot do that: its `configuration_table` gives level 1 and level 2 different `max_lazy`,
`nice_length` and `max_chain`, so a zlib at those two levels emits two different streams. Two
levels collapsing to one output is the signature of a backend that is not zlib and does not
subdivide there.

## What this changes for H.4

**The default path is not igzip.** htsjdk's BGZF default is level 5, which lands in
`deflate_medium`. So the component that decides the bytes of every BAM GATK writes is Intel's
patched zlib 1.2.13, not ISA-L.

That moves the milestone, and toward the easier side: reproducing a known zlib version plus one
extra deflate strategy is a smaller and far better understood problem than reproducing igzip.
igzip remains in scope, but only for levels 1 and 2, which nothing writes by default.

Licences are unaffected: Intel's zlib fork is under the zlib licence and ISA-L is BSD-3-Clause, so
0028's conclusion that this milestone is not licence-blocked stands.

## What still stands from 0028

The level table itself, every hash, and the portability result: 101 rows compared between Rosetta
and an `AMD EPYC 7763`, 0 differing. That experiment did not depend on which backend produced the
bytes, only on whether the bytes moved with the CPU. They do not, on either backend.
