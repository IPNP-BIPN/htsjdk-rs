# 0028. GKL is igzip below level 7 and zlib above it, and none of it is licence-blocked

**Status:** accepted; H.4 resized from "port a compression algorithm" to "reproduce igzip at levels 1 to 6"
**Date:** 2026-08-05
**Follows:** [0001](0001-deflate-backend.md), [0003](0003-deflate-fallback-is-a-status-not-a-length.md)

## The question

Milestone H.4 asks for a GKL-exact deflate, and until it exists every byte claim over BGZF has to
name the deflater it is a claim about. The entry has sat unstarted, sized as "port ISA-L", which is
a compression algorithm and a large piece of work.

Before starting it, two things were worth measuring rather than assuming: **what is actually in
`libgkl_compression.so`**, and **which part of it produces the bytes GATK writes by default**.

## What is in the library

Both, as it turns out. Extracted from the pinned container's own jar:

```text
igzip_base.c   encode_deflate_icf   IGZIP_DIST_TABLE_SIZE   igzip/adler32_avx2_4.s
deflate_fast   deflate_medium       deflate_slow            deflate_stored
"deflate 1.2.13 Copyright 1995-2022 Jean-loup Gailly and Mark Adler"
```

ISA-L's igzip **and** zlib 1.2.13, in one shared object.

## Which one runs

Measured by deflating the same 60 KB buffer through `IntelDeflater` and through
`java.util.zip.Deflater` at every level, in the pinned container:

| level | GKL | JDK zlib | |
|---|---|---|---|
| 1 | 19044 | 20041 | different |
| 2 | 19044 | 19148 | different |
| 3 | 18723 | 18411 | different |
| 4 | 18059 | 18397 | different |
| **5** | **17952** | **18102** | **different** |
| 6 | 17897 | 17854 | different |
| 7 | 17975 | 17975 | **identical** |
| 8 | 18104 | 18104 | **identical** |
| 9 | 18104 | 18104 | **identical** |

The line to notice is that **htsjdk's default BGZF level is 5**, so every BAM GATK writes without
`--use-jdk-deflater` goes through igzip, not zlib.

## What that changes

**Levels 7 to 9 need nothing.** GKL is byte-identical to the JDK's zlib there, which this port
already reproduces. A byte claim at those levels is not deflater-dependent at all, and the warning
this milestone carries — "name the deflater it is a claim about" — does not apply to them.

**Levels 1 to 6 are igzip**, and level 5 is the one that matters because it is the default.

**None of it is licence-blocked.** ISA-L is BSD-3-Clause, GKL is MIT, zlib is the zlib licence.
This is the first place in the programme where the reference implementation of a byte-deciding
component is *permissively* licensed — decisions 0013 and 0014 are about the opposite situation.
So both routes are open: linking ISA-L, or porting it. That is a materially different position
from `Math.exp`, where the only exact implementation is GPL2 and the milestone is stuck by law
rather than by effort.

## Decision

**H.4 is resized, not started.** Its statement becomes:

> reproduce ISA-L igzip's deflate output at levels 1 to 6, with level 5 first because it is the
> default; levels 7 to 9 are already exact through zlib.

And the accompanying warning is narrowed: a byte claim over BGZF at levels 7 to 9 need not name a
deflater, because there is nothing to choose between there.

## What was not measured

Whether igzip's output depends on the CPU's instruction set. The library ships AVX2 kernels
(`igzip/adler32_avx2_4.s`), and a deflate implementation that dispatches on CPU features could
produce different bytes on different silicon — which would put igzip in the same category as
`Math.pow` in decision 0007 rather than in the reproducible one. **That question has to be settled
before any igzip byte claim, and it is settled the same way 0007 settled `pow`: regenerate on a
second machine and diff.**
