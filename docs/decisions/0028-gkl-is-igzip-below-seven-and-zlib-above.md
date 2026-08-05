# 0028. GKL is igzip below level 7 and zlib above it, and none of it is licence-blocked

**Status:** partly superseded by [0029](0029-only-levels-one-and-two-are-igzip.md). The level
table and the portability result stand; the claim that levels 1 to 6 are igzip does not.
Only levels 1 and 2 are igzip, and 3 to 9 are Intel's patched zlib 1.2.13, which the JDK's
zlib 1.3.2 disagrees with below level 7. Read 0029 first.
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
this milestone carries ("name the deflater it is a claim about") does not apply to them.

**Levels 1 to 6 are igzip**, and level 5 is the one that matters because it is the default.
*(Wrong: see 0029. Levels 3 to 6 are Intel's patched zlib, not igzip. Level 5 is still the
one that matters, and it is still not the JDK's, so the rest of this section holds.)*

**None of it is licence-blocked.** ISA-L is BSD-3-Clause, GKL is MIT, zlib is the zlib licence.
This is the first place in the programme where the reference implementation of a byte-deciding
component is *permissively* licensed, where decisions 0013 and 0014 are about the opposite
situation.
So both routes are open: linking ISA-L, or porting it. That is a materially different position
from `Math.exp`, where the only exact implementation is GPL2 and the milestone is stuck by law
rather than by effort.

## Decision

**H.4 is resized, not started.** Its statement becomes:

> reproduce ISA-L igzip's deflate output at levels 1 to 6, with level 5 first because it is the
> default; levels 7 to 9 are already exact through zlib.

And the accompanying warning is narrowed: a byte claim over BGZF at levels 7 to 9 need not name a
deflater, because there is nothing to choose between there.

## The second question, and how it is answered

Whether igzip's output depends on the CPU's instruction set. The library ships AVX2 kernels
(`igzip/adler32_avx2_4.s`) and dispatches on CPU features at load time, and a deflater whose
match-finder differs per kernel would produce different bytes on different silicon, which would
put igzip in the same category as `Math.pow` in decision 0007 rather than in the reproducible one.
There would then be no fixed target to port to, and the same BAM written on two machines would
differ. That is a larger result than any amount of porting, so it is settled first, and the same
way 0007 settled `pow`: regenerate on a second machine and diff.

`tools/gkl-probe/emulated.txt` is the first machine. It was produced on Apple Silicon, where
Docker translates `linux/amd64` through Rosetta. Its `cpu` line reads `VirtualApple @ 2.50GHz`,
and Rosetta implements no AVX, so that column is igzip's SSE path. The `igzip-portability` CI job
is the second machine: it reruns the probe in the same pinned image on a real x86-64 GitHub
runner, prints that host's CPU and AVX flags, and diffs every hash.

The probe hashes rather than measures lengths, because two deflate streams of equal length are not
equal streams: the `random` fixture at level 5 produces 60074 bytes through both deflaters and
two different hashes, so a length comparison would have reported an agreement it never checked.
It also compares whole BGZF streams and not only raw deflate output, since a BAM is a sequence of
gzip members with CRCs and block lengths, and a differing block-size decision would not show up in
the deflate comparison at all.

The job carries decision 0028's own first half as an assertion, so a GKL version bump that quietly
stopped routing through igzip fails rather than silently invalidating this record. Stated per level
rather than per row: at levels 7 to 9 **every** row must match the JDK, but below 7 the claim is
asserted on the `acgt` fixture alone, because the `runs` fixture compresses identically at levels 3
and 4, measured rather than assumed. Two match-finders can agree on data that simple.

## The answer

**101 rows compared, 0 differing.**

| | first machine | second machine |
|---|---|---|
| host | Apple Silicon under Rosetta | GitHub `ubuntu-latest` |
| `cpu` line | `VirtualApple @ 2.50GHz` | `AMD EPYC 7763 64-Core Processor` |
| AVX flags | none | `avx avx2 sse4_2` (no `avx512f`) |

Every hash agrees: four fixtures times nine levels of raw deflate through both deflaters, plus the
BGZF streams at levels 1, 5 and 9, plus the input hashes and the default level.

**igzip's output is a property of the algorithm, not of the kernel that ran it.** So H.4 has a
fixed target after all: a port can be checked against a byte sequence rather than against a
machine, and a BAM written on one host is the BAM written on another. This is the opposite of what
decision 0007 found for `Math.pow`, and it is the answer that made the porting work worth sizing.

Two limits, stated rather than implied. Only two microarchitectures have been compared, and
neither of them exercises igzip's AVX512 kernel, since the EPYC 7763 reports no `avx512f`. The
result stands for the SSE and AVX2 paths; a host with AVX512 would be a third column and has not
been run.
