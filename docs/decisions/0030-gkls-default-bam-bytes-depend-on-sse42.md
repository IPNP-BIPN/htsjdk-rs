# 0030. GKL's default BAM bytes depend on SSE4.2, and now there is a port that says so

**Status:** accepted; answers the question [0028](0028-gkl-is-igzip-below-seven-and-zlib-above.md) left open
**Date:** 2026-08-05
**Follows:** [0029](0029-only-levels-one-and-two-are-igzip.md)

## What was open

Decision 0028 asked whether GKL's output is a property of the algorithm or of the CPU, and
answered "the algorithm" on the strength of 101 hashes agreeing between two hosts. Then it had to
be corrected: the two hosts advertise identical CPUID flags, so they took the same path through the
library and the experiment tested nothing about dispatch.

The question stayed open for want of a host whose flags differ. It turns out not to need one.

## The answer, from Intel's own source

`libgkl_compression.so` is jtkukunas/zlib, and that source is public. Built for `linux/amd64` and
run on the same four fixtures, it reproduces the recorded GKL column **exactly at levels 3 to 9,
28 rows of 28**. So the source is the right source, and it can be read for what the library does
rather than guessed at.

What it does is choose its hash at run time:

```c
#define UPDATE_HASH(s,h,c) ( \
    x86_cpu_has_sse42 ? UPDATE_HASH_CRC(s,h,c) : UPDATE_HASH_C(s,h,c) \
)
```

With SSE4.2 the hash is a CRC-32C of the four bytes at the position, or of three when the level is
6 or above. Without it, zlib's multiplicative rolling hash over three bytes. These are not two
spellings of one function: the chains fill in a different order, so the match finder reaches
different candidates and takes different matches.

Rebuilding the same source with that one line forced to zero and diffing:

| level | fixtures whose bytes are unchanged |
|---|---|
| 3 | 1 of 4 |
| 4 | 0 of 4 |
| **5** | **0 of 4** |
| 6 | 3 of 4 |
| 7 to 9 | 4 of 4 |

**Level 5 is htsjdk's BGZF default.** A BAM written by GATK is therefore not one byte sequence: it
is one per CPU class. Levels 7 to 9 are unaffected because their chains are long enough that the
search converges on the same match whatever order it walks; level 6 nearly so, because its hash
covers the same three bytes the multiplicative one does.

## What this does not undo

0028's 101-row result stands as measured: GKL's output does not move between an emulated Apple
Silicon host and an AMD EPYC 7763. Both report SSE4.2, so both take the CRC branch, and that is
now the explanation rather than a caveat. **A GKL byte claim is a claim about a CPU class**, and
the class every oracle run has been in is "reports SSE4.2".

## Decision

`crates/gkl-deflate` carries the branch as a parameter, `Flavour::Gkl { sse42 }`, rather than
compiling in the majority case:

- `sse42: true` is the default and is checked against the recorded GKL column, 28 of 28;
- `sse42: false` is checked against `tools/gkl-probe/no-sse42.txt`, built from Intel's source with
  the flag forced off. That file's header says plainly that it is not an oracle column, because no
  host available here lacks SSE4.2 and inventing provenance for it would be worse than the gap.

A test asserts that the two branches still differ somewhere. If they ever stop differing, the
parameter has stopped meaning anything and should go.

## The consequence for H.4

The milestone's remaining unknown is no longer "is there a fixed target". There are two, and which
one applies is decided by the machine the oracle ran on. Anything downstream that asserts BAM bytes
at the default level has to say which.
