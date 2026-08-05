# 0034. igzip's reference C does not produce igzip's bytes, and that changes what porting it means

**Status:** accepted; sizes the one piece of H.4 still open
**Date:** 2026-08-05
**Follows:** [0031](0031-gkl-levels-one-and-two-are-isal-level-one.md), [0033](0033-the-whole-of-gkl-cuts-at-sse42.md)

## The assumption worth checking first

Every port in this programme so far has been checked against readable reference code. zlib's
`deflate.c`, Intel's `deflate_medium.c`, FDLIBM's `e_pow.c`, commons-math3's `FastMathCalc`: in each
case the source that ships is the source that decides the bytes, and a port is a translation.

Decision 0031 fixed igzip's configuration, and the obvious next step was to translate
`igzip/igzip_base.c` and `igzip/huff_codes.c` the same way. Before doing that, ISA-L was built two
ways and both were run on the same fixtures.

## What the two builds do

| build | acgt | runs | random | acgt-2blocks |
|---|---|---|---|---|
| default (assembly kernels) | **19044** | **327** | 60005 | **63311** |
| `make arch=noarch` (base C only) | 19749 | 652 | 60005 | 65523 |
| **GKL** | **19044** | **327** | **60005** | **63311** |

The C in ISA-L's own source tree is not what GKL ships. It is a second implementation of the same
idea that finds different matches, and the difference is not marginal: 327 bytes against 652 on the
`runs` fixture is a factor of two.

Decision 0033 measured the same boundary from the other side. A `core2duo` host, which lacks SSE4.2,
produces exactly the `noarch` numbers, because that is where ISA-L's dispatch lands when no kernel
qualifies. So the C path is real and reachable, and it is not the path any host GATK runs on takes.

## What porting igzip therefore means

Not translating `igzip_base.c`. Reproducing the match choices of the hand-written kernels:

| file | lines |
|---|---|
| `igzip_icf_body_h1_gr_bt.asm` | 906 |
| `igzip_body.asm` | 792 |
| `igzip_icf_finish.asm` | 327 |
| `igzip_finish.asm` | 330 |

plus `igzip_gen_icf_map_lh1_04.asm`, `igzip_set_long_icf_fg_04.asm` and `encode_df_04.asm`, and
`huff_codes.c` on top of them to build the Huffman tables from the token histogram.

The C remains useful as a description of the *shape* — the hash, the window, the intermediate
compressed format — but not as the definition of any single output byte. That has to come from the
assembly, or from a differential search against the library itself.

## Why this is worth a decision rather than a note

**It is the first component in the programme where the vendor's own reference implementation does
not produce the shipped bytes.** Every previous licence and provenance question has been "may this
source be translated"; this one is "which of the vendor's two implementations is the reference at
all", and the answer is the one that is harder to read.

It also puts a floor under the estimate. `deflate_medium` was 315 lines of C and took one pass with
a symbol-stream diff to find the one function that had been misread. This is roughly 2,400 lines of
SIMD assembly with no equivalent to read against, because the readable version is known to be wrong.

## What is not decided here

Whether to do it. Three routes exist and they are not equivalent:

1. **Port the kernels.** Exact, pure Rust, and the largest piece of work left in H.4.
2. **Link ISA-L.** BSD-3-Clause permits it, and it is the trade this repository already accepts for
   the JDK deflater: decision 0001 pins `flate2` to a *vendored C zlib*, not to a Rust
   reimplementation. Levels 1 and 2 would become exact immediately.
3. **Leave levels 1 and 2 unsupported**, which is the current state. htsjdk's BGZF default is 5, so
   nothing reaches them unless a caller asks for them by name.

This decision records that the choice exists and that it is a real choice, not a formality. It is
not made here.
