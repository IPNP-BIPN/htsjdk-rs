# 0050. igzip's levels 1 and 2 are ported from the x86 kernels, so GKL's bytes come out on any CPU

**Status:** accepted
**Date:** 2026-11-13
**Follows:** [0031](0031-gkl-levels-one-and-two-are-isal-level-one.md), [0033](0033-the-whole-of-gkl-cuts-at-sse42.md), [0034](0034-igzips-reference-c-does-not-produce-igzips-bytes.md), [0043](0043-gkl-passes-levels-one-and-two-through.md)

## The decision

`gkl_deflate::deflate_gkl` answers levels 1 and 2 with `src/igzip.rs`, a Rust port of ISA-L 2.31's
stateless deflate. It no longer links ISA-L for them. Levels 1 and 2 therefore answer on every CPU,
Apple Silicon included, where they used to be refused unless an x86-64 CPU with SSE4.2 ran a linked
ISA-L. `igzip_available()` now always returns true, and is kept for callers that ask.

## What was ported

The target is the x86 assembly, not the readable C. Decision 0034 measured that the C makes 19749
bytes where GKL makes 19044. Decision 0033 measured that every kernel from SSE4.2 to AVX-512 agrees
with GKL. So one scalar transcription of the kernels is the target:

- `isal_deflate_icf_body_hash_hist_01`, from `igzip_icf_body_h1_gr_bt.asm`. It is shared by both
  levels. The kernel:
  - looks at two positions per step;
  - hashes with the `crc32` instruction (CRC32C, initial value 0, no final xor);
  - skips blind after 2 KB without a match;
  - emits matches longer than 258 as 258s, with its own off-by-one histogram writes;
  - leaves the hashes stale after a long skip.
  All of this is reproduced, the slop and padding rules included.
- `isal_deflate_icf_finish_hash_hist_01`, from `igzip_icf_finish.asm`: the last bytes of the
  input, one position at a time.
- The C around the kernels, translated as it is, from `igzip.c`, `huff_codes.c`, `flatten_ll.c`,
  `encode_df.c` and `proc_heap_base.c`. That covers:
  - the block header and the stored-block test;
  - the Huffman tables: heap, length limiting, run-length header, and the static-or-dynamic choice;
  - the repeated-byte prefix for 0x00 and 0xFF runs;
  - the stored-block fallback when the output would not fit.

The configuration is GKL's (decisions 0031 and 0043):

- the Java level is passed through as the ISA-L level;
- `ISAL_DEF_LVL1_DEFAULT` or `ISAL_DEF_LVL2_DEFAULT` bytes of level buffer;
- `end_of_stream = 1`, no flush, raw deflate.

The level buffer's size decides where blocks split. Both levels leave 261,524 bytes for tokens:
282624 - (87800 - 83088 + 16384) - 4 at level 1, and 331776 - (87800 - 83088 + 65536) - 4 at
level 2. The levels therefore differ only in the hash mask (8191 against 32767), which the C
shrinks for short inputs.

## The one bucket that cannot be reproduced

On a fresh stream, the body kernel's first-byte path hashes a register that still holds the address
of the `isal_zstream`. It uses that bucket as the third byte's. GKL therefore reads and overwrites
one hash bucket chosen by a pointer. The port uses a slot outside the table for it, holding the
table's initial value, which is what the real bucket holds in almost every run. The exception is an
address that hashes into one of the first two bytes' buckets. In that run GKL's own output depends on
where the allocator put the stream, so no fixed answer is GKL's. No fixture or corpus file has hit it.

## How it was checked

ISA-L 2.31.0 was built from source for linux/amd64 in Docker. Under Rosetta it reports SSE4.2 and
AVX2 and takes the kernel path. `tools/igzip-oracle` builds the same oracle. It reproduces GKL's
canary, so it is ISA-L as GKL ships it. Against that oracle, on Apple M4 Max:

| check | result |
|---|---|
| `igzip_canary` (GKL's own bytes) | identical, levels 1 and 2 |
| `tests/against_gkl.rs`, every level-1 and level-2 row (8 pairs) | identical |
| `tests/igzip_paths.rs`: 11 inputs aimed at each path (empty, tiny, skip, long matches, 0x00 and 0xFF prefixes, noise, stored fallback) | identical; every output inflates to its input |
| 568 generated files, 0 to 300,000 bytes, of ACGT, text, random, runs, zeros, 0xFF and mixtures | identical, both levels (1,136 streams) |

The hashes in `igzip_paths.rs` came from a developer machine under Rosetta. By decision 0008 that
makes them a pin, not a golden. The golden-grade check is on CI. On the x86-64 runner, the `isal`
feature links ISA-L next to the port, and `agrees_with_linked_isal` compares the two on 400 random
mixtures per level. On that runner a failing canary is a failure, not a skip. If CI and the pin
ever disagree, the pin is the suspect.

## What changes

- `isal-sys` is optional and off by default. The `isal` feature exists only for that cross-check.
  No code path that writes bytes uses it.
- `htsjdk-bgzf`'s `gkl-igzip` feature is a no-op. It is kept so manifests that enable it still build.
  `Deflater::Gkl` no longer refuses levels 1 and 2.
- CI keeps nasm, now only to build ISA-L for the cross-check step.

Speed, one thread, the corpus concatenated, 65,498-byte blocks, Apple M4 Max:

| level | output | throughput |
|---:|---:|---:|
| 1 | 4,190,509 bytes | 393 MB/s |
| 2 | 4,183,169 bytes | 392 MB/s |
| 5 (zlib port, for scale) | 3,571,536 bytes | 100 MB/s |

## What is not decided here

`gatk-rs` still pins the previous `htsjdk-rs`. Moving the pin is that repository's follow-up.
