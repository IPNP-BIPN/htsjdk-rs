# 0049. zlib's bytes and zlib's refusals, without zlib's cycles

**Status:** accepted
**Date:** 2026-10-31
**Refines:** [0001](0001-deflate-backend.md), which said the backend must be zlib. What it needs is zlib's bytes.

## The decision

A BGZF block's deflate goes through `gkl_deflate`, and its inflate through `htsjdk_bgzf::inflate`,
both in Rust. The C zlib stays linked. It answers every block the inflater will not, and the tests
compare the port with it.

## Deflate

`Flavour::Jdk` was already compared with C zlib byte for byte at every level, and it was as fast as
zlib while rebuilding its state for every block. It now keeps one state per thread, and it holds the
position and hash in locals, because a store through `head` or `prev` otherwise makes the compiler
reload every field of `self`. Bits leave eight bytes at a time.

Cycles on one thread, 32 MB in 65,498-byte blocks, Apple M4 Max:

| input | level | zlib | port |
|---|---:|---:|---:|
| BAM records | 5 | 2.67 G | 1.94 G |
| SAM text | 5 | 3.34 G | 1.85 G |
| VCF text | 5 | 1.89 G | 1.18 G |
| BAM records | 1 | 1.72 G | 1.53 G |
| repeating text | 5 | 0.21 G | 0.13 G |
| random bytes (16 MB) | 5 | 0.83 G | 0.67 G |
| random bytes (16 MB) | 1 | 0.79 G | 0.82 G |

Level 1 on random bytes is the one row the port does not win.

`examples/jdk_speed check` compares every block of those BAM, SAM and VCF inputs with zlib at all
nine levels, through one reused state; all 27 runs are identical. A reused state must answer as a new
one does, which `a_reused_deflater_answers_as_a_new_one` checks for all three flavours. Level 0 and
the no-compression fallback are `gkl_deflate::stored`, which `stored_is_zlibs_no_compression` checks.

## Inflate, and why it is not libdeflate

The first version put libdeflate in front of zlib and accepted a block only when its CRC matched,
reasoning that any stream the two disagree on would decode to something else. The differential test
written for it refuted that on its first run: with one bit flipped, a static-code 285 becomes 287.
libdeflate's table reads 286 and 287 as a length of 258, so the payload and the CRC are both right.
zlib says `invalid literal/length code`, and so would an htsjdk running zlib. The reader's contract
is acceptance, so a CRC cannot be the guard.

So the Rust inflater accepts a strict subset of what zlib accepts, and declining is not an answer:
the caller runs zlib from the start of the block, and only zlib ever produces an error.

| zlib's refusal | here |
|---|---|
| `invalid block type` | `inflate_raw`, type 3 |
| `invalid stored block lengths` | `stored` |
| `too many length or distance symbols` | `read_dynamic_tables` |
| `invalid code lengths set` | `build`: the code must fill its codespace |
| `invalid bit length repeat` | `read_dynamic_tables`: no previous length, or past the end |
| `invalid code -- missing end-of-block` | `read_dynamic_tables` |
| `invalid literal/lengths set`, `invalid distances set` | `build` |
| `invalid literal/length code`, `invalid distance code` | `KIND_INVALID` in `codes` |
| `invalid distance too far back` | `codes` |
| a stream that ends early, or a size that disagrees | `codes`, `stored`, `pos == size` |

`build` is stricter than zlib on purpose: zlib allows one code of length 1 or no distance code at
all, and those blocks simply go to zlib. `it_accepts_nothing_zlib_refuses_and_answers_as_zlib_does`
runs over 100,000 damaged streams from all three encoders, and
`a_damaged_block_answers_as_zlib_alone_did` flips every bit of three blocks through the reader.

Cycles on one thread: 1.39 G to 1.09 G for a 79 MB BAM, 0.18 G to 0.13 G for a 7.7 MB `.vcf.gz`.

## What it does not settle

- Every figure above is from one Apple M4 Max, on a machine other work was running on, which is why
  they are cycles and not seconds. The README's tables are from the x86-64 `benchmark` job and were
  not re-run for this.
- The inflater is not libdeflate's speed. It has no multi-literal or vector copy paths yet, and the
  gain could grow. Any such change has to keep the subset property, not only the bytes.
- CRAM's range coder is still a linear walk over each model, half of CRAM's active samples. Cumulative
  frequencies with a branchless search or a vector count took 30% more cycles, separate arrays walked
  eight at a time 3% more, and one arena for all a stream's models 4% more instructions for 3% less
  peak memory. None was kept.
- A BAM record holds 607 bytes measured, 4 allocations, for 319 on disk, down from 717 and 7: a tag
  block the list would write back is kept as it is, and freed when the list is made, so a record
  holds one or the other. `Tags::find` and `Tags::decoded` answer without leaving a list behind, and
  the SAM and CRAM readers compact what they parsed. `get` and `iter` still make the list, and a
  record then costs a little more than before: `Tags` is 32 bytes larger inline.
- A parsed VCF record of four samples holds 1,724 bytes, 24 allocations, down from 1,916; the fields
  are public, so going further is an API decision. It also holds its genotypes twice, parsed and as
  the text a writer copies, which is what opt-in deferred genotypes (PR #274) stop: they measure 960.
- An exact `reserve` of the INFO attributes as well gave 1,612 bytes, and 30% more instructions in
  `vcf-bench` under mimalloc, which either of the two reservations alone does not. The cause was not
  found; that reservation was dropped. Any capacity change has to be re-measured with the allocator
  the tools ship.
