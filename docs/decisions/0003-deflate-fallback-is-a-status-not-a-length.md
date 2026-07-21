# 0003. The BGZF no-compression fallback is a status test, not a length test

**Status:** accepted
**Date:** 2026-07-21
**Source:** htsjdk 4.2.0, `BlockCompressedOutputStream.deflateBlock`
**Found by:** differential test against the official htsjdk 4.2.0 jar

## The bug

The first port of `deflateBlock` triggered the no-compression fallback on

```rust
if compressed.len() > COMPRESSED_BUFFER_SIZE   // WRONG
```

which reads as a faithful translation of htsjdk's "did it fit in `compressedBuffer`". It is
not. 33 of 35 golden vectors passed; the incompressible payload failed, producing 200197 bytes
against htsjdk's 200152, and the resulting file was structurally corrupt.

## Why it is wrong

htsjdk deflates into a fixed array of `COMPRESSED_BUFFER_SIZE` (65518) bytes and then tests
`deflater.finished()`:

```java
int compressedSize = deflater.deflate(compressedBuffer, 0, compressedBuffer.length);
if (!deflater.finished()) { /* redo at NO_COMPRESSION */ }
```

When the deflated output *exactly* fills the buffer, zlib returns with `avail_out == 0`. It has
emitted every byte, but it cannot signal end-of-stream without room to write the final marker,
so `finished()` is **false** and the fallback is taken even though the data did fit.

That boundary is not hypothetical. It is exactly where incompressible data lands: 65498 bytes
of random input, which is precisely one full uncompressed block, deflates to **exactly 65518
bytes** at level 1, the buffer size to the byte. htsjdk therefore stores the block uncompressed
at 65503 bytes, while the length test kept the 65518-byte compressed form.

The corruption followed from that. `BSIZE` is a `u16` holding `total - 1`; with a 65518-byte
payload plus 26 bytes of framing, `total - 1` is 65543, which overflows and wraps to 7. Every
subsequent block offset was then wrong.

## Decision

Mirror Java's control flow, not its prose. Bound the output by the `Vec`'s capacity, which is
the analogue of Java's fixed array, and branch on the returned `Status`:

```rust
let mut compressed = Vec::with_capacity(COMPRESSED_BUFFER_SIZE);
let status = c.compress_vec(&self.buffer, &mut compressed, FlushCompress::Finish)?;
if status != Status::StreamEnd {
    // no-compression fallback
}
```

`flate2`'s `compress_vec` writes only into spare capacity and returns `Status::Ok` rather than
`StreamEnd` when it runs out, which reproduces `finished() == false` exactly, boundary included.

## What this says about the method

The length test is the reading a careful developer produces from the Java *semantics*. The
status test is what the Java *code* does. They agree on all but one input, and that input is
reachable in practice: any BAM containing a run of incompressible data (already-compressed
embedded content, high-entropy tags) hits it.

Two working rules follow:

1. **Port the control flow, not the intent.** Where the reference branches on an API's return
   state, find the Rust API whose state machine matches, rather than reconstructing the
   condition from what the state is believed to mean.
2. **Adversarial payloads earn their place in the corpus.** This was caught only because the
   corpus deliberately includes an incompressible payload sized to land on a block boundary.
   The other 33 vectors were happy to agree with a wrong implementation.

This is the third instance of the same shape, after
[0001](0001-deflate-backend.md) (miniz_oxide emits valid but different bytes) and
[0002](0002-sorting-collection-tie-order.md) (`sort_unstable_by` reorders ties). In each case
the wrong choice produces output that is valid, plausible, and accepted by other tools. Only
byte comparison against the reference finds it.
