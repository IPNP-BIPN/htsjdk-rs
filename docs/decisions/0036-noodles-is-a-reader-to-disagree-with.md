# 0036. noodles is a reader to disagree with, not a source of bytes

**Status:** accepted; sets the rule for reusing third-party format crates
**Date:** 2026-08-05
**Follows:** [0001](0001-deflate-backend.md), [0034](0034-igzips-reference-c-does-not-produce-igzips-bytes.md)

## The question

The obvious thing to ask, before hand-porting 169 Java files of CRAM, is whether an existing Rust
crate already does it. Two do: [noodles](https://github.com/zaeleus/noodles) (MIT, `noodles-cram`
0.98.0 on 2026-08-03, CRAM 3.0 and 3.1) and `rust-htslib`, bindings to the C htslib. Both are
permissively licensed, both are maintained, and both are far ahead of where this port is.

So: reuse rather than write, where that is possible. What is possible was measured rather than
assumed.

## The primitives cannot be reused, and that is not a judgement call

`noodles-cram`'s public surface is whole-file: `fs`, `io::Reader`, `Record`, `FileDefinition`,
`crai`. Everything a port would want to reuse is crate-private:

```rust
// noodles-cram 0.98.0, src/io/reader.rs
pub(crate) mod num;                                  // read_itf8, read_ltf8

// noodles-cram 0.98.0, src/codecs/rans_4x8.rs
pub use self::order::Order;
pub(crate) use self::{decode::decode, encode::encode};
```

`codecs::rans_4x8` is a `pub mod` whose `encode` and `decode` are not. Compiling
`rans_4x8::encode(...)` from outside gives `error[E0603]: function 'encode' is private`. The same
holds for `rans_nx16` and for the ITF8 and LTF8 readers. There is no version of "prefer the crate"
that reaches them.

## And where it is reachable, it disagrees on purpose

This matters more than the visibility, because it would still hold if noodles exported everything
tomorrow. noodles is written from the specification; this port is written from htsjdk. Where the
two differ, the specification is not the target.

Measured on the `cram-varint` suite's own cases:

| stream | htsjdk | noodles |
|---|---|---|
| `f0 00 00 01 12` | 18 | 18 |
| `f0 00 00 01 f2` | 18 | 18 |
| `80` (truncated) | **-1** | **`UnexpectedEof`** |

The first two agree, and noodles' own unit tests assert that pair explicitly, which is worth
having: an unrelated implementation reached the same reading of the five-byte ITF8's redundant
nibble. That is the third opinion this decision is about.

The third row is the disagreement. noodles reads with `read_exact`, so a truncated stream is an
error; htsjdk's `InputStream.read()` returns `-1` and the arithmetic carries on. Both are
defensible. Only one of them is what this programme reproduces, and it is not the defensible-looking
one.

## The same shape as decision 0034

0034 found that ISA-L ships two implementations of its own compressor and they disagree, so
translating the readable one would have been a confident wrong answer. This is that lesson applied
to a *third-party* implementation rather than a second one inside the same project: agreeing with
the format is not the same as agreeing with htsjdk, and only the second is the goal.

The counter-example is decision 0001, where `flate2` **is** a dependency — because it is pinned to
the *same* vendored zlib the JDK links, so it is not a second implementation at all.

## Decision

**Third-party format crates are dev dependencies, never dependencies.** The rule has two halves.

**Hand-port the primitives, against the oracle.** ITF8, LTF8, rANS, the container model. Not because
writing them is better, but because the alternative is unreachable and would be wrong if it were
reachable.

**Use noodles as a second reader over what this port writes.** That is a claim no golden makes: a
golden says "these bytes are the bytes htsjdk produced", and an independent reader says "these bytes
are a BAM". The two coincide where a suite compares whole files and come apart wherever one compares
a rendering. `crates/htsjdk-bam/tests/noodles_cross_check.rs` is the first of these: this crate
writes a BAM, `noodles-bam` reads it back, and the header, the record count, the names, the
positions, the mapping qualities, the sequence lengths and the tag count all have to agree.

**A disagreement is a question, not a failure.** When noodles and a golden disagree, the golden
wins by construction; what the disagreement buys is knowing *where* htsjdk departs from the
specification, which is exactly the class of finding this programme exists to record.

## What this does not close

`rust-htslib` was not measured. It wraps the C htslib, a third implementation again, and its value
would be the same as noodles': a reader, not a source. Worth adding when there is a CRAM file to
read; nothing to check yet.

And the crates that produce bytes stay refused for the reason 0001 gives: `libdeflater` is a
different deflate, and htsjdk's bzip2 and LZMA come from Commons Compress's **pure-Java**
implementations rather than from libbzip2 or liblzma. Measured in the pinned container, both of
those are not even reachable — `NoClassDefFoundError`, because Commons Compress is not on the
oracle's classpath at all, so the oracle as it stands cannot produce a bzip2 or LZMA CRAM block to
compare against.
