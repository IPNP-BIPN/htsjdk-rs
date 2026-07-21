# 0001. Deflate backend must be zlib, never miniz_oxide

**Status:** accepted
**Date:** 2026-07-21
**Risk addressed:** R2 (zlib version skew), rated Critical in the program plan

## Context

BGZF blocks are raw DEFLATE streams. htsjdk's `BlockCompressedOutputStream` compresses them
with `new Deflater(level, true)` (`nowrap = true`), which routes to `java.util.zip.Deflater`
and therefore to the zlib bundled with the JDK.

For byte-identical BAM output, the Rust side must produce the *same compressed bytes*, not
merely a valid DEFLATE stream that decompresses to the same payload. `flate2`, the obvious
Rust choice, can be backed by several implementations, and they are not interchangeable for
this purpose.

This blocks every tool that writes a BAM, so it was settled before any porting began.

## Experiment

Identical 64 KiB deterministic input (LCG-generated, `md5=b986e8a303205123db89e21b302d03b4`
verified equal on both sides), compressed with `nowrap` at levels 1, 5, 6, 9.

Reference: OpenJDK 17.0.19 (Homebrew), `java.util.zip.Deflater`.

| Level | JDK 17 (len / md5) | `flate2` miniz_oxide | `flate2` zlib |
|---|---|---|---|
| 1 | 50520 / `ba6902264075ef2c7b61e9f7ee0c283b` | 49356 / `5d76f7d9...` **differs** | 50520 / `ba690226...` **matches** |
| 5 | 49603 / `e8158c4d339fb75758482f9f90789aba` | 49783 / `2d9cfd89...` **differs** | 49603 / `e8158c4d...` **matches** |
| 6 | 49603 / `e8158c4d339fb75758482f9f90789aba` | 49783 / `2d9cfd89...` **differs** | 49603 / `e8158c4d...` **matches** |
| 9 | 49603 / `e8158c4d339fb75758482f9f90789aba` | 49783 / `2d9cfd89...` **differs** | 49603 / `e8158c4d...` **matches** |

## Decision

`flate2` is configured with the zlib backend and the default pure-Rust backend is disabled:

```toml
flate2 = { version = "1", default-features = false, features = ["zlib"] }
```

**The default configuration is wrong for this project and fails silently.** miniz_oxide
produces valid DEFLATE that decompresses correctly, so every test that checks round-tripping
passes while every byte-comparison against a reference BAM fails. Anyone adding `flate2 = "1"`
to a `Cargo.toml` here has introduced a latent bug.

## Consequences and remaining work

The match above used the macOS **system** zlib, which reports **1.2.12**, against a JDK whose
`libzip` reports **1.2.11**. They agreed on this vector, but agreement across zlib versions is
not guaranteed in general: it is a property of these versions and this input, not a promise.

Two follow-ups are therefore required before the claim is solid:

1. **Pin zlib, do not inherit it.** Depending on whatever zlib the host provides makes the
   build non-reproducible across machines. Use `libz-sys` with the `static` feature so a known
   zlib version is vendored and compiled in, and record that version in provenance.
2. **Broaden the conformance corpus.** Four levels on one 64 KiB vector is a smoke test, not
   proof. The real gate is levels 0 through 9 over a corpus of representative payloads
   (real BAM record blocks, highly compressible runs, incompressible random data, empty and
   single-byte inputs), executed against the JDK inside the pinned `linux/amd64` container
   rather than against a local macOS JDK.

Until (2) passes in the container, this decision is validated but not confirmed.

## Reproduction

See `tools/zlib-conformance/`.
