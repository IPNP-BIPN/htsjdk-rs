# htsjdk-rs

Pure-Rust port of [htsjdk](https://github.com/samtools/htsjdk), targeting **byte-identical**
output against a pinned reference build. Work in progress.

> **This is not the official htsjdk.** It is an independent reimplementation, not affiliated
> with or endorsed by the Broad Institute or the samtools organization.

## Why this exists, and why not `rust-htslib`

[`rust-htslib`](https://github.com/rust-bio/rust-htslib) provides FFI bindings to **HTSlib**,
the **C** library from samtools. htsjdk is the **Java** library, an independent implementation
of the same file formats. They are not two names for one codebase.

That distinction is the whole reason this crate exists. BGZF leaves block size, deflate level,
and block-fill policy to the implementation, so byte equality between HTSlib and htsjdk output
is guaranteed by nothing. GATK and Picard are built on htsjdk, so reproducing their output
byte-for-byte requires porting **htsjdk specifically**.

`rust-htslib` remains excellent for reading, and this project uses it as an independent
cross-check oracle.

## Reference version

Ported from htsjdk `4.2.0`, the version pinned by GATK 4.6.2.0's `build.gradle`.

## Scope

| Module | Status |
|---|---|
| BGZF reader/writer (`BlockCompressedOutputStream` semantics) | **byte-identical**, 35 goldens |
| BAM record codec (byte level, incl. `bin` and tag type promotion) | **byte-identical**, 136 goldens |
| SAM text header (`SAMTextHeaderCodec.encode`) | **byte-identical**, 9 goldens |
| BAM file writer (header + dictionary + framing) | **byte-identical**, 4 whole files |
| SAM text records | planned |
| BAM index | planned |
| VCF / tribble index | planned |
| CRAM | planned, later phase |

## Bit-identity contract

Output is compared byte-for-byte against goldens produced by the pinned reference running in
a digest-pinned `linux/amd64` container on JDK 17, on real x86-64 CI. Fields that are
legitimately allowed to vary (timestamps, version strings, command lines in headers) are
canonicalized under explicitly declared rules, and every comparison records which fields were
compared raw and which canonicalized.

Where a value cannot be matched exactly, it is quarantined and reported with its measured
divergence rate, and the affected output is described as **bio-identical** rather than
**bit-identical**. That vocabulary comes from
[broadinstitute/gatk#9384](https://github.com/broadinstitute/gatk/pull/9384), which
established that Java `Math.log/exp/pow` differ by roughly 1 ULP across CPU architectures
while only `StrictMath` is portable.

## Part of a three-repository program

| Repo | Ports | Depends on |
|---|---|---|
| `htsjdk-rs` | htsjdk 4.2.0 | (none) |
| `picard-rs` | Picard 3.4.0 | `htsjdk-rs` |
| `gatk-rs` | GATK 4.6.2.0 | `picard-rs`, `htsjdk-rs` |

The topology mirrors upstream, including the direction of the dependencies.

## License

MIT, matching the htsjdk sources this ports. See `LICENSE`.

One caveat, taken from htsjdk's own README: htsjdk is **not uniformly MIT**. Licensing is
per-file, and notably the CRAM code is Apache License 2.0. Since CRAM is in scope for this
port, code derived from it will carry Apache 2.0 rather than MIT, and will say so in the file
header. Check the notice on the reference file before porting anything new.
