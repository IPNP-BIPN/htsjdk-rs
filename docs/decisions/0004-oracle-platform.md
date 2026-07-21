# 0004. The emulated linux/amd64 container is a valid oracle; the plan's AVX assumption was wrong

**Status:** accepted
**Date:** 2026-07-21
**Corrects:** program plan, "Oracle and the bit-identity contract" and risk R4

## What the plan assumed

> Docker's amd64 emulation on Apple Silicon does not expose AVX, so GKL native paths can
> silently fail to load inside the container and produce goldens that match no real machine.

That reasoning drove the decision to produce all authoritative goldens on real x86-64 CI. It is
wrong, and it was never measured.

## What is actually true

Measured on this machine, Docker Desktop 29.6.1, macOS 25.5.0, Apple Silicon.

CPU flags visible inside `--platform linux/amd64`:

```
avx avx2 bmi2 pclmulqdq sse4_2      (model: VirtualApple @ 2.50GHz)
```

Intel GKL 0.8.11 native loading, on the pinned base image
`eclipse-temurin@sha256:068a8f9ae4b74d9a20de3ecca771ba1a6437f1d7a8a8ad6deaf9dbdd2274397a`:

| Platform | `os.arch` | `usingIntelDeflater` | `usingIntelInflater` |
|---|---|---|---|
| emulated `linux/amd64` | `amd64` | **true** | **true** |
| native `linux/arm64` | `aarch64` | false | false |

Rosetta 2 translates AVX and AVX2, and GKL's x86 natives load and run. The arm64 result is the
expected one: GKL ships x86-only shared objects, which is precisely the gap
[broadinstitute/gatk#9384](https://github.com/broadinstitute/gatk/pull/9384) exists to close.

Goldens confirmed identical between the local macOS JDK and the container JDK: all 70 zlib
vectors and all 35 BGZF vectors, zero mismatches, despite the two JDKs bundling different zlib
versions (macOS 1.2.11, container 1.3.2).

## Decision

**The pinned `linux/amd64` container is the authoritative oracle, and it runs locally.** Real
x86-64 CI is no longer a prerequisite for producing goldens. This removes a slow remote
dependency from the inner loop.

## A trap worth recording

The first run of this probe reported `usingIntelDeflater = false` on **both** platforms, which
would have "confirmed" the plan's assumption. The cause was a missing `commons-io` on the test
classpath: GKL extracts its native library from the jar using
`org.apache.commons.io.FilenameUtils`, and without it the load fails with a warning and falls
back silently.

The failure mode is the point. GKL logs a `WARNING` and degrades to the JDK deflater rather
than throwing. An oracle that does not *assert* provider state will happily emit goldens from a
degraded configuration, and they will look entirely normal. This is exactly why the runner
asserts rather than logs.

## Limits of this result

1. **No AVX-512.** Only `avx` and `avx2` are exposed. On a real AVX-512 host,
   `FASTEST_AVAILABLE` would select a different PairHMM implementation. The oracle contract must
   therefore pin the PairHMM implementation explicitly rather than rely on the default. This is
   already an open question in the plan and remains one.
2. **GKL loading is not GKL byte-equality.** That GKL runs under translation does not prove its
   output matches a native x86 host. Integer SIMD should translate exactly, but "should" is not
   "measured". Any future work that depends on GKL-produced bytes, in particular the
   ISA-L/igzip exact-deflate module, must be cross-checked on real x86-64 hardware.
3. The current oracle contract pins `--use_jdk_deflater`, so GKL is not on the golden path
   today. Its availability is still asserted, because a change in that state must break the
   build rather than silently change the meaning of every golden.
