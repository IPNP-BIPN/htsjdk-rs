# 0033. The whole of GKL cuts at SSE4.2, and nowhere else

**Status:** accepted; sharpens [0030](0030-gkls-default-bam-bytes-depend-on-sse42.md) and closes the question [0031](0031-gkl-levels-one-and-two-are-isal-level-one.md) left open
**Date:** 2026-08-05

## Two questions, one answer

Decision 0030 established that GKL's zlib emits different bytes with and without SSE4.2, by forcing
the flag in Intel's source. That says a boundary exists; it does not say where the boundary is, and
"depends on the CPU" is a much larger claim than the measurement supports. 0031 then left the same
question open on the igzip side.

Both are answered by running the real libraries under emulated CPUs of each generation rather than
by patching either of them.

## igzip

`isal_deflate_stateless` at level 1, the configuration 0031 identified, on the four fixtures:

| CPU | acgt | verdict |
|---|---|---|
| `core2duo` (SSE4.1, no SSE4.2) | 19749 | differs from GKL on 3 of 4 fixtures |
| `Nehalem` (SSE4.2, no AVX) | 19044 | **identical to GKL, all four** |
| `SandyBridge` (AVX) | 19044 | identical |
| `Haswell` (AVX2) | 19044 | identical |
| `Skylake-Server` | 19044 | identical |
| a `noarch` build, base C only | 19749 | differs, and equals `core2duo` |

So the pre-SSE4.2 answer is not a fourth kernel: it is the pure C fallback, the same bytes a build
with no assembly at all produces. Every assembly kernel from SSE4.2 upward agrees with every other.

## Intel's zlib

The same sweep, on the same fixture, at the levels that matter:

| CPU | L3 | L4 | L5 | L6 |
|---|---|---|---|---|
| `core2duo` | 18411 | 18451 | 18184 | 17897 |
| `Nehalem` | 18723 | 18059 | **17952** | 17897 |
| `Haswell` | 18723 | 18059 | **17952** | 17897 |
| native | 18723 | 18059 | **17952** | 17897 |

17952 is what GKL recorded at level 5. The cut is in the same place and for the same reason: the
CRC-32C hash is gated on `x86_cpu_has_sse42`. Level 6 is unaffected because its hash covers three
bytes, which is what the multiplicative one covers too.

## What this changes

**"GKL's output depends on the CPU" is true and imprecise.** It depends on one bit of CPUID.
Above SSE4.2 there is one behaviour, and AVX, AVX2 and AVX512-capable hosts all produce it; below
it there is a second, which both backends reach by falling back rather than by dispatching. Intel
shipped SSE4.2 in Nehalem in 2008 and AMD in Bulldozer in 2011, so the second class is historical
rather than hypothetical-but-current.

For the port that means a single target, stated once:

> `crates/gkl-deflate` reproduces GKL as it behaves on a CPU with SSE4.2. On a CPU without it, GKL
> produces different bytes and so does `Flavour::Gkl { sse42: false }`, which is checked against
> Intel's source rather than against the library.

**And it retires a worry rather than confirming it.** 0028 asked whether igzip's AVX2 kernels emit
different bytes from its SSE ones, and treated that as the thing that could put this milestone
beside `Math.pow`. They do not. The kernels agree; only the assembly-versus-C boundary moves the
output, and no x86-64 host GATK runs on is on the far side of it.

## What was not measured

**AVX512 was requested and not delivered.** `Skylake-Server` under QEMU's TCG drops `avx512f` and
five more bits with a warning, so that row is an AVX2 run wearing a Skylake name. It agrees with
the others, which is evidence about AVX2 and nothing about AVX512. A host with real AVX512 remains
untested, and the `igzip-portability` job will keep saying so as long as its two columns advertise
the same flags.
