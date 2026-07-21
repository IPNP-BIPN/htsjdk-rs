# 0012. NaN sign bits are chosen by the FPU, so the port is not bit-identical to itself across architectures

**Status:** accepted; scoped, measured, and gated on x86-64
**Date:** 2026-07-21
**Follows:** [0005](0005-java-math-has-three-implementations.md)

## Finding

Porting `Histogram` produced 338 statistics compared against htsjdk as raw double bits. 337
matched. The one that did not:

```
single/standardDeviation: ours=NaN (7ff8000000000000) htsjdk=NaN (fff8000000000000)
```

A standard deviation over a single observation divides by `count - 1`, which is `0.0 / 0.0`.
Both sides produce a quiet NaN. They disagree on its **sign bit**.

Measured directly, in both languages, on both architectures:

| | `0.0 / 0.0` | `sqrt(0.0 / 0.0)` |
|---|---|---|
| Java on `amd64` | `fff8000000000000` | `fff8000000000000` |
| Rust on `aarch64` | `7ff8000000000000` | `7ff8000000000000` |

x86-64 returns the "real indefinite" QNaN, which has the sign bit set. AArch64 returns the
default QNaN with the sign bit clear. Neither language chooses this; both return what the FPU
gives. IEEE 754 leaves the sign of a NaN produced by an invalid operation unspecified.

## Why this is not a porting defect, and what it is instead

**The same Rust source compiled for x86-64 matches the golden exactly.** There is nothing in
the port to fix.

What it establishes is something else, and it is worth stating plainly: **htsjdk-rs is not
bit-identical to itself across architectures.** That is the same shape as the finding this
whole program grew out of — PR #9384 showed GATK's floating-point results differ by ~1 ULP
between architectures — except one level lower. This one is not in an algorithm or a library.
It is in the silicon's choice of a bit that IEEE 754 declines to specify.

## Scope: where it can and cannot be observed

It only arises from an **invalid operation** producing a fresh NaN: `0/0`, `inf - inf`,
`inf * 0`, `sqrt(negative)`. A NaN that is merely propagated keeps its payload, so the divide is
where the sign is decided.

Whether it reaches an artefact depends on the output path:

| path | observable? |
|---|---|
| Picard metrics file | **no** — `FormatUtil` renders every NaN as `?`, sign discarded (decision 0011) |
| GATK report | **no** — same formatter |
| BAM `f` tag | **yes** — the 4 raw bytes are written, sign included |
| VCF `Float` field | text, so no |

So for the metrics archetype currently being ported, the difference is unobservable in output.
For a tool that writes a NaN into a float tag, it is not.

## Decision

**Gate on x86-64 and exempt elsewhere, explicitly and countably.**

The oracle platform is x86-64 (decision 0004), so x86-64 is where the claim is made. The
conformance test compares raw bits and:

- on `x86_64`, permits **no** exemption at all, and asserts the exemption list is empty;
- on any other architecture, permits a difference that is *only* the sign bit of a NaN, and
  asserts the exempted set equals exactly one named case.

A second test pins the cause rather than asserting it, checking `0.0 / 0.0` against the
expected constant for the architecture it is compiled for. If a future toolchain or target
changed that, the cause test fails first and names the reason, instead of the conformance test
failing with an unexplained bit.

The alternative was to force the sign to the x86 value inside the port. That was rejected: it
would make an ARM build disagree with what an ARM build actually computes, hiding the finding
inside the code that discovered it, and it would have to be applied at every site that can
generate a NaN rather than at the one the corpus happened to reach.

## Consequence for the bit-identity claim

The claim's wording needs one more qualifier, alongside the JDK deflater and the pinned locale:
bit-identity is asserted **on x86-64**. On other architectures the port is bit-identical except
for the sign of NaNs arising from invalid operations, which is a strictly smaller exemption
than PR #9384's ~1 ULP and is unobservable in every text output format.

Worth noting for the ARM port that motivated all of this: an arm64 GATK has the same property,
and nothing in PR #9384 addresses it, because Java's own `Double.toString` renders every NaN as
`"NaN"` and hides it too.
