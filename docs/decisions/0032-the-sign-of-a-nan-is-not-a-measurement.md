# 0032. The sign of a NaN is not a measurement, and 160 of the 577 were only that

**Status:** accepted; corrects the count in [0024](0024-the-references-two-branches-disagree.md) and the reasoning in [0012](0012-nan-sign.md)
**Date:** 2026-08-05

## How it surfaced

A CI run failed on a branch that touched neither jmath nor its harness:

```text
FAIL fastmath-exp/FastMathTablesDump: compared=14405
  first diff at line 5051
    real     =table  RECOMPUTED_EXP_INT_TABLE_A  1  9221120237041090560
    committed=table  RECOMPUTED_EXP_INT_TABLE_A  1  -2251799813685248
```

Both are NaN. `0x7FF8000000000000` and `0xFFF8000000000000` differ in the sign bit and nowhere
else.

## Why the existing reasoning did not cover it

Decision 0012 knew about the NaN sign and attributed it to the architecture. The Rust test carried
the conclusion directly:

```rust
if cfg!(target_arch = "x86_64") {
    assert_eq!(nan_sign_exemptions, 0,
        "on x86-64 there is nothing to exempt; the FPU produces the same NaN as the oracle");
}
```

That is a claim about x86-64, and x86-64 has now produced both. The entries come from
`FastMathCalc`'s reciprocal path, which divides an overflowed exponential by itself. `inf / inf`
computed by the FPU gives x86's default quiet NaN, which is the negative one; folded by the JIT it
gives Java's canonical NaN, which is positive. Same architecture, same container, same code. Which
one appears depends on how far the JIT got, and that is not a property of the machine either.

**The sign of a NaN is not a property of anything a port could reproduce.** Not of the arithmetic,
not of the CPU, not of the target.

## The correction to 0024

0024 says the reference's two branches disagree on **577 of 5,050** entries, and that number is
asserted in a test. Recomputed with the NaN sign out of the comparison:

| | entries |
|---|---|
| differ, raw bits | 577 |
| differ, NaN canonicalised | **417** |
| difference that was only a NaN sign | 160 |

So 160 of the 577 were the same NaN written two ways. The reference's branches disagree on 417
entries, and 0024's conclusion stands: the port carries literals rather than recomputing. Only the
size of the disagreement was inflated, by a comparison that asserted a bit neither branch chose.

## Decision

**The golden is not regenerated.** It is a faithful record of what the oracle held on the host that
produced it, and replacing it would trade a true record for a differently-true one while spending a
CI round trip. What changes is the comparison: `tools/conformance/manifest.json` declares
`canonicalise_nan` on the `fastmath-exp` suite, with the reason, and `compare.py` is the only code
that applies it — the same route every other declared correction takes.

The rule collapses a field only when its exponent is all ones **and** its mantissa is non-zero, so
both infinities and every finite value still travel bit for bit. No index or count in any dump comes
within nine quintillion of the smallest NaN pattern.

On the Rust side the architecture-conditional exemption is gone, along with the constant that
counted it. The comparison canonicalises and asserts everything else.

## What this says about the other goldens

This one was caught because a suite regenerates its golden on every CI run and two runs landed on
different hosts. A suite whose golden is never re-derived would carry the same defect silently.
That is the argument for `oracle-backed` over `unchecked`, made by a failure rather than by an
assertion, and it is worth one more look at the one suite still marked `unchecked`.
