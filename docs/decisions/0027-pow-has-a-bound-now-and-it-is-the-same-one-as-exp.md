# 0027. `pow` has a bound now, and it is the same one as `exp`

**Status:** accepted; `StrictMath.pow` exact, `Math.pow` bounded at 1 ulp
**Date:** 2026-08-04
**Follows:** [0005](0005-java-math-has-three-implementations.md),
[0007](0007-pow-may-not-be-portable-across-x86-cpus.md),
[0025](0025-fdlibm-is-portable-and-is-the-worse-stand-in-for-the-intrinsic.md)

## What was missing

Decision 0007 deferred `Math.pow` for a specific reason: HotSpot's intrinsic uses `rcpps`, the
packed approximate reciprocal, at six sites and — unlike `log` — does not refine the approximation
away. It then tested the hazard that raised, that `pow`'s bits might depend on the silicon, and
found it did not materialise: 404,883 points regenerated on real AMD EPYC matched the emulated
corpus exactly.

What no record held was a **bound**. What existed was a rate — 99.9378% agreement with the host
libm — and a rate says how often, not how far. Every call site reaching `pow` therefore had nothing
to reason with. gatk-rs's `AllelePseudoDepth` is the concrete case: `calculateWeights` reaches
`Math.pow` twice, its conformance suite is clean, and nobody could say whether that was a property
or a coincidence of the corpus.

## What was done

`fdlibm/e_pow.c` ported into `crates/jmath/src/strict_pow.rs`, exactly as decision 0025 ported
`e_exp.c`: same source, same Sun notice preserved, same reason it is allowed at all. `StrictMath`
is *specified* to be fdlibm, so the first question is exactness rather than closeness.

## Measured

Over the existing corpus, **404,964 `pow` points**, on the pinned oracle:

| claim | result |
|---|---|
| `strict_pow` == `StrictMath.pow` | **404,964 / 404,964** |
| `Math.pow` vs the host libm | 99.9378% |
| `Math.pow` vs fdlibm | 98.5317% |
| **worst fdlibm divergence** | **1 ulp**, at `pow(2, -0.5)` |
| divergences that are not a last-bit difference | **0** |

The last row is asserted, not just reported: a sign flip or a non-finite mismatch would mean no ulp
bound covers the function, and the suite names any such point rather than folding it into a rate.

## The lesson from 0025, repeated exactly

Decision 0025's finding about `exp` was that the permissively licensed implementation is the
**worse** stand-in for the intrinsic by rate — 98.6443% against the system libm's 99.9711% — and
that the useful result was the bound rather than the win. `pow` reproduces that shape line for line:

| | libm | fdlibm | worst fdlibm divergence |
|---|---|---|---|
| `exp` | 99.9711% | 98.6443% | 1 ulp |
| `pow` | 99.9378% | 98.5317% | 1 ulp |

So the agreement rate was the wrong statistic for both. It measures how often two implementations
happen to round the same way, which is a fact about the corpus. The bound is a fact about the
function.

**Which implementation a call site should use does not follow from the rate.** It follows from
host-independence: fdlibm is fixed and the system libm is whatever the host ships, so a port that
must produce the same bytes on every machine takes fdlibm and its 1 ulp, exactly as
`NaturalLogUtils` does for `exp` in gatk-rs. A call site that only needs to be close can keep the
libm and be closer more often. Both are defensible; picking by rate alone is not.

## What this changes

- **`StrictMath.pow` call sites are now exactly reproducible.** They were not before;
- **`Math.pow` call sites are bounded** at 1 ulp against a fixed, permissively licensed
  implementation, which is what they had for `exp` and did not have here;
- **decision 0007's "deferred" status is superseded** for the part that mattered. `pow` is still
  not reproducible bit-for-bit against `Math.pow` — that would need HotSpot's GPL2 source, and
  gatk-rs #71's argument applies unchanged — but it is no longer an unmeasured gap.

## What it does not change

`pow(x, 2)` is still exactly `x * x`, which is why `Histogram.getStandardDeviation` was never
blocked (0007, addendum 2). Nothing about that measurement is affected, and the fast path in this
port returns `x * x` for `y == 2` because fdlibm does.

## Verification

* `cargo test -p jmath` — `strict_pow_is_strictmath` over the whole corpus, and
  `the_pow_gap_is_measured_for_both_stand_ins`, which prints both rates and the worst distance.
* The gap test asserts fdlibm does **not** reach 100% against `Math.pow`. If it ever does, the gap
  these records describe has closed and they are wrong rather than merely stale.
* `tools/audit/provenance.py` resolves the new `Ported from fdlibm/e_pow.c` claim to a
  licence-compatible source.
