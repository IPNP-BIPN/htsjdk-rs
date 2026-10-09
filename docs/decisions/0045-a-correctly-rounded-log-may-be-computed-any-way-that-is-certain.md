# 0045. A correctly rounded `log` may be computed any way that is certain

**Status:** accepted
**Date:** 2026-10-09
**Refines:** [0006](0006-correct-rounding-is-the-target-for-log-and-log10.md)
**Result:** `jmath::math::log` and `log10` return the same bits 24 times faster

## The question

Decision 0006 made the target "round the true result once", and `ln_dd` reaches it with a 22-term
atanh series in double-double: 22 double-double divisions per call, about 277 ns. gatk-rs's speed
baseline (Milestone S there, #110) then found that series under 94% of the slowest tool's run
(`ModelSegments`, 18 times a warmed-up JVM) and most of the second (`VariantRecalibrator`). The
reference pays a hardware intrinsic for the same call.

The question is whether the series can be avoided without moving a bit.

## The answer is in 0006's own premise

A correctly rounded function has exactly one answer per input. Any algorithm that provably
returns it returns what `ln_dd` returns; the algorithm is not part of the contract, the value is.
The standard way to be fast and provably correct is Ziv's: compute quickly with an error bound,
answer only when the whole interval rounds to one `f64`, and fall back to the accurate method when
it does not.

## What was built (`crates/jmath/src/log.rs`)

- **Fast phase** (`ln_fast`). `x = 2^e * m`; a 128-entry table gives `r` near `1/m` and `-ln(r)`
  in double-double, the latter computed once by `ln_dd` itself, so no constant was typed in.
  `z = m * r - 1` is formed exactly (`two_prod`, then a Sterbenz subtraction), `|z| < 2^-7`, and
  `ln(1 + z)` is `z - z^2/2` in double-double plus an `f64` tail to `z^10`. Within 2^-7 of one,
  `z = x - 1` exactly and nothing else, because below one the decomposition would cancel two
  terms near `ln(2)`. The claimed bound is `2^-60 |z| + 2^-90 (|head| + |z|)`, several orders wider
  than the error the analysis in the code finds.
- **Test** (`round_if_certain`). The answer is returned only if `value - 2 * bound` and
  `value + 2 * bound` round to the same `f64`. Twice the bound is what makes it `ln_dd`'s answer
  and not only the true value's rounding: the true value is within `bound` of `value`, and
  `ln_dd`'s is within about 2^-100 of the true one, so both lie in the tested interval.
- **Fallback.** Anything undecided, and every subnormal, goes to `ln_dd` unchanged.

## Evidence

- 200 million inputs, half within 2^-7 of one: 199,880,214 decided by the fast phase (99.94%),
  **none differing** from `ln_dd` (`long_tests::two_hundred_million_points_agree`, run with
  `--ignored`, kept in the tree).
- On every CI run: two million inputs from every region the fast phase treats differently, the
  same agreement for `log` and `log10`, and a test that the fast phase's double-double lies within
  its claimed bound of `ln_dd`'s.
- The `jmath` suite, including the corpus of hard-to-round points 0006 was built on: unchanged.
- 277 ns per call before, 11.6 ns after, on an M-series core. In gatk-rs, `ModelSegments`' row
  goes from 5.0 s to 0.68 s and `VariantRecalibrator`'s from 0.34 s to 0.10 s natively, with
  every output identical.

## What this does not license

A faster function that is merely *close*: the fast phase never answers on its own authority, only
where its answer is provably the one the slow phase would give. `exp` and `pow` are not correctly
rounded in the JVM (0005, 0007) and cannot be treated this way.
