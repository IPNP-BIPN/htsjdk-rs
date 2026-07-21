# 0006. `log` and `log10` are correctly rounded, so round rather than port the intrinsic

**Status:** accepted
**Date:** 2026-07-21
**Refines:** [0005](0005-java-math-has-three-implementations.md)
**Result:** `Math.log` and `Math.log10` are now bit-identical over the whole corpus

## The question

Decision 0005 established that Rust's libm differs from `java.lang.Math` by up to 1 ULP and
that a hand port is required. The obvious reading is that HotSpot's x86 intrinsic must be
reproduced instruction for instruction, which is a large and unpleasant job.

Before committing to that, a cheaper hypothesis was worth testing: **what if the intrinsic is
simply correctly rounded?** If so, any correctly-rounded implementation matches it, and no
algorithm needs porting at all.

## The experiment

For every point where the host libm disagreed with the JVM, the true value was computed at
300-bit precision with MPFR and rounded once. That tells us which of the two is right.

| fn | divergent points | `Math` correct | host libm correct |
|---|---:|---:|---:|
| `log10` | 2 | **2** | 0 |
| `log` | 4 | **4** | 0 |
| `exp` | 13 | 6 | 7 |
| `pow` | 252 | 202 | 50 |

Then the stronger claim, on a random sample of *all* corpus points rather than just the
divergent ones:

| fn | sampled | `Math` == correctly rounded |
|---|---:|---:|
| `log10` | 3,358 | **100.0000%** |
| `log` | 3,293 | **100.0000%** |
| `exp` | 6,000 | 99.9833% |

## Decision

**The strategy splits by function, and the split is worth a lot of work.**

- `log`, `log10`: the target is **correct rounding**. Implemented with double-double
  arithmetic, no algorithm port. **Now 44,996/44,996 exact against `java.lang.Math`.**
- `exp`, `pow`: measurably **not** correctly rounded, so correct rounding is the wrong target
  and HotSpot's algorithm does have to be ported. **`exp` is now done and exact**
  (44,996/44,996); see the note below.
- `sqrt`: already exact everywhere, mandated by IEEE-754.

Testing the cheap hypothesis first turned the two highest-priority functions from an assembly
port into a numerics exercise.

## Implementation

`ln(x)` reduces `x = m * 2^e` with `m` centred on 1 so that `s = (m-1)/(m+1)` satisfies
`|s| <= 0.1716`, then sums the atanh series in double-double. `log10` multiplies by a
double-double `log10(e)`.

Exact powers of ten are special-cased. The series lands within an ulp of the integer but not
necessarily *on* it, and `log10(1000)` rendering as `2.9999999999999996` is the kind of value
that survives all the way into a report column.

## The bug that mattered, and why it was invisible

The first implementation reached 99.9956%, better than the libm baseline but not exact. Four
points failed, each by exactly 1 ULP.

The cause was one line: series coefficients were applied as `term * (1.0 / (2k+1))`. `1.0/3.0`
is not representable in `f64`, so that reciprocal carries ~2^-53 relative error, which lands
near **2^-60** once weighted by the largest term.

MPFR showed the failing points needed only **62 to 64 bits** to resolve, not the 106 the
double-double provides. So the failures were not a precision ceiling; they sat exactly in the
band the rounded reciprocal had corrupted. Dividing by the exactly-representable odd integer
instead fixed all four and reached 100%.

The shape is worth noting because it recurs. A rounded constant costs accuracy that is
invisible on ordinary inputs and fatal on precisely the hard-to-round ones, which is to say on
exactly the points a correctness claim depends on. "99.9956% agreement" reads like a near miss
caused by fundamental limits; it was a fixable defect, and the diagnostic that distinguished
the two was measuring how many bits the failures actually needed.


## Addendum: `exp` ported, 44,996/44,996 exact

The other branch of the split has been walked. `Math.exp` reproduces HotSpot's Intel LIBM
intrinsic from `macroAssembler_x86_exp.cpp` at `jdk-17-ga`, and matches on every corpus point.

Two things made it work, and both are the opposite of how one would naturally write it:

**The SIMD lane structure is load-bearing.** The original is packed-double code in which the
two lanes of each register carry *different* polynomial terms, combined at specific points:
lane 0 accumulates `r^5 * (c2 + c4*r)`, lane 1 accumulates `r^3 * (c3 + c5*r)`, and the `r^2/2`
term arrives separately. Rewriting that as a single Horner evaluation is algebraically
identical and produces different bits, because every intermediate rounds to `f64`.

**The final step is a multiply then an add, not an FMA.** The original emits `mulsd` followed
by `addsd`, so the product rounds before the sum. Using `mul_add` would be more accurate and
therefore wrong.

The table was generated from the source rather than transcribed: 64 entries of two values each
is precisely the kind of copying where one wrong nibble yields a plausible wrong answer.

`pow` remains, and has its own intrinsic.
