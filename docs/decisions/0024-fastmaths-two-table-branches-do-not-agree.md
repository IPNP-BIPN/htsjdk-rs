# 0024. `FastMath` has two sources for its tables, and they produce different numbers

**Status:** accepted
**Date:** 2026-07-31
**Follows:** [0023](0023-commons-math3-is-portable-where-the-jdk-is-not.md),
[0012](0012-nan-sign-bits-are-chosen-by-the-fpu.md)

## Finding

`FastMath.exp` is table-driven, and commons-math3 3.5 obtains the four tables in one of two ways,
chosen by a compile-time constant:

```java
private static final boolean RECOMPUTE_TABLES_AT_RUNTIME = false;
...
EXP_INT_TABLE_A = FastMathLiteralArrays.loadExpIntA();   // 6,175 lines of literals
// or
FastMathCalc.expint(i, tmp); FastMathCalc.splitReciprocal(tmp, recip);
```

The port was written on the computing branch, to avoid transcribing 3,550 doubles by hand, and the
oracle was asked whether the two branches agree.

**They do not.** Of the 5,050 entries in the four tables, **577 differ**, and every one of them is
in `EXP_INT_TABLE_A` or `EXP_INT_TABLE_B` — the integer table, which is the only one built through
`splitReciprocal`. The fractional tables, built by `slowexp`, agree entry for entry.

The dump reads both columns out of the same JVM: the literals from the private nested classes'
static fields, and the recomputation by invoking `FastMathCalc`'s package-private statics. So this
is not a Rust-versus-Java comparison. It is the reference disagreeing with itself.

## Why it matters

`RECOMPUTE_TABLES_AT_RUNTIME` reads like a build-time convenience: same numbers, obtained two
ways. It is not. Flipping that constant changes what `FastMath.exp` computes, and therefore what
`Gamma`, `Erf`, `NormalDistribution`, GATK's `MannWhitneyU` and every rank-sum annotation compute.
An upstream build that flipped it would move results without touching a line of algorithm.

For this port the consequence is direct: **the tables must be the literals**, because the literals
are the branch that runs upstream. A port that computed them would be defensible, self-consistent,
and wrong.

## Decision

`crates/jmath/src/fast_math_tables.rs` carries the 5,050 literals as raw bit patterns, **generated
from the oracle's own arrays** rather than typed in, and the conformance suite compares every entry
against the golden so a drift in either fails.

`fast_math_exp::recomputed_tables()` keeps the `FastMathCalc` port beside them, and the suite
asserts the disagreement's exact size in both directions. If a future commons-math3 makes the
branches agree, that assertion fails and this record can be retired; if the recomputation drifts,
it fails too, and the golden names the entries that moved.

## A third sighting of decision 0012

160 of the recomputed entries are `NaN`. The reciprocal path divides by an exponential that has
already overflowed, so the far end of the table is `inf / inf`, and IEEE 754 does not specify the
sign of the resulting NaN. x86-64 sets it, AArch64 does not.

The rule decision 0012 set applies unchanged: exempt a NaN-sign-only difference off x86-64, count
it, and assert zero exemptions on the oracle's own architecture. That this keeps arriving by new
routes — a standard deviation over one observation, a percentile over infinities, and now a
reciprocal table — is itself the argument for having written the rule down once.

## What this does not say

It does not say the literals are more correct than the recomputation, or the reverse. Neither was
checked against a higher-precision exponential here, and neither needs to be: the target is the
bits upstream produces, and upstream ships the literals.
