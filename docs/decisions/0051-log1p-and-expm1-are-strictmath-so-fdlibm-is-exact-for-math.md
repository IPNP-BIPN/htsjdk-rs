# 0051. `log1p` and `expm1` are `StrictMath`'s, so FDLIBM is exact for `Math` too

**Status:** accepted; `Math.log1p`, `Math.expm1`, `StrictMath.log1p` and `StrictMath.expm1` ported
**Date:** 2026-11-25
**Follows:** [0025](0025-fdlibm-is-portable-and-is-the-worse-stand-in-for-the-intrinsic.md),
[0027](0027-pow-has-a-bound-now-and-it-is-the-same-one-as-exp.md)

## The question

Decisions 0025 and 0027 ported FDLIBM's `exp` and `pow` and found them exact against `StrictMath`
and one ulp from `Math`, because `Math.exp` and `Math.pow` are HotSpot intrinsics. gatk-rs still
called the host's `ln_1p` and `exp_m1` wherever GATK calls `Math.log1p` and `Math.expm1`
(`NaturalLogUtils.log1mexp`, the somatic clustering model, the PairHMM's `log10OneMinusPow10`,
DRAGstr), at 99.5889% agreement on the corpus. Is `log1p` in the same position as `exp`?

It is not. In the JDK, `Math.log1p` is `return StrictMath.log1p(x);` and is not an intrinsic
candidate, and the same holds for `expm1`. So on every JVM and every architecture `Math.log1p` is
FDLIBM, and an FDLIBM port is the port.

## What was done

`fdlibm/s_log1p.c` (FDLIBM 5.3, 1995) and `fdlibm/s_expm1.c` (2004) ported into
`crates/jmath/src/strict_log1p.rs` and `strict_expm1.rs`, from netlib's copy, with Sun's notice
preserved. Not from OpenJDK's tree, whose copies are GPL2 (decision 0014). Exposed as
`strict_math::{log1p, expm1}` and, because `Math` delegates, as `math::{log1p, expm1}`, so a call
site names the class its Java names.

## Measured

| claim | result |
|---|---|
| the corpus's `Math` and `StrictMath` columns agree | **44,996 / 44,996** for each function |
| `math::log1p` == `Math.log1p`, `math::expm1` == `Math.expm1` (corpus) | **44,996 / 44,996** each |
| the same against `StrictMath` (corpus) | **44,996 / 44,996** each |
| 4,000,000 seeded random inputs (uniform in (-1, 1), log-uniform over every exponent, random bit patterns, [-800, 800], [-60, 60]), against `Math` and `StrictMath`, on `gatk-rs-oracle:4.7.0.0` (17.0.19, amd64) | **0 differ** |
| the same inputs on `eclipse-temurin:17-jdk` (17.0.21, native aarch64), `-Xbatch` and `-Xint` | **0 differ** |

The 4,000,000-point check is a one-off harness (`Check.java` writing raw bits, a Rust reader
comparing them), not a committed corpus: the committed corpus is what CI re-checks.

The same harness answers a question gatk-rs needed for `Math.pow`: over 4,000,000 pairs, arm64
HotSpot's `Math.pow` equals `StrictMath.pow` on **every** pair, compiled and interpreted, so aarch64
has no `pow` stub, as it has no `exp` stub. On the amd64 oracle the two differ on 225,671 of the
same pairs, which is the intrinsic of decision 0007.

## What this changes

- `unported_functions_are_not_yet_exact` loses `log1p` and `expm1`; they join `sqrt`, `log` and
  `log10` in the exact-match test, and a new test pins them against the `StrictMath` column too.
- Every gatk-rs call site that reaches `Math.log1p` or `Math.expm1` can be exact on every host,
  which no `exp` or `pow` call site can be on x86-64.
