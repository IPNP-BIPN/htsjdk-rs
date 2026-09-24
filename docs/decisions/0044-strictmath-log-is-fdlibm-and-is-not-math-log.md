# 0044. StrictMath.log is FDLIBM, and it is not Math.log

**Status:** accepted
**Date:** 2026-09-24
**Extends:** [0006](0006-correct-rounding-is-the-target-for-log-and-log10.md), [0025](0025-fdlibm-is-portable-and-is-the-worse-stand-in-for-the-intrinsic.md)

## Why a second log

`jmath::math::log` is correctly rounded, because `Math.log` measured correctly rounded on every
point of the corpus (0006). `StrictMath.log` is a different function: it is FDLIBM's `e_log.c`,
within one ulp, and on the corpus it differs from `Math.log` on 186 of 44,996 points. A call site
that reaches `StrictMath.log` cannot use the correctly rounded one.

gatk-rs needs it for `java.util.Random.nextGaussian`, whose polar method calls `StrictMath.log`.
`QualByDepth` adds a Gaussian jitter to every QD above 35, so a GenotypeGVCFs QD is reproducible
only through the fdlibm logarithm.

## The answer

`jmath::strict_log::log`, exposed as `jmath::strict_math::log`, is a transcription of
`__ieee754_log` from FDLIBM 5.3. The licence is the one 0025 already relied on for `exp`: Sun's
notice grants copying and modification provided it is preserved, and the module preserves it. The
operations keep the C function's order and parentheses.

`strict_log_is_strictmath` asserts it bit-identical to the corpus's StrictMath column on all
44,996 points. It also asserts that the StrictMath and Math columns differ somewhere, so the test
cannot pass by the port calling the correctly rounded `log`.
