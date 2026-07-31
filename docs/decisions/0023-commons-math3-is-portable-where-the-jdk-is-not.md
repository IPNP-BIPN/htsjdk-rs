# 0023. commons-math3 may be transcribed; the JDK may not

**Status:** accepted
**Date:** 2026-07-31
**Follows:** [0013](0013-the-last-divergences-are-blocked-by-a-licence-not-by-difficulty.md),
[0014](0014-math-exp-was-withdrawn-it-was-a-gpl2-transcription.md)

## Finding

Decision 0014 withdrew a working, bit-exact `Math.exp` because it was a transcription of
GPL2-only HotSpot source and this crate ships under MIT. Decision 0013 generalised that: what
remains between this crate and a complete `jmath` is a licence, not difficulty.

That conclusion is correct for `java.lang.Math` and wrong as a statement about "Java math".
Three implementations sit behind the functions GATK and Picard call, and they do not share a
licence:

| implementation | licence | transcribable under MIT |
|---|---|---|
| `java.lang.Math` / `StrictMath` (JDK, HotSpot intrinsics) | GPL2, no Classpath Exception on the intrinsic sources | **no** |
| commons-math3 `FastMath`, `Percentile`, `Gamma`, … | Apache 2.0 | **yes** |

Apache 2.0 grants the right to reproduce and distribute derivative works "in any medium, with
or without modifications", under terms of the recipient's choosing, provided the notice and
attribution requirements are met. MIT is such a choice. The GPL2 in the JDK's intrinsic sources
is not, which is the whole content of decision 0014.

## Why this matters now

The roadmap records "the jmath corpus must reach 100%, and it blocks most of GATK's 54
annotations". Both halves needed correcting, and this record corrects the first:

- **100% of the current corpus is unreachable by construction.** Its columns are
  `java.lang.Math`, and the functions still divergent there (`exp`, `pow`, `sin`, `cos`, `cbrt`,
  `log1p`, `expm1`) are exactly the ones whose exactness would require transcribing GPL2 source.
  A target that can only be met by violating a licence is not a target.
- **Much of what a call site actually reaches is Apache 2.0.** `MathUtils.median` goes through
  commons-math3 `Percentile`, and `MathUtils` finishes several of its own functions with
  `FastMath`. Those can be ported exactly, today, with no licence problem.

## Decision

`org.apache.commons.math3` joins the allowed sources in `tools/audit/provenance.py`, beside
htsjdk, Picard and GATK. The JDK entries stay forbidden.

The measure of `jmath` stops being "the corpus reaches 100%" and becomes **"every function a
ported call site reaches is exact, and every function that cannot be is named at the call site
that reaches it"**. The corpus keeps its role: it measures the gap the licence costs, which is
what decision 0014 asked it to do.

## Consequence, and one thing this does not fix

The first port under this decision is `Percentile` / `Median` plus `FastMath.round`, which is
what four GATK annotations reach through `MathUtils.median`.

It does not fix `Math.pow`, `Math.exp` or the trigonometric functions: those call sites name
`java.lang.Math`, and routing them to a commons-math3 equivalent would be a different function
producing different bits. Where a call site names `Math`, the port must name `Math`, and the
divergence stays measured and declared. Decision 0013 still holds for exactly those.

## Attribution

Ported commons-math3 code carries its source in the module header, as every ported symbol does,
and the crate's licence file records the Apache 2.0 origin. The version is pinned to **3.5**:
GATK's `build.gradle` declares `strictly '3.5'`, so 3.5 is the version whose numbers reach a
golden downstream. The `jmath` corpus job downloads 3.6.1 for its `FastMath` columns, which is a
disagreement between two pins in this repository and is recorded in the manifest until one of
them moves.
