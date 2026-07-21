# 0007. `Math.pow` is deferred: its intrinsic depends on an approximate hardware instruction

**Status:** accepted; hazard partially cleared, see the addendum
**Date:** 2026-07-21
**Follows:** [0006](0006-correct-rounding-is-the-target-for-log-and-log10.md)

## Finding

HotSpot's `pow` intrinsic uses `rcpps` / `rcpss`, the packed approximate reciprocal:

| intrinsic | uses `rcpps`/`rcpss` | correctly rounded | status |
|---|---|---|---|
| `exp` | no, 0 sites | no | **ported, exact** |
| `log` | yes, 5 sites | **yes** | **ported, exact** |
| `pow` | yes, 6 sites | no | **deferred, this decision** |

`rcpps` is specified by an *error bound*, not by exact results: the Intel SDM guarantees a
relative error below `1.5 * 2^-12` and states the values are implementation-dependent. Nothing
architecturally fixes the bits.

For `log` this is harmless. The instruction only seeds a refinement, and the final result was
measured to be correctly rounded, so whatever `rcpps` returned is washed out. That is why a
correctly-rounded double-double implementation matched it on every point without reproducing
the intrinsic at all.

For `pow` there is no such absorption: the result is measurably not correctly rounded
(99.9378% agreement with the host libm, and MPFR shows `Math` is the correctly-rounded answer
in only 202 of 252 divergent cases). So the approximation can reach the output.

## Why that is a problem for this project, not just for `pow`

If `Math.pow`'s exact bits depend on the CPU's `rcpps` table, then **`pow` is not reproducible
across x86-64 implementations**, and a golden containing a `pow`-derived value is only valid
for the silicon that produced it. That would undercut the oracle contract for every tool whose
output passes through `pow`.

There is a second, sharper problem specific to this setup: the corpus was generated in an
**emulated** amd64 container. Rosetta 2 must emulate `rcpps` somehow, and its choice of values
is its own. So the committed `pow` column may already describe Rosetta rather than any real
processor.

This does **not** retroactively threaten `exp` (no `rcpps` at all) or `log`, `log10`, `sqrt`
(correctly rounded, therefore hardware-independent by construction). The four ported functions
are safe. It is scoped precisely to `pow`, and to any future intrinsic port that turns out to
depend on an approximate instruction without refining it away.

## Decision

**Defer `pow` until the hazard is measured.** Porting 2,220 lines of assembly against a target
that may not be stable is the wrong order of work.

The test is cheap and decisive: regenerate the `pow` column of the corpus on real x86-64
hardware, ideally both Intel and AMD, and diff it against the emulated one. Three outcomes:

1. **Identical everywhere.** `rcpps` is de facto uniform; port the intrinsic normally.
2. **Emulated differs from real silicon.** The oracle must move to real hardware for `pow`,
   and the current corpus column is discarded.
3. **Intel differs from AMD.** `Math.pow` is genuinely not portable across x86-64, which is a
   finding worth reporting upstream, and `pow`-derived outputs must be quarantined as
   bio-identical rather than bit-identical.

## Note on decision 0004

Decision 0004 corrected the plan's claim that emulation forces goldens onto real x86-64 CI,
because the stated reason (no AVX, GKL degrades) was measured and found false. This is a
*different* reason to want real hardware in the loop, and it is specific rather than general:
not "emulation is untrustworthy", but "one instruction in one intrinsic is specified loosely,
so verify that one thing on real silicon".

That distinction matters. The first version of the argument would have slowed every golden for
a reason that did not hold; this one gates exactly the part that is actually at risk.

## Priority

`pow` is not on the current tool ladder. `RecalDatum` uses `Math.log10`, and
`DuplicationMetrics.estimateLibrarySize` uses `Math.exp`; both are ported and exact.
`MathUtils.log10SumLog10` reaches `Math.pow`, so this becomes blocking at BQSR and beyond, not
before.


## Addendum: measured, and the hazard did not materialise on AMD

CI ran the experiment. The whole corpus was regenerated inside the same pinned oracle image on
a real **AMD EPYC 7763** (`AuthenticAMD`, GitHub Actions `ubuntu-latest`) and diffed against
the emulated corpus produced under Rosetta on Apple Silicon.

| fn | compared | differing |
|---|---:|---:|
| `exp` | 44,987 | 0 |
| `log` | 44,987 | 0 |
| `log10` | 44,987 | 0 |
| `sqrt` | 44,987 | 0 |
| **`pow`** | **404,883** | **0** |

Zero differences anywhere, including every `pow` point. So Rosetta's emulation of `rcpps` and
AMD's hardware agree, at least across the sampled inputs.

Two caveats keep this from closing the question outright:

1. **Intel is still untested.** The experiment compared emulated-on-ARM against AMD. The
   original worry was vendor divergence, and only one vendor has been sampled. An Intel runner
   would complete it.
2. **Agreement on a sample is not a proof of agreement everywhere.** 404,883 points is a good
   sample of a space of `2^128` argument pairs.

What this does establish is that the risk is smaller than feared and that `pow` is no longer
blocked on the *portability* question. It still needs its 2,220-line intrinsic ported, which is
a large job, and it still first becomes relevant at BQSR.

The generated goldens were also confirmed on real silicon: the BGZF corpus regenerated in the
container on the AMD host matches the committed goldens exactly, which independently validates
every claim in decisions 0001 and 0003 outside the emulated environment.
