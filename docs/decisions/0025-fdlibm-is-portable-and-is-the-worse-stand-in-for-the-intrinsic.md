# 0025. FDLIBM is portable, and it is the *worse* stand-in for `Math.exp`

**Status:** accepted; `StrictMath.exp` ported, `Math.exp` still unported
**Date:** 2026-08-04
**Amends:** [0014](0014-math-exp-was-withdrawn-it-was-a-gpl2-transcription.md)
**Follows:** [0005](0005-java-math-has-three-implementations.md), [0013](0013-the-last-divergences-are-blocked-by-a-licence-not-by-difficulty.md)

## The question

Decision 0014 removed `jmath::math::exp` because it was a line-by-line transcription of HotSpot's
x86 intrinsic, whose source is GPL2 **only**, and recorded the cost: the corpus `exp` points were
routed to the system libm, which agrees with `java.lang.Math` on 99.9711% of them.

That left an obvious question unanswered. HotSpot's intrinsic is not the only implementation of
`exp` in the world, and it is not even the only one Java ships. `java.lang.StrictMath` is
*specified* to be FDLIBM, and FDLIBM is Sun's freely-distributable library whose notice grants the
right to use, copy, modify and distribute provided the notice is preserved. So: port FDLIBM, and
see whether it closes the gap.

It does not. It is worse.

## What was measured

`crates/jmath/src/strict_exp.rs` is FDLIBM's `e_exp.c`, ported with the notice preserved. Over the
44,996-point `exp` corpus:

| implementation | agrees with `java.lang.Math.exp` | agrees with `java.lang.StrictMath.exp` |
|---|---:|---:|
| system libm | 99.9711% | not measured |
| **FDLIBM (this port)** | **98.6443%** | **100%** |

No FDLIBM divergence from `Math.exp` exceeds **1 ulp**.

Two results, and the second is the one that would have been assumed wrong.

**FDLIBM is exactly `StrictMath.exp`.** Not approximately: on every point. The specification says
so, and the corpus now checks it rather than citing it. Any call site that reaches `StrictMath.exp`
is portable from today, bit for bit, under a permissive licence.

**FDLIBM is a worse approximation of `Math.exp` than the host's libm.** By 1.3 percentage points,
which is roughly 45 times as many divergent points. This is not paradoxical once stated: modern
libm implementations and HotSpot's intrinsic are both descended from newer, more accurate work than
FDLIBM's 1993 algorithm, so they land on the correctly-rounded result more often and therefore on
each other more often. FDLIBM is faithful to a specification that `Math` does not follow.

## What this settles

**The licence is not the only reason `Math.exp` is unported, and it never was the interesting one.**
0014 framed the gap as a cost imposed by copyright. That framing is right about why the
transcription had to go, and it quietly implied that a permissively-licensed algorithm would close
the gap. This measurement says it would not. `Math.exp` is unportable because it is *that specific
intrinsic*, not because it is an exponential.

**A permissive re-implementation is the right move anyway, for a different reason.** It bounds the
error. Before this, "the port gets `exp` wrong somewhere" had no size attached to it: the system
libm's divergence was measured in a *rate*, not a magnitude, and a rate says nothing about whether
a wrong answer is one ulp out or catastrophically out. Now it is bounded at 1 ulp by an
implementation we control and can reason about, on a host-independent basis. The system libm's
magnitude remains a property of whatever machine the port runs on, which is not a property at all.

**Call sites can now be classified rather than lumped together.** A GATK path through
`StrictMath.exp` is portable. A path through `Math.exp` is not, and the honest claim for it is
bio-identical with a named quarantine bounded at 1 ulp, not bit-identical. `AllelePseudoDepth`
(gatk-rs G1.7) is in the second group: `SomaticLikelihoodsEngine.alleleFractionsPosterior` reaches
`Math.exp` through `NaturalLogUtils`, so it stays unported for byte-identity, with the difference
that the gap now has a measured size instead of an open end.

## The rule this generalises to

For every function blocked by a licence rather than by difficulty:

1. **Do not transcribe the blocked source.** Translation into another language is a derivative
   work; changing the language changes nothing about the copyright. This is what 0014 learned and
   it is not up for revisiting.
2. **Port the best permissively-licensed implementation instead**, cite it, and preserve whatever
   notice it requires.
3. **Measure it against the reference on the corpus**, both as an agreement rate and as a worst-case
   magnitude. Report both.
4. **Record which specification it is exact for**, if any. FDLIBM is exact for `StrictMath` and
   inexact for `Math`; that distinction is the whole value of step 2, and it is invisible without
   step 3.

The outcome of step 3 can be that the permissive implementation is worse, as it is here. That is a
result, not a failure: it converts "blocked by a licence" into "blocked by this intrinsic
specifically, and by at most one ulp".
