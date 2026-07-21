# 0013. The last formatting divergences are blocked by a licence, not by difficulty

**Status:** accepted; the port is blocked, and the block is not technical
**Date:** 2026-07-21
**Follows:** [0011](0011-metrics-number-formatting-depends-on-the-jvm-locale.md)

## The situation

Two decisions name the same outstanding work.

- **0011**: `FormatUtil` agrees with htsjdk on 41,566 of 41,678 values (99.73%). The 112 that
  disagree have two diagnosed causes, and both require Java 17's `FloatingDecimal` — the exact
  decimal expansion of a double — and `DigitList.shouldRoundUp`.
- **0011, addendum**: the SAM text writer's one divergence in 67, `Float.MIN_VALUE` rendering
  as `1.4E-45` rather than `1E-45`, is the same cause for `float`.

Porting `FloatingDecimal` closes both at once, which made it the most valuable remaining item.
It was scoped, and it stopped before a line was written.

## Why it stopped

`FloatingDecimal` and `DigitList` are `java.base` classes. Checked in the pinned oracle image
rather than recalled:

```
$JAVA_HOME/legal/java.base/LICENSE      -> "The GNU General Public License (GPL), Version 2"
$JAVA_HOME/legal/java.base/ASSEMBLY_EXCEPTION
```

The exception is the OpenJDK Assembly Exception, and its scope is explicit:

> Linking this OpenJDK Code statically or dynamically with other code is making a combined work
> based on this library. Thus the terms and conditions of GPL2 cover the whole combination.
>
> As a special exception, Oracle gives you permission to **link** this OpenJDK Code with certain
> code licensed by Oracle […]

It grants permission to **link**, to a named set of Oracle modules. It does not grant permission
to translate the source into another language and redistribute the translation under a different
licence. A line-by-line port of `FloatingDecimal` is a derivative work of GPL2 code, and
htsjdk-rs is MIT because htsjdk is MIT.

This is the first place in the program where the reference implementation is **not** under a
licence compatible with the port. htsjdk is MIT, Picard is MIT, GATK is Apache 2.0. The JVM
underneath all three is not.

## Why the obvious ways round are all worse

**Port it anyway and relicense the crate GPL2.** Correct legally and wrong for the project: it
would make the number formatter, and anything linking it, GPL2. That is every Picard metrics
tool and every GATK report, which is most of the program.

**Reconstruct the algorithm from the literature.** The pre-JDK19 behaviour derives from Steele
and White's dragon4. But this violates the project's founding rule, stated in the plan and
followed in every decision so far: *do not reconstruct behaviour from documentation, papers, or
memory*. The whole reason that rule exists is findings like decision 0009, where the correct
behaviour is written down nowhere and is an emergent property of a data-structure choice. A
paper describes an algorithm; it does not describe what Java 17 does, and Java 17's known
deviation from shortest-repr is precisely the thing that has to be reproduced.

**Fit an implementation to the corpus.** Decision 0011 already rejected this explicitly, and
rejected a rule that scored *better* on the corpus for the same reason. Reversing that here
would discard the principle at the exact moment it costs something.

**Call the JVM at run time.** Defeats the purpose.

## Decision

**The 112 double divergences and the 1 float divergence stay open, and they are reclassified.**

They are not "not yet ported". They are **licence-blocked**, which is a different status with a
different resolution path, and the difference matters to anyone reading the coverage numbers.
The decision records and the pinned divergence lists are updated to say so.

The affected values are quarantined in the vocabulary PR #9384 established: any output that
passes through one of them is **bio-identical**, not bit-identical, and the list of exactly
which values is committed rather than described.

Practically the exposure is small and measurable: 112 of 41,678 sampled doubles, 0.27%, all at
tie boundaries or in the region where Java 17 emits non-shortest digits. But it is not zero, and
a metrics file containing one of those values will differ in one digit.

## What could actually resolve it

1. **Ask.** The values are a tiny, self-contained function. A request to Oracle or to the
   OpenJDK project for permission to redistribute a translation under MIT is cheap to make and
   might succeed.
2. **JDK 19+.** From JDK 19 `Double.toString` is Raffaello Giulietti's Schubfach and *is*
   shortest-round-trip, which Rust already produces. Cause A disappears entirely against a
   JDK 19+ oracle. That does not help while the contract pins JDK 17, but it means the problem
   has a shelf life, and it reframes the choice of oracle JDK as a licensing decision as well as
   a fidelity one.
3. **An independent implementation with its own provenance.** If a permissively licensed
   implementation of the pre-JDK19 algorithm exists, using it is clean. This has not been
   searched.

## The general lesson, which is bigger than this function

The plan's risk register lists ten risks and none of them is "the reference implementation's
licence is incompatible with the port". That gap is worth closing now rather than at the next
occurrence, because the next occurrence is predictable: **every place the port needs JVM
behaviour rather than htsjdk behaviour hits the same wall.** `String.hashCode` ordering,
`Arrays.sort` tie-breaking, `Math` intrinsics, `DecimalFormat` — the ported code sits on a GPL2
runtime, and every time its exact behaviour leaks into an output byte, the same question
returns.

Decision 0005 already found that `Math` has three incompatible implementations and that the port
must track which one each call site used. It did not notice that two of the three are GPL2. The
`jmath` crate was built against an empirically measured corpus and reproduces `Math.log` by
being *correctly rounded* rather than by transcription (decision 0006), so it is clean. That was
luck rather than design, and this record makes it design.
