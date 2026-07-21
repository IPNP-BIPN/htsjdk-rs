# 0017. The JVM's `%f` rounds a short decimal, not the value

**Status:** accepted; modelled independently, residue measured at 99.8537%
**Date:** 2026-07-21
**Follows:** [0011](0011-metrics-number-formatting-depends-on-the-jvm-locale.md),
[0013](0013-the-last-divergences-are-blocked-by-a-licence-not-by-difficulty.md),
[0014](0014-math-exp-was-withdrawn-it-was-a-gpl2-transcription.md)

## The problem

`VCFEncoder` puts every double it writes through `String.format(Locale.US, ...)`:

```java
private static String formatQualValue(final double qual) {
    String s = String.format(Locale.US, QUAL_FORMAT_STRING, qual);   // "%.2f"
```

Decision 0013 blocked exactly this path. `String.format` resolves into
`java.util.Formatter`, `FormattedFloatingDecimal` and `FloatingDecimal`, all in `java.base`,
which decision 0014 established is GPL2 with the Classpath Exception, and the exception permits
**linking**, not translation. Those symbols cannot be transcribed into an MIT crate.

The obvious fallback is Rust's own `format!("{:.2}", d)`. It is wrong, and it is wrong in a way
that produces valid numbers:

| value | JVM | Rust `{:.2}` |
|---|---|---|
| `0.125` | `0.13` | `0.12` |
| `2.675` | `2.68` | `2.67` |
| `1e300` | `1` and 300 zeros | the 301 digits of its exact value |

## The finding

The two disagree because they round **different things**. Rust's formatter rounds the exact
binary value, half-to-even. The JVM's formatter never sees the exact value: it takes the digit
string `FloatingDecimal` produces, which is a short decimal that round-trips, and rounds *that*
half-up.

That is why `2.675` behaves as if it were a tie when it is really `2.67499999999999982…`, and
why `1e300` prints as padding rather than as its exact expansion of
`1000000000000000052504760255204420248704…`. The JVM is rounding `1e300`'s three-character
short form.

## Decision

The contract is stated and implemented independently, in `crates/htsjdk-vcf/src/jformat.rs`,
which says at the top that it is **not a port**:

> Format the shortest decimal representation that round-trips to the double, then round that
> digit string half-up at the requested precision, padding with zeros when the request outruns
> the digits.

Nothing is transcribed. Rust's `{:e}` supplies a shortest round-tripping representation, and the
rounding and rendering are ordinary arithmetic on a digit string. The provenance audit is
satisfied because no `Ported from` claim is made against `java.base`.

## What it is worth, measured

`tools/vcf-conformance/JavaFormatSweep.java` dumps `%.2f`, `%.3f` and `%.3e` for 127,803 doubles
from the pinned oracle: exact ties at two and three decimals, the phred and allele-frequency
ranges VCF actually carries, uniform random bit patterns, and every power of ten from 1e-300 to
1e300 with both neighbours.

| format | agreement |
|---|---|
| `%.2f` | 99.8537% |
| `%.3f` | 99.8537% |
| `%.3e` | **100.0000%** |

All 187 divergences are the same root cause, and it is one already on the record. Java 17's
`FloatingDecimal` predates the Schubfach rewrite that landed in JDK 19, so its "shortest"
representation sometimes is not:

```
Double.toString(-6.451108167902569e16)  =  -6.4511081679025688E16    (17 digits)
Rust's shortest                          =  -6.451108167902569e16     (16 digits)
```

That extra digit changes `%f` only when it falls to the left of the requested decimal place. So:

* `%.3e` is unaffected at **every** magnitude, because four significant digits never reach it;
* the smallest diverging magnitude in the whole sample is **6.9e14**;
* below 1e15 the agreement is **99.996947%**.

VCF carries phred qualities and allele frequencies. Nothing it encodes comes within eleven
orders of magnitude of the divergence, and the VCF record conformance suite is byte-equal on
all 29 `formatVCFDouble` cases and all 54 record cases.

## Why this is better than the quarantine it replaces

Decision 0013 quarantined this class of formatting outright, on the grounds that the licence
blocked it. That was right about the licence and too pessimistic about the consequence. The
licence blocks **transcription**; it does not block observing behaviour, stating it, and
implementing it. What survives is not "double formatting is unavailable" but a bounded,
measured, and located residue: the doubles for which Java 17's own shortest representation is
not shortest, above 1e15, in fixed-point formats only.

The 112 `FormatUtil` divergences recorded in 0013 are a different mechanism (`DecimalFormat`
and `HALF_DOWN`, not `Formatter`) and remain where they were.

## Verification

* `cargo test -p htsjdk-vcf` covers the model directly and through the encoder.
* CI regenerates the sweep in the pinned oracle and fails if the agreement rate falls.
* Sabotage: replacing `format_fixed` with Rust's native `{:.*}` fails 3 sweep cases and 2
  record cases, so the suite detects precisely the mistake this module exists to avoid.
