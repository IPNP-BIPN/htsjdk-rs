# 0011. Metrics number formatting depends on the JVM's locale, and nothing pins it

**Status:** accepted; locale hazard measured, formatting port deliberately incomplete
**Date:** 2026-07-21
**Follows:** [0005](0005-java-math-has-three-implementations.md)

## Finding 1: the locale, which is the important one

`FormatUtil` is what every Picard metrics file and every GATK report is written through. Its
constructor is four lines and each one is a decision:

```java
this.floatFormat = NumberFormat.getNumberInstance();
this.floatFormat.setGroupingUsed(false);
this.floatFormat.setMaximumFractionDigits(6);
this.floatFormat.setRoundingMode(RoundingMode.HALF_DOWN);
decimalFormatSymbols.setNaN("?");
decimalFormatSymbols.setInfinity("?");
```

`NumberFormat.getNumberInstance()` takes the **default locale**. A grep of htsjdk, Picard's
`cmdline` package and GATK's `Main` finds no `Locale.setDefault`, no `-Duser.language`, no
pinning of any kind.

Measured in the pinned oracle container, formatting the same doubles under five locales:

| locale | `1/3` | `2/3` | `1234567.891234567` |
|---|---|---|---|
| `en-US` | `0.333333` | `0.666667` | `1234567.891235` |
| `fr-FR` | `0,333333` | `0,666667` | `1234567,891235` |
| `de-DE` | `0,333333` | `0,666667` | `1234567,891235` |
| `hi-IN` | `0.333333` | `0.666667` | `1234567.891235` |
| `ar-EG` | `٠٫٣٣٣٣٣٣` | `٠٫٦٦٦٦٦٧` | `١٢٣٤٥٦٧٫٨٩١٢٣٥` |

So a Picard metrics file written on a machine with `LANG=fr_FR.UTF-8` has commas for decimal
points, and one written under `ar-EG` uses Eastern Arabic-Indic digits and inserts a
right-to-left mark before negative signs. Integers go through a separate formatter and are
localised too: `1234567` becomes `١٢٣٤٥٦٧`.

This is worth stating plainly. **This project began from PR #9384's finding that GATK's
floating-point results differ by roughly 1 ULP across CPU architectures.** That is a real
reproducibility problem, and it is a smaller one than this. A 1-ULP difference in the last
decimal place of a metric is invisible in practice; a metrics file whose numbers are written
in a different script is not even parseable by the tool that expects to read it back. And
unlike the ULP problem, this one needs no exotic hardware: it needs a laptop configured in
French.

### Consequence

The oracle contract pins `Locale.US`, and the claim this project makes about metrics output is
scoped to it. `tools/oracle/OracleProbe.java` should assert the locale the same way it asserts
`os.arch` and the GKL provider state, so a golden cannot be generated under a locale that would
silently change it. That is the same fail-closed principle as decision 0004.

## Finding 2: `HALF_DOWN` is not a decimal tie rule, and `Double.toString` is not shortest

Two further properties, both diagnosed against the JVM rather than assumed, and both blocking a
complete port.

**Java 17's `Double.toString` does not produce the shortest round-trip decimal.** It still uses
the pre-JDK19 algorithm (JDK-4511638, fixed in JDK 19). `DecimalFormat` rounds *those* digits,
via `DigitList.set` → `FloatingDecimal`, so the difference reaches the output:

| bits | Java 17 `toString` | Rust shortest |
|---|---|---|
| `0x438f67ea69ed3795` | `2.82879384806159008E17` | `2.82879384806159e17` |
| `0x45300c520a43f0af` | `1.9400994884341944E25` | `1.9400994884341945e25` |

The first emits three digits more than needed; the second picks a different final digit.

**`RoundingMode.HALF_DOWN` in `DecimalFormat` is not "ties toward zero".**
`DigitList.shouldRoundUp` consults whether the digit string represents the value *exactly*
(`valueExactAsDecimal`) and whether `FloatingDecimal` already rounded up to produce it
(`alreadyRounded`). A decimal that merely looks like a tie is resolved against the true binary
value. `0.9999995` prints as `1` because the true value is `0.99999950000000004…`, above the
tie; `0.1234565` prints as `0.123456` because its true value is below.

Measured over the corpus, restricted to apparent ties:

| `sign(true value − decimal)` | htsjdk rounds up |
|---|---|
| positive | 39 of 39 |
| negative | 7 of 46 |

Which is to say the direction is mostly explained by comparing against the true value, and not
entirely. The residue is `alreadyRounded`.

## Decision: stop at a measured partial port rather than fit a rule

The port implements the locale-pinned layout, the `?` symbols, the six-digit limit, and the
`valueExactAsDecimal` half of the tie rule. It agrees with htsjdk on **41,566 of 41,678 values,
99.73%**. The 112 that disagree are listed by bit pattern in
`crates/htsjdk-metrics/tests/data/known_divergences.tsv`, with the answer each side gives.

Completing it means porting Java 17's `FloatingDecimal` and `DigitList` — the exact decimal
expansion of a double — not adjusting a rounding condition.

Three intermediate rules were tried and measured, and the best-scoring one was **not** kept:

| rule | divergences of 41,678 |
|---|---:|
| plain decimal HALF_DOWN | 109 |
| plain HALF_UP | 115 |
| HALF_DOWN with `valueExactAsDecimal` (kept) | 112 |

The kept rule scores worse than the plain one. It is kept because it encodes a real clause of
the documented Java algorithm, while the plain rule happens to score better on this corpus by
accident. Choosing the lower number here would be fitting to the corpus, which passes the test
at hand and diverges on the next input — the precise failure mode this project exists to
prevent, and the same lesson as decision 0006, where chasing the last 0.0044% of `log` by
tweaking was the wrong move and finding the cause was the right one.

The failing values are pinned rather than deleted. A conformance suite trimmed to what already
passes reports 100% and means nothing.

## Bearing on the calibration gate

`FormatUtil` is shared by all 57 metrics collectors and by every GATK report, so it is the
first thing the metrics archetype needs and it was attempted first for that reason. It is
**not** representative of archetype cost: it is shared infrastructure paid once, and the
marginal cost of collector number two is still unmeasured.


## Addendum: cause A confirmed for `float`, not only `double`

Porting the SAM text writer produced one divergence in 67 lines, and it is the same cause.
`TextTagCodec` renders a float tag with `Float.toString`, and Java 17's is not the shortest
round-trip decimal there either:

| value | Java 17 | Rust shortest |
|---|---|---|
| `Float.MIN_VALUE` (`0x00000001`) | `1.4E-45` | `1E-45` |

Both parse back to the same subnormal, so both are valid round trips; they are simply different
decimals. The `FloatingDecimal` port this decision already lists as outstanding covers both
widths, so the SAM case is pinned in
`crates/htsjdk-bam/tests/sam_text_conformance.rs` rather than patched, on the same reasoning:
patching one rendering to match while the underlying algorithm differs would pass that test and
fail the next value.
