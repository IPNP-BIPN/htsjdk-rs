# 0026. The formatting tie rule was never licence-blocked

**Status:** accepted; 112 divergences reduced to 68, and the remainder is one cause
**Date:** 2026-08-04
**Follows:** [0011](0011-metrics-number-formatting-depends-on-the-jvm-locale.md),
[0013](0013-the-last-divergences-are-blocked-by-a-licence-not-by-difficulty.md),
[0017](0017-the-jvm-formats-a-short-decimal-not-the-value.md)

## What was believed

Decision 0013 quarantined 112 of `FormatUtil`'s 41,678 sampled values as **licence-blocked**, and
the reasoning was explicit:

> `DigitList.shouldRoundUp` consults whether that string is exact and whether it was already
> rounded up. Both require the exact decimal expansion of the double.

The conclusion drawn was that both required `FloatingDecimal`, which is GPL2 and may be linked but
not translated. The port's tie rule was left with a comment saying so, and an approximation in
place of the missing facts:

```rust
// INCOMPLETE. Java also consults `alreadyRounded` […] Reproducing it needs the exact
// decimal expansion of the double, which is most of `FloatingDecimal` itself.
let round_up = first_dropped > b'5' || (first_dropped == b'5' && (past_half || !exact));
```

Read that last line as a decision procedure: when the digit string is not the double's exact
value, round up. The true value is below the halfway point about as often as it is above, so this
rule is wrong on roughly half the cases it decides.

## The mistake in the reasoning

The premise is true and the inference does not follow. Both facts do need the exact decimal
expansion of the double. Getting that expansion does not need `FloatingDecimal`.

Every finite double is a finite decimal, exactly:

```text
m * 2^e   with e >= 0  is  m doubled e times
m * 2^e   with e <  0  is  m * 5^-e, with the point moved -e places left
```

because `2^e = 5^-e / 10^-e`. Multiplying a decimal digit string by five is one pass over its
digits. That is thirty lines of schoolbook arithmetic on the bits of an IEEE 754 double, and it is
in no sense a translation of anything in `java.base`. What `FloatingDecimal` is needed for is the
*shortest* representation, which is a genuinely hard problem and a genuinely different one — and
Rust's `{:e}` already supplies it.

So the port never needed the blocked code. It needed a comparison nobody had written:

```rust
fn compare_to_exact(digits: &str, exp10: i32, magnitude: f64) -> Ordering
```

`Equal` is `valueExactAsDecimal`. `Less` is `alreadyRounded` — the shortest form was rounded up to
reach it, so the true value is below the halfway point and the apparent tie is not one.

## What it was worth, measured

The corpus is unchanged: 41,678 values from the pinned oracle, `tools/metrics-conformance/FormatDump.java`.

| | divergences | agreement |
|---|---|---|
| before | 112 | 99.73% |
| after | **68** | **99.84%** |

44 fixed, 0 introduced. The values that moved are exactly the shape the old rule could not see:

| value | the double | old | htsjdk |
|---|---|---|---|
| `0.1234565` | `0.1234564999999999967…` | `0.123457` | `0.123456` |
| `0.1234575` | `0.1234574999999999977…` | `0.123458` | `0.123457` |

Both print a `5` at the seventh fraction digit with nothing after it. No rule reading only the
digit string separates them.

## What the remaining 68 are

One cause, and it is already on the record. All 68 are Java 17 emitting digits that the shortest
form does not have: **66 are above 2^53**, and **2 need sixteen significant digits**. That is the
pre-Schubfach `FloatingDecimal` behaviour decision 0017 located and measured for `String.format`,
reached here through `DecimalFormat` instead.

Not one remaining divergence is a rounding decision, and the suite asserts that rather than
stating it: `nothing_left_is_a_rounding_decision` fails if a divergence appears below 2^53 with
fewer than sixteen significant digits.

## One place the old rule was right, and why

The branch where the whole value falls below the last place the pattern can show is left exactly
as it was. There a lone `5` rounds **up** unless the digit string is the exact value — the
opposite of the rule one digit to the right. `5e-7` is `4.99999999999999977…e-7`, below the
halfway point, and htsjdk still prints `0.000001`.

Changing that branch to match the other one was tried, and it broke that value. The asymmetry is
real, it is measured, and it is not understood. It is recorded here rather than tidied away.

## The general lesson

Decision 0013's closing section is about a risk the plan had not registered: the reference
implementation's licence being incompatible with the port. This record adds the failure mode that
comes with it. **A licence blocks the code, not the fact.** Once a piece of work is labelled
licence-blocked it stops being examined, and the label covers whatever was nearby when it was
applied. Here it covered a comparison that needed no privileged access at all.

Decision 0017 had already drawn the same distinction — "the licence blocks transcription; it does
not block observing behaviour, stating it, and implementing it" — and drew it for `Formatter`
while leaving `DecimalFormat` where it was, on the explicit grounds that it was a different
mechanism. It is a different mechanism. It was not a different licence question.

The remaining item worth re-examining in the same spirit is `Math.exp` (0014, and gatk-rs #71),
where the argument for out-of-scope is not the licence but the absence of a specification to
implement against. That argument is unaffected by this one, and this record does not weaken it.

## Verification

* `cargo test -p htsjdk-metrics` — 68 declared divergences, each still exactly as recorded.
* `nothing_left_is_a_rounding_decision` — the claim of this record, as an assertion.
* Sabotage: restoring `!exact` in place of the comparison reintroduces 44 divergences and the
  suite names all of them.
