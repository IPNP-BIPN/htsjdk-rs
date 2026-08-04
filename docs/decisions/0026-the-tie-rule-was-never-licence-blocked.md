# 0026. The formatting tie rule was never licence-blocked

**Status:** accepted; 112 divergences reduced to 66, and the remainder is one cause
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
| after | **66** | **99.84%** |

46 fixed, 0 introduced. The values that moved are exactly the shape the old rule could not see:

| value | the double | old | htsjdk |
|---|---|---|---|
| `0.1234565` | `0.1234564999999999967…` | `0.123457` | `0.123456` |
| `0.1234575` | `0.1234574999999999977…` | `0.123458` | `0.123457` |

Both print a `5` at the seventh fraction digit with nothing after it. No rule reading only the
digit string separates them.

## What the remaining 66 are

One cause, and it is already on the record. **Every one is above 2^53**, where Java 17 stops
printing the shortest decimal form. That is the pre-Schubfach `FloatingDecimal` behaviour decision
0017 located and measured for `String.format`, reached here through `DecimalFormat` instead.

Below that line the port and the reference agree on every value in the corpus, and the suite
asserts that rather than stating it: `nothing_left_is_below_two_to_the_fifty_three` fails if any
divergence appears under it.

## Addendum, 2026-08-04: two of the 68 were ours

The first pass left 68 and attributed all of them to Java. Two were not Java's doing.

`6.985838094673373e14` is exactly `698583809467337.25`, so the two sixteen-digit forms `…337.2` and
`…337.3` are **equidistant** and both round-trip. "Shortest, then nearest, then even" is the
specification Ryu and Schubfach implement, and therefore what Java 19 and later do — and what
Java 17 does here too. Rust's formatter gets the length right and picks the odd form.

`shortest_decimal` now corrects it: among the decimals of the length Rust chose, take the one
nearest the double, ties to even. That is implementing to the specification, not to Java 17.

The cost had to be managed. At seventeen digits some twenty neighbouring decimals round-trip, so a
naive check ran the exact expansion on nearly every value and the suite went from 0.12s to 6s. A
tie needs the double to *be* a short decimal, which is testable without expanding it: written as
`odd * 2^power`, the exact decimal has at least `digits(odd) + floor(0.699 * -power)` significant
digits when `power` is negative and `digits(odd) + floor(0.301 * power)` when it is positive, and a
tie between forms of at most eighteen digits cannot happen once either exceeds nineteen. With that
gate the suite is back to 0.15s.

The same correction, and the same gate, are in gatk-rs's `DecimalFormat` port, where they took its
divergences from 4 to 2 — likewise leaving only values above 2^53.

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
