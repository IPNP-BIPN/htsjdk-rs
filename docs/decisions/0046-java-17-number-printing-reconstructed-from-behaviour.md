# 0046. Java 17's number printing, reconstructed from its behaviour

**Status:** accepted
**Date:** 2026-10-13
**Closes:** the licence block of [0013](0013-the-last-divergences-are-blocked-by-a-licence-not-by-difficulty.md), and the residues of [0011](0011-metrics-number-formatting-depends-on-the-jvm-locale.md) and [0017](0017-the-jvm-formats-a-short-decimal-not-the-value.md)
**Result:** `FormatUtil` 41,678 of 41,678, the SAM text writer 67 of 67, `%.2f` / `%.3f` / `%.3e` 127,803 of 127,803

## What was blocked, and what was not

Decision 0013 blocked the last formatting divergences on a licence. Java 17's `FloatingDecimal` is
GPL-2.0, and a translation of it cannot ship in an MIT crate. That conclusion stands: nothing here
translates it.

0013 also rejected "reconstruct the algorithm from the literature", because a paper describes an
algorithm and not what Java 17 does. That objection is answered by measuring rather than by
reading: every rule below was found as a cluster of disagreements between a candidate model and the
oracle, and each is closed only when the cluster is. The source and the bytecode of `java.base`
were not consulted. The model that did the fitting may have seen OpenJDK's source in training, so
this is not a clean room in the strict sense; what the rules rest on is the measurement, and the
code is shaped by Steele and White's and Gay's published descriptions, not by `FloatingDecimal`'s.

## Method

`Double.toString` and `Float.toString` were queried in the pinned oracle image
(Temurin 17.0.19, `linux/amd64`) on stratified samples:

| family | values |
|---|---|
| uniform random bit patterns | 2,000,000 |
| short decimals (1 to 17 digits, every exponent) and their ±3-ulp neighbours | 4,100,000 |
| every subnormal significand from 1 to 2,000,000 | 2,000,000 |
| every odd significand up to 8191 at every binary exponent | 8,548,351 |
| uniform significands at binary exponents 63 to 200 | 3,000,000 |
| dense sample at binary exponent 84 | 400,000 |
| random `float` bit patterns, and dense samples at `float` exponents 82 and 85 | 4,000,000 |

The starting model was Steele and White's free-format generation (PLDI 1990), which is what the
pre-Schubfach behaviour is documented to descend from. Each disagreement cluster was then isolated
by binary exponent, by digit count and by the generator's internal state, and fitted.

## What Java 17 does

1. **Bounds are exclusive.** A decimal exactly half an ulp from the double is rejected even where
   round-half-even reading would land back on it. This is `1e23` printing `9.999999999999999E22`,
   and it alone accounted for 1,020 of the first 1,091 disagreements above 2^63.
2. **A power of two gets the narrow interval on both sides.** Its lower gap is half its upper gap;
   Java 17 uses the smaller for both. 345 of the 2,046 normal powers of two print a digit more
   because of it, and every one of the 2,035 checked agrees. Subnormal powers of two get the same
   treatment although their gaps are equal.
3. **An exact tie at the last digit goes to the even digit.** 46 ties rounded 7 to 8 and 43 kept 2.
4. **Integers between 2^52 and 2^63 are printed from their exact value** (2^23 to 2^63 for
   `float`), rounded half-up after dropping `digits(2^(s - 2)) - 1` low digits, `s` being the shift
   of the significand. The same formula gives the double table (0, 1, 2 for shifts up to 5, 8, 10)
   and the float one (0 to 11). A tie cannot occur there.
5. **The high-side test overflows a signed 64-bit integer.** At `double` exponent 84 the last
   digit is kept where it should round up once the remainder passes 0.94737 of a unit. The
   predicted threshold, `(r + m) * 2 * 10^k >= 2^(63 + k)`, is 947,366,450 in units of the value;
   the measured one lies between 947,359,232 and 947,367,424. At `float` exponent 85 the prediction
   is 0.789007089 and the measurement lies between 0.7890069 and 0.7890083. `k` is Gay's estimate
   of the decimal exponent (1990), and it is observable: at `float` exponent 82 it is one too high
   from `9.5826e24` up, which makes the overflow certain there, and that boundary is where the
   estimate crosses 25. The overflow happens only up to exponent 84 for `double` and 85 for
   `float`; above that, applying the rule mispredicts.
6. **Six subnormals of each type** print a digit more, or a different digit, than the rules above
   give: `double` significands 1, 10, 12, 14, 16 and 18 (`4.9E-324`, then `4.9E-323` to
   `8.9E-323`), and `float` significands 2, 3, 6, 7, 21 and 29 (`2.8E-45` to `4.1E-44`). The
   estimate is exact for all of them, so rule 5 is not the cause. The `float` six were missed by
   random sampling and found by the CI sweep's enumeration of the smallest subnormals, which is
   why that enumeration is in it. For `double` no seventh exists among the subnormal significands
   up to 2,000,000 or the subnormals with up to 13 significant bits; for `float` there is none, by the exhaustive check below. They
   are a table.

## Verification

- The 25,039,046 sampled strings above, compared as strings: no difference.
- **Every `float`**: `Float.toString` of all 2^31 non-negative bit patterns (NaNs skipped), CRC32
  per block of 2^20, computed in the oracle and in Rust: all 2,048 blocks agree. A negative float
  differs only by its sign, so `float` is verified exhaustively.
- The repository's own goldens, which were produced independently of this work: `FormatUtil`'s
  41,678 values (the 66 quarantined ones now match), the SAM text writer's 67 records
  (`Float.MIN_VALUE` is `1.4E-45`), and the `%f` / `%e` sweep of 127,803 doubles (was 99.85%).

## Where it lives

`jmath::to_string`: `double_to_string`, `float_to_string`, and `double_digits` / `float_digits`
for the callers that round the digit string themselves (`htsjdk-metrics::format`,
`htsjdk-vcf::jformat`). The arithmetic needs integers up to about 1,200 bits, held in a small
stack-allocated big integer rather than a dependency; values whose scale fits in 2^122, which is
nearly every value a metrics file or a VCF carries, run the same loop in `u128`. It costs about
twice Rust's own `{:e}` on such values.

`tools/jmath-conformance/ToStringSweep.java` regenerates 3.1 million lines over the families above
in CI, and `crates/jmath/tests/to_string_sweep.rs` requires every one to match.

## What this changes elsewhere

- 0013's quarantine is empty. "Licence-blocked" no longer describes any formatting output.
- The divergence lists in `format_conformance.rs` and `sam_text_conformance.rs` are gone; a
  divergence there is now a regression.
- gatk-rs's `UnrenderableDouble` was the same gap and can close on the next pin.
