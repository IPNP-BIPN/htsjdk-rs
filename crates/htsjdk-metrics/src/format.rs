//! Number formatting for metrics files.
//!
//! Ported from `htsjdk.samtools.util.FormatUtil`, which is what every Picard metrics file and
//! every GATK report is written through.
//!
//! `FormatUtil`'s constructor is four lines long and each of them is a decision:
//!
//! ```java
//! this.floatFormat = NumberFormat.getNumberInstance();
//! this.floatFormat.setGroupingUsed(false);
//! this.floatFormat.setMaximumFractionDigits(6);
//! this.floatFormat.setRoundingMode(RoundingMode.HALF_DOWN);
//! decimalFormatSymbols.setNaN("?");
//! decimalFormatSymbols.setInfinity("?");
//! ```
//!
//! - `getNumberInstance()` takes the **default locale**. Nothing in htsjdk, Picard or GATK
//!   pins it, so the bytes of a metrics file depend on the JVM's locale. See decision 0011:
//!   measured, and it changes both the decimal separator and, for some locales, the digits
//!   themselves.
//! - `HALF_DOWN` is not `DecimalFormat`'s default, which is `HALF_EVEN`. The two disagree on
//!   exact ties at the sixth fraction digit.
//! - NaN and infinity share the symbol `?`, so a metrics file cannot distinguish `NaN` from
//!   `+Infinity`. `-Infinity` becomes `-?`, because the sign is applied separately.

use std::cmp::Ordering;

/// `FormatUtil.DECIMAL_DIGITS_TO_PRINT`.
pub const DECIMAL_DIGITS_TO_PRINT: usize = 6;

/// `FormatUtil.format(long)`. Grouping is disabled, so this is the plain decimal form.
pub fn format_long(value: i64) -> String {
    value.to_string()
}

/// `FormatUtil.format(boolean)`.
pub fn format_bool(value: bool) -> &'static str {
    if value {
        "Y"
    } else {
        "N"
    }
}

/// `FormatUtil.format(double)`, under the pinned `en-US` locale.
///
/// Never uses scientific notation: `Double.MAX_VALUE` comes out as its full 309-digit
/// expansion, which is what htsjdk does and therefore what this does.
pub fn format_double(value: f64) -> String {
    if value.is_nan() {
        // `setNaN("?")`. The sign of a NaN is not printed.
        return "?".to_string();
    }
    let negative = value.is_sign_negative();
    if value.is_infinite() {
        // `setInfinity("?")`, with the sign applied by the format's negative pattern.
        return if negative { "-?" } else { "?" }.to_string();
    }

    let (digits, exp10) = shortest_decimal(value.abs());
    let body = round_to_fraction_digits(&digits, exp10, DECIMAL_DIGITS_TO_PRINT, value.abs());
    // Negative zero keeps its sign: htsjdk prints `-0`.
    if negative {
        format!("-{body}")
    } else {
        body
    }
}

/// The shortest decimal that round-trips, as `(digits, exp10)` where the value is
/// `0.<digits> * 10^exp10`.
///
/// `DecimalFormat` rounds this shortest representation rather than the exact binary value,
/// because `DigitList.set` goes through `FloatingDecimal`, the same code as `Double.toString`.
/// That distinction is observable: `0.1` is exactly
/// `0.1000000000000000055511151231257827…`, and formatting it to six digits gives `0.1` rather
/// than a rounding of the true value.
fn shortest_decimal(value: f64) -> (String, i32) {
    if value == 0.0 {
        return ("0".to_string(), 1);
    }
    // Rust's `{:e}` gives the shortest round-trip form, `d.ddddde±xx`.
    let s = format!("{value:e}");
    let (mantissa, exponent) = s.split_once('e').expect("scientific form");
    let exponent: i32 = exponent.parse().expect("exponent");
    let digits: String = mantissa.chars().filter(|c| c.is_ascii_digit()).collect();
    let digits = digits.trim_end_matches('0');
    let digits = if digits.is_empty() { "0" } else { digits };
    // `{:e}` normalises to one digit before the point, so exp10 (digits before the point when
    // written as 0.<digits>) is exponent + 1.
    (digits.to_string(), exponent + 1)
}

/// Where the double's exact value sits relative to its own shortest decimal form.
///
/// This is the pair of facts `DigitList.shouldRoundUp` consults and that this port used to be
/// missing: `Equal` is `valueExactAsDecimal`, and `Less` is `alreadyRounded` — the shortest form
/// was rounded **up** to reach it, so an apparent tie is really below the halfway point.
///
/// The previous version approximated both with "is the decimal exact", and rounded up whenever it
/// was not. That is right half the time by construction, which is what most of the quarantined
/// divergences of decision 0011 were.
///
/// Every double is a finite decimal, so this is an exact comparison and not an estimate. Nothing
/// here is a transcription of `FloatingDecimal`: `m * 2^e` with `e` negative is `m * 5^-e / 10^-e`,
/// and multiplying a decimal by five is one pass over its digits. That is arithmetic on the
/// double's own bits, which is why decision 0013's licence blocker does not reach it.
fn compare_to_exact(digits: &str, exp10: i32, magnitude: f64) -> Ordering {
    let (exact_digits, exact_exp10) = exact_decimal(magnitude);
    // Both are `0.<digits> * 10^exp` with no leading zero, so a differing exponent settles it. It
    // does differ for a value like 1e23, whose shortest form is `1` at exponent 24 while the double
    // is 0.99999999999999991611392 at exponent 23.
    if exp10 != exact_exp10 {
        return exact_exp10.cmp(&exp10);
    }
    let shortest = digits.as_bytes();
    for index in 0..shortest.len().max(exact_digits.len()) {
        let ours = exact_digits.get(index).copied().unwrap_or(0);
        let theirs = shortest.get(index).map_or(0, |b| b - b'0');
        match ours.cmp(&theirs) {
            Ordering::Equal => {}
            other => return other,
        }
    }
    Ordering::Equal
}

/// The double's exact decimal expansion, as `(digits, exp10)` with the value `0.<digits> * 10^exp10`.
fn exact_decimal(value: f64) -> (Vec<u8>, i32) {
    let bits = value.to_bits();
    let biased = ((bits >> 52) & 0x7ff) as i32;
    let fraction = bits & ((1u64 << 52) - 1);
    // Subnormals carry no implicit leading one and have a fixed exponent.
    let (mantissa, exponent) = if biased == 0 {
        (fraction, -1074)
    } else {
        (fraction | (1u64 << 52), biased - 1075)
    };

    // Least significant digit first, which is the direction a carry travels.
    let mut digits: Vec<u8> = Vec::with_capacity(800);
    let mut rest = mantissa;
    if rest == 0 {
        digits.push(0);
    }
    while rest > 0 {
        digits.push((rest % 10) as u8);
        rest /= 10;
    }

    // `m * 2^e` for a non-negative exponent is repeated doubling; for a negative one it is
    // `m * 5^-e` with the point moved -e places left, because `2^e = 5^-e / 10^-e`.
    let (factor, shift) = if exponent >= 0 {
        (2u8, 0)
    } else {
        (5u8, -exponent)
    };
    for _ in 0..exponent.abs() {
        let mut carry = 0u8;
        for digit in digits.iter_mut() {
            let product = *digit * factor + carry;
            *digit = product % 10;
            carry = product / 10;
        }
        while carry > 0 {
            digits.push(carry % 10);
            carry /= 10;
        }
    }

    // Back to most significant first, and to the `0.<digits>` convention.
    digits.reverse();
    let exp10 = digits.len() as i32 - shift;
    (digits, exp10)
}

/// Rounds `0.<digits> * 10^exp10` to at most `max_fraction` fraction digits, HALF_DOWN, and
/// renders it without grouping, without trailing zeros, and with at least one integer digit.
fn round_to_fraction_digits(digits: &str, exp10: i32, max_fraction: usize, value: f64) -> String {
    // Position, counted in digits after the decimal point, of the last digit we keep.
    // A digit at index i in `digits` sits at fraction position i + 1 - exp10.
    let keep = max_fraction as i32 + exp10;

    let mut kept: Vec<u8>;
    let mut exp10 = exp10;

    if keep <= 0 {
        // Everything is beyond the last kept place, so the answer is either zero or one unit in
        // that place. This branch is left as it was, because measurement says it was already right
        // and the reason is worth writing down: here a lone '5' rounds **up** unless the digit
        // string is the exact value, which is the opposite of what happens one place to the right.
        // `5e-7` is `4.99999999999999977…e-7`, below the halfway point, and htsjdk still prints
        // `0.000001`. Only a genuine tie goes down, which is `HALF_DOWN` doing what it says.
        let first = digits.as_bytes()[0];
        let lone_five = first == b'5' && digits.len() == 1;
        let round_up = keep == 0
            && first >= b'5'
            && !(lone_five && compare_to_exact(digits, exp10, value) == Ordering::Equal);
        if !round_up {
            return "0".to_string();
        }
        kept = vec![1];
        exp10 += 1;
    } else {
        let keep = keep as usize;
        if keep >= digits.len() {
            kept = digits.bytes().map(|b| b - b'0').collect();
        } else {
            kept = digits.bytes().take(keep).map(|b| b - b'0').collect();
            let rest = &digits[keep..];
            let first_dropped = rest.as_bytes()[0];
            let past_half = rest.bytes().skip(1).any(|b| b != b'0');
            // `DigitList.shouldRoundUp`, HALF_DOWN branch. A leading '5' with nothing after it
            // looks like a tie, but it only *is* one when the digit string is the double's exact
            // value. Otherwise the true value sits to one side of the halfway point and decides,
            // and which side it is cannot be read off the digits: `0.1234565` and `0.1234575`
            // print the same shape and go opposite ways.
            //
            // A real tie goes down, because `FormatUtil` sets `HALF_DOWN` rather than the
            // `DecimalFormat` default.
            let round_up = first_dropped > b'5'
                || (first_dropped == b'5'
                    && (past_half || compare_to_exact(digits, exp10, value) == Ordering::Greater));
            if round_up {
                let mut i = kept.len();
                loop {
                    if i == 0 {
                        kept.insert(0, 1);
                        exp10 += 1;
                        break;
                    }
                    i -= 1;
                    if kept[i] == 9 {
                        kept[i] = 0;
                    } else {
                        kept[i] += 1;
                        break;
                    }
                }
            }
        }
    }

    // Drop trailing zeros: minimumFractionDigits is 0.
    while kept.len() > 1 && *kept.last().unwrap() == 0 && kept.len() as i32 > exp10 {
        kept.pop();
    }

    let mut out = String::new();
    if exp10 <= 0 {
        // minimumIntegerDigits is 1, so a pure fraction gets a leading "0".
        out.push('0');
        out.push('.');
        for _ in 0..(-exp10) {
            out.push('0');
        }
        for d in &kept {
            out.push((b'0' + d) as char);
        }
    } else {
        let int_len = exp10 as usize;
        for i in 0..int_len {
            out.push(match kept.get(i) {
                Some(d) => (b'0' + d) as char,
                None => '0',
            });
        }
        if kept.len() > int_len {
            out.push('.');
            for d in &kept[int_len..] {
                out.push((b'0' + d) as char);
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nan_and_positive_infinity_are_indistinguishable() {
        assert_eq!(format_double(f64::NAN), "?");
        assert_eq!(format_double(f64::INFINITY), "?");
        assert_eq!(
            format_double(f64::NEG_INFINITY),
            "-?",
            "only the negative infinity is distinguishable, and only by its sign"
        );
    }

    #[test]
    fn negative_zero_keeps_its_sign() {
        assert_eq!(format_double(0.0), "0");
        assert_eq!(format_double(-0.0), "-0");
    }

    #[test]
    fn trailing_zeros_are_dropped_but_a_leading_zero_is_kept() {
        assert_eq!(format_double(1.0), "1");
        assert_eq!(format_double(0.5), "0.5");
        assert_eq!(format_double(1.5), "1.5");
        assert_eq!(format_double(-0.5), "-0.5");
    }

    #[test]
    fn six_fraction_digits_is_the_maximum() {
        assert_eq!(format_double(1.0 / 3.0), "0.333333");
        assert_eq!(format_double(2.0 / 3.0), "0.666667");
        assert_eq!(format_double(1.0 / 7.0), "0.142857");
        assert_eq!(format_double(std::f64::consts::PI), "3.141593");
    }

    /// The tie rule, which used to round up whenever the digit string was not the exact value and
    /// so was wrong on about half of the cases it decided.
    ///
    /// Every one of these prints a `5` at the seventh fraction digit and nothing after it, so no
    /// rule that looks only at the digit string can tell them apart. What separates them is which
    /// side of the halfway point the double sits on.
    #[test]
    fn an_apparent_tie_is_settled_by_the_true_value() {
        // 0.1234564999999999967…, below the halfway point, so down.
        assert_eq!(format_double(0.123_456_5), "0.123456");
        assert_eq!(format_double(0.123_457_5), "0.123457");
        assert_eq!(format_double(-0.123_456_5), "-0.123456");
        // 0.1234525000000000066…, above it, so up — and note the digit before the 5 is even, so
        // a half-even reading of the digit string would have sent this one the other way too.
        assert_eq!(format_double(0.123_452_5), "0.123453");
    }

    /// A genuine tie, where the digit string *is* the value, is the only place the rounding mode
    /// itself shows. `HALF_DOWN` is not `DecimalFormat`'s default, and this is what it buys.
    #[test]
    fn a_real_tie_goes_down_because_the_mode_is_half_down() {
        // 2^-21 is exactly 0.0000004768371582031250, a tie at the seventh digit.
        assert_eq!(format_double(4.768_371_582_031_25e-7), "0");
        // And a tie that is not at the underflow boundary.
        assert_eq!(format_double(0.062_500_05), "0.0625");
    }

    #[test]
    fn very_small_values_round_away_entirely() {
        assert_eq!(format_double(f64::MIN_POSITIVE), "0");
        assert_eq!(format_double(1e-7), "0");
    }

    /// No scientific notation, ever. `Double.MAX_VALUE` is written out in full.
    #[test]
    fn large_values_are_written_in_full() {
        let s = format_double(f64::MAX);
        assert!(!s.contains('e') && !s.contains('E'));
        assert_eq!(s.len(), 309);
        assert!(s.starts_with("17976931348623157"));
    }

    #[test]
    fn integers_have_no_grouping_separators() {
        assert_eq!(format_long(1_234_567), "1234567");
        assert_eq!(format_long(-1_234_567), "-1234567");
        assert_eq!(format_long(i64::MIN), "-9223372036854775808");
    }

    #[test]
    fn booleans_are_y_and_n() {
        assert_eq!(format_bool(true), "Y");
        assert_eq!(format_bool(false), "N");
    }
}
