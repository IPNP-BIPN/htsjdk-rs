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
    let body = round_to_fraction_digits(&digits, exp10, DECIMAL_DIGITS_TO_PRINT);
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

/// Whether `0.<digits> * 10^exp10` is the double's exact value, `valueExactAsDecimal` in
/// `DigitList`.
///
/// A decimal `n / 10^k` is exactly representable in binary only when it is a dyadic rational,
/// which happens exactly when `5^k` divides `n`. Since the digit string is already the shortest
/// decimal that round-trips, the double equals it exactly precisely in that case.
fn is_exact_decimal(digits: &str, exp10: i32) -> bool {
    let k = digits.len() as i32 - exp10;
    if k <= 0 {
        return true; // an integer, exactly representable if it round-trips at all
    }
    // 5^k outgrows any 17-digit mantissa well before k = 28, so it cannot divide.
    if k > 27 {
        return false;
    }
    let Ok(n) = digits.parse::<u128>() else {
        return false;
    };
    let pow = 5u128.pow(k as u32);
    n % pow == 0
}

/// Rounds `0.<digits> * 10^exp10` to at most `max_fraction` fraction digits, HALF_DOWN, and
/// renders it without grouping, without trailing zeros, and with at least one integer digit.
fn round_to_fraction_digits(digits: &str, exp10: i32, max_fraction: usize) -> String {
    // Position, counted in digits after the decimal point, of the last digit we keep.
    // A digit at index i in `digits` sits at fraction position i + 1 - exp10.
    let keep = max_fraction as i32 + exp10;

    let mut kept: Vec<u8>;
    let mut exp10 = exp10;

    // Whether the digit string represents the double *exactly*, which is what
    // `DigitList.shouldRoundUp` calls `valueExactAsDecimal` and what decides a tie.
    let exact = is_exact_decimal(digits, exp10);

    if keep <= 0 {
        // Everything is beyond the last kept place.
        let lone_five = keep == 0 && digits.as_bytes()[0] == b'5' && digits.len() == 1;
        if keep == 0 && digits.as_bytes()[0] >= b'5' && !(lone_five && exact) {
            kept = vec![1];
            exp10 += 1;
        } else {
            return "0".to_string();
        }
    } else {
        let keep = keep as usize;
        if keep >= digits.len() {
            kept = digits.bytes().map(|b| b - b'0').collect();
        } else {
            kept = digits.bytes().take(keep).map(|b| b - b'0').collect();
            let rest = &digits[keep..];
            let first_dropped = rest.as_bytes()[0];
            let past_half = rest.bytes().skip(1).any(|b| b != b'0');
            // `DigitList.shouldRoundUp`, HALF_DOWN branch, partially ported. A leading '5'
            // with nothing after it looks like a tie, but it only *is* one when the digits
            // represent the value exactly; otherwise the true binary value sits off the tie.
            //
            // INCOMPLETE. Java also consults `alreadyRounded`, which says whether
            // `FloatingDecimal` rounded *up* to produce the digit string, and returns false in
            // that case. Reproducing it needs the exact decimal expansion of the double, which
            // is most of `FloatingDecimal` itself. See decision 0011: the residual divergence
            // is measured and pinned rather than fitted away.
            let round_up = first_dropped > b'5' || (first_dropped == b'5' && (past_half || !exact));
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
