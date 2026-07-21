//! `%.Nf` and `%.Ne` as the JVM produces them.
//!
//! **Not a port.** This module is written against *measured behaviour*, not against source, and
//! the distinction is the whole reason it exists. `VCFEncoder` formats every double through
//! `String.format(Locale.US, "%.2f", d)` and friends, which resolve into `java.util.Formatter`
//! and `FloatingDecimal` in `java.base`. Decision 0014 established that `java.base` is GPL2
//! (with the Classpath Exception, which permits linking and not translation), so those symbols
//! cannot be transcribed. The provenance audit forbids claiming them.
//!
//! What *can* be done is state the contract the goldens exhibit and implement it independently:
//!
//! > The JVM formats the **shortest decimal representation that round-trips** to the double,
//! > then rounds *that digit string* half-up (away from zero) at the requested precision,
//! > padding with zeros when the request outruns the digits.
//!
//! Rust's `{:e}` produces a shortest round-tripping representation too, so the two agree
//! wherever the JVM's own shortest representation is genuinely shortest. It is not always:
//! decision 0011 measured that Java 17's `FloatingDecimal` predates the Schubfach rewrite that
//! landed in JDK 19 and emits a non-shortest digit string for a small fraction of doubles. That
//! fraction is this module's divergence, it is measured rather than assumed, and it is the same
//! residue already quarantined in decision 0013.
//!
//! The consequence is worth stating plainly, because it is not what a reader expects. Rust's
//! `format!("{:.2}", d)` rounds the **exact binary value**; the JVM rounds a **short decimal
//! approximation of it**. Those differ for any tie that the shortest representation invents:
//! `0.125` prints `0.13` on the JVM and `0.12` in Rust, and `1e300` prints 301 digits of zeros
//! on the JVM against the 301 digits of its exact value in Rust. So the ergonomic Rust default
//! is wrong here in both the small and the large, and neither case announces itself.

/// The digit string and decimal exponent of a finite, non-zero magnitude.
///
/// `value = 0.<digits> * 10^(exp10 + 1)`, or equivalently `digits[0].digits[1..] * 10^exp10`.
struct Shortest {
    digits: Vec<u8>,
    exp10: i32,
}

fn shortest(magnitude: f64) -> Shortest {
    // `{:e}` renders `d.dddde±X`, shortest round-trip, which is exactly the pair we need.
    let s = format!("{:e}", magnitude);
    let (mantissa, exp) = s.split_once('e').expect("{:e} always has an exponent");
    let digits: Vec<u8> = mantissa.bytes().filter(|b| b.is_ascii_digit()).collect();
    Shortest {
        digits,
        exp10: exp.parse().expect("{:e} exponent is an integer"),
    }
}

/// Rounds `s` half-up (away from zero) to `keep` leading digits.
///
/// Returns the kept digits; `s.exp10` is incremented in place when the carry adds a digit, since
/// `999` rounding to two digits is `10` at one greater exponent, not `100`.
fn round_half_up(s: &mut Shortest, keep: usize) -> Vec<u8> {
    if keep >= s.digits.len() {
        return s.digits.clone();
    }
    let mut kept = s.digits[..keep].to_vec();
    // Half-up on the *digit string*: only the first dropped digit is consulted, so 0.1249999
    // stays 0.12 while 0.125 becomes 0.13. A rule that looked at the remaining digits would
    // agree on every case but the exact tie, which is the case that matters.
    if s.digits[keep] < b'5' {
        return kept;
    }
    for i in (0..kept.len()).rev() {
        if kept[i] == b'9' {
            kept[i] = b'0';
        } else {
            kept[i] += 1;
            return kept;
        }
    }
    // Every digit carried: the value grew an order of magnitude.
    kept.insert(0, b'1');
    kept.pop();
    s.exp10 += 1;
    kept
}

/// The sign prefix. `-0.0` keeps its sign, as the JVM does.
fn sign_of(d: f64) -> &'static str {
    if d.is_sign_negative() {
        "-"
    } else {
        ""
    }
}

/// `String.format(Locale.US, "%.<prec>f", d)`.
pub fn format_fixed(d: f64, prec: usize) -> String {
    if let Some(s) = non_finite(d) {
        return s;
    }
    let sign = sign_of(d);
    let magnitude = d.abs();
    if magnitude == 0.0 {
        return format!("{sign}{}", zero_with(prec));
    }

    let mut s = shortest(magnitude);
    // The digit at index `i` has place value 10^(exp10 - i), so digits down to 10^-prec are the
    // ones that survive: `exp10 - i >= -prec`.
    let keep = s.exp10 + prec as i32 + 1;
    if keep <= 0 {
        // No digit of the value reaches the requested precision. It still rounds *up* when the
        // first digit is 5 or more and it sits immediately below the last place, which is how
        // 0.005 at two decimals becomes 0.01 rather than 0.00.
        return if keep == 0 && s.digits[0] >= b'5' {
            format!("{sign}{}", one_in_last_place(prec))
        } else {
            format!("{sign}{}", zero_with(prec))
        };
    }
    let digits = round_half_up(&mut s, keep as usize);

    let mut out = String::from(sign);
    if s.exp10 >= 0 {
        let int_len = s.exp10 as usize + 1;
        for i in 0..int_len {
            out.push(*digits.get(i).unwrap_or(&b'0') as char);
        }
        if prec > 0 {
            out.push('.');
            for i in 0..prec {
                out.push(*digits.get(int_len + i).unwrap_or(&b'0') as char);
            }
        }
    } else {
        out.push('0');
        if prec > 0 {
            out.push('.');
            let leading_zeros = (-s.exp10 - 1) as usize;
            for i in 0..prec {
                if i < leading_zeros {
                    out.push('0');
                } else {
                    out.push(*digits.get(i - leading_zeros).unwrap_or(&b'0') as char);
                }
            }
        }
    }
    out
}

/// `String.format(Locale.US, "%.<prec>e", d)`.
pub fn format_scientific(d: f64, prec: usize) -> String {
    if let Some(s) = non_finite(d) {
        return s;
    }
    let sign = sign_of(d);
    let magnitude = d.abs();
    if magnitude == 0.0 {
        return format!("{sign}{}e+00", zero_with(prec));
    }

    let mut s = shortest(magnitude);
    let digits = round_half_up(&mut s, prec + 1);

    let mut out = String::from(sign);
    out.push(digits[0] as char);
    if prec > 0 {
        out.push('.');
        for i in 0..prec {
            out.push(*digits.get(1 + i).unwrap_or(&b'0') as char);
        }
    }
    // The exponent is always signed and at least two digits wide.
    out.push('e');
    out.push(if s.exp10 < 0 { '-' } else { '+' });
    out.push_str(&format!("{:02}", s.exp10.abs()));
    out
}

fn non_finite(d: f64) -> Option<String> {
    if d.is_nan() {
        Some("NaN".to_string())
    } else if d.is_infinite() {
        Some(if d < 0.0 { "-Infinity" } else { "Infinity" }.to_string())
    } else {
        None
    }
}

fn one_in_last_place(prec: usize) -> String {
    if prec == 0 {
        "1".to_string()
    } else {
        format!("0.{}1", "0".repeat(prec - 1))
    }
}

fn zero_with(prec: usize) -> String {
    if prec == 0 {
        "0".to_string()
    } else {
        format!("0.{}", "0".repeat(prec))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The two cases where rounding the shortest decimal differs from rounding the exact value.
    /// Both are what the JVM prints; both are what Rust's own `{:.2}` does not.
    #[test]
    fn ties_round_away_from_zero_on_the_short_decimal() {
        assert_eq!(format_fixed(0.125, 2), "0.13");
        assert_eq!(format!("{:.2}", 0.125), "0.12", "Rust rounds the exact value");

        // 2.675 is really 2.67499999999999982236431605997495353221893310546875.
        assert_eq!(format_fixed(2.675, 2), "2.68");
        assert_eq!(format!("{:.2}", 2.675), "2.67", "Rust rounds the exact value");
    }

    /// A large double prints its *short* digits padded with zeros, not its exact expansion.
    #[test]
    fn a_large_double_is_padded_not_expanded() {
        let ours = format_fixed(1e300, 2);
        assert!(ours.starts_with('1'));
        assert_eq!(ours.len(), 1 + 300 + 3, "one digit, 300 zeros, '.00'");
        assert!(ours[1..301].bytes().all(|b| b == b'0'));
        assert_ne!(format!("{:.2}", 1e300), ours, "Rust expands the exact value");
    }

    #[test]
    fn the_carry_can_grow_the_number() {
        assert_eq!(format_fixed(9.999, 2), "10.00");
        assert_eq!(format_fixed(0.999, 2), "1.00");
        assert_eq!(format_scientific(9.9999, 3), "1.000e+01");
    }

    #[test]
    fn a_value_below_the_precision_rounds_to_zero_or_one_ulp_of_it() {
        assert_eq!(format_fixed(0.001, 2), "0.00");
        assert_eq!(format_fixed(0.005, 2), "0.01");
        assert_eq!(format_fixed(-0.001, 2), "-0.00", "the sign survives");
    }

    #[test]
    fn negative_zero_keeps_its_sign() {
        assert_eq!(format_fixed(-0.0, 2), "-0.00");
        assert_eq!(format_fixed(0.0, 2), "0.00");
    }

    #[test]
    fn scientific_pads_the_exponent_to_two_digits() {
        assert_eq!(format_scientific(0.001, 3), "1.000e-03");
        assert_eq!(format_scientific(1e10, 3), "1.000e+10");
        assert_eq!(format_scientific(-0.5, 3), "-5.000e-01");
        assert_eq!(format_scientific(1e-100, 3), "1.000e-100");
    }
}
