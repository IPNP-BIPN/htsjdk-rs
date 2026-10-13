//! `Double.toString` and `Float.toString` as the pinned Java 17 runtime prints them.
//!
//! Java 17 does **not** print the shortest decimal that round-trips. Rust does, and so does every
//! JDK from 19 on, so the two agree on most values and disagree on a measurable set. This module
//! reproduces the Java 17 output exactly, including where it is longer than necessary or not the
//! nearest of two candidates.
//!
//! # Provenance
//!
//! The JDK's own converter is GPL-2.0 and this crate is MIT, so nothing here is a translation of
//! it (decision 0013). The behaviour was reconstructed from the outside: the oracle image's
//! `Double.toString` and `Float.toString` were queried on stratified samples, and an algorithm from
//! the public literature was fitted to what they print. That algorithm is Steele and White's
//! free-format digit generation ("How to Print Floating-Point Numbers Accurately", PLDI 1990),
//! with the decimal-exponent estimate from Gay ("Correctly Rounded Binary-Decimal and
//! Decimal-Binary Conversions", 1990). Each deviation from the textbook form below was found as a
//! cluster of disagreements, and each is pinned by the measurement that found it. Decision 0046
//! records the method and the counts.
//!
//! What Java 17 does, as measured:
//!
//! 1. **Bounds are exclusive.** A candidate exactly half an ulp away is rejected even when the
//!    round-half-even read would land back on the double. This is the whole of `1e23` printing as
//!    `9.999999999999999E22`.
//! 2. **A power of two gets the narrow interval on both sides.** The gap below a power of two is
//!    half the gap above it; Java 17 uses the smaller one for both, so it can print a digit more
//!    than needed. Subnormal powers of two get the same treatment although their gaps are equal.
//! 3. **An exact tie at the last digit goes to the even digit.**
//! 4. **Integers between 2^52 (2^23 for `float`) and 2^63 are printed from their exact value**,
//!    rounded half-up after dropping `digits(2^(s-2)) - 1` low digits, where `s` is the shift of
//!    the significand. The result is often 17 or 18 digits where 16 would round-trip.
//! 5. **The high-side test overflows a signed 64-bit integer** for some values near `10^25`:
//!    binary exponent 84 for `double`, up to 85 for `float`. Where it does, the last digit is kept
//!    instead of rounded up. The threshold is `(r + m) * 2 * 10^k >= 2^(63 + k)` with `k` the
//!    estimated decimal exponent, and it was located to within one representable value.
//! 6. **Six subnormals of each type** print a digit more, or a different digit, than the rules
//!    above give: `double` significands 1, 10, 12, 14, 16 and 18, and `float` significands 2, 3,
//!    6, 7, 21 and 29. For `double`, every subnormal significand up to 2,000,000 and every
//!    subnormal with up to 13 significant bits was checked; for `float`, every value.

use std::cmp::Ordering;

/// The decimal digits Java 17 prints for `|value|`, as `(digits, exp10)` where the magnitude is
/// `0.<digits> * 10^exp10`. `digits` has no leading or trailing zero.
///
/// `value` must be finite and non-zero.
pub fn double_digits(value: f64) -> (Vec<u8>, i32) {
    debug_assert!(value.is_finite() && value != 0.0);
    let bits = value.to_bits() & !(1u64 << 63);
    if let Some(exception) = SUBNORMAL_EXCEPTIONS.iter().find(|e| e.0 == bits) {
        return (exception.1.to_vec(), exception.2);
    }
    let biased = (bits >> 52) as i32;
    let fraction = bits & ((1u64 << 52) - 1);
    let (significand, exponent, subnormal) = if biased == 0 {
        (fraction, -1074, true)
    } else {
        (fraction | (1u64 << 52), biased - 1075, false)
    };
    digits_of(Binary {
        significand,
        exponent,
        subnormal,
        precision: 52,
        unbiased: biased - 1023,
        overflow_limit: 84,
    })
}

/// The decimal digits Java 17 prints for `|value|`, as [`double_digits`] does for a double.
pub fn float_digits(value: f32) -> (Vec<u8>, i32) {
    debug_assert!(value.is_finite() && value != 0.0);
    let bits = value.to_bits() & !(1u32 << 31);
    if let Some(exception) = FLOAT_SUBNORMAL_EXCEPTIONS.iter().find(|e| e.0 == bits) {
        return (exception.1.to_vec(), exception.2);
    }
    let biased = (bits >> 23) as i32;
    let fraction = u64::from(bits & ((1u32 << 23) - 1));
    let (significand, exponent, subnormal) = if biased == 0 {
        (fraction, -149, true)
    } else {
        (fraction | (1u64 << 23), biased - 150, false)
    };
    digits_of(Binary {
        significand,
        exponent,
        subnormal,
        precision: 23,
        unbiased: biased - 127,
        overflow_limit: 85,
    })
}

/// `Double.toString(value)`.
pub fn double_to_string(value: f64) -> String {
    if value.is_nan() {
        return "NaN".to_string();
    }
    let sign = if value.is_sign_negative() { "-" } else { "" };
    if value.is_infinite() {
        return format!("{sign}Infinity");
    }
    if value == 0.0 {
        return format!("{sign}0.0");
    }
    let (digits, exp10) = double_digits(value);
    format!("{sign}{}", layout(&digits, exp10))
}

/// `Float.toString(value)`.
pub fn float_to_string(value: f32) -> String {
    if value.is_nan() {
        return "NaN".to_string();
    }
    let sign = if value.is_sign_negative() { "-" } else { "" };
    if value.is_infinite() {
        return format!("{sign}Infinity");
    }
    if value == 0.0 {
        return format!("{sign}0.0");
    }
    let (digits, exp10) = float_digits(value);
    format!("{sign}{}", layout(&digits, exp10))
}

/// The `Double.toString` layout, which the javadoc specifies: plain notation for magnitudes in
/// `[10^-3, 10^7)`, computerized scientific notation otherwise, and always at least one digit after
/// the point.
fn layout(digits: &[u8], exp10: i32) -> String {
    let text: String = digits.iter().map(|d| char::from(b'0' + d)).collect();
    let length = text.len() as i32;
    if exp10 > 0 && exp10 <= 7 {
        if length <= exp10 {
            format!("{text}{}.0", "0".repeat((exp10 - length) as usize))
        } else {
            let (whole, fraction) = text.split_at(exp10 as usize);
            format!("{whole}.{fraction}")
        }
    } else if exp10 > -3 && exp10 <= 0 {
        format!("0.{}{text}", "0".repeat((-exp10) as usize))
    } else {
        let (lead, rest) = text.split_at(1);
        let rest = if rest.is_empty() { "0" } else { rest };
        format!("{lead}.{rest}E{}", exp10 - 1)
    }
}

/// The six measured `double` subnormal exceptions: raw bits, digits, `exp10`.
const SUBNORMAL_EXCEPTIONS: [(u64, &[u8], i32); 6] = [
    (1, &[4, 9], -323),
    (10, &[4, 9], -322),
    (12, &[5, 9], -322),
    (14, &[6, 9], -322),
    (16, &[7, 9], -322),
    (18, &[8, 9], -322),
];

/// The six measured `float` subnormal exceptions: raw bits, digits, `exp10`.
const FLOAT_SUBNORMAL_EXCEPTIONS: [(u32, &[u8], i32); 6] = [
    (2, &[2, 8], -44),
    (3, &[4, 2], -44),
    (6, &[8, 4], -44),
    (7, &[9, 8], -44),
    (21, &[2, 9], -43),
    (29, &[4, 1], -43),
];

struct Binary {
    /// The integer significand, with the implicit bit for a normal value.
    significand: u64,
    /// The power of two the significand is multiplied by.
    exponent: i32,
    subnormal: bool,
    /// Explicit fraction bits: 52 for `double`, 23 for `float`.
    precision: i32,
    /// The biased exponent minus the bias, so `-1023` (or `-127`) for every subnormal.
    unbiased: i32,
    /// The largest unbiased exponent at which the high-side test can overflow.
    overflow_limit: i32,
}

fn digits_of(value: Binary) -> (Vec<u8>, i32) {
    if (value.precision..=62).contains(&value.unbiased) {
        return exact_integer_digits(value.significand << (value.unbiased - value.precision), {
            let shift = value.unbiased - value.precision;
            if shift >= 2 {
                decimal_length(1u64 << (shift - 2)) - 1
            } else {
                0
            }
        });
    }
    free_format(&value)
}

/// Point 4: an integer below 2^63, printed from its exact value with `drop` low digits rounded
/// off half-up. A tie cannot arise, since the integer is a multiple of a power of two larger than
/// the dropped part.
fn exact_integer_digits(integer: u64, drop: u32) -> (Vec<u8>, i32) {
    let text = integer.to_string();
    let length = text.len() as i32;
    let keep = text.len() - drop as usize;
    let mut head: u64 = text[..keep].parse().expect("decimal digits");
    let tail = &text[keep..];
    if !tail.is_empty() && tail.as_bytes()[0] >= b'5' {
        head += 1;
    }
    let rounded = head.to_string();
    let exp10 = length + i32::from(rounded.len() > keep);
    (
        strip_trailing_zeros(rounded.bytes().map(|b| b - b'0').collect()),
        exp10,
    )
}

fn decimal_length(value: u64) -> u32 {
    value.to_string().len() as u32
}

/// Points 1, 2, 3 and 5: Steele and White's free-format generation with the measured choices.
fn free_format(value: &Binary) -> (Vec<u8>, i32) {
    let f = value.significand;
    let e = value.exponent;
    let narrow = f.is_power_of_two() && (f == 1u64 << value.precision || value.subnormal);
    let k = estimate_exp10(f, e);
    let overflow_applies = value.unbiased <= value.overflow_limit && k >= 0;

    // Every quantity stays below `11 * s` in the loop, so a scale under 2^122 runs in `u128`.
    let (mut digits, mut exp10) = if let Some((r, s, m)) = scaled_u128(f, e, narrow, k) {
        let overflow = overflow_applies.then(|| (s.to_big().shl(63 + k as u32), k as u32));
        generate(r, s, m, overflow.as_ref(), k)
    } else {
        let (r, s, m) = scaled_big(f, e, narrow, k);
        let overflow = overflow_applies.then(|| (s.shl(63 + k as u32), k as u32));
        generate(r, s, m, overflow.as_ref(), k)
    };

    // Carry a rounded-up 10 leftwards.
    let mut index = digits.len() - 1;
    while digits[index] == 10 {
        digits[index] = 0;
        if index == 0 {
            digits.insert(0, 1);
            exp10 += 1;
            break;
        }
        digits[index - 1] += 1;
        index -= 1;
    }
    let leading = digits.iter().take_while(|&&d| d == 0).count();
    digits.drain(..leading);
    exp10 -= leading as i32;
    (strip_trailing_zeros(digits), exp10)
}

/// The starting `(r, s, m)`: `r / s` is the value and `m / s` half the gap to a neighbour, all
/// scaled by two so that half an ulp is an integer, and by `10^k` on whichever side keeps them
/// integers. Point 2 doubles the value and the scale, which halves the interval on both sides.
fn scaled_big(f: u64, e: i32, narrow: bool, k: i32) -> (Big, Big, Big) {
    let (mut r, mut s, mut m) = if e >= 0 {
        (
            Big::from(f).shl(e as u32 + 1),
            Big::from(2),
            Big::from(1).shl(e as u32),
        )
    } else {
        (
            Big::from(f).shl(1),
            Big::from(1).shl((1 - e) as u32),
            Big::from(1),
        )
    };
    if narrow {
        r = r.shl(1);
        s = s.shl(1);
    }
    if k >= 0 {
        s = s.mul_pow10(k as u32);
    } else {
        r = r.mul_pow10((-k) as u32);
        m = m.mul_pow10((-k) as u32);
    }
    (r, s, m)
}

/// [`scaled_big`] in `u128`, when all three stay below 2^122.
fn scaled_u128(f: u64, e: i32, narrow: bool, k: i32) -> Option<(u128, u128, u128)> {
    const LIMIT: u32 = 122;
    let shift = |value: u128, by: i32| -> Option<u128> {
        (by >= 0 && 128 - value.leading_zeros() + by as u32 <= LIMIT).then(|| value << by)
    };
    let (mut r, mut s, mut m) = if e >= 0 {
        (shift(u128::from(f), e + 1)?, 2u128, shift(1, e)?)
    } else {
        (u128::from(f) << 1, shift(1, 1 - e)?, 1u128)
    };
    if narrow {
        r = shift(r, 1)?;
        s = shift(s, 1)?;
    }
    let power = 10u128.checked_pow(k.unsigned_abs())?;
    let fits = |value: u128| (128 - value.leading_zeros() <= LIMIT).then_some(value);
    if k >= 0 {
        s = fits(s.checked_mul(power)?)?;
    } else {
        r = fits(r.checked_mul(power)?)?;
        m = fits(m.checked_mul(power)?)?;
    }
    Some((r, s, m))
}

/// The digit loop, from the scaled value `r / s` and half-gap `m / s`. Returns the digits, last
/// one possibly 10 before the carry, and the decimal exponent.
///
/// `overflow` is point 5's threshold and `k` with it. Generic so that the common case avoids the
/// big integer; both instantiations run the same decisions.
fn generate<T: Scaled>(
    mut r: T,
    s: T,
    mut m: T,
    overflow: Option<&(Big, u32)>,
    k: i32,
) -> (Vec<u8>, i32) {
    let mut exp10 = k;
    let mut digits: Vec<u8> = Vec::with_capacity(20);
    let mut first = true;
    let (digit, low, high) = loop {
        r = r.times_ten();
        m = m.times_ten();
        let digit = r.take_digit(s);
        let low = r < m;
        let mut high = r.plus(m) > s;
        if high {
            if let Some((threshold, k)) = overflow {
                // Point 5: `(r + m) * 2 * 10^k >= s * 2^(63 + k)`.
                if r.plus(m).to_big().shl(1).mul_pow10(*k) >= *threshold {
                    high = false;
                }
            }
        }
        if first && digit == 0 && !(low || high) {
            // The estimate was one too high: drop the leading zero.
            exp10 -= 1;
            first = false;
            continue;
        }
        first = false;
        if low || high {
            break (digit, low, high);
        }
        digits.push(digit as u8);
    };

    let last = match (low, high) {
        (true, false) => digit,
        (false, true) => digit + 1,
        _ => match r.plus(r).cmp(&s) {
            Ordering::Less => digit,
            Ordering::Greater => digit + 1,
            // Point 3.
            Ordering::Equal => digit + (digit & 1),
        },
    };
    digits.push(last as u8);
    (digits, exp10)
}

/// The arithmetic [`generate`] needs.
trait Scaled: Copy + Ord {
    fn times_ten(self) -> Self;
    fn plus(self, other: Self) -> Self;
    /// Replaces `self` with `self % divisor` and returns the quotient, below ten here.
    fn take_digit(&mut self, divisor: Self) -> u32;
    fn to_big(self) -> Big;
}

impl Scaled for u128 {
    fn times_ten(self) -> Self {
        self * 10
    }

    fn plus(self, other: Self) -> Self {
        self + other
    }

    fn take_digit(&mut self, divisor: Self) -> u32 {
        // At most nine subtractions, which is cheaper than a 128-bit division.
        let mut quotient = 0;
        while *self >= divisor {
            *self -= divisor;
            quotient += 1;
        }
        quotient
    }

    fn to_big(self) -> Big {
        Big::from(self as u64).add(&Big::from((self >> 64) as u64).shl(64))
    }
}

impl Scaled for Big {
    fn times_ten(self) -> Self {
        self.mul_small(10)
    }

    fn plus(self, other: Self) -> Self {
        self.add(&other)
    }

    fn take_digit(&mut self, divisor: Self) -> u32 {
        self.div_rem_small_quotient(&divisor)
    }

    fn to_big(self) -> Big {
        self
    }
}

/// Gay's estimate of the decimal exponent, in the `0.d * 10^k` convention: `floor(log10 v) + 1`
/// or one more. Evaluated in `f64` in the published order, because where it rounds is observable
/// through point 5.
// `0.301029995663981` is Gay's published constant, not `LOG10_2`, which is a different double;
// where the estimate rounds is observable.
#[allow(clippy::approx_constant)]
fn estimate_exp10(significand: u64, exponent: i32) -> i32 {
    let bits = 64 - significand.leading_zeros() as i32;
    let binary_exponent = bits - 1 + exponent;
    let normalized = (significand << (53 - bits)) & ((1u64 << 52) - 1);
    let x = f64::from_bits((1023u64 << 52) | normalized);
    let estimate = (x - 1.5) * 0.289529654602168
        + 0.1760912590558
        + f64::from(binary_exponent) * 0.301029995663981;
    estimate.floor() as i32 + 1
}

fn strip_trailing_zeros(mut digits: Vec<u8>) -> Vec<u8> {
    while digits.len() > 1 && digits.last() == Some(&0) {
        digits.pop();
    }
    digits
}

/// The few operations the generator needs on unsigned integers wider than 64 bits: little-endian
/// base-2^32 limbs on the stack, never negative.
///
/// 48 limbs is 1,536 bits. The largest value the generator forms is a subnormal's remainder,
/// scaled by `10^323` and then by ten once per digit, under 1,200 bits.
#[derive(Clone, Copy, Debug)]
struct Big {
    len: usize,
    limbs: [u32; LIMBS],
}

const LIMBS: usize = 48;

impl Big {
    fn from(value: u64) -> Self {
        let mut limbs = [0u32; LIMBS];
        limbs[0] = value as u32;
        limbs[1] = (value >> 32) as u32;
        let mut big = Big { len: 2, limbs };
        big.trim();
        big
    }

    fn trim(&mut self) {
        while self.len > 1 && self.limbs[self.len - 1] == 0 {
            self.len -= 1;
        }
    }

    fn shl(&self, bits: u32) -> Self {
        let words = (bits / 32) as usize;
        let bits = bits % 32;
        let mut limbs = [0u32; LIMBS];
        let mut carry = 0u32;
        for index in 0..self.len {
            let limb = self.limbs[index];
            if bits == 0 {
                limbs[index + words] = limb;
            } else {
                limbs[index + words] = (limb << bits) | carry;
                carry = limb >> (32 - bits);
            }
        }
        limbs[self.len + words] = carry;
        let mut big = Big {
            len: self.len + words + 1,
            limbs,
        };
        big.trim();
        big
    }

    fn mul_small(&self, factor: u32) -> Self {
        let mut big = *self;
        let mut carry = 0u64;
        for limb in &mut big.limbs[..self.len] {
            let product = u64::from(*limb) * u64::from(factor) + carry;
            *limb = product as u32;
            carry = product >> 32;
        }
        big.limbs[self.len] = carry as u32;
        big.len += 1;
        big.trim();
        big
    }

    fn mul_pow10(&self, mut power: u32) -> Self {
        let mut big = *self;
        while power >= 9 {
            big = big.mul_small(1_000_000_000);
            power -= 9;
        }
        big.mul_small(10u32.pow(power))
    }

    fn add(&self, other: &Big) -> Self {
        let length = self.len.max(other.len);
        let mut big = Big {
            len: length + 1,
            limbs: [0u32; LIMBS],
        };
        let mut carry = 0u64;
        for index in 0..length {
            let sum = u64::from(self.limbs[index]) + u64::from(other.limbs[index]) + carry;
            big.limbs[index] = sum as u32;
            carry = sum >> 32;
        }
        big.limbs[length] = carry as u32;
        big.trim();
        big
    }

    /// `self - other`, which must not be negative.
    fn sub_assign(&mut self, other: &Big) {
        let mut borrow = 0i64;
        for index in 0..self.len {
            let difference = i64::from(self.limbs[index]) - i64::from(other.limbs[index]) - borrow;
            if difference < 0 {
                self.limbs[index] = (difference + (1i64 << 32)) as u32;
                borrow = 1;
            } else {
                self.limbs[index] = difference as u32;
                borrow = 0;
            }
        }
        debug_assert_eq!(borrow, 0);
        self.trim();
    }

    /// Replaces `self` with `self % divisor` and returns the quotient, which the generator keeps
    /// below ten.
    fn div_rem_small_quotient(&mut self, divisor: &Big) -> u32 {
        let mut quotient = 0;
        while *self >= *divisor {
            self.sub_assign(divisor);
            quotient += 1;
        }
        quotient
    }
}

impl PartialEq for Big {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for Big {}

impl PartialOrd for Big {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Big {
    fn cmp(&self, other: &Self) -> Ordering {
        self.len.cmp(&other.len).then_with(|| {
            self.limbs[..self.len]
                .iter()
                .rev()
                .cmp(other.limbs[..other.len].iter().rev())
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn values_named_in_the_module_documentation() {
        assert_eq!(double_to_string(1e23), "9.999999999999999E22");
        assert_eq!(double_to_string(f64::from_bits(1)), "4.9E-324");
        assert_eq!(double_to_string(0.1), "0.1");
        assert_eq!(double_to_string(1.0), "1.0");
        assert_eq!(double_to_string(-0.0), "-0.0");
        assert_eq!(double_to_string(1e7), "1.0E7");
        assert_eq!(double_to_string(0.001), "0.001");
        assert_eq!(double_to_string(1e-4), "1.0E-4");
        assert_eq!(float_to_string(f32::from_bits(1)), "1.4E-45");
        assert_eq!(float_to_string(0.1), "0.1");
        assert_eq!(float_to_string(f32::MAX), "3.4028235E38");
    }
}
