//! Correctly-rounded `log` and `log10`, matching `java.lang.Math`.
//!
//! `Math.log` and `Math.log10` were measured to be **correctly rounded** on every point of the
//! conformance corpus, so the target here is not "reproduce HotSpot's algorithm" but simply
//! "round the true result once". That is a much smaller job than porting the intrinsic, and it
//! is why these two land before `exp` and `pow`, which are *not* correctly rounded and do need
//! the algorithm ported.
//!
//! See `docs/decisions/0006-correct-rounding-is-the-target-for-log-and-log10.md`.
//!
//! `clippy::approx_constant` is allowed module-wide. It flags the `hi` halves of the
//! double-double constants as approximations of `std::f64::consts::LN_2` and `LOG10_E`, and
//! taking that advice would defeat the entire point: these are pairs precisely so they carry
//! ~53 bits *more* than a single `f64`, and `hi` alone is deliberately only the leading half.
//! Substituting the std constant silently discards the precision that makes correct rounding
//! possible, and the corpus would fail on exactly the hard-to-round points.
//!
//! `clippy::excessive_precision` is allowed for the same reason: the pairs were emitted at
//! 400-bit precision and the digits record what was generated, not what `f64` can hold.
#![allow(clippy::approx_constant, clippy::excessive_precision)]

use crate::dd::{self, DoubleDouble};

// ln(2) and log10(e) to ~106 bits, generated at 400-bit precision.
const LN2: DoubleDouble = DoubleDouble::new(6.93147180559945286e-01, 2.31904681384629956e-17);
const LOG10_E: DoubleDouble = DoubleDouble::new(4.34294481903251817e-01, 1.09831965021676507e-17);

/// `ln(x)` in double-double, for finite positive normal-or-subnormal `x`.
///
/// Argument reduction puts the mantissa in `[sqrt(1/2), sqrt(2))` so that `s = (m-1)/(m+1)`
/// satisfies `|s| <= 0.1716`. The atanh series in `s^2` then converges by a factor of ~0.029
/// per term, so 22 terms carry it past 106 bits.
fn ln_dd(x: f64) -> DoubleDouble {
    // Decompose x = m * 2^e with m in [1, 2).
    let mut e = ((x.to_bits() >> 52) & 0x7ff) as i64 - 1023;
    let mut m = f64::from_bits((x.to_bits() & 0x000f_ffff_ffff_ffff) | 0x3ff0_0000_0000_0000);
    if x.to_bits() >> 52 & 0x7ff == 0 {
        // Subnormal: scale into the normal range and correct the exponent afterwards.
        let scaled = x * f64::from_bits(0x4350_0000_0000_0000); // 2^54
        e = ((scaled.to_bits() >> 52) & 0x7ff) as i64 - 1023 - 54;
        m = f64::from_bits((scaled.to_bits() & 0x000f_ffff_ffff_ffff) | 0x3ff0_0000_0000_0000);
    }

    // Centre the mantissa on 1 to shrink |s|.
    if m > std::f64::consts::SQRT_2 {
        m *= 0.5;
        e += 1;
    }

    // s = (m - 1) / (m + 1), computed in double-double because the subtraction cancels.
    let num = dd::add_f64(DoubleDouble::from_f64(m), -1.0);
    let den = dd::add_f64(DoubleDouble::from_f64(m), 1.0);
    let s = dd::div(num, den);
    let s2 = dd::mul(s, s);

    // atanh(s) = s + s^3/3 + s^5/5 + ...  and  ln(m) = 2 * atanh(s)
    //
    // The coefficients are divided, not multiplied by a precomputed reciprocal. `1.0/3.0` is
    // not representable, so multiplying by it injects a relative error of about 2^-53 into the
    // largest term, which lands near 2^-60 overall. That is precisely the range where the
    // hard-to-round corpus cases live (they need 62 to 64 bits to resolve), so the reciprocal
    // form fails exactly the points this function exists to get right, and passes everywhere
    // else. Dividing by the exactly-representable odd integer keeps full precision.
    let mut sum = DoubleDouble::from_f64(1.0);
    let mut term = DoubleDouble::from_f64(1.0);
    for k in 1..=22 {
        term = dd::mul(term, s2);
        let denom = DoubleDouble::from_f64((2 * k + 1) as f64);
        sum = dd::add(sum, dd::div(term, denom));
    }
    let ln_m = dd::mul_f64(dd::mul(s, sum), 2.0);

    dd::add(ln_m, dd::mul_f64(LN2, e as f64))
}

/// `java.lang.Math.log`, correctly rounded.
pub fn log(x: f64) -> f64 {
    if x.is_nan() || x < 0.0 {
        return f64::NAN;
    }
    if x == 0.0 {
        return f64::NEG_INFINITY;
    }
    if x.is_infinite() {
        return f64::INFINITY;
    }
    if x == 1.0 {
        return 0.0;
    }
    ln_dd(x).to_f64()
}

/// `java.lang.Math.log10`, correctly rounded.
pub fn log10(x: f64) -> f64 {
    if x.is_nan() || x < 0.0 {
        return f64::NAN;
    }
    if x == 0.0 {
        return f64::NEG_INFINITY;
    }
    if x.is_infinite() {
        return f64::INFINITY;
    }
    if x == 1.0 {
        return 0.0;
    }

    // Exact powers of ten must return exact integers. The series route would land within one
    // ulp but not necessarily *on* the integer, and log10 of a power of ten showing up as
    // 2.9999999999999996 is the kind of thing that survives all the way into a report column.
    if x.fract() == 0.0 && x > 0.0 && x <= 1e22 {
        let mut p = 1.0f64;
        for k in 0..=22 {
            if p == x {
                return k as f64;
            }
            p *= 10.0;
        }
    }

    dd::mul(ln_dd(x), LOG10_E).to_f64()
}
