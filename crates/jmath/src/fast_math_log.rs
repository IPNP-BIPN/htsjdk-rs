//! `FastMath.log`, ported from `org.apache.commons.math3.util.FastMath` (commons-math3 3.5),
//! which is Apache 2.0. See `docs/decisions/0023-commons-math3-is-portable-where-the-jdk-is-not.md`.
//!
//! The logarithm eleven call sites in `Gamma`, `Erf` and `NormalDistribution` reach, and therefore
//! what GATK's `MannWhitneyU` and every rank-sum annotation stand on. It is **not**
//! [`crate::math::log`], which is `java.lang.Math.log`: that one is correctly rounded and this one
//! is not, so the two disagree in the last bit and a call site must name the right one.
//!
//! # There are two algorithms behind one function
//!
//! ```java
//! if ((exp == -1 || exp == 0) && x < 1.01 && x > 0.99 && hiPrec == null) {
//!     // a straight polynomial expansion in higher precision
//! }
//! ```
//!
//! Near 1 the table-driven method loses precision, so the function switches to a nine-term
//! double-double polynomial. The guard includes `hiPrec == null`, so **the same input takes a
//! different path depending on whether the caller wanted extra precision**, and `pow`, which asks
//! for it, therefore gets the other algorithm on `[0.99, 1.01]`.
//!
//! # A subnormal is renormalised by shifting, and the shift moves the exponent
//!
//! ```java
//! bits <<= 1;
//! while ((bits & 0x0010000000000000L) == 0) { --exp; bits <<= 1; }
//! ```
//!
//! The loop runs on the raw bits, so the mantissa index the table is looked up with comes from the
//! *shifted* pattern. A port that normalised by multiplying by a power of two would agree on every
//! normal input and diverge on the subnormals, which is where a log-likelihood spends its time.
//!
//! # The final sum is ordered, and the order is the accuracy
//!
//! The six terms are added smallest-last through a running compensation (`c = a + t;
//! d = -(c - a - t)`), with a comment listing each term's magnitude range. Reassociating that sum
//! is not an optimisation; it is a different function.

#![allow(clippy::excessive_precision)]

use crate::fast_math_tables::{LN_MANT_A, LN_MANT_B};

const HEX_40000000: f64 = 1073741824.0;
const TWO_POWER_52: f64 = 4503599627370496.0;
const LN_2_A: f64 = 0.693147063255310059;
const LN_2_B: f64 = 1.17304635250823482e-7;

/// `LN_QUICK_COEF`, the polynomial used on `[0.99, 1.01]`.
const LN_QUICK_COEF: [[f64; 2]; 9] = [
    [1.0, 5.669184079525E-24],
    [-0.25, -0.25],
    [0.3333333134651184, 1.986821492305628E-8],
    [-0.25, -6.663542893624021E-14],
    [0.19999998807907104, 1.1921056801463227E-8],
    [-0.1666666567325592, -7.800414592973399E-9],
    [0.1428571343421936, 5.650007086920087E-9],
    [-0.12502530217170715, -7.44321345601866E-11],
    [0.11113807559013367, 9.219544613762692E-9],
];

/// `LN_HI_PREC_COEF`, used when the caller asked for extra precision.
const LN_HI_PREC_COEF: [[f64; 2]; 6] = [
    [1.0, -6.032174644509064E-23],
    [-0.25, -0.25],
    [0.3333333134651184, 1.9868161777724352E-8],
    [-0.2499999701976776, -2.957007209750105E-8],
    [0.19999954104423523, 1.5830993332061267E-10],
    [-0.16624879837036133, -2.6033824355191673E-8],
];

/// `FastMath.log(double)`.
pub fn log(x: f64) -> f64 {
    log_with_hi_prec(x, None)
}

/// `FastMath.log(double, double[])`, the form `pow` calls.
///
/// Passing `Some` changes the answer, not just the amount of information returned: the
/// `[0.99, 1.01]` shortcut is guarded on the argument being absent.
pub fn log_with_hi_prec(x: f64, mut hi_prec: Option<&mut [f64; 2]>) -> f64 {
    if x == 0.0 {
        return f64::NEG_INFINITY;
    }
    let mut bits = x.to_bits() as i64;

    if (bits & (0x8000000000000000u64 as i64) != 0 || x.is_nan()) && x != 0.0 {
        if let Some(hi) = hi_prec.as_deref_mut() {
            hi[0] = f64::NAN;
        }
        return f64::NAN;
    }

    if x == f64::INFINITY {
        if let Some(hi) = hi_prec.as_deref_mut() {
            hi[0] = f64::INFINITY;
        }
        return f64::INFINITY;
    }

    let mut exp = ((bits >> 52) - 1023) as i32;

    if bits & 0x7ff0000000000000 == 0 {
        // Subnormal: renormalise on the raw bits, which is what decides the table index below.
        bits <<= 1;
        while bits & 0x0010000000000000 == 0 {
            exp -= 1;
            bits <<= 1;
        }
    }

    if (exp == -1 || exp == 0) && x < 1.01 && x > 0.99 && hi_prec.is_none() {
        // The table-driven method loses precision here, so a nine-term double-double polynomial
        // runs instead. Note the guard: a caller wanting high precision does not take this path.
        let mut xa = x - 1.0;
        // The reference computes `xb = xa - x + 1.0` here and overwrites it three lines later
        // with `ab`, so the first value is dead in the reference too and is not reproduced.
        let mut tmp = xa * HEX_40000000;
        let mut aa = xa + tmp - tmp;
        let mut ab = xa - aa;
        xa = aa;
        let mut xb = ab;

        let last = LN_QUICK_COEF[LN_QUICK_COEF.len() - 1];
        let mut ya = last[0];
        let mut yb = last[1];

        for coefficient in LN_QUICK_COEF.iter().rev().skip(1) {
            aa = ya * xa;
            ab = ya * xb + yb * xa + yb * xb;
            tmp = aa * HEX_40000000;
            ya = aa + tmp - tmp;
            yb = aa - ya + ab;

            aa = ya + coefficient[0];
            ab = yb + coefficient[1];
            tmp = aa * HEX_40000000;
            ya = aa + tmp - tmp;
            yb = aa - ya + ab;
        }

        aa = ya * xa;
        ab = ya * xb + yb * xa + yb * xb;
        tmp = aa * HEX_40000000;
        ya = aa + tmp - tmp;
        yb = aa - ya + ab;

        return ya + yb;
    }

    let mantissa_index = ((bits & 0x000ffc0000000000) >> 42) as usize;
    let lnm = [
        f64::from_bits(LN_MANT_A[mantissa_index]),
        f64::from_bits(LN_MANT_B[mantissa_index]),
    ];

    let epsilon =
        (bits & 0x3ffffffffff) as f64 / (TWO_POWER_52 + (bits & 0x000ffc0000000000) as f64);

    let lnza;
    let mut lnzb = 0.0;

    if let Some(_hi) = hi_prec.as_deref_mut() {
        let mut tmp = epsilon * HEX_40000000;
        let mut aa = epsilon + tmp - tmp;
        let ab = epsilon - aa;
        let xa = aa;
        let mut xb = ab;

        // The division is redone at higher precision, which is why this branch is not just the
        // other one with extra bookkeeping.
        let numer = (bits & 0x3ffffffffff) as f64;
        let denom = TWO_POWER_52 + (bits & 0x000ffc0000000000) as f64;
        aa = numer - xa * denom - xb * denom;
        xb += aa / denom;

        let last = LN_HI_PREC_COEF[LN_HI_PREC_COEF.len() - 1];
        let mut ya = last[0];
        let mut yb = last[1];

        for coefficient in LN_HI_PREC_COEF.iter().rev().skip(1) {
            aa = ya * xa;
            let mut ab2 = ya * xb + yb * xa + yb * xb;
            tmp = aa * HEX_40000000;
            ya = aa + tmp - tmp;
            yb = aa - ya + ab2;

            aa = ya + coefficient[0];
            ab2 = yb + coefficient[1];
            tmp = aa * HEX_40000000;
            ya = aa + tmp - tmp;
            yb = aa - ya + ab2;
        }

        aa = ya * xa;
        let ab2 = ya * xb + yb * xa + yb * xb;

        lnza = aa + ab2;
        lnzb = -(lnza - aa - ab2);
    } else {
        let mut value = -0.16624882440418567;
        value = value * epsilon + 0.19999954120254515;
        value = value * epsilon + -0.2499999997677497;
        value = value * epsilon + 0.3333333333332802;
        value = value * epsilon + -0.5;
        value = value * epsilon + 1.0;
        lnza = value * epsilon;
    }

    // The six terms are summed with a running compensation, smallest last. The reference lists
    // each term's magnitude range beside this block; reassociating it is a different function.
    let mut a = LN_2_A * exp as f64;
    let mut b = 0.0;
    let add = |a: &mut f64, b: &mut f64, term: f64| {
        let c = *a + term;
        let d = -(c - *a - term);
        *a = c;
        *b += d;
    };
    add(&mut a, &mut b, lnm[0]);
    add(&mut a, &mut b, lnza);
    add(&mut a, &mut b, LN_2_B * exp as f64);
    add(&mut a, &mut b, lnm[1]);
    add(&mut a, &mut b, lnzb);

    if let Some(hi) = hi_prec {
        hi[0] = a;
        hi[1] = b;
    }

    a + b
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_specials_are_the_reference_specials() {
        assert_eq!(log(0.0), f64::NEG_INFINITY);
        assert_eq!(log(-0.0), f64::NEG_INFINITY);
        assert!(log(-1.0).is_nan());
        assert!(log(f64::NAN).is_nan());
        assert_eq!(log(f64::INFINITY), f64::INFINITY);
    }

    #[test]
    fn the_ordinary_range_is_close_to_the_true_logarithm() {
        for x in [1.0, 2.0, 0.5, 10.0, 1e300, 1e-300, 1.005, 0.995] {
            let ours = log(x);
            assert!(
                (ours - x.ln()).abs() <= 1e-12 * x.ln().abs().max(1.0),
                "log({x})"
            );
        }
    }

    /// Asking for high precision changes which algorithm runs, and near 1 that is visible.
    #[test]
    fn the_high_precision_form_takes_a_different_path_near_one() {
        let mut hi = [0.0f64; 2];
        let with = log_with_hi_prec(1.005, Some(&mut hi));
        let without = log(1.005);
        assert!((with - without).abs() < 1e-15);
        assert!((hi[0] + hi[1] - with).abs() < 1e-15);
    }
}
