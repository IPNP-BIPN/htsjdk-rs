//! `SaddlePointExpansion` and the hypergeometric log-probability, ported from
//! `org.apache.commons.math3.distribution` (commons-math3 3.5), which is Apache 2.0.
//!
//! GATK's `FisherExactTest` builds a `HypergeometricDistribution` and asks it for the log
//! probability of every point of the support, so this is the arithmetic under `FS`, the
//! Fisher-strand annotation. Three of its decisions are not what a probability formula suggests.
//!
//! # A table of thirty-one exact values, and a series after that
//!
//! ```java
//! if (z < 15.0) {
//!     double z2 = 2.0 * z;
//!     if (FastMath.floor(z2) == z2) { ret = EXACT_STIRLING_ERRORS[(int) z2]; }
//!     else { ret = Gamma.logGamma(z + 1.0) - (z + 0.5) * FastMath.log(z) + z - HALF_LOG_2_PI; }
//! }
//! ```
//!
//! Below 15 and on a half-integer the answer is read out of a table of literals; below 15 and off
//! a half-integer it goes through `logGamma`; at 15 and above it is a five-term asymptotic series.
//! Three algorithms, and the boundaries are exact comparisons on a double.
//!
//! The table's **last entry is dead**: the guard is `z < 15.0`, so `z = 15.0` takes the series and
//! the tabulated value for 15.0 is never read. The two disagree in the thirteenth digit, which is
//! the size of the discontinuity the table exists to avoid.
//!
//! # The deviance part iterates until a float stops moving
//!
//! ```java
//! while (s1 != s) { s = s1; ej *= v; s1 = s + ej / ((j * 2) + 1); ++j; }
//! ```
//!
//! The loop's termination is `s1 != s` on doubles: it stops when the addition stops changing the
//! value, which is a property of the rounding rather than of a tolerance. A port with an epsilon
//! would stop at a different term.
//!
//! # `x == 0` and `x == n` are different formulas, and which one runs depends on `p`
//!
//! `logBinomialProbability` has four branches at the ends of its range, chosen by comparing `p`
//! or `q` against `0.1`, and only the middle one uses the Stirling error at all.

#![allow(clippy::excessive_precision)]

use crate::{fast_math, gamma};

const EXACT_STIRLING_ERRORS: [f64; 31] = [
    0.0,
    0.1534264097200273452913848,
    0.0810614667953272582196702,
    0.0548141210519176538961390,
    0.0413406959554092940938221,
    0.03316287351993628748511048,
    0.02767792568499833914878929,
    0.02374616365629749597132920,
    0.02079067210376509311152277,
    0.01848845053267318523077934,
    0.01664469118982119216319487,
    0.01513497322191737887351255,
    0.01387612882307074799874573,
    0.01281046524292022692424986,
    0.01189670994589177009505572,
    0.01110455975820691732662991,
    0.010411265261972096497478567,
    0.009799416126158803298389475,
    0.009255462182712732917728637,
    0.008768700134139385462952823,
    0.008330563433362871256469318,
    0.007934114564314020547248100,
    0.007573675487951840794972024,
    0.007244554301320383179543912,
    0.006942840107209529865664152,
    0.006665247032707682442354394,
    0.006408994188004207068439631,
    0.006171712263039457647532867,
    0.005951370112758847735624416,
    0.005746216513010115682023589,
    0.005554733551962801371038690,
];

/// `MathUtils.TWO_PI`, which is `2 * FastMath.PI` and therefore exactly `2 * std::f64::consts::PI`.
const TWO_PI: f64 = 2.0 * std::f64::consts::PI;

/// `SaddlePointExpansion.HALF_LOG_2_PI`, computed with `FastMath.log` as the reference computes it.
fn half_log_2_pi() -> f64 {
    0.5 * fast_math::log(TWO_PI)
}

/// `SaddlePointExpansion.getStirlingError`.
pub fn stirling_error(z: f64) -> f64 {
    if z < 15.0 {
        let z2 = 2.0 * z;
        if z2.floor() == z2 {
            return EXACT_STIRLING_ERRORS[z2 as usize];
        }
        return gamma::log_gamma(z + 1.0) - (z + 0.5) * fast_math::log(z) + z - half_log_2_pi();
    }
    let z2 = z * z;
    (0.083333333333333333333
        - (0.00277777777777777777778
            - (0.00079365079365079365079365
                - (0.000595238095238095238095238 - 0.0008417508417508417508417508 / z2) / z2)
                / z2)
            / z2)
        / z
}

/// `SaddlePointExpansion.getDeviancePart`.
pub fn deviance_part(x: f64, mu: f64) -> f64 {
    if (x - mu).abs() < 0.1 * (x + mu) {
        let d = x - mu;
        let mut v = d / (x + mu);
        let mut s1 = v * d;
        let mut s = f64::NAN;
        let mut ej = 2.0 * x * v;
        v *= v;
        let mut j = 1i32;
        // The loop stops when the addition stops moving the double, not at a tolerance.
        while s1 != s {
            s = s1;
            ej *= v;
            s1 = s + ej / ((j * 2) + 1) as f64;
            j += 1;
        }
        return s1;
    }
    x * fast_math::log(x / mu) + mu - x
}

/// `SaddlePointExpansion.logBinomialProbability`.
pub fn log_binomial_probability(x: i32, n: i32, p: f64, q: f64) -> f64 {
    let (xf, nf) = (x as f64, n as f64);
    if x == 0 {
        return if p < 0.1 {
            -deviance_part(nf, nf * q) - nf * p
        } else {
            nf * fast_math::log(q)
        };
    }
    if x == n {
        return if q < 0.1 {
            -deviance_part(nf, nf * p) - nf * q
        } else {
            nf * fast_math::log(p)
        };
    }
    let ret = stirling_error(nf)
        - stirling_error(xf)
        - stirling_error(nf - xf)
        - deviance_part(xf, nf * p)
        - deviance_part(nf - xf, nf * q);
    let f = (TWO_PI * xf * (nf - xf)) / nf;
    -0.5 * fast_math::log(f) + ret
}

/// `HypergeometricDistribution.getDomain`, which is where the support of the test comes from.
pub fn hypergeometric_domain(population: i32, successes: i32, sample: i32) -> (i32, i32) {
    (
        (sample - (population - successes)).max(0),
        successes.min(sample),
    )
}

/// `HypergeometricDistribution.logProbability(x)`.
///
/// Three binomial log-probabilities combined as `p1 + p2 - p3`, which is not the same double as
/// the logarithm of the hypergeometric formula written directly.
pub fn hypergeometric_log_probability(population: i32, successes: i32, sample: i32, x: i32) -> f64 {
    let (lo, hi) = hypergeometric_domain(population, successes, sample);
    if x < lo || x > hi {
        return f64::NEG_INFINITY;
    }
    let p = sample as f64 / population as f64;
    let q = (population - sample) as f64 / population as f64;
    let p1 = log_binomial_probability(x, successes, p, q);
    let p2 = log_binomial_probability(sample - x, population - successes, p, q);
    let p3 = log_binomial_probability(sample, population, p, q);
    p1 + p2 - p3
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_table_covers_every_half_integer_below_fifteen() {
        assert_eq!(EXACT_STIRLING_ERRORS.len(), 31);
        assert_eq!(stirling_error(0.0), 0.0);
        // The last entry of the table is dead code: the guard is `z < 15.0`, so 15.0 itself
        // takes the asymptotic series and the tabulated value for it is never read. The two
        // disagree in the thirteenth digit, which is the size of the discontinuity the table was
        // there to avoid.
        assert_ne!(stirling_error(15.0), EXACT_STIRLING_ERRORS[30]);
        assert!((stirling_error(15.0) - EXACT_STIRLING_ERRORS[30]).abs() < 1e-15);
        assert_eq!(stirling_error(14.5), EXACT_STIRLING_ERRORS[29]);
    }

    #[test]
    fn the_domain_is_the_overlap_of_the_two_ways_of_running_out() {
        assert_eq!(hypergeometric_domain(10, 5, 5), (0, 5));
        assert_eq!(hypergeometric_domain(10, 8, 5), (3, 5));
    }

    #[test]
    fn a_point_outside_the_support_is_minus_infinity() {
        assert_eq!(
            hypergeometric_log_probability(10, 8, 5, 0),
            f64::NEG_INFINITY
        );
    }
}
