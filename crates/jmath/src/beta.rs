//! `Beta.logBeta` and `Beta.regularizedBeta`, ported from `org.apache.commons.math3.special.Beta`
//! (commons-math 3.5).
//!
//! The log of the beta function, which is what every beta-binomial likelihood is built from. GATK
//! reaches it through `BetaBinomialDistribution.logProbability`, and Mutect's somatic clustering
//! model sums that over its clusters.
//!
//! # It is not `logGamma(p) + logGamma(q) - logGamma(p + q)`
//!
//! That identity is true of the mathematics and false of the doubles. commons-math implements the
//! NSWC `DBETLN` routine, which is four branches on `min(p, q)`, each a different arrangement of
//! the same terms:
//!
//! ```java
//! final double a = FastMath.min(p, q);
//! final double b = FastMath.max(p, q);
//! if (a >= 10.0) { ... } else if (a > 2.0) { ... } else if (a >= 1.0) { ... } else { ... }
//! ```
//!
//! Only the last branch, and only when `b < 10`, uses gammas directly — and even there it is
//! `log(gamma(a) * gamma(b) / gamma(a + b))` rather than a sum of three `logGamma`s, with a comment
//! in the reference saying the original NSWC form was less accurate.
//!
//! # `logBeta(1, 1)` is negative zero
//!
//! The `a >= 1, b <= 2` branch is `logGamma(a) + logGamma(b) - logGammaSum(a, b)`, which at `(1, 1)`
//! is `0 + 0 - 0` — and comes out `-0.0`, because `logGamma(1.0)` is itself `-0.0` through
//! `logGamma1p(0.0)`. A caller comparing with `==` cannot see the difference; one formatting the
//! value can.
//!
//! # The delta table is Didonato and Morris, not Stirling
//!
//! The `a >= 10` branch corrects a Stirling-like approximation with `Δ(p) + Δ(q) - Δ(p + q)`,
//! evaluated from fifteen transcribed constants over a Horner scheme in `(10/a)^2`. The constants
//! are written here exactly as the reference writes them, at a precision Rust's parser accepts and
//! `clippy::excessive_precision` objects to.

#![allow(clippy::excessive_precision)]

use crate::continued_fraction::{self, ContinuedFractionError};
use crate::fast_math;
use crate::gamma::{gamma, log_gamma, log_gamma1p, GammaError};

/// `HALF_LOG_TWO_PI`, as the reference spells it rather than as a computed constant.
const HALF_LOG_TWO_PI: f64 = 0.9189385332046727;

/// `DELTA`, the coefficients of the Didonato and Morris correction.
const DELTA: [f64; 15] = [
    0.833333333333333333333333333333E-01,
    -0.277777777777777777777777752282E-04,
    0.793650793650793650791732130419E-07,
    -0.595238095238095232389839236182E-09,
    0.841750841750832853294451671990E-11,
    -0.191752691751854612334149171243E-12,
    0.641025640510325475730918472625E-14,
    -0.295506514125338232839867823991E-15,
    0.179643716359402238723287696452E-16,
    -0.139228964661627791231203060395E-17,
    0.133802855014020915603275339093E-18,
    -0.154246009867966094273710216533E-19,
    0.197701992980957427278370133333E-20,
    -0.234065664793997056856992426667E-21,
    0.171348014966398575409015466667E-22,
];

/// What the helpers refuse, each of them a range the reference checks before computing.
#[derive(Debug, Clone, PartialEq)]
pub enum BetaError {
    /// `OutOfRangeException` from `logGammaSum` or `deltaMinusDeltaSum`.
    OutOfRange { value: f64, low: f64, high: f64 },
    /// `NumberIsTooSmallException` from the three helpers that need `b >= 10`.
    TooSmall { value: f64, bound: f64 },
    /// A gamma the reference could compute and this port has not measured.
    Gamma(GammaError),
    /// The continued fraction of `regularizedBeta` refused to converge.
    ContinuedFraction(ContinuedFractionError),
}

impl From<GammaError> for BetaError {
    fn from(error: GammaError) -> Self {
        BetaError::Gamma(error)
    }
}

impl From<ContinuedFractionError> for BetaError {
    fn from(error: ContinuedFractionError) -> Self {
        BetaError::ContinuedFraction(error)
    }
}

/// `logGammaSum(a, b)`: `log(Gamma(a + b))` for `1 <= a, b <= 2`.
pub fn log_gamma_sum(a: f64, b: f64) -> Result<f64, BetaError> {
    if !(1.0..=2.0).contains(&a) {
        return Err(BetaError::OutOfRange {
            value: a,
            low: 1.0,
            high: 2.0,
        });
    }
    if !(1.0..=2.0).contains(&b) {
        return Err(BetaError::OutOfRange {
            value: b,
            low: 1.0,
            high: 2.0,
        });
    }
    // `(a - 1) + (b - 1)`, which is not `a + b - 2`.
    let x = (a - 1.0) + (b - 1.0);
    if x <= 0.5 {
        Ok(log_gamma1p(1.0 + x)?)
    } else if x <= 1.5 {
        Ok(log_gamma1p(x)? + fast_math::log1p(x))
    } else {
        Ok(log_gamma1p(x - 1.0)? + fast_math::log(x * (1.0 + x)))
    }
}

/// `deltaMinusDeltaSum(a, b)`: `Δ(b) - Δ(a + b)` for `0 <= a <= b` and `b >= 10`.
pub fn delta_minus_delta_sum(a: f64, b: f64) -> Result<f64, BetaError> {
    if a < 0.0 || a > b {
        return Err(BetaError::OutOfRange {
            value: a,
            low: 0.0,
            high: b,
        });
    }
    if b < 10.0 {
        return Err(BetaError::TooSmall {
            value: b,
            bound: 10.0,
        });
    }
    let h = a / b;
    let p = h / (1.0 + h);
    let q = 1.0 / (1.0 + h);
    let q2 = q * q;
    // `s[i] = 1 + q + ... - q**(2 * i)`, built forwards and consumed backwards.
    let mut s = [0.0f64; DELTA.len()];
    s[0] = 1.0;
    for i in 1..s.len() {
        s[i] = 1.0 + (q + q2 * s[i - 1]);
    }
    let sqrt_t = 10.0 / b;
    let t = sqrt_t * sqrt_t;
    let mut w = DELTA[DELTA.len() - 1] * s[s.len() - 1];
    for i in (0..DELTA.len() - 1).rev() {
        w = t * w + DELTA[i] * s[i];
    }
    Ok(w * p / b)
}

/// `sumDeltaMinusDeltaSum(p, q)`: `Δ(p) + Δ(q) - Δ(p + q)` for `p, q >= 10`.
pub fn sum_delta_minus_delta_sum(p: f64, q: f64) -> Result<f64, BetaError> {
    if p < 10.0 {
        return Err(BetaError::TooSmall {
            value: p,
            bound: 10.0,
        });
    }
    if q < 10.0 {
        return Err(BetaError::TooSmall {
            value: q,
            bound: 10.0,
        });
    }
    let a = p.min(q);
    let b = p.max(q);
    let sqrt_t = 10.0 / a;
    let t = sqrt_t * sqrt_t;
    let mut z = DELTA[DELTA.len() - 1];
    for i in (0..DELTA.len() - 1).rev() {
        z = t * z + DELTA[i];
    }
    Ok(z / a + delta_minus_delta_sum(a, b)?)
}

/// `logGammaMinusLogGammaSum(a, b)`: `log(Gamma(b) / Gamma(a + b))` for `a >= 0` and `b >= 10`.
pub fn log_gamma_minus_log_gamma_sum(a: f64, b: f64) -> Result<f64, BetaError> {
    if a < 0.0 {
        return Err(BetaError::TooSmall {
            value: a,
            bound: 0.0,
        });
    }
    if b < 10.0 {
        return Err(BetaError::TooSmall {
            value: b,
            bound: 10.0,
        });
    }
    // `d = a + b - 0.5`, written the way the smaller argument decides.
    let (d, w) = if a <= b {
        (b + (a - 0.5), delta_minus_delta_sum(a, b)?)
    } else {
        (a + (b - 0.5), delta_minus_delta_sum(b, a)?)
    };
    let u = d * fast_math::log1p(a / b);
    let v = a * (fast_math::log(b) - 1.0);
    // The larger term is subtracted last, which is not the same double as the other order.
    Ok(if u <= v { (w - u) - v } else { (w - v) - u })
}

/// `logBeta(p, q)`, all four branches.
///
/// `NaN` for a non-positive or NaN argument, as the reference answers rather than refusing.
pub fn log_beta(p: f64, q: f64) -> Result<f64, BetaError> {
    if p.is_nan() || q.is_nan() || p <= 0.0 || q <= 0.0 {
        return Ok(f64::NAN);
    }
    let a = p.min(q);
    let b = p.max(q);
    if a >= 10.0 {
        let w = sum_delta_minus_delta_sum(a, b)?;
        let h = a / b;
        let c = h / (1.0 + h);
        let u = -(a - 0.5) * fast_math::log(c);
        let v = b * fast_math::log1p(h);
        return Ok(if u <= v {
            (((-0.5 * fast_math::log(b) + HALF_LOG_TWO_PI) + w) - u) - v
        } else {
            (((-0.5 * fast_math::log(b) + HALF_LOG_TWO_PI) + w) - v) - u
        });
    }
    if a > 2.0 {
        if b > 1000.0 {
            let n = (a - 1.0).floor() as i32;
            let mut prod = 1.0;
            let mut ared = a;
            for _ in 0..n {
                ared -= 1.0;
                prod *= ared / (1.0 + ared / b);
            }
            return Ok((fast_math::log(prod) - n as f64 * fast_math::log(b))
                + (log_gamma(ared) + log_gamma_minus_log_gamma_sum(ared, b)?));
        }
        let mut prod1 = 1.0;
        let mut ared = a;
        while ared > 2.0 {
            ared -= 1.0;
            let h = ared / b;
            prod1 *= h / (1.0 + h);
        }
        if b < 10.0 {
            let mut prod2 = 1.0;
            let mut bred = b;
            while bred > 2.0 {
                bred -= 1.0;
                prod2 *= bred / (ared + bred);
            }
            return Ok(fast_math::log(prod1)
                + fast_math::log(prod2)
                + (log_gamma(ared) + (log_gamma(bred) - log_gamma_sum(ared, bred)?)));
        }
        return Ok(fast_math::log(prod1)
            + log_gamma(ared)
            + log_gamma_minus_log_gamma_sum(ared, b)?);
    }
    if a >= 1.0 {
        if b > 2.0 {
            if b < 10.0 {
                let mut prod = 1.0;
                let mut bred = b;
                while bred > 2.0 {
                    bred -= 1.0;
                    prod *= bred / (a + bred);
                }
                return Ok(fast_math::log(prod)
                    + (log_gamma(a) + (log_gamma(bred) - log_gamma_sum(a, bred)?)));
            }
            return Ok(log_gamma(a) + log_gamma_minus_log_gamma_sum(a, b)?);
        }
        // Where `logBeta(1, 1)` becomes negative zero.
        return Ok(log_gamma(a) + log_gamma(b) - log_gamma_sum(a, b)?);
    }
    if b >= 10.0 {
        return Ok(log_gamma(a) + log_gamma_minus_log_gamma_sum(a, b)?);
    }
    // The reference's own comment: the NSWC form was `logGamma(a) + (logGamma(b) - logGamma(a + b))`
    // and this one "turns out to be more accurate".
    Ok(fast_math::log(gamma(a)? * gamma(b)? / gamma(a + b)?))
}

/// `Beta.DEFAULT_EPSILON`, which is **not** `ContinuedFraction`'s own `10e-9`.
pub const DEFAULT_EPSILON: f64 = 1E-14;

/// `regularizedBeta(x, a, b)`, the regularized incomplete beta function `I(x, a, b)`.
///
/// `BinomialDistribution.cumulativeProbability(k)` is `1 - regularizedBeta(p, k + 1, n - k)`, which
/// is how Mutect's normal-artifact filter gets its p-value.
pub fn regularized_beta(x: f64, a: f64, b: f64) -> Result<f64, BetaError> {
    regularized_beta_with(x, a, b, DEFAULT_EPSILON, i32::MAX)
}

/// `regularizedBeta(x, a, b, epsilon, maxIterations)`.
///
/// # The symmetry test's second conjunct is implied by its first
///
/// `1 - x <= (b + 1) / (2 + b + a)` rearranges to `x >= (a + 1) / (2 + b + a)`, which the strict
/// `x > (a + 1) / (2 + b + a)` beside it already asserts. In the reals the `&&` is redundant; in
/// doubles the two sides are computed differently, so they can disagree by a rounding at the
/// boundary alone. Both are transcribed rather than simplified, because which of the two
/// arrangements a call takes decides its last digits: at `a = 2, b = 5`, one ULP either side of
/// `1/3` answers the same double by two different routes, and `1/3` itself answers another.
// The seven guards are one `||` chain in the reference, in this order. `!(0.0..=1.0).contains(&x)`
// would be the same predicate and a different program to read against the Java.
#[allow(clippy::manual_range_contains)]
pub fn regularized_beta_with(
    x: f64,
    a: f64,
    b: f64,
    epsilon: f64,
    max_iterations: i32,
) -> Result<f64, BetaError> {
    if x.is_nan() || a.is_nan() || b.is_nan() || x < 0.0 || x > 1.0 || a <= 0.0 || b <= 0.0 {
        return Ok(f64::NAN);
    }
    if x > (a + 1.0) / (2.0 + b + a) && 1.0 - x <= (b + 1.0) / (2.0 + b + a) {
        // The recursion swaps `a` and `b` as well as reflecting `x`.
        return Ok(1.0 - regularized_beta_with(1.0 - x, b, a, epsilon, max_iterations)?);
    }
    let fraction = continued_fraction::evaluate(
        |_, _| 1.0,
        |n, x| {
            if n % 2 == 0 {
                let m = f64::from(n) / 2.0;
                (m * (b - m) * x) / ((a + (2.0 * m) - 1.0) * (a + (2.0 * m)))
            } else {
                let m = (f64::from(n) - 1.0) / 2.0;
                -((a + m) * (a + b + m) * x) / ((a + (2.0 * m)) * (a + (2.0 * m) + 1.0))
            }
        },
        x,
        epsilon,
        max_iterations,
    )?;
    // `FastMath.exp(...) * 1.0 / fraction`, which is `(exp(...) * 1.0) / fraction`: the reference
    // writes the numerator of the fraction's reciprocal out, and it is exactly one.
    let numerator = fast_math::exp(
        (a * fast_math::log(x)) + (b * fast_math::log1p(-x)) - fast_math::log(a) - log_beta(a, b)?,
    );
    Ok(numerator * 1.0 / fraction)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn log_beta_of_one_and_one_is_negative_zero() {
        let value = log_beta(1.0, 1.0).expect("in range");
        assert_eq!(value, 0.0, "which is what a careless test would assert");
        assert!(value.is_sign_negative(), "and this is what it would miss");
    }

    #[test]
    fn the_four_branches_answer_what_the_reference_answers() {
        // a < 1, b < 10: the gamma-product branch.
        assert_eq!(log_beta(0.5, 0.5).expect("in range"), 1.1447298858494004);
        assert_eq!(log_beta(0.01, 0.01).expect("in range"), 5.298155239985667);
        // a < 1, b >= 10.
        assert_eq!(log_beta(0.1, 100.0).expect("in range"), 1.7926462324527925);
        // 1 <= a <= 2.
        assert_eq!(
            log_beta(1.0, 2.0).expect("in range"),
            -std::f64::consts::LN_2
        );
        // a > 2, b < 10.
        assert_eq!(log_beta(10.0, 2.0).expect("in range"), -4.7004803657924175);
        // a >= 10.
        assert_eq!(log_beta(10.0, 10.0).expect("in range"), -13.736229227036555);
        assert_eq!(
            log_beta(100.0, 100.0).expect("in range"),
            -139.66525908670664
        );
    }

    #[test]
    fn a_non_positive_argument_is_not_a_number_rather_than_a_refusal() {
        assert!(log_beta(0.0, 1.0).expect("answered").is_nan());
        assert!(log_beta(-1.0, 1.0).expect("answered").is_nan());
        assert!(log_beta(f64::NAN, 1.0).expect("answered").is_nan());
    }

    #[test]
    fn the_helpers_check_their_own_ranges() {
        assert!(matches!(
            log_gamma_sum(0.5, 1.0),
            Err(BetaError::OutOfRange { .. })
        ));
        assert!(matches!(
            delta_minus_delta_sum(1.0, 9.0),
            Err(BetaError::TooSmall { .. })
        ));
        assert!(matches!(
            sum_delta_minus_delta_sum(9.0, 10.0),
            Err(BetaError::TooSmall { .. })
        ));
        // And `a > b` is out of range even when both are large.
        assert!(matches!(
            delta_minus_delta_sum(20.0, 10.0),
            Err(BetaError::OutOfRange { .. })
        ));
    }

    /// Every value here is the reference's, from gatk-rs's `normal-artifact-filter` dump.
    #[test]
    fn the_two_branches_answer_what_the_reference_answers() {
        // The direct branch: `x` at or below `(a + 1) / (2 + b + a)`.
        assert_eq!(
            regularized_beta(0.01, 2.0, 5.0).unwrap(),
            0.001460447605000002
        );
        assert_eq!(regularized_beta(0.2, 2.0, 5.0).unwrap(), 0.3446400000000001);
        // The symmetry branch, which recurses on `1 - x` with `a` and `b` swapped.
        assert_eq!(
            regularized_beta(0.5, 2.0, 2.0).unwrap(),
            0.49999999999999994
        );
        assert_eq!(regularized_beta(0.99, 2.0, 5.0).unwrap(), 0.999999999405);
        assert_eq!(regularized_beta(0.4, 2.0, 5.0).unwrap(), 0.7667200000000001);
        // The shapes a binomial CDF asks for.
        assert_eq!(
            regularized_beta(0.001, 1.0, 100.0).unwrap(),
            0.09520785288629108
        );
        assert_eq!(
            regularized_beta(0.001, 11.0, 90.0).unwrap(),
            1.3053208102366178E-19
        );
        assert_eq!(
            regularized_beta(0.25, 4.0, 1.0).unwrap(),
            0.003906250000000001
        );
        assert_eq!(
            regularized_beta(0.25, 1.0, 4.0).unwrap(),
            0.6835937500000003
        );
        assert_eq!(
            regularized_beta(0.5, 500.0, 500.0).unwrap(),
            0.5000000000000142
        );
        assert_eq!(
            regularized_beta(0.5, 0.001, 0.001).unwrap(),
            0.49999999999999944
        );
        // The endpoints, where the logarithms go to infinity and the guards do not catch them.
        assert_eq!(regularized_beta(0.0, 2.0, 5.0).unwrap(), 0.0);
        assert_eq!(regularized_beta(1.0, 2.0, 5.0).unwrap(), 1.0);
    }

    /// One ULP either side of the branch boundary, which is where the redundant conjunct lives.
    #[test]
    fn the_boundary_is_the_branch() {
        let boundary = 3.0 / 9.0;
        // At the boundary itself `x > (a + 1) / (2 + b + a)` is false: the direct branch.
        assert_eq!(
            regularized_beta(boundary, 2.0, 5.0).unwrap(),
            0.6488340192043897
        );
        // One ULP below, still direct, and one ULP above, which reflects -- and the two answer the
        // same double, which is not the one the boundary answers.
        assert_eq!(
            regularized_beta(f64::from_bits(boundary.to_bits() - 1), 2.0, 5.0).unwrap(),
            0.6488340192043895
        );
        assert_eq!(
            regularized_beta(f64::from_bits(boundary.to_bits() + 1), 2.0, 5.0).unwrap(),
            0.6488340192043895
        );
    }

    /// Every guard answers NaN rather than refusing, as `logBeta`'s do.
    #[test]
    fn the_guards_answer_not_a_number() {
        for (x, a, b) in [
            (-0.5, 2.0, 5.0),
            (1.5, 2.0, 5.0),
            (0.5, 0.0, 5.0),
            (0.5, 2.0, -1.0),
            (f64::NAN, 2.0, 5.0),
        ] {
            assert!(regularized_beta(x, a, b).unwrap().is_nan(), "{x}, {a}, {b}");
        }
    }

    /// `BinomialDistribution.cumulativeProbability`, which is the one caller that matters here.
    #[test]
    fn the_binomial_cumulative_probability_is_one_minus_this() {
        // 100 trials at Q30, nine or fewer successes.
        let p = 10.0f64.powf(-3.0);
        assert_eq!(1.0 - regularized_beta(p, 10.0, 91.0).unwrap(), 1.0);
        // Ten trials, four or fewer.
        assert_eq!(
            1.0 - regularized_beta(p, 5.0, 6.0).unwrap(),
            0.9999999999997491
        );
        // Twenty fair trials, nine or fewer.
        assert_eq!(
            1.0 - regularized_beta(0.5, 10.0, 11.0).unwrap(),
            0.4119014739990241
        );
    }
}
