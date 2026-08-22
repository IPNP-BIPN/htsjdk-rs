//! `PoissonDistribution`, ported from
//! `org.apache.commons.math3.distribution.PoissonDistribution` (commons-math 3.5).
//!
//! GATK's `NuMTFilterTool` asks it one question: at what alternate depth does a call stop being
//! explicable as a nuclear insertion of mitochondrial DNA? The answer is a quantile,
//! `new PoissonDistribution(coverage * copies / 2).inverseCumulativeProbability(0.99)`, and it is
//! the only arithmetic in that tool.
//!
//! # The CDF is the regularized gamma, and its epsilon is not the gamma's default
//!
//! ```java
//! public double cumulativeProbability(int x) {
//!     if (x < 0) { return 0; }
//!     if (x == Integer.MAX_VALUE) { return 1; }
//!     return Gamma.regularizedGammaQ((double) x + 1, mean, epsilon, maxIterations);
//! }
//! ```
//!
//! `epsilon` is the distribution's own `DEFAULT_EPSILON`, `1e-12`, and `maxIterations` is
//! `10_000_000`. `Gamma`'s own defaults are `1e-14` and `Integer.MAX_VALUE`, so passing those
//! instead is a plausible port and a different function.
//!
//! # The quantile is a bracket and a bisection, not a formula
//!
//! `AbstractIntegerDistribution.inverseCumulativeProbability` narrows a bracket with a one-sided
//! Chebyshev inequality and then bisects it, exactly as [`crate::binomial`] does for its own
//! distribution. Three details decide the answer:
//!
//!  * **the lower bound is decremented before the search**, so `cdf(lower) < p <= cdf(upper)` holds
//!    from the first iteration;
//!  * **the narrowing takes `ceil` and then subtracts one** at both ends;
//!  * **and the bisection answers the UPPER end**, so the result is the smallest `x` whose CDF
//!    reaches `p` rather than the largest that falls short.
//!
//! For a Poisson the mean and the variance are the same number, so the bracket is
//! `mu -/+ k * sqrt(mu)`, and the variance is zero only for a mean of zero, which the constructor
//! refuses.

use crate::gamma::{regularized_gamma_q, GammaError};

/// `PoissonDistribution.DEFAULT_EPSILON`. Not `Gamma`'s `1e-14`.
pub const DEFAULT_EPSILON: f64 = 1e-12;

/// `PoissonDistribution.DEFAULT_MAX_ITERATIONS`. Not `Gamma`'s `Integer.MAX_VALUE`.
pub const DEFAULT_MAX_ITERATIONS: i32 = 10_000_000;

/// What the distribution refuses.
#[derive(Debug, Clone, PartialEq)]
pub enum PoissonError {
    /// `NotStrictlyPositiveException(MEAN, p)`.
    MeanNotStrictlyPositive(f64),
    /// `OutOfRangeException(p, 0, 1)`.
    ProbabilityOutOfRange(f64),
    /// The regularized gamma underneath refused.
    Gamma(GammaError),
}

impl From<GammaError> for PoissonError {
    fn from(error: GammaError) -> Self {
        PoissonError::Gamma(error)
    }
}

/// The constructor's one check.
fn check_mean(mean: f64) -> Result<(), PoissonError> {
    if !(mean > 0.0) {
        return Err(PoissonError::MeanNotStrictlyPositive(mean));
    }
    Ok(())
}

/// `cumulativeProbability(x)`.
pub fn cumulative_probability(mean: f64, x: i32) -> Result<f64, PoissonError> {
    check_mean(mean)?;
    if x < 0 {
        return Ok(0.0);
    }
    if x == i32::MAX {
        return Ok(1.0);
    }
    Ok(regularized_gamma_q(
        f64::from(x) + 1.0,
        mean,
        DEFAULT_EPSILON,
        DEFAULT_MAX_ITERATIONS,
    )?)
}

/// `inverseCumulativeProbability(p)`: the smallest `x` whose CDF reaches `p`.
///
/// A `p` outside the unit interval is refused. A NaN is NOT, for the same reason as in
/// [`crate::binomial`]: `NaN < 0` and `NaN > 1` are both false, the equality tests fail, the
/// narrowing produces NaNs that narrow nothing, and the bisection walks the whole support.
pub fn inverse_cumulative_probability(mean: f64, p: f64) -> Result<i32, PoissonError> {
    check_mean(mean)?;
    // Two comparisons rather than a `RangeInclusive::contains`, because a NaN has to pass BOTH of
    // them: `contains` would refuse it, and the reference does not.
    #[allow(clippy::manual_range_contains)]
    if p < 0.0 || p > 1.0 {
        return Err(PoissonError::ProbabilityOutOfRange(p));
    }

    // `getSupportLowerBound()` is 0 and `getSupportUpperBound()` is `Integer.MAX_VALUE`. Neither
    // is `Integer.MIN_VALUE`, so the branch that checks the CDF at the lower bound is unreachable
    // and the decrement always happens.
    if p == 0.0 {
        return Ok(0);
    }
    let mut lower: i32 = -1;
    let mut upper: i32 = i32::MAX;
    if p == 1.0 {
        return Ok(upper);
    }

    // For a Poisson the variance IS the mean, so sigma is its square root and is zero only for a
    // mean of zero, which the constructor already refused.
    let sigma = mean.sqrt();
    let chebyshev_applies = !(mean.is_infinite()
        || mean.is_nan()
        || sigma.is_infinite()
        || sigma.is_nan()
        || sigma == 0.0);
    if chebyshev_applies {
        let mut k = ((1.0 - p) / p).sqrt();
        let mut tmp = mean - k * sigma;
        if tmp > f64::from(lower) {
            lower = (tmp.ceil() as i32) - 1;
        }
        k = 1.0 / k;
        tmp = mean + k * sigma;
        if tmp < f64::from(upper) {
            upper = (tmp.ceil() as i32) - 1;
        }
    }

    // `solveInverseCumulativeProbability`: bisect, and answer the UPPER end.
    //
    // The midpoint is a WRAPPING add. Java's `(lower + upper) / 2` overflows silently when the
    // bracket still reaches `Integer.MAX_VALUE`, which is why the reference re-computes it as
    // `lower + (upper - lower) / 2` whenever the result falls outside the bracket. That guard is
    // reachable: a NaN quantile narrows nothing, so the first midpoint is taken over the whole
    // support. A port that adds without wrapping panics there in debug and answers the same number
    // in release, which is the worst of both.
    while lower + 1 < upper {
        let mut middle = lower.wrapping_add(upper) / 2;
        if middle < lower || middle > upper {
            middle = lower + (upper - lower) / 2;
        }
        if cumulative_probability(mean, middle)? >= p {
            upper = middle;
        } else {
            lower = middle;
        }
    }
    Ok(upper)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every value here is the reference's, from gatk-rs's `numt-filter` measurement: twenty-one
    /// cutoffs over seven autosomal coverages and three copy counts, taken in the pinned container
    /// on a real x86-64 runner. `NuMTFilterTool` asks for the mean `coverage * copies / 2` and the
    /// quantile `1 - 0.01`.
    #[test]
    fn the_cutoffs_are_the_references() {
        for (mean, cutoff) in [
            (0.125, 1),
            (0.25, 2),
            (0.5, 3),
            (1.0, 4),
            (1.25, 4),
            (2.0, 6),
            (2.5, 7),
            (5.0, 11),
            (7.5, 15),
            (10.0, 18),
            (15.0, 25),
            (20.0, 31),
            (25.0, 37),
            (50.0, 67),
            (60.0, 79),
            (200.0, 234),
            (250.0, 288),
            (500.0, 553),
            (2000.0, 2105),
        ] {
            assert_eq!(
                inverse_cumulative_probability(mean, 1.0 - 0.01).unwrap(),
                cutoff,
                "mean {mean}"
            );
        }
    }

    /// The answer is the smallest x whose CDF reaches p, so the cutoff's own CDF is at or above
    /// the quantile and its predecessor's is below it.
    #[test]
    fn the_bisection_answers_the_upper_end() {
        let p = 1.0 - 0.01;
        assert!(cumulative_probability(60.0, 79).unwrap() >= p);
        assert!(cumulative_probability(60.0, 78).unwrap() < p);
    }

    /// The two guards of the CDF, neither of which is a refusal.
    #[test]
    fn the_cdf_has_two_guards() {
        assert_eq!(cumulative_probability(60.0, -1).unwrap(), 0.0);
        assert_eq!(cumulative_probability(60.0, i32::MAX).unwrap(), 1.0);
    }

    /// The degenerate quantiles answer the support's bounds without any search.
    #[test]
    fn the_degenerate_quantiles_answer_the_bounds() {
        assert_eq!(inverse_cumulative_probability(60.0, 0.0).unwrap(), 0);
        assert_eq!(inverse_cumulative_probability(60.0, 1.0).unwrap(), i32::MAX);
    }

    #[test]
    fn a_mean_that_is_not_strictly_positive_is_refused() {
        assert_eq!(
            cumulative_probability(0.0, 1),
            Err(PoissonError::MeanNotStrictlyPositive(0.0))
        );
        assert_eq!(
            inverse_cumulative_probability(-1.0, 0.5),
            Err(PoissonError::MeanNotStrictlyPositive(-1.0))
        );
    }

    #[test]
    fn a_probability_outside_the_unit_interval_is_refused() {
        assert_eq!(
            inverse_cumulative_probability(60.0, 1.5),
            Err(PoissonError::ProbabilityOutOfRange(1.5))
        );
        assert_eq!(
            inverse_cumulative_probability(60.0, -0.5),
            Err(PoissonError::ProbabilityOutOfRange(-0.5))
        );
    }

    /// A NaN is not refused: it passes both comparisons, narrows nothing, and walks the support.
    ///
    /// It is also the only case that reaches the overflow guard, since the bracket still spans to
    /// `Integer.MAX_VALUE` when the first midpoint is taken.
    #[test]
    fn a_nan_quantile_is_not_refused() {
        assert_eq!(
            inverse_cumulative_probability(60.0, f64::NAN).unwrap(),
            i32::MAX
        );
    }
}
