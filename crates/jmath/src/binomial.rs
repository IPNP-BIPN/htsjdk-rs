//! `BinomialDistribution.cumulativeProbability`, ported from
//! `org.apache.commons.math3.distribution.BinomialDistribution` (commons-math 3.5).
//!
//! Mutect's normal-artifact filter asks this for a p-value: "how surprising is it that a normal of
//! this depth carries this many alternate reads, if every one of them were a base-calling error?"
//!
//! # It is three arms, and the middle one is the only arithmetic
//!
//! ```java
//! if (x < 0) { ret = 0.0; }
//! else if (x >= numberOfTrials) { ret = 1.0; }
//! else { ret = 1.0 - Beta.regularizedBeta(probabilityOfSuccess, x + 1.0, numberOfTrials - x); }
//! ```
//!
//! Both guards matter to the caller. `x < 0` is reached whenever the normal has no alternate read
//! at all, because the filter asks for `cumulativeProbability(normalAltDepth - 1)`: the answer is
//! `0.0`, so the p-value is exactly `1` and no filtering follows. `x >= trials` is reached when
//! every read supports the allele: the answer is `1.0`, the p-value is `0`, below any threshold,
//! and the filter returns a hard `1.0`. Neither is a refusal.
//!
//! A distribution of zero trials answers `1.0` to every non-negative `x`, since `x >= 0` is
//! `x >= numberOfTrials`.

use crate::beta::{regularized_beta, BetaError};

/// What the constructor refuses, before any probability is asked for.
#[derive(Debug, Clone, PartialEq)]
pub enum BinomialError {
    /// `NotPositiveException(NUMBER_OF_TRIALS, trials)`.
    NegativeTrials(i32),
    /// `OutOfRangeException(p, 0, 1)`.
    ProbabilityOutOfRange(f64),
    /// The regularized beta underneath refused.
    Beta(BetaError),
}

impl From<BetaError> for BinomialError {
    fn from(error: BetaError) -> Self {
        BinomialError::Beta(error)
    }
}

/// `new BinomialDistribution(null, trials, p).cumulativeProbability(x)`.
///
/// The `null` is the reference's own: the filter passes no random generator, and the constructor
/// accepts it because nothing samples from the distribution.
// clippy offers `!(0.0..=1.0).contains(&p)` for the range check below. That is a different
// program: `contains` is false for NaN, so the negation refuses NaN, where `p < 0.0 || p > 1.0` is
// false for NaN and lets it through to the beta. The reference lets it through.
#[allow(clippy::manual_range_contains)]
pub fn cumulative_probability(trials: i32, p: f64, x: i32) -> Result<f64, BinomialError> {
    if trials < 0 {
        return Err(BinomialError::NegativeTrials(trials));
    }
    if p < 0.0 || p > 1.0 {
        return Err(BinomialError::ProbabilityOutOfRange(p));
    }
    if x < 0 {
        Ok(0.0)
    } else if x >= trials {
        Ok(1.0)
    } else {
        // `numberOfTrials - x` is int arithmetic widened on the way in, and `x + 1.0` is not.
        Ok(1.0 - regularized_beta(p, f64::from(x) + 1.0, f64::from(trials - x))?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every value here is the reference's, from gatk-rs's `normal-artifact-filter` golden.
    #[test]
    fn the_middle_arm_answers_what_the_reference_answers() {
        // A hundred reads at Q30.
        assert_eq!(
            cumulative_probability(100, 0.001, 0).unwrap(),
            0.9047921471137089
        );
        assert_eq!(cumulative_probability(100, 0.001, 9).unwrap(), 1.0);
        assert_eq!(cumulative_probability(100, 0.001, 99).unwrap(), 1.0);
        // Ten reads, and twenty fair ones.
        assert_eq!(
            cumulative_probability(10, 0.001, 4).unwrap(),
            0.9999999999997491
        );
        assert_eq!(
            cumulative_probability(20, 0.5, 9).unwrap(),
            0.4119014739990241
        );
        // A quality so low that nine successes out of a hundred is no longer surprising at all,
        // and one so high that it is certain.
        assert_eq!(
            cumulative_probability(100, 0.6309573444801932, 9).unwrap(),
            0.0
        );
        assert_eq!(cumulative_probability(100, 1.0E-6, 9).unwrap(), 1.0);
        // The degenerate probabilities, which the guards let through to the beta.
        assert_eq!(cumulative_probability(20, 0.0, 9).unwrap(), 1.0);
        assert_eq!(cumulative_probability(20, 1.0, 9).unwrap(), 0.0);
    }

    /// The two arms the filter reaches by arithmetic rather than by choice.
    #[test]
    fn the_guards_are_answers_and_not_refusals() {
        // `normalAltDepth - 1` when the normal is clean.
        assert_eq!(cumulative_probability(100, 0.001, -1).unwrap(), 0.0);
        // Every read supports the allele.
        assert_eq!(cumulative_probability(100, 0.001, 100).unwrap(), 1.0);
        // And a normal of no depth at all answers 1.0 to x = 0, since 0 >= 0.
        assert_eq!(cumulative_probability(0, 0.001, 0).unwrap(), 1.0);
    }

    #[test]
    fn the_constructor_checks_its_own_arguments() {
        assert_eq!(
            cumulative_probability(-1, 0.5, 0),
            Err(BinomialError::NegativeTrials(-1))
        );
        assert_eq!(
            cumulative_probability(10, 1.5, 0),
            Err(BinomialError::ProbabilityOutOfRange(1.5))
        );
        // A NaN probability is NOT out of range: `p < 0` and `p > 1` are both false for NaN, so it
        // passes the constructor untouched and reaches the beta, which answers NaN.
        assert!(cumulative_probability(10, f64::NAN, 5).unwrap().is_nan());
    }
}
