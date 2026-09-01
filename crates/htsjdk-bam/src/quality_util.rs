//! `htsjdk.samtools.util.QualityUtil`.
//!
//! Two conversions and a sum, and all three are arithmetic before they are utilities: which
//! `log10`, which `round`, and in which order. `picard-rs` ported the same functions three times,
//! in three tools, each with `f64::log10` and `f64::round`, and neither is the function Java calls:
//!
//! * `Math.log10` is correctly rounded (decision 0006), so [`jmath::math::log10`] is the port and
//!   `f64::log10` is a different function that usually agrees.
//! * `Math.round` is **not** `floor(x + 0.5)` and Rust's `f64::round` is "half away from zero".
//!   They disagree on `0.49999999999999994` and on every negative half, and `getPhredScoreFrom
//!   ErrorProbability` takes a negative argument whenever the probability is above 1.
//!
//! The error-probability table is built by the class's static initialiser and is a table, not a
//! function: `getErrorProbabilityFromPhredScore` indexes it, so a score outside `0..=100` throws
//! rather than extrapolating, and the port returns `None` where Java throws
//! `ArrayIndexOutOfBoundsException`.
//!
//! One value in that table is **not** reproduced by construction. The initialiser calls
//! `Math.pow`, which decision 0007 deferred and decision 0027 bounded at 1 ulp, so the port
//! computes the table with [`jmath::strict_math::pow`] and the two can differ in the last bit. The
//! `quality-util` suite is what settles it: the reference's own 101 entries are dumped as bit
//! patterns, and the difference is measured rather than assumed. Until that golden is committed
//! the table here is `StrictMath`'s answer, which is stated rather than hidden.

/// `QualityUtil.errorProbabilityByPhredScore`, the static table of 101 entries.
///
/// `1d / Math.pow(10d, i / 10d)`, in that order: a reciprocal of a power rather than
/// `Math.pow(10d, -i / 10d)`, which is a different double for most `i`.
fn error_probability_table() -> [f64; 101] {
    let mut table = [0.0f64; 101];
    for (i, slot) in table.iter_mut().enumerate() {
        *slot = 1.0 / jmath::strict_math::pow(10.0, i as f64 / 10.0);
    }
    table
}

/// `QualityUtil.getErrorProbabilityFromPhredScore(i)`.
///
/// `None` is the Java `ArrayIndexOutOfBoundsException`: the method indexes a 101-entry table and
/// makes no promise for anything else.
pub fn error_probability_from_phred_score(score: i32) -> Option<f64> {
    if !(0..=100).contains(&score) {
        return None;
    }
    Some(error_probability_table()[score as usize])
}

/// `QualityUtil.getPhredScoreFromErrorProbability(probability)`.
///
/// `(int) Math.round(-10 * Math.log10(probability))`, with Java's `round` and Java's `log10`, and
/// with Java's narrowing cast: `Math.round` answers a `long`, and `(int)` of a long keeps its low
/// 32 bits rather than saturating. A probability of zero therefore does not answer `Integer.MAX_
/// VALUE`; `-10 * log10(0)` is `Infinity`, `Math.round(Infinity)` is `Long.MAX_VALUE`, and `(int)
/// Long.MAX_VALUE` is `-1`.
pub fn phred_score_from_error_probability(probability: f64) -> i32 {
    let rounded = jmath::math::round(-10.0 * jmath::math::log10(probability));
    rounded as i32
}

/// `QualityUtil.getPhredScoreFromObsAndErrors(observations, errors)`.
pub fn phred_score_from_obs_and_errors(observations: f64, errors: f64) -> i32 {
    phred_score_from_error_probability(errors / observations)
}

/// `QualityUtil.sumOfErrorProbabilities(bases, quals)`.
///
/// A no-call base contributes a whole error rather than its quality's probability, which is the
/// one line of this function that is not arithmetic.
pub fn sum_of_error_probabilities(bases: &[u8], quals: &[u8]) -> f64 {
    let table = error_probability_table();
    let mut sum = 0.0;
    for (i, &base) in bases.iter().enumerate() {
        if crate::sequence::is_no_call(base) {
            sum += 1.0;
        } else {
            sum += table[quals[i] as usize];
        }
    }
    sum
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_table_is_a_reciprocal_of_a_power_not_a_negative_power() {
        // 1/pow(10, 0.7) and pow(10, -0.7) are different doubles, and the class computes the first.
        let table = error_probability_table();
        assert_eq!(table[7], 1.0 / jmath::strict_math::pow(10.0, 0.7));
        assert_eq!(table[0], 1.0);
    }

    #[test]
    fn a_score_outside_the_table_is_not_extrapolated() {
        assert!(error_probability_from_phred_score(-1).is_none());
        assert!(error_probability_from_phred_score(101).is_none());
        assert!(error_probability_from_phred_score(100).is_some());
    }

    #[test]
    fn the_round_is_javas_and_not_rusts() {
        // Math.round is half UP, toward positive infinity; f64::round is half AWAY FROM ZERO. They
        // agree on every positive half and disagree on every negative one, and this function takes
        // a negative argument whenever the probability is above 1, which an observed error rate
        // above one produces.
        assert_eq!(jmath::math::round(-2.5), -2);
        assert_eq!((-2.5f64).round() as i64, -3);
        // And the javadoc's floor(x + 0.5), which neither of them is, on the double below a half.
        let below_half = 0.499_999_999_999_999_94_f64;
        assert_eq!(jmath::math::round(below_half), 0);
        assert_eq!((below_half + 0.5).floor() as i64, 1);
    }

    #[test]
    fn a_zero_probability_wraps_rather_than_saturating() {
        // Math.round(Infinity) is Long.MAX_VALUE and (int) of it is -1. A port that saturated to
        // i32::MAX would write a different number into a metrics file.
        assert_eq!(phred_score_from_error_probability(0.0), -1);
    }

    #[test]
    fn a_no_call_costs_a_whole_error() {
        let sum = sum_of_error_probabilities(b"AN", &[10, 10]);
        assert!((sum - (0.1 + 1.0)).abs() < 1e-12, "{sum}");
    }
}
