//! `CombinatoricsUtils.binomialCoefficientLog`, ported from commons-math 3.5.
//!
//! The other half of a beta-binomial likelihood. Three routes to the same number, chosen by `n`:
//!
//! ```java
//! if (n < 67)   return FastMath.log(binomialCoefficient(n, k));        // exact long
//! if (n < 1030) return FastMath.log(binomialCoefficientDouble(n, k));  // rounded product
//! // and past that, a sum of logs
//! ```
//!
//! # The routes do not agree with each other
//!
//! A port that took any single one of them would be right for some `n` and wrong for the rest. The
//! exact route computes an integer and logs it; the double route multiplies fractions and rounds
//! the product with `floor(x + 0.5)` before logging; the third never forms the coefficient at all.
//! `binomialCoefficientLog(1000, 500)` is `689.4672615678512` through the double route and a
//! different last digit through the sum of logs.
//!
//! # What this port refuses
//!
//! Two arms of `binomialCoefficient` beyond `n <= 61`: a gcd-splitting one for `n <= 66` and an
//! overflow-checking one past that. Nothing measured reaches either. They are refused rather than
//! guessed, on the same rule as `gamma` past 20.
//!
//! # The sum-of-logs route is measured now
//!
//! It was refused with them until gatk-rs's `contamination-filter` golden reached it: a site of
//! depth 1050 goes through `BetaBinomialDistribution.logProbability`, which asks for
//! `binomialCoefficientLog(1050, 200)`. The route never forms the coefficient at all, and its
//! `k > n / 2` reflection comes **before** the loops, so a large `k` is answered by a shorter sum
//! rather than a longer one.

use crate::fast_math;

/// What the coefficients refuse.
#[derive(Debug, Clone, PartialEq)]
pub enum CombinatoricsError {
    /// `NumberIsTooLargeException`: `k` past `n`.
    TooLarge { n: i64, k: i64 },
    /// `NotPositiveException`: a negative `n`.
    Negative { n: i64 },
    /// A total the reference could compute through an arm this port has not measured.
    Unmeasured { n: i64 },
}

/// `checkBinomial(n, k)`.
fn check_binomial(n: i64, k: i64) -> Result<(), CombinatoricsError> {
    if n < k {
        return Err(CombinatoricsError::TooLarge { n, k });
    }
    if n < 0 {
        return Err(CombinatoricsError::Negative { n });
    }
    Ok(())
}

/// `binomialCoefficient(n, k)`, the exact `long` route, for `n <= 61`.
pub fn binomial_coefficient(n: i64, k: i64) -> Result<i64, CombinatoricsError> {
    check_binomial(n, k)?;
    if n == k || k == 0 {
        return Ok(1);
    }
    if k == 1 || k == n - 1 {
        return Ok(n);
    }
    if k > n / 2 {
        return binomial_coefficient(n, n - k);
    }
    if n > 61 {
        return Err(CombinatoricsError::Unmeasured { n });
    }
    // `(n choose k) == ((n - k + 1) * ... * n) / (1 * ... * k)`, in an order where every partial
    // product is already an integer, so the truncating division never loses anything.
    let mut result: i64 = 1;
    for j in 1..=k {
        result = result * (n - k + j) / j;
    }
    Ok(result)
}

/// `binomialCoefficientDouble(n, k)`, the rounded-product route.
pub fn binomial_coefficient_double(n: i64, k: i64) -> Result<f64, CombinatoricsError> {
    check_binomial(n, k)?;
    if n == k || k == 0 {
        return Ok(1.0);
    }
    if k == 1 || k == n - 1 {
        return Ok(n as f64);
    }
    if k > n / 2 {
        return binomial_coefficient_double(n, n - k);
    }
    if n < 67 {
        return Ok(binomial_coefficient(n, k)? as f64);
    }
    let mut result = 1.0f64;
    for i in 1..=k {
        result *= (n - k + i) as f64 / i as f64;
    }
    // The product is an integer up to rounding, and this is the reference's way of saying so.
    Ok((result + 0.5).floor())
}

/// `binomialCoefficientLog(n, k)`, all three routes.
pub fn binomial_coefficient_log(n: i64, k: i64) -> Result<f64, CombinatoricsError> {
    check_binomial(n, k)?;
    if n == k || k == 0 {
        return Ok(0.0);
    }
    if k == 1 || k == n - 1 {
        return Ok(fast_math::log(n as f64));
    }
    if n < 67 {
        return Ok(fast_math::log(binomial_coefficient(n, k)? as f64));
    }
    if n < 1030 {
        return Ok(fast_math::log(binomial_coefficient_double(n, k)?));
    }
    // Past 1030 the coefficient itself would overflow a double, so it is never formed: the ratio
    // is summed in log space instead. The reflection comes FIRST, so `k > n / 2` is answered by
    // the same route with a smaller `k` rather than by a longer loop.
    if k > n / 2 {
        return binomial_coefficient_log(n, n - k);
    }
    let mut log_sum = 0.0;
    // `n! / (n - k)!`
    for i in (n - k + 1)..=n {
        log_sum += fast_math::log(i as f64);
    }
    // divided by `k!`
    for i in 2..=k {
        log_sum -= fast_math::log(i as f64);
    }
    Ok(log_sum)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_three_routes_answer_what_the_reference_answers() {
        assert_eq!(binomial_coefficient_log(0, 0).expect("in range"), 0.0);
        assert_eq!(binomial_coefficient_log(1, 1).expect("in range"), 0.0);
        assert_eq!(binomial_coefficient_log(10, 0).expect("in range"), 0.0);
        // `k == 1`, which is `log(n)` rather than any coefficient.
        assert_eq!(
            binomial_coefficient_log(10, 1).expect("in range"),
            std::f64::consts::LN_10
        );
        assert_eq!(
            binomial_coefficient_log(1000, 1).expect("in range"),
            6.907755278982137
        );
        // The exact route.
        assert_eq!(
            binomial_coefficient_log(10, 5).expect("in range"),
            5.529429087511423
        );
        // The rounded-product route, twice.
        assert_eq!(
            binomial_coefficient_log(100, 50).expect("in range"),
            66.78384165201743
        );
        assert_eq!(
            binomial_coefficient_log(1000, 500).expect("in range"),
            689.4672615678512
        );
    }

    #[test]
    fn the_exact_route_is_exact() {
        assert_eq!(binomial_coefficient(10, 5).expect("in range"), 252);
        // Where the naive product would overflow if it multiplied before dividing.
        assert_eq!(
            binomial_coefficient(61, 30).expect("in range"),
            232714176627630544
        );
    }

    #[test]
    fn what_it_refuses() {
        assert_eq!(
            binomial_coefficient_log(10, 11),
            Err(CombinatoricsError::TooLarge { n: 10, k: 11 })
        );
        assert_eq!(
            binomial_coefficient(-1, -1),
            Err(CombinatoricsError::Negative { n: -1 })
        );
        // The two unported arms of the exact route. The sum-of-logs route used to be refused
        // beside them and is measured now.
        assert_eq!(
            binomial_coefficient(62, 31),
            Err(CombinatoricsError::Unmeasured { n: 62 })
        );
        assert_eq!(
            binomial_coefficient_double(62, 31),
            Err(CombinatoricsError::Unmeasured { n: 62 })
        );
    }

    /// The sum-of-logs route, past 1030, where the coefficient itself would overflow a double.
    ///
    /// Every value is the reference's, from the pinned oracle container; gatk-rs's
    /// `contamination-filter` golden is what exercises the route end to end, through a site of
    /// depth 1050.
    #[test]
    fn the_sum_of_logs_route_answers_what_the_reference_answers() {
        assert_eq!(
            binomial_coefficient_log(1050, 200).expect("in range"),
            507.79546692773414
        );
        assert_eq!(
            binomial_coefficient_log(1050, 50).expect("in range"),
            198.16406048906853
        );
        // At the boundary itself: 1030 is the first `n` the route takes.
        assert_eq!(
            binomial_coefficient_log(1030, 2).expect("in range"),
            13.180509636497787
        );
        // `k > n / 2`, which reflects before the loops rather than summing more terms.
        assert_eq!(
            binomial_coefficient_log(1030, 515).expect("in range"),
            710.2469048650747
        );
        assert_eq!(
            binomial_coefficient_log(2000, 1900).expect("in range"),
            393.83377418920224
        );
        assert_eq!(
            binomial_coefficient_log(5000, 2500).expect("in range"),
            3461.2514648514098
        );
    }

    /// One below the boundary is the double route, and `k == n - 1` never reaches either.
    #[test]
    fn the_boundary_is_the_route() {
        assert_eq!(
            binomial_coefficient_log(1029, 200).expect("in range"),
            503.3102211260805
        );
        // `k == n - 1` returns `log(n)` before any route is chosen.
        assert_eq!(
            binomial_coefficient_log(1030, 1029).expect("in range"),
            6.937314081223682
        );
        // And `n == k` is zero, however large.
        assert_eq!(binomial_coefficient_log(1050, 1050).expect("in range"), 0.0);
    }
}
