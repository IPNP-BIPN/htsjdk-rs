//! `ContinuedFraction.evaluate`, ported from `org.apache.commons.math3.util.ContinuedFraction`
//! (commons-math 3.5).
//!
//! The modified Lentz algorithm, which is how commons-math evaluates the regularized beta and the
//! regularized gamma. GATK reaches it through `BinomialDistribution.cumulativeProbability`, which
//! Mutect's normal-artifact filter asks for a p-value.
//!
//! # The caller supplies the epsilon, and it is not this module's own
//!
//! [`DEFAULT_EPSILON`] here is `10e-9`, and nothing in the ported call sites uses it:
//! `Beta.regularizedBeta` passes its own `1E-14`. A port that reached for the algorithm's default
//! would stop the loop five digits early. The constant is here because the reference has it, not
//! because anything should call [`evaluate`] without an epsilon.
//!
//! # The zero test is a tolerance, not an equality
//!
//! `Precision.equals(v, 0.0, small)` with `small = 1e-50` is
//! `equals(v, 0.0, 1) || abs(0.0 - v) <= small`: within one ULP of zero, or within `1e-50` of it.
//! Every value the first arm accepts (`0.0`, `-0.0`, `±Double.MIN_VALUE`) the second accepts too,
//! and both reject NaN, so the whole test is `v.abs() <= 1e-50` -- which is emphatically not
//! `v == 0.0`. A convergent that merely underflows towards zero is replaced by `1e-50` rather than
//! dividing.
//!
//! # An infinite or NaN convergent is refused
//!
//! Not returned. The reference throws `ConvergenceException`, and the two messages differ, so the
//! refusal names which one it was.

/// `ContinuedFraction.DEFAULT_EPSILON`. See the module note: no ported call site uses it.
pub const DEFAULT_EPSILON: f64 = 10e-9;

/// `small`, the value that stands in for a convergent too close to zero, and the tolerance of the
/// test that decides "too close".
const SMALL: f64 = 1e-50;

/// Which of the two `ConvergenceException`s the reference throws.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Divergence {
    /// `CONTINUED_FRACTION_INFINITY_DIVERGENCE`.
    Infinite,
    /// `CONTINUED_FRACTION_NAN_DIVERGENCE`.
    NotANumber,
}

/// What [`evaluate`] refuses, each of them an exception the reference throws.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ContinuedFractionError {
    /// `ConvergenceException`: the convergent left the finite numbers.
    Diverged { x: f64, kind: Divergence },
    /// `MaxCountExceededException`: the loop ran out of iterations before `|delta - 1| < epsilon`.
    NotConvergent { max_iterations: i32, x: f64 },
}

/// `Precision.equals(value, 0.0, SMALL)`, which is a tolerance and rejects NaN.
fn is_zero(value: f64) -> bool {
    value.abs() <= SMALL
}

/// `evaluate(x, epsilon, maxIterations)`, with the abstract `getA` and `getB` passed in.
///
/// `get_a(0, x)` is the seed, and the loop runs from `n = 1`.
pub fn evaluate(
    get_a: impl Fn(i32, f64) -> f64,
    get_b: impl Fn(i32, f64) -> f64,
    x: f64,
    epsilon: f64,
    max_iterations: i32,
) -> Result<f64, ContinuedFractionError> {
    let mut h_prev = get_a(0, x);

    // The reference's comment: "use the value of small as epsilon criteria for zero checks".
    if is_zero(h_prev) {
        h_prev = SMALL;
    }

    let mut n = 1;
    let mut d_prev = 0.0;
    let mut c_prev = h_prev;
    let mut h_n = h_prev;

    while n < max_iterations {
        let a = get_a(n, x);
        let b = get_b(n, x);

        let mut d_n = a + b * d_prev;
        if is_zero(d_n) {
            d_n = SMALL;
        }
        let mut c_n = a + b / c_prev;
        if is_zero(c_n) {
            c_n = SMALL;
        }

        d_n = 1.0 / d_n;
        let delta_n = c_n * d_n;
        h_n = h_prev * delta_n;

        if h_n.is_infinite() {
            return Err(ContinuedFractionError::Diverged {
                x,
                kind: Divergence::Infinite,
            });
        }
        if h_n.is_nan() {
            return Err(ContinuedFractionError::Diverged {
                x,
                kind: Divergence::NotANumber,
            });
        }

        // `<`, not `<=`: a delta exactly one epsilon away from one runs another iteration.
        if (delta_n - 1.0).abs() < epsilon {
            break;
        }

        d_prev = d_n;
        c_prev = c_n;
        h_prev = h_n;
        n += 1;
    }

    if n >= max_iterations {
        return Err(ContinuedFractionError::NotConvergent { max_iterations, x });
    }

    Ok(h_n)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The golden ratio, `1 + 1/(1 + 1/(1 + ...))`, whose `a` and `b` are both one.
    #[test]
    fn the_all_ones_fraction_is_the_golden_ratio() {
        let value = evaluate(|_, _| 1.0, |_, _| 1.0, 0.0, 1e-14, i32::MAX).expect("converges");
        assert!((value - (1.0 + 5.0f64.sqrt()) / 2.0).abs() < 1e-13);
    }

    /// The loop's exit is `n >= max_iterations` **after** the `break`, so a fraction that has
    /// converged at the last permitted iteration is still a refusal if it never broke out.
    #[test]
    fn running_out_of_iterations_is_a_refusal() {
        assert_eq!(
            evaluate(|_, _| 1.0, |_, _| 1.0, 0.0, 0.0, 5),
            Err(ContinuedFractionError::NotConvergent {
                max_iterations: 5,
                x: 0.0
            })
        );
    }

    /// A convergent that overflows is refused, and the overflow check comes **before** the
    /// convergence check: a seed at `Double.MAX_VALUE` against a delta a hair above one leaves the
    /// finite numbers on the first iteration, and the epsilon is never consulted.
    #[test]
    fn an_infinite_convergent_is_refused() {
        let value = evaluate(
            |n, _| if n == 0 { f64::MAX } else { 1.0 },
            |_, _| 1e300,
            0.5,
            1e-14,
            i32::MAX,
        );
        assert_eq!(
            value,
            Err(ContinuedFractionError::Diverged {
                x: 0.5,
                kind: Divergence::Infinite
            })
        );
    }

    /// And a NaN convergent is a different refusal, with a different message in the reference.
    #[test]
    fn a_nan_convergent_is_refused() {
        assert_eq!(
            evaluate(|_, _| 1.0, |_, _| f64::NAN, 0.5, 1e-14, i32::MAX),
            Err(ContinuedFractionError::Diverged {
                x: 0.5,
                kind: Divergence::NotANumber
            })
        );
    }

    /// Zero is not tested with `==`: a seed of `1e-60` is "zero" and is replaced by `1e-50`.
    #[test]
    fn the_zero_test_is_a_tolerance() {
        assert!(is_zero(1e-60));
        assert!(is_zero(-0.0));
        assert!(!is_zero(1e-40));
        assert!(!is_zero(f64::NAN));
    }
}
