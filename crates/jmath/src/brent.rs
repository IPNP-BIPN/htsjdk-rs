//! `BrentOptimizer`, ported from `org.apache.commons.math3.optim.univariate.BrentOptimizer`
//! (commons-math 3.5), with the `SearchInterval` and evaluation budget around it.
//!
//! Brent's golden-section search with parabolic interpolation, over one variable. GATK reaches it
//! through `OptimizationUtils.max`, which Mutect's strand-artifact filter uses to re-estimate its
//! beta shape between passes.
//!
//! # Maximising is negating inside the search, not negating the objective
//!
//! ```java
//! double fx = computeObjectiveValue(x);
//! if (!isMinim) { fx = -fx; }
//! ...
//! current = new UnivariatePointValuePair(u, isMinim ? fu : -fu);
//! ```
//!
//! Every comparison in the search is on the negated value; the pair that comes back carries the
//! objective's own. A port that negated the objective once at the top would agree on the point and
//! on the value, and would still be a different program: the negation is applied to each evaluation
//! separately, and `-(-0.0)` is `0.0`.
//!
//! # It returns the best point ever seen, not the last one
//!
//! `best` is folded over the previous and current pairs at every iteration, starting from the
//! initial guess. So a search that wanders away from a maximum still reports the maximum it passed,
//! and a flat objective reports the guess.
//!
//! # It does not reach the bounds
//!
//! The stopping criterion is `|x - m| <= 2 * tol - 0.5 * (b - a)`, and the step is bounded below by
//! `tol1`. A strictly increasing objective over `[0.01, 100]` with both tolerances at `0.01` stops
//! at `98.88158923714171`: the bound is an interval end, never a candidate.
//!
//! # The tolerances and the budget are refusals
//!
//! A relative tolerance below `2 * ulp(1)` and an absolute tolerance of zero or less are refused at
//! construction rather than clamped, and the evaluation budget is refused when it is exceeded, not
//! when it is reached: the call that would be the `(max + 1)`-th throws.

/// `FastMath.abs(double)`, which is `Math.abs`: the sign bit cleared, so `abs(NaN)` is NaN.
fn java_abs(x: f64) -> f64 {
    x.abs()
}

/// `GOLDEN_SECTION`, `0.5 * (3 - FastMath.sqrt(5))`.
///
/// `sqrt` is const-unavailable in Rust, so the root is written out; a test asserts it is the same
/// double `5.0f64.sqrt()` produces, which IEEE-754 mandates to be correctly rounded.
const GOLDEN_SECTION: f64 = 0.5 * (3.0 - 2.23606797749979);

/// `MIN_RELATIVE_TOLERANCE`, `2 * ulp(1)`.
pub const MIN_RELATIVE_TOLERANCE: f64 = 2.0 * f64::EPSILON;

/// `GoalType`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GoalType {
    Minimize,
    Maximize,
}

/// `UnivariatePointValuePair`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PointValuePair {
    pub point: f64,
    pub value: f64,
}

/// What the optimiser refuses, each of them an exception the reference throws.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BrentError {
    /// `NumberIsTooSmallException(rel, MIN_RELATIVE_TOLERANCE, true)`.
    RelativeToleranceTooSmall { value: f64, minimum: f64 },
    /// `NotStrictlyPositiveException(abs)`.
    AbsoluteToleranceNotPositive { value: f64 },
    /// `NumberIsTooLargeException(lo, hi, false)` out of `SearchInterval`.
    IntervalNotIncreasing { lower: f64, upper: f64 },
    /// `OutOfRangeException(init, lo, hi)` out of `SearchInterval`.
    StartOutOfRange { value: f64, lower: f64, upper: f64 },
    /// `TooManyEvaluationsException(max)`.
    TooManyEvaluations { maximum: i32 },
}

impl BrentError {
    /// The exception class the reference throws.
    pub fn class(&self) -> &'static str {
        match self {
            BrentError::RelativeToleranceTooSmall { .. } => {
                "org.apache.commons.math3.exception.NumberIsTooSmallException"
            }
            BrentError::AbsoluteToleranceNotPositive { .. } => {
                "org.apache.commons.math3.exception.NotStrictlyPositiveException"
            }
            BrentError::IntervalNotIncreasing { .. } => {
                "org.apache.commons.math3.exception.NumberIsTooLargeException"
            }
            BrentError::StartOutOfRange { .. } => {
                "org.apache.commons.math3.exception.OutOfRangeException"
            }
            BrentError::TooManyEvaluations { .. } => {
                "org.apache.commons.math3.exception.TooManyEvaluationsException"
            }
        }
    }

    /// The message, with the numbers rendered as `MessageFormat` renders them.
    ///
    /// That is why a relative tolerance of `1e-17` is refused with "0 is smaller than the minimum
    /// (0)": both numbers are below the format's three fraction digits.
    pub fn message(&self) -> String {
        match self {
            BrentError::RelativeToleranceTooSmall { value, minimum } => format!(
                "{} is smaller than the minimum ({})",
                format_number(*value),
                format_number(*minimum)
            ),
            BrentError::AbsoluteToleranceNotPositive { value } => format!(
                "{} is smaller than, or equal to, the minimum (0)",
                format_number(*value)
            ),
            BrentError::IntervalNotIncreasing { lower, upper } => format!(
                "{} is larger than, or equal to, the maximum ({})",
                format_number(*lower),
                format_number(*upper)
            ),
            BrentError::StartOutOfRange {
                value,
                lower,
                upper,
            } => format!(
                "{} out of [{}, {}] range",
                format_number(*value),
                format_number(*lower),
                format_number(*upper)
            ),
            BrentError::TooManyEvaluations { maximum } => {
                format!("illegal state: maximal count ({maximum}) exceeded: evaluations")
            }
        }
    }
}

/// `NumberFormat.getInstance(Locale.US)` as `MessageFormat` uses it: grouping by threes, at most
/// three fraction digits, half-even.
///
/// The reference rounds the double's exact decimal value; this rounds the value scaled by a
/// thousand. The two agree on everything the goldens carry, and can differ where the scaling lands
/// exactly on a tie that the exact value is not on -- `0.0005` is the nearest such double, and
/// nothing measured formats one. Ties are not claimed.
fn format_number(value: f64) -> String {
    if value.is_nan() {
        return "NaN".to_string();
    }
    if value.is_infinite() {
        return if value > 0.0 {
            "∞".to_string()
        } else {
            "-∞".to_string()
        };
    }
    // Three fraction digits, half-even, then trailing zeros dropped.
    let scaled = value * 1000.0;
    let rounded = round_half_even(scaled) / 1000.0;
    let negative = rounded < 0.0 || (rounded == 0.0 && value.is_sign_negative() && rounded != 0.0);
    let magnitude = rounded.abs();
    let whole = magnitude.trunc();
    let fraction = ((magnitude - whole) * 1000.0).round() as i64;
    let mut text = group(whole);
    if fraction != 0 {
        let mut digits = format!("{fraction:03}");
        while digits.ends_with('0') {
            digits.pop();
        }
        text.push('.');
        text.push_str(&digits);
    }
    if negative {
        format!("-{text}")
    } else {
        text
    }
}

/// `RoundingMode.HALF_EVEN`, which is what `NumberFormat` uses by default.
fn round_half_even(value: f64) -> f64 {
    let floor = value.floor();
    let difference = value - floor;
    // Above a half rounds up, below rounds down, and exactly a half goes to the even neighbour --
    // which is `floor` when `floor` is even. The two "round down" arms are the same expression and
    // different reasons.
    #[allow(clippy::if_same_then_else)]
    if difference > 0.5 {
        floor + 1.0
    } else if difference < 0.5 {
        floor
    } else if (floor / 2.0).fract() == 0.0 {
        floor
    } else {
        floor + 1.0
    }
}

/// The grouping separator every three digits, which is why a thousand renders as `1,000`.
fn group(whole: f64) -> String {
    let digits = format!("{}", whole as i64);
    let mut grouped = String::new();
    for (index, character) in digits.chars().enumerate() {
        if index > 0 && (digits.len() - index) % 3 == 0 {
            grouped.push(',');
        }
        grouped.push(character);
    }
    grouped
}

/// `OptimizationUtils.max(function, min, max, guess, relative, absolute, maxEvaluations)`.
///
/// The objective is evaluated at most `max_evaluations` times; the call that would be one past that
/// is [`BrentError::TooManyEvaluations`].
pub fn maximize(
    objective: impl Fn(f64) -> f64,
    min: f64,
    max: f64,
    guess: f64,
    relative_tolerance: f64,
    absolute_tolerance: f64,
    max_evaluations: i32,
) -> Result<PointValuePair, BrentError> {
    optimize(
        objective,
        GoalType::Maximize,
        min,
        max,
        guess,
        relative_tolerance,
        absolute_tolerance,
        max_evaluations,
    )
}

/// `BrentOptimizer(rel, abs).optimize(objective, goal, new SearchInterval(min, max, guess), maxEval)`.
#[allow(clippy::too_many_arguments)]
pub fn optimize(
    objective: impl Fn(f64) -> f64,
    goal: GoalType,
    min: f64,
    max: f64,
    guess: f64,
    relative_threshold: f64,
    absolute_threshold: f64,
    max_evaluations: i32,
) -> Result<PointValuePair, BrentError> {
    // The optimiser is constructed before the interval is built, so its refusals come first.
    if relative_threshold < MIN_RELATIVE_TOLERANCE {
        return Err(BrentError::RelativeToleranceTooSmall {
            value: relative_threshold,
            minimum: MIN_RELATIVE_TOLERANCE,
        });
    }
    if absolute_threshold <= 0.0 {
        return Err(BrentError::AbsoluteToleranceNotPositive {
            value: absolute_threshold,
        });
    }
    if min >= max {
        return Err(BrentError::IntervalNotIncreasing {
            lower: min,
            upper: max,
        });
    }
    if guess < min || guess > max {
        return Err(BrentError::StartOutOfRange {
            value: guess,
            lower: min,
            upper: max,
        });
    }

    let is_minimum = goal == GoalType::Minimize;
    let mut evaluations = 0;
    // `computeObjectiveValue` increments first and refuses when the count passes the maximum.
    let evaluate = |x: f64, evaluations: &mut i32| -> Result<f64, BrentError> {
        *evaluations += 1;
        if *evaluations > max_evaluations {
            return Err(BrentError::TooManyEvaluations {
                maximum: max_evaluations,
            });
        }
        Ok(objective(x))
    };

    // `lo < hi` is guaranteed above, so `a` is the lower bound; the reference sorts anyway.
    let (mut a, mut b) = if min < max { (min, max) } else { (max, min) };

    let mut x = guess;
    let mut v = x;
    let mut w = x;
    let mut d = 0.0;
    let mut e = 0.0;
    let mut fx = evaluate(x, &mut evaluations)?;
    if !is_minimum {
        fx = -fx;
    }
    let mut fv = fx;
    let mut fw = fx;

    let mut previous: Option<PointValuePair> = None;
    let mut current = PointValuePair {
        point: x,
        value: if is_minimum { fx } else { -fx },
    };
    let mut best = current;

    loop {
        let m = 0.5 * (a + b);
        let tol1 = relative_threshold * java_abs(x) + absolute_threshold;
        let tol2 = 2.0 * tol1;

        let stop = java_abs(x - m) <= tol2 - 0.5 * (b - a);
        if stop {
            return Ok(pick(
                best,
                pick_option(previous, Some(current), is_minimum),
                is_minimum,
            ));
        }

        let mut p;
        let mut q;
        let r;

        if java_abs(e) > tol1 {
            // Fit a parabola.
            let r_first = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r_first;
            q = 2.0 * (q - r_first);

            if q > 0.0 {
                p = -p;
            } else {
                q = -q;
            }

            r = e;
            e = d;

            if p > q * (a - x) && p < q * (b - x) && java_abs(p) < java_abs(0.5 * q * r) {
                d = p / q;
                let candidate = x + d;

                // `f` must not be evaluated too close to `a` or `b`.
                if candidate - a < tol2 || b - candidate < tol2 {
                    d = if x <= m { tol1 } else { -tol1 };
                }
            } else {
                // Golden section step.
                e = if x < m { b - x } else { a - x };
                d = GOLDEN_SECTION * e;
            }
        } else {
            // Golden section step.
            e = if x < m { b - x } else { a - x };
            d = GOLDEN_SECTION * e;
        }

        // Update by at least `tol1`. The reference declares `u` at zero and assigns it in the
        // parabolic branch too; that assignment is dead, since this one is unconditional.
        let u = if java_abs(d) < tol1 {
            if d >= 0.0 {
                x + tol1
            } else {
                x - tol1
            }
        } else {
            x + d
        };

        let mut fu = evaluate(u, &mut evaluations)?;
        if !is_minimum {
            fu = -fu;
        }

        previous = Some(current);
        current = PointValuePair {
            point: u,
            value: if is_minimum { fu } else { -fu },
        };
        best = pick(
            best,
            pick_option(previous, Some(current), is_minimum),
            is_minimum,
        );

        // No convergence checker is supplied by `OptimizationUtils`, so there is no early return.

        if fu <= fx {
            if u < x {
                b = x;
            } else {
                a = x;
            }
            v = w;
            fv = fw;
            w = x;
            fw = fx;
            x = u;
            fx = fu;
        } else {
            if u < x {
                a = u;
            } else {
                b = u;
            }
            if fu <= fw || equals_within_one_ulp(w, x) {
                v = w;
                fv = fw;
                w = u;
                fw = fu;
            } else if fu <= fv || equals_within_one_ulp(v, x) || equals_within_one_ulp(v, w) {
                v = u;
                fv = fu;
            }
        }
    }
}

/// `best(a, b, isMinim)` with both present.
fn pick(a: PointValuePair, b: PointValuePair, is_minimum: bool) -> PointValuePair {
    if is_minimum {
        if a.value <= b.value {
            a
        } else {
            b
        }
    } else if a.value >= b.value {
        a
    } else {
        b
    }
}

/// The same, where either may be absent, as `previous` is on the first iteration.
fn pick_option(
    a: Option<PointValuePair>,
    b: Option<PointValuePair>,
    is_minimum: bool,
) -> PointValuePair {
    match (a, b) {
        (None, Some(b)) => b,
        (Some(a), None) => a,
        (Some(a), Some(b)) => pick(a, b, is_minimum),
        (None, None) => unreachable!("the current pair is always present"),
    }
}

/// `Precision.equals(x, y)`, which is `equals(x, y, 1)`: within one ulp, and false for NaN.
fn equals_within_one_ulp(x: f64, y: f64) -> bool {
    const SIGN_MASK: i64 = i64::MIN;
    const POSITIVE_ZERO_BITS: i64 = 0;
    const NEGATIVE_ZERO_BITS: i64 = i64::MIN;
    let x_bits = x.to_bits() as i64;
    let y_bits = y.to_bits() as i64;
    let is_equal = if ((x_bits ^ y_bits) & SIGN_MASK) == 0 {
        // Same sign, so the subtraction cannot overflow.
        (x_bits - y_bits).abs() <= 1
    } else {
        let (delta_plus, delta_minus) = if x_bits < y_bits {
            (
                y_bits - POSITIVE_ZERO_BITS,
                x_bits.wrapping_sub(NEGATIVE_ZERO_BITS),
            )
        } else {
            (
                x_bits - POSITIVE_ZERO_BITS,
                y_bits.wrapping_sub(NEGATIVE_ZERO_BITS),
            )
        };
        if delta_plus > 1 {
            false
        } else {
            delta_minus <= 1 - delta_plus
        }
    };
    is_equal && !x.is_nan() && !y.is_nan()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every value here is the reference's, from gatk-rs's `strand-artifact-mstep` golden.
    #[test]
    fn the_analytic_objectives_answer_what_the_reference_answers() {
        let quadratic = |x: f64| -(x - 3.0) * (x - 3.0);
        let pair = maximize(quadratic, 0.01, 100.0, 1.0, 0.01, 0.01, 100).expect("optimised");
        assert_eq!(pair.point, 3.0000000000000013);
        assert_eq!(pair.value, -1.7749370367472766E-30);

        let pair = maximize(
            |x: f64| -(x - 50.0) * (x - 50.0),
            0.01,
            100.0,
            1.0,
            0.01,
            0.01,
            100,
        )
        .expect("optimised");
        assert_eq!(pair.point, 50.0);
        assert_eq!(pair.value, -0.0, "which is negative zero");
        assert!(pair.value.is_sign_negative());

        // The bounds are interval ends, never candidates.
        let pair = maximize(|x: f64| -x, 0.01, 100.0, 1.0, 0.01, 0.01, 100).expect("optimised");
        assert_eq!(pair.point, 0.020762640150789262);
        let pair = maximize(|x: f64| x, 0.01, 100.0, 1.0, 0.01, 0.01, 100).expect("optimised");
        assert_eq!(pair.point, 98.88158923714171);

        // A flat objective answers the guess.
        let pair = maximize(|_| 1.0, 0.01, 100.0, 1.0, 0.01, 0.01, 100).expect("optimised");
        assert_eq!(pair.point, 1.0);
        assert_eq!(pair.value, 1.0);

        // Two maxima in one interval: which is found is the trajectory's business.
        let pair = maximize(f64::sin, 0.01, 20.0, 1.0, 0.01, 0.01, 100).expect("optimised");
        assert_eq!(pair.point, 7.868010716446599);
        assert_eq!(pair.value, 0.9999015940364933);

        // The guess on each bound.
        let pair = maximize(quadratic, 0.01, 100.0, 0.01, 0.01, 0.01, 100).expect("optimised");
        assert_eq!(pair.point, 3.0000000000000004);
        let pair = maximize(quadratic, 0.01, 100.0, 100.0, 0.01, 0.01, 100).expect("optimised");
        assert_eq!(pair.point, 3.0000000000000018);

        // A tolerance ten million times tighter walks to the same double.
        let pair = maximize(quadratic, 0.01, 100.0, 1.0, 1e-10, 1e-10, 1000).expect("optimised");
        assert_eq!(pair.point, 3.0000000000000013);

        // A NaN objective stops at the second point and reports NaN.
        let pair = maximize(|_| f64::NAN, 0.01, 100.0, 1.0, 0.01, 0.01, 100).expect("optimised");
        assert_eq!(pair.point, 1.02);
        assert!(pair.value.is_nan());
    }

    /// The five refusals, and the messages the reference formats.
    #[test]
    fn what_it_refuses_and_how_it_says_so() {
        let quadratic = |x: f64| -(x - 3.0) * (x - 3.0);

        let error = maximize(quadratic, 100.0, 0.01, 1.0, 0.01, 0.01, 100).expect_err("refused");
        assert_eq!(
            error.message(),
            "100 is larger than, or equal to, the maximum (0.01)"
        );

        let error = maximize(quadratic, 0.01, 100.0, 1.0, 0.01, 0.01, 3).expect_err("refused");
        assert_eq!(
            error.message(),
            "illegal state: maximal count (3) exceeded: evaluations"
        );

        // Both numbers are below three fraction digits, so both render as zero.
        let error = maximize(quadratic, 0.01, 100.0, 1.0, 1e-17, 0.01, 100).expect_err("refused");
        assert_eq!(error.message(), "0 is smaller than the minimum (0)");

        let error = maximize(quadratic, 0.01, 100.0, 1.0, 0.01, 0.0, 100).expect_err("refused");
        assert_eq!(
            error.message(),
            "0 is smaller than, or equal to, the minimum (0)"
        );

        let error = maximize(quadratic, 0.01, 100.0, 200.0, 0.01, 0.01, 100).expect_err("refused");
        assert_eq!(error.message(), "200 out of [0.01, 100] range");
    }

    /// The number format's own corners: grouping, three fraction digits, half-even.
    #[test]
    fn the_message_numbers_are_formatted_not_printed() {
        assert_eq!(format_number(100.0), "100");
        assert_eq!(format_number(0.01), "0.01");
        assert_eq!(format_number(1e-17), "0");
        assert_eq!(format_number(1000.0), "1,000");
        assert_eq!(format_number(1234567.0), "1,234,567");
        // Every number the measured messages carry, and its rendering.
        assert_eq!(format_number(200.0), "200");
        assert_eq!(format_number(MIN_RELATIVE_TOLERANCE), "0");
        assert_eq!(format_number(0.0), "0");
    }

    /// The written-out root is the one `sqrt` produces.
    #[test]
    fn the_golden_section_is_computed_from_the_correctly_rounded_root() {
        assert_eq!(2.23606797749979_f64.to_bits(), 5.0_f64.sqrt().to_bits());
        assert_eq!(GOLDEN_SECTION, 0.5 * (3.0 - 5.0_f64.sqrt()));
    }

    /// `Precision.equals(x, y)` is one ulp, not equality, and NaN is not equal to itself.
    #[test]
    fn the_ulp_comparison_is_the_reference_one() {
        assert!(equals_within_one_ulp(1.0, 1.0));
        assert!(equals_within_one_ulp(
            1.0,
            f64::from_bits(1.0f64.to_bits() + 1)
        ));
        assert!(!equals_within_one_ulp(
            1.0,
            f64::from_bits(1.0f64.to_bits() + 2)
        ));
        assert!(equals_within_one_ulp(0.0, -0.0));
        assert!(!equals_within_one_ulp(f64::NAN, f64::NAN));
    }
}
