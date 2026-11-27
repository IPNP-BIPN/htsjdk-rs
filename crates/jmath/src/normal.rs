//! `NormalDistribution`, ported from `org.apache.commons.math3.distribution` (commons-math3 3.5),
//! which is Apache 2.0.
//!
//! Two methods, and GATK's `MannWhitneyU` uses both: the cumulative probability turns a Z score
//! into a p-value, and the inverse turns a p-value back into a Z score. They are **not** inverses
//! of each other in code. The first goes through `Erf.erf`, an incomplete gamma function evaluated
//! by series; the second goes through `Erf.erfInv`, a rational approximation that touches no gamma
//! function at all. A round trip need not return its input bit for bit, so neither can be derived
//! from the other.
//!
//! `MannWhitneyU` (GATK 4.7.0.0) holds one standard normal,
//! `new NormalDistribution(NORMAL_MEAN, NORMAL_SD)` with `NORMAL_MEAN = 0` and `NORMAL_SD = 1`,
//! and uses the two methods on different branches, never on the same value:
//!
//! ```java
//! p = 2 * NORMAL.cumulativeProbability(NORMAL_MEAN + z * NORMAL_SD);
//! ```
//!
//! when either sample reaches `MINIMUM_NORMAL_N`, and otherwise, after the exact permutation test,
//!
//! ```java
//! z = NORMAL.inverseCumulativeProbability(p);
//! ```
//!
//! # The 40-sigma shortcut
//!
//! ```java
//! if (FastMath.abs(dev) > 40 * standardDeviation) { return dev < 0 ? 0.0d : 1.0d; }
//! ```
//!
//! Beyond forty standard deviations the answer is exactly 0 or 1 without computing anything, so
//! the function is discontinuous at a value the author chose rather than at one the arithmetic
//! forced.
//!
//! # Reading this file without reading Rust
//!
//! - `struct` with `pub` fields is a Java class with public fields; its methods are in the
//!   `impl NormalDistribution` block, with `&self` as Java's `this`.
//! - `impl Default` is the no-argument constructor, here the standard normal.
//! - `Option<f64>` is "a double or nothing": `Some(x)` is an answer and `None` stands for the
//!   reference's `OutOfRangeException`.

use crate::gamma;

/// `AbstractRealDistribution.SQRT2`.
///
/// The correctly rounded square root of two. Both methods use it, so the
/// `gamma_erf_normal` golden's 230 `cumulativeProbability` and 109 `inverseCumulativeProbability`
/// rows are what check that it is the reference's double.
const SQRT2: f64 = std::f64::consts::SQRT_2;

/// The distribution, which `MannWhitneyU` builds as the standard normal.
#[derive(Debug, Clone, Copy)]
pub struct NormalDistribution {
    pub mean: f64,
    pub standard_deviation: f64,
}

impl Default for NormalDistribution {
    fn default() -> Self {
        NormalDistribution {
            mean: 0.0,
            standard_deviation: 1.0,
        }
    }
}

impl NormalDistribution {
    /// `new NormalDistribution(mean, sd)`. Nothing here checks the arguments: a standard
    /// deviation that is not positive is accepted.
    pub fn new(mean: f64, standard_deviation: f64) -> Self {
        NormalDistribution {
            mean,
            standard_deviation,
        }
    }

    /// `cumulativeProbability(x)`: the 40-sigma shortcut of the module note, then
    /// `0.5 * (1 + erf(dev / (sd * SQRT2)))` in that order, with commons-math's `Erf.erf`, not a
    /// libm `erf`.
    pub fn cumulative_probability(&self, x: f64) -> f64 {
        let dev = x - self.mean;
        if dev.abs() > 40.0 * self.standard_deviation {
            return if dev < 0.0 { 0.0 } else { 1.0 };
        }
        0.5 * (1.0 + gamma::erf(dev / (self.standard_deviation * SQRT2)))
    }

    /// `inverseCumulativeProbability(p)`. Outside `[0, 1]` the reference throws
    /// `OutOfRangeException`; here that is `None`.
    pub fn inverse_cumulative_probability(&self, p: f64) -> Option<f64> {
        // `p < 0.0 || p > 1.0`, and not a range check: both comparisons are false for NaN, so a
        // NaN probability is not out of range to the reference. It flows through and comes back
        // as NaN.
        #[allow(clippy::manual_range_contains)]
        if p < 0.0 || p > 1.0 {
            return None;
        }
        Some(self.mean + self.standard_deviation * SQRT2 * gamma::erf_inv(2.0 * p - 1.0))
    }
}
