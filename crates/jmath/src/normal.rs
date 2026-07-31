//! `NormalDistribution`, ported from `org.apache.commons.math3.distribution` (commons-math3 3.5),
//! which is Apache 2.0.
//!
//! Two methods, and GATK's `MannWhitneyU` uses both: the cumulative probability turns a Z score
//! into a p-value, and the inverse turns a p-value back into a Z score. They are **not** inverses
//! of each other in code. The first goes through `Erf.erf`, an incomplete gamma function evaluated
//! by series; the second goes through `Erf.erfInv`, a rational approximation that touches no gamma
//! function at all. A round trip does not return its input bit for bit, and `MannWhitneyU` makes
//! exactly that round trip.
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

use crate::gamma;

/// `AbstractRealDistribution.SQRT2`.
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
    pub fn new(mean: f64, standard_deviation: f64) -> Self {
        NormalDistribution {
            mean,
            standard_deviation,
        }
    }

    /// `cumulativeProbability(x)`.
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
        if !(0.0..=1.0).contains(&p) {
            return None;
        }
        Some(self.mean + self.standard_deviation * SQRT2 * gamma::erf_inv(2.0 * p - 1.0))
    }
}
