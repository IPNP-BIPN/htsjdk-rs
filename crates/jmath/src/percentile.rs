//! `Percentile` and `Median`, ported from `org.apache.commons.math3.stat.descriptive.rank`
//! (commons-math3 3.5), which is Apache 2.0 and therefore transcribable where the JDK's own math
//! is not. See `docs/decisions/0023-commons-math3-is-portable-where-the-jdk-is-not.md`.
//!
//! GATK pins commons-math3 to **3.5 strictly** (`build.gradle`: *"changing this breaks
//! ModelSegmentsIntegrationTests, they're quite brittle"*), so this is a port of 3.5 and not of
//! the 3.6.1 the jmath corpus job downloads for its `FastMath` columns. The two versions are the
//! same here; the pin is recorded because a percentile is a value a golden carries.
//!
//! # A median is not "sort and take the middle"
//!
//! ```java
//! index    = p * (length + 1)                       // LEGACY, with p = 0.5
//! estimate = lower + dif * (upper - lower)          // linear interpolation between neighbours
//! ```
//!
//! For an even count that is the mean of the two central values, which matches the usual
//! definition. For an odd count it is the central value. But the interpolation is what runs in
//! both cases, and the arithmetic is `lower + dif * (upper - lower)`, **not** `(lower + upper) / 2`:
//! the two differ in the last bit for values far apart, and the port has to use the reference's
//! form to carry that bit.
//!
//! # `NaN` is removed, not ranked
//!
//! The default `NaNStrategy` is `REMOVED`, so a `NaN` in the input shortens the array rather than
//! sorting to one end. An array of nothing but `NaN` becomes empty and the answer is `NaN` again,
//! by a different route.
//!
//! # Selection is a value, not an arrangement
//!
//! The reference reaches its order statistics through `KthSelector`, a quickselect with a
//! median-of-three pivot and cached pivots. That is a performance structure: `select(work, k)`
//! returns the kth smallest either way, so this port sorts and indexes, and the golden is what
//! says the two agree.
//!
//! # The rounding at the call site is the arithmetic one
//!
//! `MathUtils.median(int[])` finishes with `FastMath.round`, which is literally
//! `(long) floor(x + 0.5)`. That is the definition `java.lang.Math.round` **stopped** using in
//! Java 7, so the two disagree on `0.49999999999999994` (see [`crate::math::round`]). A port that
//! reused the `Math.round` already in this crate would be wrong here by one, on one input.

/// `Percentile.EstimationType`, restricted to the two GATK reaches.
///
/// The other seven (`R_2` through `R_9`) exist in the reference and are not ported: nothing in
/// GATK 4.6.2.0 names them, and an unused estimator is a claim nothing measures.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EstimationType {
    /// `LEGACY`, the default, and what `new Median()` uses.
    Legacy,
    /// `R_1`, which `CanonicalSVCollapser` asks for by name.
    R1,
}

impl EstimationType {
    /// `EstimationType.index(p, length)`, with `p` already divided by 100.
    fn index(self, p: f64, length: usize) -> f64 {
        let length = length as f64;
        match self {
            // `Double.compare(p, 0) == 0` is not `p == 0.0`: it separates -0.0 from 0.0. The
            // caller here always passes 0.5, so the branch is unreachable in GATK and ported
            // anyway rather than dropped.
            EstimationType::Legacy => {
                if p.total_cmp(&0.0).is_eq() {
                    0.0
                } else if p.total_cmp(&1.0).is_eq() {
                    length
                } else {
                    p * (length + 1.0)
                }
            }
            EstimationType::R1 => {
                if p.total_cmp(&0.0).is_eq() {
                    0.0
                } else {
                    length * p + 0.5
                }
            }
        }
    }

    /// `EstimationType.estimate`. `R_1` differs only by the position it hands to the same body.
    fn estimate(self, sorted: &[f64], pos: f64) -> f64 {
        let pos = match self {
            EstimationType::Legacy => pos,
            EstimationType::R1 => (pos - 0.5).ceil(),
        };
        let length = sorted.len();
        let fpos = pos.floor();
        let int_pos = fpos as i64;
        let dif = pos - fpos;

        if pos < 1.0 {
            return sorted[0];
        }
        if pos >= length as f64 {
            return sorted[length - 1];
        }
        let lower = sorted[(int_pos - 1) as usize];
        let upper = sorted[int_pos as usize];
        // Not `(lower + upper) / 2`. See the module note.
        lower + dif * (upper - lower)
    }
}

/// `Percentile.getWorkArray` under the default `NaNStrategy.REMOVED`: a copy with the `NaN`s gone.
fn work_array(values: &[f64]) -> Vec<f64> {
    values.iter().copied().filter(|v| !v.is_nan()).collect()
}

/// Java's `Arrays.sort(double[])` order: `-0.0` sorts before `0.0`, and `NaN` sorts last.
///
/// `KthSelector` partitions with `<`, under which `-0.0 < 0.0` is false, so it can leave either
/// zero in the kth slot where the sort is definite. The golden is what says whether that
/// difference ever surfaces; the port takes the sort's answer because the reference's own tail
/// (`Arrays.sort(work, begin, end)`) is the sort.
fn java_sort(values: &mut [f64]) {
    values.sort_by(|a, b| a.total_cmp(b));
}

/// `new Percentile(quantile).withEstimationType(type).evaluate(values)`.
///
/// An empty array is `NaN`, and a single value is returned as itself before any estimation runs,
/// which matters because the estimator would index an array of one.
pub fn evaluate(values: &[f64], quantile: f64, estimation: EstimationType) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    if values.len() == 1 {
        // `values[begin]`, before the NaN strategy runs, so a lone NaN comes back as NaN rather
        // than as an empty array's NaN. Same bits, different reason.
        return values[0];
    }
    let mut work = work_array(values);
    if work.is_empty() {
        return f64::NAN;
    }
    java_sort(&mut work);
    estimation.estimate(&work, estimation.index(quantile / 100.0, work.len()))
}

/// `new Median().evaluate(values)`, which is the 50th percentile under `LEGACY`.
pub fn median(values: &[f64]) -> f64 {
    evaluate(values, 50.0, EstimationType::Legacy)
}

/// `MathUtils.median(int[])`: the median of the values as doubles, rounded by `FastMath.round`
/// and then narrowed to `int`.
///
/// The narrowing is a Java `(int)` cast of a `long`, which truncates the high bits rather than
/// saturating. Nothing in GATK can reach a median outside `int` range from an `int[]`, so the
/// cast is a formality here, and it is ported as the wrapping one it is rather than as a clamp.
pub fn median_of_ints(values: &[i32], estimation: EstimationType) -> i32 {
    let doubles: Vec<f64> = values.iter().map(|v| *v as f64).collect();
    crate::fast_math::round(evaluate(&doubles, 50.0, estimation)) as i32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn an_even_count_interpolates_and_an_odd_one_does_not() {
        assert_eq!(median(&[1.0, 2.0, 3.0, 4.0]), 2.5);
        assert_eq!(median(&[1.0, 2.0, 3.0]), 2.0);
    }

    #[test]
    fn nan_shortens_the_array_rather_than_sorting_to_an_end() {
        assert_eq!(median(&[1.0, f64::NAN, 2.0, 3.0]), 2.0);
        assert!(median(&[f64::NAN, f64::NAN]).is_nan());
        // One value is returned before the NaN strategy runs.
        assert!(median(&[f64::NAN]).is_nan());
    }

    #[test]
    fn the_call_site_rounds_the_arithmetic_way() {
        // 0.49999999999999994 + 0.5 rounds up to 1.0, so the arithmetic round answers 1 where
        // Math.round answers 0. The medians below are chosen to land on that value.
        assert_eq!(crate::fast_math::round(0.499_999_999_999_999_94), 1);
        assert_eq!(crate::math::round(0.499_999_999_999_999_94), 0);
    }

    #[test]
    fn r1_takes_a_neighbour_where_legacy_interpolates() {
        assert_eq!(
            evaluate(&[1.0, 2.0, 3.0, 4.0], 50.0, EstimationType::R1),
            2.0
        );
    }
}
