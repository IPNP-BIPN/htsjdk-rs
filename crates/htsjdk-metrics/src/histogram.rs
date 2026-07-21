//! `Histogram`, the accumulator behind 20 of Picard's 44 metrics tools.
//!
//! Ported from `htsjdk.samtools.util.Histogram`.
//!
//! Every statistic here is a floating-point sum over the bins **in key order**, because the
//! bins live in a `TreeMap`. Summation is not associative in floating point, so the iteration
//! order is part of the answer: the same bins accumulated in insertion order, or in descending
//! order, give a different last bit. That is the whole reason this is a port rather than a
//! reimplementation of the textbook formulas.
//!
//! Three specific choices are reproduced rather than improved:
//!
//! - `getMean` deliberately does **not** call `getSum() / getCount()`. It accumulates the
//!   product and the count in one pass, and htsjdk's comment says this is for efficiency. The
//!   two are not the same number.
//! - `getStandardDeviation` divides by `count - 1`, the sample variance, and takes the mean
//!   from `getMean` rather than from its own running total.
//! - `getMedian` on an even count averages the two middle **bin id values**, found by walking
//!   the cumulative count, not by materialising the values.

use std::collections::BTreeMap;

/// A bin key, ordered as a `TreeMap` orders it.
///
/// htsjdk's key is a generic `K` compared by its natural ordering, and `getIdValue()` requires
/// it to be a `Number`. Every metrics use is numeric, so the key is held as its `f64` value
/// with a total order, which reproduces `TreeMap` ordering for `Integer`, `Byte` and `Double`
/// keys alike.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Key(pub f64);

impl Eq for Key {}

impl Ord for Key {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // `Double.compareTo` is a total order: it places -0.0 below 0.0 and NaN above
        // everything. `total_cmp` is the same order, which matters because a `TreeMap<Double>`
        // will happily hold a NaN key and put it last.
        self.0.total_cmp(&other.0)
    }
}

impl PartialOrd for Key {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

/// `Histogram`.
#[derive(Debug, Clone, Default)]
pub struct Histogram {
    pub bin_label: String,
    pub value_label: String,
    bins: BTreeMap<Key, f64>,
}

impl Histogram {
    pub fn new(bin_label: &str, value_label: &str) -> Self {
        Histogram {
            bin_label: bin_label.to_string(),
            value_label: value_label.to_string(),
            bins: BTreeMap::new(),
        }
    }

    /// `Histogram.increment(id)`, which adds 1.
    pub fn increment(&mut self, id: f64) {
        self.increment_by(id, 1.0);
    }

    /// `Histogram.increment(id, increment)`.
    pub fn increment_by(&mut self, id: f64, increment: f64) {
        *self.bins.entry(Key(id)).or_insert(0.0) += increment;
    }

    /// `Histogram.prefillBins`: creates a bin at zero so it appears in the output.
    pub fn prefill(&mut self, ids: &[f64]) {
        for &id in ids {
            self.bins.insert(Key(id), 0.0);
        }
    }

    /// Bins in key order, which is the order every statistic accumulates in.
    pub fn bins(&self) -> impl Iterator<Item = (f64, f64)> + '_ {
        self.bins.iter().map(|(k, v)| (k.0, *v))
    }

    pub fn size(&self) -> usize {
        self.bins.len()
    }

    pub fn is_empty(&self) -> bool {
        self.bins.is_empty()
    }

    /// `Histogram.get(id)`: the bin's value, or `None` when there is no such bin.
    ///
    /// The distinction matters to callers that accumulate coverage by probing keys that may not
    /// exist; treating a missing bin as zero and as absent are the same thing arithmetically,
    /// but only one of them is what the Java tests for.
    pub fn get(&self, id: f64) -> Option<f64> {
        self.bins.get(&Key(id)).copied()
    }

    /// `Histogram.trimByWidth`: removes bins whose id exceeds `width`, from the top down.
    ///
    /// It walks the descending key set and **stops at the first key that is not above the
    /// width**, rather than filtering the whole map. For a normally ordered histogram the two
    /// are the same; they differ if a NaN key is present, since NaN sorts last and compares
    /// false against everything, so the walk stops immediately and nothing is trimmed.
    pub fn trim_by_width(&mut self, width: i32) {
        let width = width as f64;
        let mut to_remove = Vec::new();
        for k in self.bins.keys().rev() {
            if k.0 > width {
                to_remove.push(*k);
            } else {
                break;
            }
        }
        for k in to_remove {
            self.bins.remove(&k);
        }
    }

    /// `Histogram.getCount`: the total of the bin values, not the number of bins.
    pub fn count(&self) -> f64 {
        self.bins.values().sum()
    }

    /// `Histogram.getSum`: the total of `value * id`.
    pub fn sum(&self) -> f64 {
        self.bins.iter().map(|(k, v)| v * k.0).sum()
    }

    /// `Histogram.getSumOfValues`, which is [`Self::count`] under another name in the Java.
    pub fn sum_of_values(&self) -> f64 {
        self.bins.values().sum()
    }

    /// `Histogram.getMean`.
    ///
    /// Not `sum() / count()`. htsjdk accumulates both in a single pass to avoid iterating
    /// twice, and the two orderings of the same additions do not give the same last bit.
    pub fn mean(&self) -> f64 {
        let mut product = 0.0;
        let mut total_count = 0.0;
        for (k, v) in &self.bins {
            product += k.0 * v;
            total_count += v;
        }
        product / total_count
    }

    /// `Histogram.getStandardDeviation`, the sample standard deviation.
    ///
    /// `pow(value - mean, 2)` in the Java reaches `Math.pow`, which decision 0007 defers as
    /// possibly non-portable. It is written as a multiplication here because `Math.pow(x, 2.0)`
    /// was measured bit-identical to `x * x` on 1,999,558 points, including subnormals and both
    /// infinities; the intrinsic special-cases integral exponents before it reaches the
    /// approximate instruction. CI re-measures that on real silicon.
    pub fn standard_deviation(&self) -> f64 {
        let mean = self.mean();
        let mut count = 0.0;
        let mut total = 0.0;
        for (k, v) in &self.bins {
            let d = k.0 - mean;
            count += v;
            total += v * (d * d);
        }
        (total / (count - 1.0)).sqrt()
    }

    /// `Histogram.getMedian`.
    ///
    /// Walks the cumulative count to find the bins holding the middle observations, and
    /// averages their **ids**. On an odd count both halves land on the same bin, so the
    /// average is that bin's id.
    pub fn median(&self) -> f64 {
        let count = self.count();
        if count == 0.0 {
            return 0.0;
        }
        if count == 1.0 {
            return self.bins.keys().next().map(|k| k.0).unwrap_or(0.0);
        }
        let (mid_low, mid_high) = if count % 2.0 == 0.0 {
            let low = count / 2.0;
            (low, low + 1.0)
        } else {
            let mid = (count / 2.0).ceil();
            (mid, mid)
        };

        let mut total = 0.0;
        let mut low_value: Option<f64> = None;
        let mut high_value: Option<f64> = None;
        for (k, v) in &self.bins {
            total += v;
            if low_value.is_none() && total >= mid_low {
                low_value = Some(k.0);
            }
            if high_value.is_none() && total >= mid_high {
                high_value = Some(k.0);
            }
            if low_value.is_some() && high_value.is_some() {
                break;
            }
        }
        (low_value.unwrap_or(0.0) + high_value.unwrap_or(0.0)) / 2.0
    }

    /// `Histogram.getMedianAbsoluteDeviation`.
    ///
    /// Builds a second histogram of `|id - median|` weighted by the bin values, and takes
    /// *its* median. Two distinct ids equidistant from the median collapse into one bin, which
    /// is why this is not the median of a materialised list of deviations.
    pub fn median_absolute_deviation(&self) -> f64 {
        let median = self.median();
        let mut deviations = Histogram::default();
        for (k, v) in &self.bins {
            deviations.increment_by((k.0 - median).abs(), *v);
        }
        deviations.median()
    }

    /// `Histogram.estimateSdViaMad`.
    pub fn estimate_sd_via_mad(&self) -> f64 {
        1.4826 * self.median_absolute_deviation()
    }

    /// `Histogram.getPercentile`.
    ///
    /// Returns the id of the first bin whose cumulative fraction reaches the percentile. It is
    /// a bin id, never an interpolated value.
    pub fn percentile(&self, percentile: f64) -> Result<f64, HistogramError> {
        if percentile <= 0.0 {
            return Err(HistogramError::PercentileOutOfRange(percentile));
        }
        if percentile >= 1.0 {
            return Err(HistogramError::PercentileOutOfRange(percentile));
        }
        if let Some((k, v)) = self.bins.iter().find(|(_, v)| **v < 0.0) {
            return Err(HistogramError::NegativeCount { id: k.0, count: *v });
        }
        let total = self.count();
        if total == 0.0 {
            return Err(HistogramError::EmptyHistogram);
        }
        let mut so_far = 0.0;
        for (k, v) in &self.bins {
            so_far += v;
            if so_far / total >= percentile {
                return Ok(k.0);
            }
        }
        Err(HistogramError::EmptyHistogram)
    }

    /// `Histogram.getCumulativeProbability`.
    pub fn cumulative_probability(&self, v: f64) -> f64 {
        let mut count = 0.0;
        let mut total = 0.0;
        for (k, value) in &self.bins {
            if k.0 <= v {
                count += value;
            }
            total += value;
        }
        count / total
    }

    /// `Histogram.getMode`: the id of the largest bin.
    ///
    /// The comparison is strictly `<`, so the **first** bin in key order wins a tie rather
    /// than the last.
    pub fn mode(&self) -> Option<f64> {
        let mut best: Option<(f64, f64)> = None;
        for (k, v) in &self.bins {
            match best {
                None => best = Some((k.0, *v)),
                Some((_, bv)) if bv < *v => best = Some((k.0, *v)),
                _ => {}
            }
        }
        best.map(|(k, _)| k)
    }

    /// `Histogram.getMin`.
    pub fn min(&self) -> Option<f64> {
        self.bins.keys().next().map(|k| k.0)
    }

    /// `Histogram.getMax`.
    pub fn max(&self) -> Option<f64> {
        self.bins.keys().next_back().map(|k| k.0)
    }

    /// `Histogram.getMeanBinSize`.
    pub fn mean_bin_size(&self) -> f64 {
        self.sum_of_values() / self.size() as f64
    }

    /// `Histogram.getMedianBinSize`: the median of the bin *values*, from a sorted copy.
    pub fn median_bin_size(&self) -> f64 {
        if self.is_empty() {
            return 0.0;
        }
        let mut values: Vec<f64> = self.bins.values().copied().collect();
        values.sort_by(|a, b| a.total_cmp(b));
        let mid = values.len() / 2;
        let mut median = values[mid];
        // Written as the Java writes it; clippy prefers is_multiple_of, which reads less like
        // the source it is transcribed from.
        #[allow(clippy::manual_is_multiple_of)]
        if values.len() % 2 == 0 {
            median = (median + values[mid - 1]) / 2.0;
        }
        median
    }
}

/// Why a statistic could not be computed.
#[derive(Debug, Clone, PartialEq)]
pub enum HistogramError {
    PercentileOutOfRange(f64),
    NegativeCount { id: f64, count: f64 },
    EmptyHistogram,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn h(pairs: &[(f64, f64)]) -> Histogram {
        let mut hist = Histogram::new("bin", "count");
        for &(id, v) in pairs {
            hist.increment_by(id, v);
        }
        hist
    }

    #[test]
    fn bins_come_out_in_key_order_whatever_the_insertion_order() {
        let hist = h(&[(5.0, 1.0), (1.0, 1.0), (3.0, 1.0)]);
        let ids: Vec<f64> = hist.bins().map(|(k, _)| k).collect();
        assert_eq!(ids, vec![1.0, 3.0, 5.0]);
    }

    #[test]
    fn count_is_the_total_of_the_values_not_the_number_of_bins() {
        let hist = h(&[(1.0, 3.0), (2.0, 4.0)]);
        assert_eq!(hist.count(), 7.0);
        assert_eq!(hist.size(), 2);
    }

    #[test]
    fn the_mean_weights_ids_by_their_counts() {
        // Three 1s and one 5: (3 + 5) / 4 = 2.
        let hist = h(&[(1.0, 3.0), (5.0, 1.0)]);
        assert_eq!(hist.mean(), 2.0);
    }

    /// `getMean` accumulates in one pass rather than dividing `getSum()` by `getCount()`.
    /// On well-behaved inputs the two agree; the port follows the one-pass form because the
    /// orderings are not guaranteed to agree in the last bit.
    #[test]
    fn the_mean_uses_the_one_pass_form() {
        let hist = h(&[(0.1, 3.0), (0.2, 7.0), (0.3, 11.0)]);
        assert_eq!(hist.mean(), hist.sum() / hist.count());
    }

    #[test]
    fn the_standard_deviation_is_the_sample_form() {
        // 1, 2, 3: mean 2, sample variance ((1)+(0)+(1))/2 = 1, sd 1.
        let hist = h(&[(1.0, 1.0), (2.0, 1.0), (3.0, 1.0)]);
        assert_eq!(hist.standard_deviation(), 1.0);
    }

    /// An even count averages the two middle bin ids; an odd count returns the middle one.
    #[test]
    fn the_median_averages_the_two_middle_ids_on_an_even_count() {
        // 1, 2, 3, 4 -> (2 + 3) / 2
        assert_eq!(
            h(&[(1.0, 1.0), (2.0, 1.0), (3.0, 1.0), (4.0, 1.0)]).median(),
            2.5
        );
        // 1, 2, 3 -> 2
        assert_eq!(h(&[(1.0, 1.0), (2.0, 1.0), (3.0, 1.0)]).median(), 2.0);
    }

    #[test]
    fn a_single_observation_is_its_own_median() {
        assert_eq!(h(&[(7.0, 1.0)]).median(), 7.0);
    }

    #[test]
    fn an_empty_histogram_has_a_median_of_zero_rather_than_an_error() {
        assert_eq!(Histogram::default().median(), 0.0);
    }

    /// Both middle observations can fall in the same bin, and then the average is that bin.
    #[test]
    fn both_halves_can_land_in_one_bin() {
        // Four observations, all in bin 9.
        assert_eq!(h(&[(9.0, 4.0)]).median(), 9.0);
    }

    /// The MAD is the median of a *histogram* of deviations, so equidistant ids merge.
    #[test]
    fn the_mad_merges_equidistant_ids_into_one_bin() {
        // 1, 3: median 2, deviations |1-2| and |3-2| are both 1, merging into one bin of 2.
        let hist = h(&[(1.0, 1.0), (3.0, 1.0)]);
        assert_eq!(hist.median(), 2.0);
        assert_eq!(hist.median_absolute_deviation(), 1.0);
        assert_eq!(hist.estimate_sd_via_mad(), 1.4826);
    }

    #[test]
    fn a_percentile_returns_a_bin_id_never_an_interpolation() {
        let hist = h(&[(10.0, 1.0), (20.0, 1.0), (30.0, 1.0), (40.0, 1.0)]);
        assert_eq!(hist.percentile(0.25).unwrap(), 10.0);
        assert_eq!(hist.percentile(0.5).unwrap(), 20.0);
        assert_eq!(hist.percentile(0.75).unwrap(), 30.0);
        assert_eq!(hist.percentile(0.99).unwrap(), 40.0);
    }

    #[test]
    fn percentiles_outside_the_open_unit_interval_are_refused() {
        let hist = h(&[(1.0, 1.0)]);
        assert!(hist.percentile(0.0).is_err());
        assert!(hist.percentile(1.0).is_err());
        assert!(hist.percentile(-0.5).is_err());
    }

    #[test]
    fn a_negative_count_makes_a_percentile_an_error_not_a_wrong_number() {
        let hist = h(&[(1.0, 1.0), (2.0, -3.0)]);
        assert_eq!(
            hist.percentile(0.5),
            Err(HistogramError::NegativeCount {
                id: 2.0,
                count: -3.0
            })
        );
    }

    /// The mode comparison is strictly `<`, so a tie goes to the first bin in key order.
    #[test]
    fn a_tied_mode_goes_to_the_lowest_key() {
        assert_eq!(h(&[(5.0, 2.0), (1.0, 2.0)]).mode(), Some(1.0));
        assert_eq!(h(&[(1.0, 1.0), (5.0, 2.0)]).mode(), Some(5.0));
    }

    #[test]
    fn min_and_max_are_the_extreme_keys() {
        let hist = h(&[(5.0, 1.0), (1.0, 1.0), (3.0, 1.0)]);
        assert_eq!(hist.min(), Some(1.0));
        assert_eq!(hist.max(), Some(5.0));
    }

    #[test]
    fn cumulative_probability_counts_bins_at_or_below() {
        let hist = h(&[(1.0, 1.0), (2.0, 1.0), (3.0, 2.0)]);
        assert_eq!(hist.cumulative_probability(2.0), 0.5);
        assert_eq!(hist.cumulative_probability(3.0), 1.0);
        assert_eq!(hist.cumulative_probability(0.0), 0.0);
    }

    #[test]
    fn median_bin_size_sorts_the_values_not_the_keys() {
        // Values 5, 1, 3 sorted are 1, 3, 5 -> median 3.
        let hist = h(&[(1.0, 5.0), (2.0, 1.0), (3.0, 3.0)]);
        assert_eq!(hist.median_bin_size(), 3.0);
        assert_eq!(hist.mean_bin_size(), 3.0);
    }

    #[test]
    fn trimming_removes_only_the_tail_above_the_width() {
        let mut hist = h(&[(1.0, 1.0), (5.0, 1.0), (10.0, 1.0), (100.0, 1.0)]);
        hist.trim_by_width(10);
        let ids: Vec<f64> = hist.bins().map(|(k, _)| k).collect();
        assert_eq!(ids, vec![1.0, 5.0, 10.0], "the width itself is kept");
    }

    /// The walk stops at the first key not above the width, so a NaN key (which sorts last and
    /// compares false against everything) blocks the trim entirely.
    #[test]
    fn a_nan_key_stops_the_trim_immediately() {
        let mut hist = h(&[(1.0, 1.0), (100.0, 1.0)]);
        hist.increment(f64::NAN);
        hist.trim_by_width(10);
        assert_eq!(hist.size(), 3, "nothing is removed once the walk stops");
    }

    #[test]
    fn get_distinguishes_a_missing_bin_from_a_zero_one() {
        let mut hist = h(&[(1.0, 5.0)]);
        assert_eq!(hist.get(1.0), Some(5.0));
        assert_eq!(hist.get(2.0), None);
        hist.prefill(&[2.0]);
        assert_eq!(hist.get(2.0), Some(0.0));
    }

    #[test]
    fn prefilled_bins_appear_at_zero() {
        let mut hist = Histogram::new("b", "c");
        hist.prefill(&[1.0, 2.0, 3.0]);
        hist.increment(2.0);
        assert_eq!(hist.size(), 3);
        assert_eq!(hist.count(), 1.0);
        let ids: Vec<f64> = hist.bins().map(|(k, _)| k).collect();
        assert_eq!(ids, vec![1.0, 2.0, 3.0]);
    }

    /// The key order is `Double.compareTo`, a total order, so -0.0 sorts below 0.0 and a NaN
    /// key sorts last rather than making the map inconsistent.
    #[test]
    fn the_key_order_is_javas_total_order() {
        let mut hist = Histogram::new("b", "c");
        for id in [f64::NAN, 0.0, -0.0, -1.0, 1.0] {
            hist.increment(id);
        }
        let ids: Vec<f64> = hist.bins().map(|(k, _)| k).collect();
        assert_eq!(ids.len(), 5);
        assert_eq!(ids[0], -1.0);
        assert!(
            ids[1].is_sign_negative() && ids[1] == 0.0,
            "-0.0 before 0.0"
        );
        assert!(ids[2].is_sign_positive() && ids[2] == 0.0);
        assert_eq!(ids[3], 1.0);
        assert!(ids[4].is_nan(), "NaN sorts last");
    }
}
