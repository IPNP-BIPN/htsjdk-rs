//! `htsjdk.samtools.ConstantMemoryDownsamplingIterator` and the statistics its base class keeps.
//!
//! The decision is a pure function of the **read name** and the seed, which is what makes it
//! constant-memory and what makes both ends of a template share a fate: a record is kept when
//! `Murmur3(seed).hashUnencodedChars(name)` is at or below a threshold derived from the proportion.
//!
//! The threshold is the part worth porting carefully:
//!
//! ```java
//! final long range = (long) Integer.MAX_VALUE - (long) Integer.MIN_VALUE;
//! this.maxHashValue = Integer.MIN_VALUE + (int) Math.round(range * proportion);
//! ```
//!
//! `range` is 2^32 - 1 as a `long`, the product is a `double`, `Math.round` answers a `long`, and
//! the `(int)` cast **truncates to the low 32 bits** rather than saturating. At a proportion of 1
//! that truncation is the whole behaviour: `Integer.MIN_VALUE + 4294967295` overflows an `int` back
//! to `Integer.MAX_VALUE`, which is the only value that keeps every record. A port that saturated
//! would answer `Integer.MAX_VALUE` too and agree here by luck; one that used 64-bit arithmetic
//! throughout would keep everything at proportions where the reference does not.
//!
//! `Math.round` is Java's, not Rust's: [`jmath::math::round`]. The proportion that reaches the
//! difference between the two is about 1.16e-10, which no corpus here carries, so this is a latent
//! difference rather than a measured one, and it is closed rather than left open.
//!
//! `picard-rs` ports the keep predicate inside `DownsampleSam`, which is the tool that drives the
//! iterator. The iterator, its threshold and its statistics are htsjdk's.

use crate::murmur3::Murmur3;

/// `ConstantMemoryDownsamplingIterator`'s keep test, and the statistics `DownsamplingIterator`
/// keeps while it runs.
pub struct ConstantMemoryDownsampler {
    hasher: Murmur3,
    max_hash_value: i32,
    target_proportion: f64,
    records_seen: u64,
    records_accepted: u64,
}

impl ConstantMemoryDownsampler {
    /// `new ConstantMemoryDownsamplingIterator(iterator, proportion, seed)`.
    ///
    /// # Panics
    ///
    /// On a proportion outside `0..=1`, which is `DownsamplingIterator`'s own
    /// `IllegalArgumentException`.
    pub fn new(proportion: f64, seed: i32) -> Self {
        assert!(proportion >= 0.0, "targetProportion must be >= 0");
        assert!(proportion <= 1.0, "targetProportion must be <= 1");
        ConstantMemoryDownsampler {
            hasher: Murmur3::new(seed),
            max_hash_value: max_hash_value(proportion),
            target_proportion: proportion,
            records_seen: 0,
            records_accepted: 0,
        }
    }

    /// The threshold a hash is compared against, exposed because it is the whole of the decision.
    pub fn max_hash_value(&self) -> i32 {
        self.max_hash_value
    }

    pub fn target_proportion(&self) -> f64 {
        self.target_proportion
    }

    /// Whether the record with this name survives. Both ends of a template hash the same name, so
    /// a pair is kept or dropped together.
    pub fn keep(&self, read_name: &str) -> bool {
        self.hasher.hash_unencoded_chars(read_name) <= self.max_hash_value
    }

    /// The same test, recording the record as seen and as accepted or discarded.
    pub fn accept(&mut self, read_name: &str) -> bool {
        let kept = self.keep(read_name);
        self.records_seen += 1;
        if kept {
            self.records_accepted += 1;
        }
        kept
    }

    pub fn seen_count(&self) -> u64 {
        self.records_seen
    }

    pub fn accepted_count(&self) -> u64 {
        self.records_accepted
    }

    pub fn discarded_count(&self) -> u64 {
        self.records_seen - self.records_accepted
    }

    /// `getAcceptedFraction`, which divides by the seen count and therefore answers `NaN` before
    /// anything has been seen rather than zero.
    pub fn accepted_fraction(&self) -> f64 {
        self.records_accepted as f64 / self.records_seen as f64
    }

    pub fn discarded_fraction(&self) -> f64 {
        self.discarded_count() as f64 / self.records_seen as f64
    }

    /// `resetStatistics`.
    pub fn reset_statistics(&mut self) {
        self.records_seen = 0;
        self.records_accepted = 0;
    }
}

/// `Integer.MIN_VALUE + (int) Math.round(range * proportion)`, in Java's arithmetic.
///
/// The `(int)` cast keeps the low 32 bits of the `long`, and the addition wraps, which is why a
/// proportion of 1 answers `Integer.MAX_VALUE` rather than overflowing.
pub fn max_hash_value(proportion: f64) -> i32 {
    let range = i32::MAX as i64 - i32::MIN as i64;
    let rounded = jmath::math::round(range as f64 * proportion);
    (i32::MIN).wrapping_add(rounded as i32)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_proportion_of_one_keeps_everything() {
        // Integer.MIN_VALUE + 4294967295 wraps to Integer.MAX_VALUE, and every hash is <= that.
        assert_eq!(max_hash_value(1.0), i32::MAX);
        let sampler = ConstantMemoryDownsampler::new(1.0, 1);
        for i in 0..200 {
            assert!(sampler.keep(&format!("read{i}")));
        }
    }

    #[test]
    fn a_proportion_of_zero_keeps_only_the_smallest_possible_hash() {
        assert_eq!(max_hash_value(0.0), i32::MIN);
        let sampler = ConstantMemoryDownsampler::new(0.0, 1);
        let kept = (0..500)
            .filter(|i| sampler.keep(&format!("read{i}")))
            .count();
        assert_eq!(kept, 0, "only a hash of exactly Integer.MIN_VALUE survives");
    }

    #[test]
    fn both_ends_of_a_template_share_a_fate() {
        let sampler = ConstantMemoryDownsampler::new(0.5, 1);
        // The hash is of the name alone, so the two ends cannot disagree.
        assert_eq!(sampler.keep("read0001"), sampler.keep("read0001"));
    }

    #[test]
    fn the_statistics_count_what_the_iterator_saw() {
        let mut sampler = ConstantMemoryDownsampler::new(0.5, 1);
        assert!(sampler.accepted_fraction().is_nan(), "0/0 before anything");
        for i in 0..100 {
            sampler.accept(&format!("read{i}"));
        }
        assert_eq!(sampler.seen_count(), 100);
        assert_eq!(
            sampler.accepted_count() + sampler.discarded_count(),
            sampler.seen_count()
        );
        sampler.reset_statistics();
        assert_eq!(sampler.seen_count(), 0);
    }

    #[test]
    fn the_seed_changes_which_records_survive() {
        let one = ConstantMemoryDownsampler::new(0.5, 1);
        let same_proportion = ConstantMemoryDownsampler::new(0.5, 42);
        assert_eq!(
            one.max_hash_value(),
            same_proportion.max_hash_value(),
            "the threshold is a function of the proportion alone"
        );
        let a = (0..200).filter(|i| one.keep(&format!("read{i}"))).count();
        let b = (0..200)
            .filter(|i| same_proportion.keep(&format!("read{i}")))
            .count();
        assert_ne!(a, b, "a different seed keeps a different set");
    }
}
