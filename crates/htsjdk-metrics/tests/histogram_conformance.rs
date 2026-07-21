//! Conformance for `Histogram` against htsjdk, over every statistic it computes.
//!
//! Goldens from `tools/histogram-conformance/HistogramDump.java` in the pinned oracle.
//!
//! Compared as **raw double bits**, not printed decimals. A decimal rendering would hide
//! exactly the last-bit differences that summation order produces, which is the whole reason
//! these statistics are ported rather than recomputed from the textbook formulas.
//!
//! One exemption exists and it is **not** a porting defect. `0.0 / 0.0` produces a NaN whose
//! sign bit is chosen by the FPU: x86-64 returns `fff8000000000000`, aarch64 returns
//! `7ff8000000000000`. Measured on both, in Java and in Rust, so it is a hardware property
//! rather than a language one. The same Rust source compiled for x86-64 matches the golden
//! exactly. See decision 0012.
//!
//! So the test exempts a NaN sign difference **only when not running on x86-64**, counts the
//! exemptions, and requires the count to be zero on the oracle's own architecture. CI runs on
//! x86-64 and is therefore the arbiter.

use std::io::Read;

use htsjdk_metrics::histogram::Histogram;

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/histogram.txt.gz");
    let f = std::fs::File::open(&p).expect("corpus");
    let mut s = String::new();
    flate2::read::GzDecoder::new(f)
        .read_to_string(&mut s)
        .expect("corpus is gzip");
    s
}

fn build(name: &str) -> Histogram {
    let mut h = Histogram::new("bin", "count");
    let mut add = |id: f64, v: f64| h.increment_by(id, v);
    match name {
        "single" => add(7.0, 1.0),
        "two" => {
            add(1.0, 1.0);
            add(3.0, 1.0);
        }
        "uniform_odd" => {
            for i in 1..=3 {
                add(i as f64, 1.0);
            }
        }
        "uniform_even" => {
            for i in 1..=4 {
                add(i as f64, 1.0);
            }
        }
        "one_bin_many" => add(9.0, 4.0),
        "weighted" => {
            add(1.0, 3.0);
            add(5.0, 1.0);
        }
        "tenths" => {
            add(0.1, 3.0);
            add(0.2, 7.0);
            add(0.3, 11.0);
        }
        "insert_size" => {
            for i in 0..600 {
                let x = 100.0 + i as f64;
                let t = (x - 350.0) / 90.0;
                add(x, (10000.0 * (-(t * t)).exp()).floor());
            }
        }
        "quality_skew" => {
            for i in 0..=45 {
                add(i as f64, (i * i) as f64);
            }
        }
        "fractional_counts" => {
            let mut rng = JavaRandom::new(20260721);
            for i in 0..50 {
                add(i as f64 * 0.37, rng.next_double() * 100.0);
            }
        }
        "wide_magnitudes" => {
            for i in 0..200 {
                add(10f64.powf((i as f64 - 100.0) / 10.0), 1.0 + (i % 7) as f64);
            }
        }
        "tied_mode" => {
            add(5.0, 2.0);
            add(1.0, 2.0);
            add(3.0, 2.0);
        }
        "negative_ids" => {
            add(-5.0, 2.0);
            add(-1.0, 3.0);
            add(0.0, 1.0);
            add(4.0, 2.0);
        }
        other => panic!("unknown case {other}"),
    }
    h
}

/// `java.util.Random`, so the fractional-count case can be reproduced exactly.
///
/// A different generator would give different counts, and the comparison would then be
/// measuring two different histograms rather than two implementations of one statistic.
struct JavaRandom {
    seed: i64,
}

impl JavaRandom {
    fn new(seed: i64) -> Self {
        JavaRandom {
            seed: (seed ^ 0x5DEECE66D) & ((1 << 48) - 1),
        }
    }
    fn next(&mut self, bits: u32) -> i32 {
        self.seed = self.seed.wrapping_mul(0x5DEECE66D).wrapping_add(0xB) & ((1 << 48) - 1);
        (self.seed >> (48 - bits)) as i32
    }
    fn next_double(&mut self) -> f64 {
        let hi = (self.next(26) as i64) << 27;
        let lo = self.next(27) as i64;
        (hi + lo) as f64 * (1.0 / ((1u64 << 53) as f64))
    }
}

fn stat(h: &Histogram, name: &str) -> Option<f64> {
    if let Some(p) = name.strip_prefix("percentile_") {
        return h.percentile(p.parse().unwrap()).ok();
    }
    if let Some(v) = name.strip_prefix("cumulativeProbability_") {
        return Some(h.cumulative_probability(v.parse().unwrap()));
    }
    Some(match name {
        "count" => h.count(),
        "sum" => h.sum(),
        "sumOfValues" => h.sum_of_values(),
        "mean" => h.mean(),
        "standardDeviation" => h.standard_deviation(),
        "median" => h.median(),
        "medianAbsoluteDeviation" => h.median_absolute_deviation(),
        "estimateSdViaMad" => h.estimate_sd_via_mad(),
        "meanBinSize" => h.mean_bin_size(),
        "medianBinSize" => h.median_bin_size(),
        "mode" => h.mode()?,
        "min" => h.min()?,
        "max" => h.max()?,
        other => panic!("unknown stat {other}"),
    })
}

/// The gate: every statistic, on every shape, bit-identical.
#[test]
fn every_statistic_matches_htsjdk_bit_for_bit() {
    let text = corpus();
    let mut checked = 0usize;
    let mut failures = Vec::new();
    let mut nan_sign_exemptions: Vec<String> = Vec::new();

    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let mut it = line.split('\t');
        let (case, name, expected) = (it.next().unwrap(), it.next().unwrap(), it.next().unwrap());
        let h = build(case);
        let ours = stat(&h, name);
        checked += 1;

        if expected == "ERROR" {
            if ours.is_some() {
                failures.push(format!("{case}/{name}: htsjdk raised, we returned a value"));
            }
            continue;
        }
        let want = f64::from_bits(u64::from_str_radix(expected, 16).unwrap());
        match ours {
            None => failures.push(format!(
                "{case}/{name}: we returned nothing, htsjdk gave {want}"
            )),
            Some(got) => {
                // Raw bits, so -0.0 and NaN payloads count. NaN == NaN is false, so compare
                // representations rather than values.
                if got.to_bits() != want.to_bits() {
                    if is_nan_sign_only(got, want) && !cfg!(target_arch = "x86_64") {
                        nan_sign_exemptions.push(format!("{case}/{name}"));
                    } else {
                        failures.push(format!(
                            "{case}/{name}: ours={got:?} ({:016x}) htsjdk={want:?} ({:016x})",
                            got.to_bits(),
                            want.to_bits()
                        ));
                    }
                }
            }
        }
    }

    assert!(checked > 300, "corpus is smaller than expected: {checked}");
    if cfg!(target_arch = "x86_64") {
        assert!(
            nan_sign_exemptions.is_empty(),
            "on x86-64 there is nothing to exempt; the FPU produces the same NaN as the oracle"
        );
    } else {
        // Reported rather than hidden, so the exemption cannot quietly grow.
        assert_eq!(
            nan_sign_exemptions,
            vec!["single/standardDeviation".to_string()],
            "the NaN-sign exemption list changed; see decision 0012"
        );
    }
    assert!(
        failures.is_empty(),
        "{} of {checked} statistics diverge:\n{}",
        failures.len(),
        failures
            .iter()
            .take(25)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
}

/// True when two values differ only in the sign bit of a NaN.
fn is_nan_sign_only(a: f64, b: f64) -> bool {
    a.is_nan() && b.is_nan() && (a.to_bits() ^ b.to_bits()) == 1 << 63
}

/// Pins the cause, so the exemption above is explained by a measurement rather than asserted.
///
/// A standard deviation over a single observation divides by `count - 1`, i.e. `0.0 / 0.0`.
/// Which NaN that produces is an FPU choice.
#[test]
fn zero_divided_by_zero_has_a_platform_dependent_nan_sign() {
    let z = std::hint::black_box(0.0f64);
    let nan = z / z;
    assert!(nan.is_nan());
    if cfg!(target_arch = "x86_64") {
        assert_eq!(
            nan.to_bits(),
            0xfff8_0000_0000_0000,
            "x86-64 returns the negative quiet NaN, which is what the oracle recorded"
        );
    } else {
        assert_eq!(
            nan.to_bits(),
            0x7ff8_0000_0000_0000,
            "aarch64 returns the positive quiet NaN"
        );
    }
}

/// The corpus must contain the shapes the port claims to handle.
#[test]
fn the_corpus_covers_the_shapes_that_matter() {
    let text = corpus();
    let cases: std::collections::BTreeSet<&str> = text
        .lines()
        .filter(|l| !l.starts_with('#') && l.contains('\t'))
        .map(|l| l.split('\t').next().unwrap())
        .collect();
    for expected in [
        "single",
        "uniform_even",
        "uniform_odd",
        "one_bin_many",
        "fractional_counts",
        "wide_magnitudes",
        "tied_mode",
        "negative_ids",
        "insert_size",
    ] {
        assert!(cases.contains(expected), "the corpus must have {expected}");
    }
}
