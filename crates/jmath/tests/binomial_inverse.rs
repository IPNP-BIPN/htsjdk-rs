//! Conformance for the binomial quantile against commons-math3 3.5.
//!
//! Golden from `tools/jmath-conformance/BinomialInverseDump.java`.
//!
//! # What this suite is for
//!
//!  * **721 quantiles**, over eight trial counts, ten probabilities and nine `p` values;
//!  * **the cumulative probabilities the search bisects over**, so a disagreement can be placed in
//!    the CDF or in the search rather than somewhere between them;
//!  * **the mean and the variance**, which decide whether the Chebyshev narrowing runs at all;
//!  * **and the NaN quantile**, which answers the number of trials rather than throwing.

use std::io::Read;

use jmath::binomial::{
    cumulative_probability, inverse_cumulative_probability, numerical_mean, numerical_variance,
    BinomialError,
};

fn golden() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/binomial_inverse.txt.gz");
    let file = std::fs::File::open(&path).expect("golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("golden is gzip");
    text
}

fn from_bits(field: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(field, 16).expect("sixteen hex digits"))
}

fn to_bits(value: f64) -> String {
    format!("{:016x}", value.to_bits())
}

#[test]
fn every_quantile_matches_the_golden() {
    let text = golden();
    let mut rows = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("inverse\t") else {
            continue;
        };
        // `unexpected-nan=10`: the NaN row, which the dump prints as an answer because the
        // reference does not refuse it.
        if rest.starts_with("unexpected-") {
            let (label, expected) = rest.split_once('=').expect("an answer");
            assert_eq!(label, "unexpected-nan");
            let answer = inverse_cumulative_probability(10, 0.5, f64::NAN)
                .expect("a NaN is not out of range");
            assert_eq!(answer.to_string(), expected, "the NaN quantile");
            rows += 1;
            continue;
        }
        let (inputs, expected) = rest.split_once('=').expect("an answer");
        let fields: Vec<&str> = inputs.split(',').collect();
        let trials: i32 = fields[0].parse().expect("a trial count");
        let probability = from_bits(fields[1]);
        let quantile = from_bits(fields[2]);
        let answer = inverse_cumulative_probability(trials, probability, quantile)
            .expect("the inputs are in range");
        assert_eq!(
            answer.to_string(),
            expected,
            "inverse({trials}, {probability}, {quantile})"
        );
        rows += 1;
    }
    assert_eq!(rows, 721, "the golden's quantile rows");
}

#[test]
fn every_cumulative_probability_matches_the_golden() {
    let text = golden();
    let mut rows = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("cumulative\t") else {
            continue;
        };
        let (inputs, expected) = rest.split_once('=').expect("a value");
        let fields: Vec<&str> = inputs.split(',').collect();
        let trials: i32 = fields[0].parse().expect("a trial count");
        let probability = from_bits(fields[1]);
        let x: i32 = fields[2].parse().expect("an x");
        let value =
            cumulative_probability(trials, probability, x).expect("the inputs are in range");
        assert_eq!(
            to_bits(value),
            expected,
            "cdf({trials}, {probability}, {x})"
        );
        rows += 1;
    }
    assert!(rows > 0, "the golden carries cumulative rows");
}

#[test]
fn the_moments_match_the_golden() {
    let text = golden();
    let mut rows = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("moments\t") else {
            continue;
        };
        let (inputs, expected) = rest.split_once('=').expect("two values");
        let fields: Vec<&str> = inputs.split(',').collect();
        let trials: i32 = fields[0].parse().expect("a trial count");
        let probability = from_bits(fields[1]);
        let (mean, variance) = expected.split_once(',').expect("two values");
        assert_eq!(to_bits(numerical_mean(trials, probability)), mean);
        assert_eq!(
            to_bits(numerical_variance(trials, probability)),
            variance,
            "variance({trials}, {probability})"
        );
        rows += 1;
    }
    assert!(rows > 0, "the golden carries moment rows");
}

/// The two refusals, and the one that is not.
#[test]
fn the_range_check_refuses_only_what_the_reference_refuses() {
    let text = golden();
    for (label, quantile) in [("below-zero", -0.1), ("above-one", 1.1)] {
        let error = inverse_cumulative_probability(10, 0.5, quantile).expect_err("out of range");
        assert_eq!(error, BinomialError::ProbabilityOutOfRange(quantile));
        assert!(
            text.lines()
                .any(|line| line.starts_with(&format!("error\t{label}\t"))),
            "the golden carries {label}"
        );
    }
    // A NaN passes the range check, because both of its comparisons are false.
    assert_eq!(
        inverse_cumulative_probability(10, 0.5, f64::NAN).expect("not refused"),
        10
    );
}

/// The degenerate probabilities skip the narrowing, because their variance is zero.
#[test]
fn a_zero_variance_skips_the_narrowing() {
    assert_eq!(numerical_variance(10, 0.0), 0.0);
    assert_eq!(numerical_variance(10, 1.0), 0.0);
    // Every mass sits at zero, so any quantile above zero answers zero.
    assert_eq!(
        inverse_cumulative_probability(10, 0.0, 0.5).expect("valid"),
        0
    );
    // And at one, every mass sits at the trial count.
    assert_eq!(
        inverse_cumulative_probability(10, 1.0, 0.5).expect("valid"),
        10
    );
}

/// `p = 0` answers the support's lower bound before any search, and `p = 1` its upper.
#[test]
fn the_two_endpoints_answer_before_the_search() {
    assert_eq!(
        inverse_cumulative_probability(30, 0.5, 0.0).expect("valid"),
        0
    );
    assert_eq!(
        inverse_cumulative_probability(30, 0.5, 1.0).expect("valid"),
        30
    );
}
