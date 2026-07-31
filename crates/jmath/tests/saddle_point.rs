//! Conformance for `SaddlePointExpansion` and the hypergeometric log probability, against the
//! oracle.
//!
//! Golden from `tools/jmath-conformance/SaddlePointDump.java`, against the pinned commons-math3
//! **3.5**, reaching the package-private class by reflection.
//!
//! This is the arithmetic under `FS`, the Fisher-strand annotation.

use std::io::Read;

use jmath::saddle_point::{
    deviance_part, hypergeometric_log_probability, log_binomial_probability, stirling_error,
};

fn golden() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/saddle_point.txt.gz");
    let file = std::fs::File::open(&path).expect("golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("golden is gzip");
    text
}

/// The rows whose answer is an invalid operation's NaN. `getDeviancePart(0, mu)` divides zero by
/// zero inside its series, and `getStirlingError(NaN)` propagates one.
const EXPECTED_NAN_SIGN_EXEMPTIONS: usize = 2;

fn from_bits(field: &str) -> f64 {
    f64::from_bits(field.parse::<i64>().expect("raw bits") as u64)
}

#[test]
fn every_answer_is_bit_identical_to_the_reference() {
    let text = golden();
    let mut counts = std::collections::BTreeMap::new();
    let mut nan_sign_exemptions = 0usize;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("saddle\t") else {
            continue;
        };
        let mut fields = rest.split('\t');
        let name = fields.next().expect("a function");
        let args: Vec<&str> = fields.next().expect("the arguments").split(',').collect();
        let expected = from_bits(fields.next().expect("a result"));
        let ours = match name {
            "getStirlingError" => stirling_error(from_bits(args[0])),
            "getDeviancePart" => deviance_part(from_bits(args[0]), from_bits(args[1])),
            "logBinomialProbability" => log_binomial_probability(
                args[0].parse().expect("x"),
                args[1].parse().expect("n"),
                from_bits(args[2]),
                from_bits(args[3]),
            ),
            "logProbability" => hypergeometric_log_probability(
                args[0].parse().expect("population"),
                args[1].parse().expect("successes"),
                args[2].parse().expect("sample"),
                args[3].parse().expect("x"),
            ),
            other => panic!("unknown function {other}"),
        };
        if ours.to_bits() != expected.to_bits() {
            // Decision 0012: an invalid operation's NaN carries the sign the FPU chose.
            if ours.is_nan() && expected.is_nan() && !cfg!(target_arch = "x86_64") {
                nan_sign_exemptions += 1;
            } else {
                panic!(
                    "{name}({}) = {ours:e}, reference {expected:e}",
                    args.join(",")
                );
            }
        }
        *counts.entry(name.to_string()).or_insert(0u32) += 1;
    }
    assert!(!counts.is_empty(), "the golden carries no rows");
    if cfg!(target_arch = "x86_64") {
        assert_eq!(
            nan_sign_exemptions, 0,
            "on x86-64 there is nothing to exempt; the FPU produces the same NaN as the oracle"
        );
    } else {
        assert_eq!(
            nan_sign_exemptions, EXPECTED_NAN_SIGN_EXEMPTIONS,
            "the NaN-sign exemption count changed; see decision 0012"
        );
    }
    for (name, count) in &counts {
        println!("{count} {name} answers identical");
    }
    println!("{nan_sign_exemptions} NaN-sign exemptions");
}

/// The table's last entry is unreachable, and the discontinuity it hides is real.
#[test]
fn fifteen_takes_the_series_rather_than_the_table() {
    let table = stirling_error(14.5);
    let series = stirling_error(15.0);
    // Both are near 0.0055; the point is that 15.0 does not read the entry tabulated for it.
    assert_ne!(
        series.to_bits(),
        jmath::saddle_point::EXACT_STIRLING_ERROR_AT_FIFTEEN.to_bits()
    );
    assert!((series - jmath::saddle_point::EXACT_STIRLING_ERROR_AT_FIFTEEN).abs() < 1e-15);
    assert!(table > series);
}
