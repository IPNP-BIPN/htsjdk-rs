//! Conformance for `Gamma`, `Erf` and `NormalDistribution`, against the oracle.
//!
//! Golden from `tools/jmath-conformance/GammaErfNormalDump.java`, against the pinned
//! commons-math3 **3.5**.
//!
//! This is the layer GATK's `MannWhitneyU` stands on, so every rank-sum annotation's last digit
//! is decided here. The rows worth naming:
//!
//! ```text
//! gamma  erf     <40.0>          exactly 1, by the shortcut
//! gamma  erf     <40.0000001>    also exactly 1, by a different route
//! gamma  erfInv  <1.0>           +inf, from a branch added to fix the polynomial's answer
//! gamma  invGamma1pm1  <-0.5>    E:NumberIsTooSmallException
//! ```

use std::io::Read;

use jmath::gamma;
use jmath::normal::NormalDistribution;

fn golden() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/gamma_erf_normal.txt.gz");
    let file = std::fs::File::open(&path).expect("golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("golden is gzip");
    text
}

/// `Gamma.DEFAULT_EPSILON`, which the two-argument form passes, written as the reference writes
/// it: `10e-15`, which is 1e-14 and not 1e-15.
const DEFAULT_EPSILON: f64 = 10e-15;

/// The rows whose answer is an invalid operation's NaN: `erfInv(NaN)` and the two normal-quantile
/// rows that go through it. Three off x86-64, none on it.
const EXPECTED_NAN_SIGN_EXEMPTIONS: usize = 3;

/// Decision 0012: two renderings that are both NaN and differ only in the sign bit.
fn is_nan_sign_only(ours: &str, expected: &str) -> bool {
    let parse = |s: &str| {
        s.parse::<i64>()
            .ok()
            .map(|bits| f64::from_bits(bits as u64))
    };
    match (parse(ours), parse(expected)) {
        (Some(a), Some(b)) => a.is_nan() && b.is_nan() && (a.to_bits() ^ b.to_bits()) == 1 << 63,
        _ => false,
    }
}

fn from_bits(field: &str) -> f64 {
    f64::from_bits(field.parse::<i64>().expect("raw bits") as u64)
}

/// The port's answer, rendered the way the dump renders the reference's.
fn call(name: &str, inputs: &[f64]) -> String {
    let render = |value: f64| (value.to_bits() as i64).to_string();
    match name {
        "logGamma" => render(gamma::log_gamma(inputs[0])),
        "lanczos" => render(gamma::lanczos(inputs[0])),
        "erf" => render(gamma::erf(inputs[0])),
        "erfc" => render(gamma::erfc(inputs[0])),
        "erfInv" => render(gamma::erf_inv(inputs[0])),
        "digamma" => render(gamma::digamma(inputs[0])),
        "trigamma" => render(gamma::trigamma(inputs[0])),
        "invGamma1pm1" => match gamma::inv_gamma1pm1(inputs[0]) {
            Ok(value) => render(value),
            Err(error) => exception(&error),
        },
        "logGamma1p" => match gamma::log_gamma1p(inputs[0]) {
            Ok(value) => render(value),
            Err(error) => exception(&error),
        },
        "regularizedGammaP" => {
            match gamma::regularized_gamma_p(inputs[0], inputs[1], DEFAULT_EPSILON, i32::MAX) {
                Ok(value) => render(value),
                Err(error) => exception(&error),
            }
        }
        "regularizedGammaQ" => {
            match gamma::regularized_gamma_q(inputs[0], inputs[1], DEFAULT_EPSILON, i32::MAX) {
                Ok(value) => render(value),
                Err(error) => exception(&error),
            }
        }
        "cumulativeProbability" => {
            render(NormalDistribution::default().cumulative_probability(inputs[0]))
        }
        "inverseCumulativeProbability" => {
            match NormalDistribution::default().inverse_cumulative_probability(inputs[0]) {
                Some(value) => render(value),
                None => "E:org.apache.commons.math3.exception.OutOfRangeException".to_string(),
            }
        }
        other => panic!("unknown function {other}"),
    }
}

/// The reference's exception class for each refusal the port models.
fn exception(error: &gamma::GammaError) -> String {
    let class = match error {
        gamma::GammaError::TooSmall { .. } => {
            "org.apache.commons.math3.exception.NumberIsTooSmallException"
        }
        gamma::GammaError::TooLarge { .. } => {
            "org.apache.commons.math3.exception.NumberIsTooLargeException"
        }
        gamma::GammaError::MaxCountExceeded { .. } => {
            "org.apache.commons.math3.exception.MaxCountExceededException"
        }
        gamma::GammaError::ContinuedFractionDiverged { .. } => {
            "org.apache.commons.math3.exception.ConvergenceException"
        }
    };
    format!("E:{class}")
}

#[test]
fn every_answer_is_bit_identical_to_the_reference() {
    let text = golden();
    let mut count = 0;
    let mut refusals = 0;
    let mut nan_sign_exemptions = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("gamma\t") else {
            continue;
        };
        let mut fields = rest.split('\t');
        let name = fields.next().expect("a function");
        let inputs: Vec<f64> = fields
            .next()
            .expect("the inputs")
            .split(',')
            .map(from_bits)
            .collect();
        let expected = fields.next().expect("a result");
        if expected.starts_with("E:") {
            refusals += 1;
        }
        let ours = call(name, &inputs);
        if ours != expected {
            // Decision 0012: an invalid operation's NaN carries the sign the FPU chose, and the
            // two architectures choose differently. Exempted only off x86-64, and counted.
            if is_nan_sign_only(&ours, expected) && !cfg!(target_arch = "x86_64") {
                nan_sign_exemptions += 1;
            } else {
                panic!("{name}{inputs:?}: ours {ours}, reference {expected}");
            }
        }
        count += 1;
    }
    assert!(count > 0, "the golden carries no rows");
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
    println!(
        "{count} answers identical, {refusals} refusals, {nan_sign_exemptions} NaN-sign exemptions"
    );
}

/// `erf` and `erfInv` are inverses in mathematics and unrelated in code.
#[test]
fn the_round_trip_through_the_two_error_functions_does_not_return_its_input() {
    let x = 0.75_f64;
    let back = gamma::erf_inv(gamma::erf(x));
    assert!((back - x).abs() < 1e-12, "they agree to within a tolerance");
    assert_ne!(back.to_bits(), x.to_bits(), "and not bit for bit");
}
