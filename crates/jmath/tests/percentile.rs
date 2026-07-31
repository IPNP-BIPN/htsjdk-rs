//! Conformance for commons-math3 `Percentile`, `Median` and `FastMath.round`, against the oracle.
//!
//! Golden from `tools/jmath-conformance/PercentileDump.java`, produced against the pinned
//! commons-math3 **3.5**, which is the version GATK declares `strictly`.
//!
//! Every double travels as raw bits, so a last-bit difference cannot be lost to formatting.
//! The rows that carry the claim:
//!
//! ```text
//! median  LEGACY  1,2      2      interpolation, then FastMath.round of 1.5
//! median  R_1     1,2      1      the neighbour, not the interpolation
//! round   <0.49999999999999994>  1  0     FastMath.round and Math.round, one apart
//! ```
//!
//! # The infinities produce a NaN, and its sign is the FPU's
//!
//! `lower + dif * (upper - lower)` over `{-inf, +inf}` is `-inf + 0.5 * inf`, an invalid
//! operation, and over `{inf, inf}` it is `inf + 0.5 * (inf - inf)`. Both are NaN, and x86-64
//! returns the NaN with the sign bit set while AArch64 returns it clear. That is decision 0012,
//! reached here by a second route: the port is exact on the oracle's architecture and exempted,
//! countably, anywhere else.

use std::io::Read;

use jmath::percentile::{evaluate, median_of_ints, EstimationType};

fn golden() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/percentile.txt.gz");
    let file = std::fs::File::open(&path).expect("golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("golden is gzip");
    text
}

/// The rows whose answer is an invalid operation's NaN, counted from the golden. Every one is
/// `{-inf, +inf}`, `{-inf, 0, +inf}` or `{+inf, +inf}` under some quantile.
const EXPECTED_NAN_SIGN_EXEMPTIONS: usize = 16;

/// Decision 0012: the two differ in the sign bit of a NaN and in nothing else.
fn is_nan_sign_only(a: f64, b: f64) -> bool {
    a.is_nan() && b.is_nan() && (a.to_bits() ^ b.to_bits()) == 1 << 63
}

fn estimation(label: &str) -> EstimationType {
    match label {
        "LEGACY" => EstimationType::Legacy,
        "R_1" => EstimationType::R1,
        other => panic!("unknown estimation type {other}"),
    }
}

/// The dump writes doubles as `Double.doubleToRawLongBits`, and an empty field is an empty array.
fn doubles(field: &str) -> Vec<f64> {
    if field.is_empty() {
        return Vec::new();
    }
    field
        .split(',')
        .map(|b| f64::from_bits(b.parse::<i64>().expect("bits") as u64))
        .collect()
}

fn ints(field: &str) -> Vec<i32> {
    if field.is_empty() {
        return Vec::new();
    }
    field
        .split(',')
        .map(|v| v.parse::<i32>().expect("an int"))
        .collect()
}

#[test]
fn every_percentile_is_bit_identical_to_the_reference() {
    let text = golden();
    let mut count = 0;
    let mut nan_sign_exemptions = 0usize;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("percentile\t") else {
            continue;
        };
        let mut fields = rest.split('\t');
        let type_label = fields.next().expect("a type");
        let quantile: f64 = fields
            .next()
            .expect("a quantile")
            .parse()
            .expect("a number");
        let input = doubles(fields.next().expect("an input"));
        let expected = fields.next().expect("a result");
        // No row of this golden is an exception: `Percentile` answers NaN where a stricter API
        // would throw, including for the empty array. If one ever throws, this says so rather
        // than parsing the message as a number.
        assert!(
            !expected.starts_with("E:"),
            "the reference threw {expected} on {input:?}, which the port does not model"
        );
        let ours = evaluate(&input, quantile, estimation(type_label));
        let want = f64::from_bits(expected.parse::<i64>().expect("bits") as u64);
        if ours.to_bits() != want.to_bits() {
            // Decision 0012: an invalid operation's NaN carries the sign the FPU chose, and the
            // two architectures choose differently. Exempted only off x86-64, and counted.
            if is_nan_sign_only(ours, want) && !cfg!(target_arch = "x86_64") {
                nan_sign_exemptions += 1;
            } else {
                panic!(
                    "{type_label} p={quantile} over {input:?}: ours={:016x} reference={:016x}",
                    ours.to_bits(),
                    want.to_bits()
                );
            }
        }
        count += 1;
    }
    assert!(count > 0, "the golden carries no percentile rows");
    if cfg!(target_arch = "x86_64") {
        assert_eq!(
            nan_sign_exemptions, 0,
            "on x86-64 there is nothing to exempt; the FPU produces the same NaN as the oracle"
        );
    } else {
        // A number rather than a list, because every one of them is the same operation on a
        // different input. It is asserted exactly so the exemption cannot quietly grow.
        assert_eq!(
            nan_sign_exemptions, EXPECTED_NAN_SIGN_EXEMPTIONS,
            "the NaN-sign exemption count changed; see decision 0012"
        );
    }
    println!("{count} percentiles compared, {nan_sign_exemptions} NaN-sign exemptions");
}

#[test]
fn every_median_of_ints_matches_the_reference() {
    let text = golden();
    let mut count = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("median\t") else {
            continue;
        };
        let mut fields = rest.split('\t');
        let type_label = fields.next().expect("a type");
        let input = ints(fields.next().expect("an input"));
        let expected: i32 = fields.next().expect("a result").parse().expect("an int");
        assert_eq!(
            median_of_ints(&input, estimation(type_label)),
            expected,
            "median({input:?}) under {type_label}"
        );
        count += 1;
    }
    assert!(count > 0, "the golden carries no median rows");
}

/// `FastMath.round` against `Math.round`, on the same inputs, in the same golden.
///
/// The two are one apart on `0.49999999999999994`, and both are ported: the call site decides
/// which one is correct, and this is the evidence that they are not interchangeable.
#[test]
fn both_roundings_match_their_own_reference_and_disagree_where_expected() {
    let text = golden();
    let mut disagreements = 0;
    let mut count = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("round\t") else {
            continue;
        };
        let mut fields = rest.split('\t');
        let x = f64::from_bits(fields.next().expect("bits").parse::<i64>().expect("an int") as u64);
        let fast: i64 = fields.next().expect("FastMath").parse().expect("a long");
        let math: i64 = fields.next().expect("Math").parse().expect("a long");
        assert_eq!(jmath::fast_math::round(x), fast, "FastMath.round({x})");
        assert_eq!(jmath::math::round(x), math, "Math.round({x})");
        if fast != math {
            disagreements += 1;
        }
        count += 1;
    }
    assert!(count > 0, "the golden carries no round rows");
    assert!(
        disagreements > 0,
        "the two roundings agreed everywhere, which means the witness input left the golden"
    );
    println!("{count} roundings compared, {disagreements} where the two definitions differ");
}
