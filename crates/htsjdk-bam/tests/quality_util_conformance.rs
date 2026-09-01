//! `QualityUtil` against the reference's own answers.
//!
//! The dump is `tools/quality-util-conformance/QualityUtilDump.java`, and every value in it is a
//! bit pattern rather than a decimal, so what is compared is the double and not its rendering.
//!
//! The golden is committed and re-derived by the `quality-util` suite on every run; the dump can
//! still be overridden with an environment variable while a harness change is being checked.

use std::io::Read;
use std::path::Path;

use htsjdk_bam::quality_util::{
    error_probability_from_phred_score, phred_score_from_error_probability,
    phred_score_from_obs_and_errors,
};

/// How far the port's error-probability table may sit from the reference's, in ulp.
///
/// Two, measured: `1 / StrictMath.pow(10, x)` and `1 / Math.pow(10, x)` differ by two units in the
/// last place at score 25 and agree everywhere else on the first dump taken. The bound is a
/// property of `Math.pow`'s intrinsic (decisions 0007 and 0027), not of this code, so it is stated
/// here and re-measured on every run rather than assumed to be zero.
const TABLE_ULP_BOUND: u64 = 2;

fn bits(value: f64) -> String {
    format!("{:016x}", value.to_bits())
}

fn from_bits(text: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(text, 16).expect("a 16-digit bit pattern"))
}

#[test]
fn every_answer_matches_the_reference() {
    // The golden was produced by the pinned container on real x86-64 and is re-derived on every
    // run; `QUALITY_UTIL_DUMP` still overrides it, which is how a local run checks a change to the
    // harness before CI does.
    let dump = match std::env::var("QUALITY_UTIL_DUMP") {
        Ok(path) => std::fs::read_to_string(path).expect("the dump named by QUALITY_UTIL_DUMP"),
        Err(_) => {
            let golden =
                Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/quality_util.txt.gz");
            let file = std::fs::File::open(&golden).expect("the committed golden");
            let mut text = String::new();
            flate2::read::GzDecoder::new(file)
                .read_to_string(&mut text)
                .expect("the golden decompresses");
            text
        }
    };

    let mut worst_table_ulp: u64 = 0;
    let mut table = 0;
    let mut phred = 0;
    let mut obs = 0;
    for line in dump.lines() {
        let fields: Vec<&str> = line.trim().split('\t').collect();
        match fields.as_slice() {
            ["table", score, expected] => {
                // The table is built by `Math.pow`, which decision 0007 deferred and decision 0027
                // bounded at 1 ulp; the port builds it with `StrictMath.pow`. So this entry is
                // MEASURED against a bound rather than asserted equal, and the distance is
                // reported: a port that silently drifted would still fail, and the one place the
                // two libraries disagree is not dressed up as agreement.
                let score: i32 = score.parse().expect("a phred score");
                let ours = error_probability_from_phred_score(score).expect("0..=100 is in range");
                let theirs = from_bits(expected);
                let distance = (ours.to_bits() as i64 - theirs.to_bits() as i64).unsigned_abs();
                assert!(
                    distance <= TABLE_ULP_BOUND,
                    "table[{score}]: ours={} theirs={expected} distance={distance} ulp",
                    bits(ours)
                );
                worst_table_ulp = worst_table_ulp.max(distance);
                table += 1;
            }
            ["phred", probability, expected] => {
                let ours = phred_score_from_error_probability(from_bits(probability));
                assert_eq!(
                    ours.to_string(),
                    *expected,
                    "phred({probability}) = {}",
                    from_bits(probability)
                );
                phred += 1;
            }
            ["obs", observations, errors, expected] => {
                let ours =
                    phred_score_from_obs_and_errors(from_bits(observations), from_bits(errors));
                assert_eq!(ours.to_string(), *expected, "obs({observations}, {errors})");
                obs += 1;
            }
            _ => panic!("unrecognized dump line: {line}"),
        }
    }
    assert_eq!(table, 101, "the table has 101 entries");
    assert!(phred > 0 && obs > 0, "the sweep and the pairs both ran");
    println!("table={table} (worst {worst_table_ulp} ulp) phred={phred} obs={obs}");
}
