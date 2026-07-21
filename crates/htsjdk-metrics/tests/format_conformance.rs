//! Conformance for `FormatUtil` against htsjdk, over 41,678 values.
//!
//! Goldens from `tools/metrics-conformance/FormatDump.java` in the pinned oracle container
//! under `Locale.US`.
//!
//! **The port is not finished, and this file says so precisely rather than passing.** 41,566 of
//! 41,678 values match, 99.73%. The remaining 112 are listed in
//! `tests/data/known_divergences.tsv` by bit pattern, with the answer each side gives, so that
//! neither an improvement nor a regression can happen silently.
//!
//! Pinning the failures rather than trimming the corpus is the point. A conformance suite that
//! contains only what already passes reports 100% and means nothing.
//!
//! The cause is diagnosed and recorded in decision 0011. It is not a rounding tweak away:
//! `DecimalFormat` rounds the digit string produced by Java 17's `FloatingDecimal`, and
//! `DigitList.shouldRoundUp` consults whether that string is exact and whether it was already
//! rounded up. Both require the exact decimal expansion of the double. Fitting a rule to this
//! corpus instead would pass here and diverge elsewhere, which is the failure mode this whole
//! project exists to avoid.

use std::io::Read;

use htsjdk_metrics::{format_double, format_long};

const KNOWN: &str = include_str!("data/known_divergences.tsv");

fn known() -> Vec<(&'static str, &'static str, &'static str)> {
    KNOWN
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let mut it = l.split('\t');
            (
                it.next().expect("bits"),
                it.next().expect("ours"),
                it.next().expect("htsjdk"),
            )
        })
        .collect()
}

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/format.txt.gz");
    let f = std::fs::File::open(&p).expect("corpus");
    let mut s = String::new();
    flate2::read::GzDecoder::new(f)
        .read_to_string(&mut s)
        .expect("corpus is gzip");
    s
}

fn ours_for(key: &str) -> String {
    if let Some(hex) = key.strip_prefix('L') {
        format_long(u64::from_str_radix(hex, 16).expect("long bits") as i64)
    } else {
        format_double(f64::from_bits(
            u64::from_str_radix(key, 16).expect("double bits"),
        ))
    }
}

/// The gate: every value matches htsjdk, except exactly the ones declared.
#[test]
fn the_only_divergences_are_the_declared_ones() {
    let text = corpus();
    let declared: std::collections::HashSet<&str> = known().iter().map(|(b, _, _)| *b).collect();

    let mut checked = 0usize;
    let mut undeclared = Vec::new();
    let mut newly_fixed = Vec::new();

    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let (key, expected) = line.split_once('\t').expect("bits TAB formatted");
        let ours = ours_for(key);
        checked += 1;

        match (ours == expected, declared.contains(key)) {
            (true, false) | (false, true) => {}
            (false, false) => {
                if undeclared.len() < 20 {
                    undeclared.push(format!("{key}: ours={ours:?} htsjdk={expected:?}"));
                }
            }
            (true, true) => newly_fixed.push(key.to_string()),
        }
    }

    assert_eq!(checked, 41_678, "corpus size changed");
    assert!(
        undeclared.is_empty(),
        "{} new divergences appeared, which is a regression:\n{}",
        undeclared.len(),
        undeclared.join("\n")
    );
    assert!(
        newly_fixed.is_empty(),
        "{} declared divergences now match, which is good news: remove them from \
         tests/data/known_divergences.tsv and update decision 0011.\n{newly_fixed:?}",
        newly_fixed.len()
    );
}

/// Each declared divergence must still produce exactly the recorded wrong answer. One that
/// drifted to a *different* wrong answer would otherwise pass unnoticed.
#[test]
fn each_declared_divergence_is_still_exactly_as_recorded() {
    for (bits, ours_expected, _) in known() {
        assert_eq!(
            ours_for(bits),
            ours_expected,
            "the port's answer for {bits} changed without the declaration being updated"
        );
    }
}

/// Keeps the rate quoted in decision 0011 tied to the code.
#[test]
fn the_measured_agreement_rate_matches_the_decision_record() {
    let total = 41_678usize;
    let diverging = known().len();
    assert_eq!(diverging, 112);
    let rate = (total - diverging) as f64 / total as f64;
    assert!(
        (0.9972..0.9974).contains(&rate),
        "decision 0011 states 99.73%, measured {:.4}%",
        rate * 100.0
    );
}

/// The corpus must contain the cases the port claims to handle, or passing proves nothing.
#[test]
fn the_corpus_covers_the_interesting_shapes() {
    let text = corpus();
    let values: Vec<&str> = text
        .lines()
        .filter(|l| !l.starts_with('#') && l.contains('\t'))
        .map(|l| l.split_once('\t').unwrap().1)
        .collect();

    assert!(values.contains(&"?"), "NaN / infinity");
    assert!(values.contains(&"-?"), "negative infinity");
    assert!(values.contains(&"-0"), "negative zero");
    assert!(
        values.iter().any(|v| v.len() > 100),
        "a full-length expansion, proving no scientific notation"
    );
    assert!(
        !values.iter().any(|v| v.contains('e') || v.contains('E')),
        "htsjdk never emits scientific notation here"
    );
    assert!(
        values
            .iter()
            .any(|v| v.split_once('.').is_some_and(|(_, f)| f.len() == 6)),
        "values that actually reach the six-digit limit"
    );
}
