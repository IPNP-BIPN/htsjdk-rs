//! Conformance for `FormatUtil` against htsjdk, over 41,678 values.
//!
//! Goldens from `tools/metrics-conformance/FormatDump.java` in the pinned oracle container
//! under `Locale.US`.
//!
//! **Every value matches.** It took three steps, and the decision records carry each:
//!
//! - 0011 and 0013 found 112 divergences with two causes, `DigitList.shouldRoundUp` and the
//!   pre-Schubfach `FloatingDecimal` digits, and blocked both on the JDK's GPL-2.0 licence.
//! - 0026 showed the first cause never needed `FloatingDecimal`: whether a digit string is the
//!   double's exact value is arithmetic on the double's own bits. That fixed 44 and left 66, every
//!   one above 2^53.
//! - 0046 reconstructed the Java 17 digits from the oracle's measured behaviour rather than from
//!   its source, in `jmath::to_string`. That fixed the last 66 and introduced none.
//!
//! There is no list of declared divergences any more. A divergence that appears is a regression.

use std::io::Read;

use htsjdk_metrics::{format_double, format_long};

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

/// The gate: every value matches htsjdk.
#[test]
fn every_value_matches() {
    let text = corpus();
    let mut checked = 0usize;
    let mut diverging = Vec::new();
    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let (key, expected) = line.split_once('\t').expect("bits TAB formatted");
        let ours = ours_for(key);
        checked += 1;
        if ours != expected && diverging.len() < 20 {
            diverging.push(format!("{key}: ours={ours:?} htsjdk={expected:?}"));
        }
    }
    assert_eq!(checked, 41_678, "corpus size changed");
    assert!(
        diverging.is_empty(),
        "divergences appeared, which is a regression:\n{}",
        diverging.join("\n")
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
