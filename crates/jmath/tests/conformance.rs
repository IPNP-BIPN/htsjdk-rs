//! Conformance against the reference JVM, from the corpus in `tests/data/jmath.csv.gz`.
//!
//! Only `sqrt` is implemented so far and it is asserted at 100%. The remaining functions are
//! measured rather than asserted, and the test *fails* if any of them silently reaches 100%
//! while still being backed by Rust's libm: that would mean the baseline moved and decision
//! 0005's numbers no longer describe reality.

use std::collections::BTreeMap;
use std::io::{BufRead, BufReader};

fn corpus() -> impl BufRead {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/jmath.csv.gz");
    let f = std::fs::File::open(&p).unwrap_or_else(|e| panic!("open {}: {e}", p.display()));
    BufReader::new(flate2::read::GzDecoder::new(f))
}

fn bits(s: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(s, 16).expect("hex bits"))
}

fn same(a: f64, b: f64) -> bool {
    a.to_bits() == b.to_bits() || (a.is_nan() && b.is_nan())
}

/// (function, agreements with `Math`, total points)
fn agreement() -> BTreeMap<String, (u64, u64)> {
    let mut acc: BTreeMap<String, (u64, u64)> = BTreeMap::new();
    for line in corpus().lines() {
        let line = line.unwrap();
        if line.starts_with('#') {
            continue;
        }
        let mut it = line.split(',');
        let (Some(f), Some(inp), Some(mb)) = (it.next(), it.next(), it.next()) else {
            continue;
        };
        let (x, y) = match inp.split_once(':') {
            Some((a, b)) => (bits(a), bits(b)),
            None => (bits(inp), f64::NAN),
        };
        let math = bits(mb);
        let ours = match f {
            "sqrt" => jmath::math::sqrt(x),
            "exp" => x.exp(),
            "log" => x.ln(),
            "log10" => x.log10(),
            "log1p" => x.ln_1p(),
            "expm1" => x.exp_m1(),
            "cbrt" => x.cbrt(),
            "sin" => x.sin(),
            "cos" => x.cos(),
            "pow" => x.powf(y),
            _ => continue,
        };
        let e = acc.entry(f.to_string()).or_insert((0, 0));
        e.1 += 1;
        if same(ours, math) {
            e.0 += 1;
        }
    }
    acc
}

#[test]
fn corpus_is_present_and_substantial() {
    let a = agreement();
    let total: u64 = a.values().map(|(_, n)| n).sum();
    assert!(
        total > 800_000,
        "corpus shrank to {total} points; decision 0005 was measured on 809,930"
    );
}

/// `sqrt` is the one function that needs no porting, because IEEE-754 mandates its rounding.
#[test]
fn sqrt_is_bit_identical_to_the_jvm() {
    let a = agreement();
    let (ok, n) = a["sqrt"];
    assert_eq!(ok, n, "sqrt must be exact over all {n} points, got {ok}");
}

/// Records which functions are not yet bit-identical, and reports the live rate.
///
/// The invariant is deliberately "not yet exact" rather than a numeric threshold: a threshold
/// would need updating every time the corpus or the host libm moves, and would fail for the
/// wrong reason. Agreement rates here are high (`log10` reaches 99.9956%) and that is exactly
/// the trap decision 0005 documents: a rate that reads like success is still millions of
/// differing values across a HaplotypeCaller run.
///
/// When one of these is ported, its entry moves to the exact-match list above and the row in
/// decision 0005 is updated. This test failing means a function became exact without anyone
/// recording it.
#[test]
fn unported_functions_are_not_yet_exact() {
    let a = agreement();
    let mut report = Vec::new();
    for f in [
        "exp", "log", "log10", "pow", "log1p", "expm1", "cbrt", "sin", "cos",
    ] {
        let (ok, n) = a[f];
        report.push(format!("{f}={:.4}%", 100.0 * ok as f64 / n as f64));
        assert!(
            ok < n,
            "`{f}` is now bit-identical to the JVM over all {n} points. If it was ported, move \
             it to the exact-match test and update decision 0005."
        );
    }
    println!("agreement with java.lang.Math: {}", report.join("  "));
}
