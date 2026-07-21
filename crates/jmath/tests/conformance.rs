//! Conformance against the reference JVM, from the corpus in `tests/data/jmath.csv.gz`.
//!
//! `sqrt`, `log` and `log10` are ported and asserted bit-identical over every point. The rest
//! delegate to Rust's libm and are measured rather than asserted; that test *fails* if one of
//! them silently reaches 100%, which would mean a function was ported without decisions 0005
//! and 0006 being updated to match.
//!
//! `exp` was in the first list and is now in the second. See decision 0014.

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
            "log" => jmath::math::log(x),
            "log10" => jmath::math::log10(x),
            // `exp` is WITHDRAWN, not missing: the port was a transcription of GPL2-only
            // HotSpot source (decision 0014). It stays in the corpus, routed to the system
            // libm, so its divergence rate is measured and reported rather than the function
            // quietly disappearing from the table.
            "exp" => x.exp(),
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

/// Functions that are bit-identical to `java.lang.Math` over the whole corpus.
///
/// Two routes to exactness, and the distinction matters:
///
/// - `sqrt` is free: IEEE-754 mandates its rounding, so every implementation already agrees.
/// - `log` and `log10` are correctly rounded in the reference, so rounding the true result
///   suffices and no algorithm port was needed.
///
/// There was a third route and it is gone. `exp` is *not* correctly rounded, so being exact
/// required reproducing HotSpot's intrinsic operation by operation, which is a transcription of
/// GPL2-only source and could not ship under this crate's MIT licence. Withdrawn in decision
/// 0014.
///
/// See decision 0006.
#[test]
fn ported_functions_are_bit_identical_to_the_jvm() {
    let a = agreement();
    for f in ["sqrt", "log", "log10"] {
        let (ok, n) = a[f];
        assert_eq!(
            ok,
            n,
            "`{f}` must match java.lang.Math on all {n} points, got {ok} ({} divergent)",
            n - ok
        );
    }
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
///
/// `exp` moved in the other direction. It was exact, by transcription of GPL2-only HotSpot
/// source, and was withdrawn in decision 0014 because that transcription could not ship under
/// this crate's MIT licence. Its rate here is the system libm's, and the gap between that rate
/// and 100% is the exact size of what the licence costs.
#[test]
fn unported_functions_are_not_yet_exact() {
    let a = agreement();
    let mut report = Vec::new();
    for f in ["exp", "pow", "log1p", "expm1", "cbrt", "sin", "cos"] {
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
