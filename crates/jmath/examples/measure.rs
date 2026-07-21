//! Measures how far Rust's libm is from the reference JVM, per function.
//! This is the R1 experiment: it decides how much of `jmath` must be hand-ported.
use std::collections::BTreeMap;
use std::io::{BufRead, BufReader};

#[derive(Default, Clone)]
struct Stat {
    n: u64,
    rust_eq_math: u64,
    rust_eq_strict: u64,
    rust_eq_fast: u64,
    math_eq_strict: u64,
    math_eq_fast: u64,
    max_ulp_vs_math: u64,
}

fn ulp_diff(a: f64, b: f64) -> u64 {
    if a.is_nan() && b.is_nan() {
        return 0;
    }
    if a == b {
        return 0;
    }
    if !a.is_finite() || !b.is_finite() {
        return u64::MAX;
    }
    let ai = a.to_bits() as i64;
    let bi = b.to_bits() as i64;
    let ord = |i: i64| if i < 0 { i64::MIN - i } else { i };
    (ord(ai) - ord(bi)).unsigned_abs()
}

fn apply(f: &str, x: f64, y: f64) -> Option<f64> {
    Some(match f {
        "exp" => x.exp(),
        "log" => x.ln(),
        "log10" => x.log10(),
        "log1p" => x.ln_1p(),
        "expm1" => x.exp_m1(),
        "sqrt" => x.sqrt(),
        "cbrt" => x.cbrt(),
        "sin" => x.sin(),
        "cos" => x.cos(),
        "pow" => x.powf(y),
        _ => return None,
    })
}

fn main() {
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| format!("{}/tests/data/jmath.csv", env!("CARGO_MANIFEST_DIR")));
    let f = std::fs::File::open(&path).expect("open corpus");
    let mut stats: BTreeMap<String, Stat> = BTreeMap::new();

    for line in BufReader::new(f).lines() {
        let line = line.unwrap();
        if line.starts_with('#') {
            continue;
        }
        let mut it = line.split(',');
        let (fname, inp, mb, sb, fb) = match (it.next(), it.next(), it.next(), it.next(), it.next())
        {
            (Some(a), Some(b), Some(c), Some(d), Some(e)) => (a, b, c, d, e),
            _ => continue,
        };
        let (x, y) = if let Some((a, b)) = inp.split_once(':') {
            (
                f64::from_bits(u64::from_str_radix(a, 16).unwrap()),
                f64::from_bits(u64::from_str_radix(b, 16).unwrap()),
            )
        } else {
            (
                f64::from_bits(u64::from_str_radix(inp, 16).unwrap()),
                f64::NAN,
            )
        };
        let math = f64::from_bits(u64::from_str_radix(mb, 16).unwrap());
        let strict = f64::from_bits(u64::from_str_radix(sb, 16).unwrap());
        let fast = f64::from_bits(u64::from_str_radix(fb, 16).unwrap());
        let Some(rust) = apply(fname, x, y) else {
            continue;
        };

        let s = stats.entry(fname.to_string()).or_default();
        s.n += 1;
        let same = |a: f64, b: f64| a.to_bits() == b.to_bits() || (a.is_nan() && b.is_nan());
        if same(rust, math) {
            s.rust_eq_math += 1;
        }
        if same(rust, strict) {
            s.rust_eq_strict += 1;
        }
        if same(rust, fast) {
            s.rust_eq_fast += 1;
        }
        if same(math, strict) {
            s.math_eq_strict += 1;
        }
        if same(math, fast) {
            s.math_eq_fast += 1;
        }
        let d = ulp_diff(rust, math);
        if d != u64::MAX && d > s.max_ulp_vs_math {
            s.max_ulp_vs_math = d;
        }
    }

    println!(
        "{:<8} {:>9} {:>11} {:>11} {:>11} {:>11} {:>11} {:>7}",
        "fn",
        "points",
        "rust=Math",
        "rust=Strict",
        "rust=Fast",
        "Math=Strict",
        "Math=Fast",
        "maxULP"
    );
    let pct = |a: u64, n: u64| format!("{:.4}%", 100.0 * a as f64 / n.max(1) as f64);
    for (k, s) in &stats {
        println!(
            "{:<8} {:>9} {:>11} {:>11} {:>11} {:>11} {:>11} {:>7}",
            k,
            s.n,
            pct(s.rust_eq_math, s.n),
            pct(s.rust_eq_strict, s.n),
            pct(s.rust_eq_fast, s.n),
            pct(s.math_eq_strict, s.n),
            pct(s.math_eq_fast, s.n),
            s.max_ulp_vs_math
        );
    }
}
