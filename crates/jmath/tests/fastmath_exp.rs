//! Conformance for `FastMath.exp` and its four tables, against the oracle.
//!
//! Golden from `tools/jmath-conformance/FastMathTablesDump.java`, against the pinned
//! commons-math3 **3.5**, whose tables are read out of the JVM by reflection.
//!
//! Two claims, and the first is what makes the second checkable:
//!
//! - every one of the 3,550 table entries the port **computes** with `FastMathCalc` equals the
//!   literal the reference **ships**. The reference can do either and picks the literals with a
//!   compile-time constant; the port picks the computation, so this is where the two branches are
//!   shown to agree;
//! - `exp` itself is bit-identical on the boundaries its own branches name, including the two
//!   shifted recursions that keep a subnormal result precise.
//!
//! The port was also run against the existing `jmath.csv` corpus's `FastMath` column while it was
//! written: 44,996 of 44,996 points bit-identical. That corpus is generated against 3.6.1 and this
//! golden against 3.5, so the two together also say the versions agree here.

use std::io::Read;

fn golden() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/fastmath_exp.txt.gz");
    let file = std::fs::File::open(&path).expect("golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("golden is gzip");
    text
}

fn from_bits(field: &str) -> f64 {
    f64::from_bits(field.parse::<i64>().expect("raw bits") as u64)
}

#[test]
fn every_table_entry_the_port_computes_equals_the_literal_the_reference_ships() {
    let text = golden();
    let tables = jmath::fast_math::exp_tables();
    let mut count = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("table\t") else {
            continue;
        };
        let mut fields = rest.split('\t');
        let name = fields.next().expect("a table name");
        let index: usize = fields.next().expect("an index").parse().expect("a number");
        let expected = from_bits(fields.next().expect("the value"));
        let ours = match name {
            "EXP_INT_TABLE_A" => tables.exp_int_a[index],
            "EXP_INT_TABLE_B" => tables.exp_int_b[index],
            "EXP_FRAC_TABLE_A" => tables.exp_frac_a[index],
            "EXP_FRAC_TABLE_B" => tables.exp_frac_b[index],
            other => panic!("unknown table {other}"),
        };
        assert_eq!(
            ours.to_bits(),
            expected.to_bits(),
            "{name}[{index}]: computed {ours:e}, reference {expected:e}"
        );
        count += 1;
    }
    assert_eq!(count, 3550, "the golden should carry every table entry");
    println!("{count} table entries identical, all computed rather than transcribed");
}

#[test]
fn every_exp_is_bit_identical_to_the_reference() {
    let text = golden();
    let mut count = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("exp\t") else {
            continue;
        };
        let (input, expected) = rest.split_once('\t').expect("an input and a result");
        let x = from_bits(input);
        let want = from_bits(expected);
        let ours = jmath::fast_math::exp(x);
        assert_eq!(ours.to_bits(), want.to_bits(), "FastMath.exp({x:e})");
        count += 1;
    }
    assert!(count > 0, "the golden carries no exp rows");
    println!("{count} exponentials bit-identical");
}

/// The whole `jmath.csv` corpus, whose `FastMath` column is what a ported call site will reach.
///
/// The corpus is generated against 3.6.1 and the golden above against 3.5, so a divergence here
/// with the golden green would mean the two versions differ, which is worth failing over.
#[test]
fn the_corpus_agrees_with_the_port_on_every_exponential() {
    use std::io::{BufRead, BufReader};
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/jmath.csv.gz");
    let file = std::fs::File::open(path).expect("corpus");
    let reader = BufReader::new(flate2::read::GzDecoder::new(file));
    let hex = |s: &str| f64::from_bits(u64::from_str_radix(s, 16).expect("hex bits"));
    let (mut matched, mut total) = (0u64, 0u64);
    for line in reader.lines() {
        let line = line.expect("a line");
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split(',').collect();
        if fields[0] != "exp" {
            continue;
        }
        total += 1;
        if jmath::fast_math::exp(hex(fields[1])).to_bits() == hex(fields[4]).to_bits() {
            matched += 1;
        }
    }
    assert!(total > 40_000, "the corpus lost its exp rows: {total}");
    assert_eq!(
        matched, total,
        "FastMath.exp diverges on {} of {total} corpus points",
        total - matched
    );
    println!("{matched}/{total} corpus exponentials bit-identical");
}
