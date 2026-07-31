//! Conformance for `FastMath.exp` and its four tables, against the oracle.
//!
//! Golden from `tools/jmath-conformance/FastMathTablesDump.java`, against the pinned
//! commons-math3 **3.5**, whose tables are read out of the JVM by reflection.
//!
//! Three claims:
//!
//! - every one of the 5,050 table entries the port carries equals the literal the reference ships;
//! - the reference's **other** branch, `RECOMPUTE_TABLES_AT_RUNTIME`, produces different tables:
//!   577 entries differ, which is decision 0024 and the reason the port carries literals;
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

/// The recomputed entries whose value is an invalid operation's NaN.
///
/// The reciprocal path divides by an overflowed exponential, so the far end of the integer table
/// is `inf / inf` and its NaN carries whichever sign the FPU chose. 160 entries off x86-64, none
/// on it.
const EXPECTED_NAN_SIGN_EXEMPTIONS: usize = 160;

/// Decision 0012: the two differ in the sign bit of a NaN and in nothing else.
fn is_nan_sign_only(a: f64, b: f64) -> bool {
    a.is_nan() && b.is_nan() && (a.to_bits() ^ b.to_bits()) == 1 << 63
}

fn from_bits(field: &str) -> f64 {
    f64::from_bits(field.parse::<i64>().expect("raw bits") as u64)
}

#[test]
fn every_table_entry_equals_the_literal_the_reference_ships() {
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
        if name.starts_with("RECOMPUTED_") {
            continue;
        }
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
            "{name}[{index}]: ours {ours:e}, reference {expected:e}"
        );
        count += 1;
    }
    assert_eq!(count, 5050, "the golden should carry every table entry");
    println!("{count} table entries identical to the shipped literals");
}

/// Decision 0024: the reference's two ways of obtaining the same tables disagree.
///
/// The count is asserted exactly, in both directions. If a future commons-math3 made the branches
/// agree, this fails and the decision can be retired; if the port's `FastMathCalc` transcription
/// drifted, it fails too, and the golden says which entries moved.
#[test]
fn the_two_branches_of_the_reference_disagree_on_exactly_577_entries() {
    let text = golden();
    let ours = jmath::fast_math::recomputed_exp_tables();
    let mut compared = 0;
    let mut differ_from_literals = 0;
    let mut nan_sign_exemptions = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("table\t") else {
            continue;
        };
        let mut fields = rest.split('\t');
        let name = fields.next().expect("a table name");
        let Some(bare) = name.strip_prefix("RECOMPUTED_") else {
            continue;
        };
        let index: usize = fields.next().expect("an index").parse().expect("a number");
        let expected = from_bits(fields.next().expect("the value"));
        let (recomputed, literal) = match bare {
            "EXP_INT_TABLE_A" => (
                ours.exp_int_a[index],
                f64::from_bits(jmath::fast_math_tables::EXP_INT_A[index]),
            ),
            "EXP_INT_TABLE_B" => (
                ours.exp_int_b[index],
                f64::from_bits(jmath::fast_math_tables::EXP_INT_B[index]),
            ),
            "EXP_FRAC_TABLE_A" => (
                ours.exp_frac_a[index],
                f64::from_bits(jmath::fast_math_tables::EXP_FRAC_A[index]),
            ),
            "EXP_FRAC_TABLE_B" => (
                ours.exp_frac_b[index],
                f64::from_bits(jmath::fast_math_tables::EXP_FRAC_B[index]),
            ),
            other => panic!("unknown table {other}"),
        };
        // The port's recomputation matches the reference's recomputation, entry for entry, up to
        // the sign of a NaN. Overflowing entries divide infinity by infinity, and which NaN that
        // produces is the FPU's choice: decision 0012, reached here by a third route.
        if recomputed.to_bits() != expected.to_bits() {
            if is_nan_sign_only(recomputed, expected) && !cfg!(target_arch = "x86_64") {
                nan_sign_exemptions += 1;
            } else {
                panic!(
                    "recomputed {bare}[{index}]: ours {:016x}, reference {:016x}",
                    recomputed.to_bits(),
                    expected.to_bits()
                );
            }
        }
        // Counted against the reference's own recomputed column, so the count is the reference's
        // disagreement with itself rather than this architecture's NaN signs.
        if expected.to_bits() != literal.to_bits() {
            differ_from_literals += 1;
        }
        compared += 1;
    }
    assert_eq!(
        compared, 5050,
        "the golden should carry both columns in full"
    );
    assert_eq!(
        differ_from_literals, 577,
        "the number of entries on which the reference's two branches disagree changed"
    );
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
        "{differ_from_literals} of {compared} entries differ between the two branches, \
         {nan_sign_exemptions} NaN-sign exemptions"
    );
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
        matched,
        total,
        "FastMath.exp diverges on {} of {total} corpus points",
        total - matched
    );
    println!("{matched}/{total} corpus exponentials bit-identical");
}
