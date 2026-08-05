//! Conformance for `FastMath.exp` and its four tables, against the oracle.
//!
//! Golden from `tools/jmath-conformance/FastMathTablesDump.java`, against the pinned
//! commons-math3 **3.5**, whose tables are read out of the JVM by reflection.
//!
//! Three claims:
//!
//! - every one of the 5,050 table entries the port carries equals the literal the reference ships;
//! - the reference's **other** branch, `RECOMPUTE_TABLES_AT_RUNTIME`, produces different tables:
//!   417 entries differ, which is decision 0024 and the reason the port carries literals;
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

/// The bits of a double, with every NaN collapsed to the canonical one.
///
/// This test used to compare raw bits and carry an architecture-conditional exemption for the
/// sign of a NaN: 160 entries off x86-64, none on it, because x86's default quiet NaN is the
/// negative one and the golden had been taken there. That reasoning was wrong in a way only a
/// second runner could show. The reference's `inf / inf` reaches the JIT, not only the FPU, and a
/// later CI host produced the **positive** NaN for the same entries on the same architecture. The
/// sign is not a property of the arithmetic, of the CPU, or of the target: it is not a property of
/// anything a port could reproduce.
///
/// So it is not compared. `f64::NAN` and Java's `Double.doubleToLongBits` agree on the canonical
/// form, and every finite value and both infinities still travel bit for bit.
fn canonical_bits(value: f64) -> u64 {
    if value.is_nan() {
        f64::NAN.to_bits()
    } else {
        value.to_bits()
    }
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
        if name.starts_with("RECOMPUTED_") || name.contains("LN_MANT") {
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
    assert_eq!(
        count, 5050,
        "the golden should carry every exponential table entry"
    );
    println!("{count} table entries identical to the shipped literals");
}

/// Decision 0024: the reference's two ways of obtaining the same tables disagree.
///
/// The count is asserted exactly, in both directions. If a future commons-math3 made the branches
/// agree, this fails and the decision can be retired; if the port's `FastMathCalc` transcription
/// drifted, it fails too, and the golden says which entries moved.
#[test]
fn the_two_branches_of_the_reference_disagree_on_exactly_417_entries() {
    let text = golden();
    let ours = jmath::fast_math::recomputed_exp_tables();
    let mut compared = 0;
    let mut differ_from_literals = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("table\t") else {
            continue;
        };
        let mut fields = rest.split('\t');
        let name = fields.next().expect("a table name");
        let Some(bare) = name.strip_prefix("RECOMPUTED_") else {
            continue;
        };
        if bare.contains("LN_MANT") {
            continue;
        }
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
        // The port's recomputation matches the reference's, entry for entry, once the sign of a
        // NaN is out of the comparison. Nothing else is exempted.
        assert_eq!(
            canonical_bits(recomputed),
            canonical_bits(expected),
            "recomputed {bare}[{index}]"
        );
        // Counted against the reference's own recomputed column, so the count is the reference's
        // disagreement with itself rather than this architecture's NaN signs.
        if canonical_bits(expected) != canonical_bits(literal) {
            differ_from_literals += 1;
        }
        compared += 1;
    }
    assert_eq!(
        compared, 5050,
        "the golden should carry both columns in full"
    );
    // 417, not the 577 this test asserted before. The other 160 were the same NaN with two
    // different signs on the two sides, counted as a disagreement by a raw-bits comparison. The
    // reference's two branches disagree on 417 entries; the remaining 160 they agree on, in the
    // only sense in which a NaN can be agreed on.
    assert_eq!(
        differ_from_literals, 417,
        "the number of entries on which the reference's two branches disagree changed"
    );
    println!("{differ_from_literals} of {compared} entries differ between the two branches");
}

#[test]
fn every_logarithm_is_bit_identical_to_the_reference() {
    let text = golden();
    let mut count = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("log\t") else {
            continue;
        };
        let (input, expected) = rest.split_once('\t').expect("an input and a result");
        let x = from_bits(input);
        let want = from_bits(expected);
        assert_eq!(
            jmath::fast_math::log(x).to_bits(),
            want.to_bits(),
            "FastMath.log({x:e})"
        );
        count += 1;
    }
    assert!(count > 0, "the golden carries no log rows");
    println!("{count} logarithms bit-identical");
}

/// `LN_MANT`, the logarithm's own table, in both of the reference's forms.
#[test]
fn the_logarithm_table_matches_and_its_two_branches_agree() {
    let text = golden();
    let mut compared = 0;
    let mut differ = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("table\t") else {
            continue;
        };
        let mut fields = rest.split('\t');
        let name = fields.next().expect("a table name");
        if !name.contains("LN_MANT") {
            continue;
        }
        let index: usize = fields.next().expect("an index").parse().expect("a number");
        let expected = from_bits(fields.next().expect("the value"));
        let ours = match name {
            "LN_MANT_A" => f64::from_bits(jmath::fast_math_tables::LN_MANT_A[index]),
            "LN_MANT_B" => f64::from_bits(jmath::fast_math_tables::LN_MANT_B[index]),
            // The recomputed column is compared against the literals rather than against the
            // port, since the port does not carry a `slowLog`: what matters is whether the
            // reference's two branches agree here as they failed to for the exponential.
            "RECOMPUTED_LN_MANT_A" => {
                let literal = jmath::fast_math_tables::LN_MANT_A[index];
                if expected.to_bits() != literal {
                    differ += 1;
                }
                compared += 1;
                continue;
            }
            "RECOMPUTED_LN_MANT_B" => {
                let literal = jmath::fast_math_tables::LN_MANT_B[index];
                if expected.to_bits() != literal {
                    differ += 1;
                }
                compared += 1;
                continue;
            }
            other => panic!("unknown table {other}"),
        };
        assert_eq!(ours.to_bits(), expected.to_bits(), "{name}[{index}]");
        compared += 1;
    }
    assert_eq!(
        compared, 4096,
        "the golden should carry LN_MANT in both forms"
    );
    // Unlike the exponential's integer table, this one agrees with its recomputation: the
    // disagreement of decision 0024 is specific to the reciprocal path.
    assert_eq!(
        differ, 0,
        "LN_MANT's two branches now disagree on {differ} entries; see decision 0024"
    );
    println!("{compared} LN_MANT entries compared, {differ} branch disagreements");
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

/// The whole `jmath.csv` corpus, whose `FastMath` column is what a ported call site will reach,
/// for both functions ported here.
///
/// The corpus is generated against 3.6.1 and the golden above against 3.5, so a divergence here
/// with the golden green would mean the two versions differ, which is worth failing over.
#[test]
fn the_corpus_agrees_with_the_port_on_every_exponential_and_logarithm() {
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
        let ours = match fields[0] {
            "exp" => jmath::fast_math::exp(hex(fields[1])),
            "log" => jmath::fast_math::log(hex(fields[1])),
            _ => continue,
        };
        total += 1;
        if ours.to_bits() == hex(fields[4]).to_bits() {
            matched += 1;
        }
    }
    assert!(total > 80_000, "the corpus lost rows: {total}");
    assert_eq!(
        matched,
        total,
        "FastMath diverges on {} of {total} corpus points",
        total - matched
    );
    println!("{matched}/{total} corpus exponentials and logarithms bit-identical");
}
