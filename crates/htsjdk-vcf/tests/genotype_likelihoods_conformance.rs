//! Conformance for `GL` to `PL` and the double parser under it, against the oracle.
//!
//! Goldens from `tools/vcf-conformance/GenotypeLikelihoodsDump.java`.
//!
//! The parsed doubles are compared as **raw bits**. A decimal rendering of a parse result hides
//! exactly the disagreements this suite exists to find: two parsers that print `1.5` may have
//! produced different doubles, and the PL that comes out of them differs by one.
//!
//! Three rows the golden corrected, all of which the port had wrong and none of which any ordinary
//! genotype would have shown:
//!
//! ```text
//! gl     0,inf,-1                2147483647,0,2147483647
//! gl     <empty>                 E:java.lang.NumberFormatException:empty String
//! round  <0.49999999999999994>   0
//! ```
//!
//! `Math.min` propagates NaN and Rust's `f64::min` returns the other operand, so the infinite
//! element, whose shifted value is `inf - inf`, comes out **zero** while the finite ones saturate.
//! `"".split(",")` is a one-element array holding the empty string, not an empty array, so an empty
//! GL field is a parse failure rather than an empty likelihood vector. And `Math.round` is not
//! `floor(x + 0.5)` despite its own javadoc saying so: on the double just below a half, the
//! addition rounds up to exactly 1.0 and the arithmetic version answers 1 where the correct half-up
//! answer is 0.

use std::io::Read;

use htsjdk_vcf::genotype_likelihoods::{
    gl_field_to_pls, java_round, parse_gl_field, parse_java_double, parse_vcf_double,
};

fn golden() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/genotype_likelihoods.txt.gz");
    let file = std::fs::File::open(&path).expect("golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("golden is gzip");
    text
}

/// The dump's `escape`: everything outside printable ASCII travels as `\uXXXX`.
fn unescape(text: &str) -> String {
    let mut out = String::new();
    let mut chars = text.chars();
    while let Some(c) = chars.next() {
        if c != '\\' {
            out.push(c);
            continue;
        }
        match chars.next() {
            Some('\\') => out.push('\\'),
            Some('u') => {
                let hex: String = (0..4).filter_map(|_| chars.next()).collect();
                let code = u32::from_str_radix(&hex, 16).expect("four hex digits");
                out.push(char::from_u32(code).expect("a character"));
            }
            other => panic!("unknown escape {other:?}"),
        }
    }
    out
}

/// Rows of one kind, as `(input, value)` with the input unescaped.
fn rows(text: &str, kind: &str) -> Vec<(String, String)> {
    text.lines()
        .filter_map(|line| line.strip_prefix(&format!("{kind}\t")))
        .map(|rest| {
            // The input may be empty, so the split has to be positional rather than on content.
            let (input, value) = rest.rsplit_once('\t').expect("an input and a value");
            (unescape(input), value.to_string())
        })
        .collect()
}

/// The golden's error rows for the two parsers are always the same exception; what is compared is
/// whether the parse succeeded and, when it did, the exact bits.
fn expected_double(value: &str) -> Option<i64> {
    if value.starts_with("E:") {
        None
    } else {
        Some(value.parse().expect("raw bits"))
    }
}

#[test]
fn every_double_parses_exactly_as_java_parses_it() {
    let text = golden();
    let cases = rows(&text, "double");
    assert!(!cases.is_empty(), "the golden carries no double rows");

    for (input, expected) in &cases {
        let ours = parse_java_double(input);
        match expected_double(expected) {
            None => assert!(
                ours.is_none(),
                "{input:?}: we parsed {ours:?} where the reference threw"
            ),
            Some(bits) => {
                let ours = ours.unwrap_or_else(|| panic!("{input:?}: we refused a valid double"));
                assert_eq!(
                    ours.to_bits() as i64,
                    bits,
                    "{input:?}: parsed to different bits"
                );
            }
        }
    }
    println!("{} doubles identical bit for bit", cases.len());
}

#[test]
fn every_vcf_double_matches_the_reference() {
    let text = golden();
    let cases = rows(&text, "vcfdouble");
    assert!(!cases.is_empty(), "the golden carries no vcfdouble rows");

    for (input, expected) in &cases {
        let ours = parse_vcf_double(input);
        match expected_double(expected) {
            None => assert!(
                ours.is_none(),
                "{input:?}: we parsed {ours:?} where the reference threw"
            ),
            Some(bits) => {
                let ours = ours.unwrap_or_else(|| panic!("{input:?}: we refused a valid double"));
                assert_eq!(ours.to_bits() as i64, bits, "{input:?}");
            }
        }
    }
    println!("{} VCF doubles identical", cases.len());
}

#[test]
fn every_gl_field_becomes_the_reference_pls() {
    let text = golden();
    let cases = rows(&text, "gl");
    assert!(!cases.is_empty(), "the golden carries no gl rows");

    for (input, expected) in &cases {
        let ours = match gl_field_to_pls(input) {
            Ok(None) => "null".to_string(),
            Ok(Some(pls)) => pls
                .iter()
                .map(|pl| pl.to_string())
                .collect::<Vec<_>>()
                .join(","),
            Err(error) => format!("E:{}:{}", error.class(), error.message()),
        };
        assert_eq!(&ours, expected, "GL field {input:?}");
    }
    println!("{} GL fields identical", cases.len());
}

/// The likelihoods themselves, before the conversion: a divergence here is a parse divergence and
/// a divergence in the PLs alone is a rounding one.
#[test]
fn every_gl_field_parses_to_the_reference_likelihoods() {
    let text = golden();
    let cases = rows(&text, "glbits");
    assert!(!cases.is_empty(), "the golden carries no glbits rows");

    for (input, expected) in &cases {
        let ours = match parse_gl_field(input) {
            Ok(None) => "null".to_string(),
            Ok(Some(values)) => values
                .iter()
                .map(|value| (value.to_bits() as i64).to_string())
                .collect::<Vec<_>>()
                .join(","),
            // A field the reference refuses has no `glbits` row at all, so reaching here means the
            // golden and the port disagree about whether it parses.
            Err(error) => panic!("{input:?}: we refused a field the reference parsed: {error:?}"),
        };
        assert_eq!(&ours, expected, "GL field {input:?}");
    }
}

#[test]
fn every_rounding_matches_the_reference() {
    let text = golden();
    let mut count = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("round\t") else {
            continue;
        };
        let (bits, expected) = rest.split_once('\t').expect("bits and a result");
        let value = f64::from_bits(bits.parse::<i64>().expect("raw bits") as u64);
        let expected: i64 = expected.parse().expect("a long");
        assert_eq!(java_round(value), expected, "Math.round({value})");
        count += 1;
    }
    assert!(count > 0, "the golden carries no round rows");
    println!("{count} roundings identical");
}

/// The rows that a port gets wrong by writing the obvious Rust, or by following the javadoc.
#[test]
fn the_rows_that_rust_gets_wrong_by_default() {
    let text = golden();
    let gl = |input: &str| -> String {
        rows(&text, "gl")
            .into_iter()
            .find(|(key, _)| key == input)
            .unwrap_or_else(|| panic!("no gl row for {input:?}"))
            .1
    };

    // `Math.min` propagates NaN; `f64::min` returns the other operand. The infinite element is the
    // one that comes out zero.
    assert_eq!(gl("0,inf,-1"), "2147483647,0,2147483647");
    // `"".split(",")` holds one empty string; `",".split(",")` holds nothing.
    assert_eq!(gl(""), "E:java.lang.NumberFormatException:empty String");
    assert_eq!(gl(","), "");
    // And a trailing empty element is dropped, so this field has two likelihoods and not three.
    assert_eq!(gl("-1,-2,"), "0,10");

    // And `Math.round` is not `floor(x + 0.5)`: the double just below a half rounds to zero.
    let just_below_a_half = 0.499_999_999_999_999_94_f64;
    assert_eq!(java_round(just_below_a_half), 0);
    assert_eq!(
        (just_below_a_half + 0.5).floor() as i64,
        1,
        "the javadoc's formula"
    );
    // The half-up rule itself survives, which is what makes the two look alike on ordinary input.
    assert_eq!(java_round(-1.5), -1);
    assert_eq!(java_round(1.5), 2);
    assert_eq!(java_round(2.5), 3);
}
