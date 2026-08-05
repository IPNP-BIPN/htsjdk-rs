//! Conformance for the CRAM substitution matrix, against
//! `htsjdk.samtools.cram.structure.SubstitutionMatrix`.
//!
//! Goldens from `tools/cram-conformance/CramSubstitutionMatrixDump.java` in the pinned oracle.
//!
//! The rows that justify the suite:
//!
//! ```text
//! vector  all-zero                              A  0,0,0,0,0            27  C=0,G=1,T=2,N=3
//! vector  difference-is-two-to-the-32           G  0,4294967296,0,0,0   27  A=0,C=1,T=2,N=3
//! vector  difference-is-two-to-the-32-plus-one  G  0,4294967297,0,0,0   75  A=1,C=0,T=2,N=3
//! lowercase  1b1b1b1b1b  N  0,0,0,0
//! ```
//!
//! A substitution seen 4294967296 times ranks *behind* one never seen, because the comparator
//! narrows a `long` difference to an `int`; one more occurrence and it ranks first. And `N` is the
//! one base whose row is not copied to its lower case twin.

use std::io::Read;

use htsjdk_cram::substitution_matrix::{
    MatrixError, SubstitutionMatrix, BASES, BASES_SIZE, CODES_PER_BASE, SYMBOL_SPACE_SIZE,
};

fn corpus() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/cram_substitution_matrix.txt.gz");
    let file = std::fs::File::open(&path).expect("corpus");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("corpus is gzip");
    text
}

fn rows<'a>(corpus: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    let prefix = format!("{kind}\t");
    corpus
        .lines()
        .filter_map(|line| line.strip_prefix(&prefix))
        .map(|rest| rest.split('\t').collect())
        .collect()
}

fn unhex(text: &str) -> [u8; BASES_SIZE] {
    let bytes: Vec<u8> = (0..text.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&text[i..i + 2], 16).expect("hex"))
        .collect();
    bytes.try_into().expect("five bytes")
}

/// The dump writes every character outside printable ASCII as `\uXXXX`, because these messages
/// carry a NUL and a U+FFFF and a golden that held them raw would have control bytes in it.
fn escape(message: &str) -> String {
    message
        .chars()
        .map(|c| {
            if (' '..='~').contains(&c) {
                c.to_string()
            } else {
                format!("\\u{:04x}", c as u32)
            }
        })
        .collect()
}

/// The three sizes, measured rather than assumed.
#[test]
fn the_matrix_is_five_bases_four_codes_and_a_symbol_space_of_one_hundred_and_twenty_eight() {
    let corpus = corpus();
    let sizes = rows(&corpus, "sizes");
    assert_eq!(sizes.len(), 1);
    assert_eq!(sizes[0][0], BASES_SIZE.to_string());
    assert_eq!(sizes[0][1], CODES_PER_BASE.to_string());
    assert_eq!(sizes[0][2], SYMBOL_SPACE_SIZE.to_string());
}

/// The packed code vector, and the forward lookup it fills on the way past.
///
/// This is the whole of the ranking: the double sort, the ordinal tie-break, and the comparator's
/// narrowing. Every one of them is visible in the byte.
#[test]
fn every_code_vector_is_the_one_the_reference_packed() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "vector") {
        let (label, reference) = (row[0], row[1].as_bytes()[0]);
        let counts: Vec<i64> = row[2]
            .split(',')
            .map(|value| value.parse().expect("count"))
            .collect();
        let mut frequencies = [0i64; SYMBOL_SPACE_SIZE];
        for (base, count) in BASES.iter().zip(&counts) {
            frequencies[*base as usize] = *count;
        }

        let mut matrix = SubstitutionMatrix::from_encoded([0u8; BASES_SIZE]);
        let vector = matrix.substitution_code_vector(reference, &frequencies);
        assert_eq!(
            vector.to_string(),
            row[3],
            "{label}/{}: code vector",
            reference as char
        );

        let mine: Vec<String> = BASES
            .iter()
            .filter(|base| **base != reference)
            .map(|base| {
                format!(
                    "{}={}",
                    *base as char,
                    matrix
                        .code(reference, *base)
                        .expect("a code for a substitute")
                )
            })
            .collect();
        assert_eq!(
            mine.join(","),
            row[4],
            "{label}/{}: forward lookup",
            reference as char
        );
        compared += 1;
    }
    assert_eq!(compared, 35, "code vectors compared");
}

/// What a serialized matrix decodes to, and the one base whose row is not copied to lower case.
#[test]
fn n_is_the_base_whose_row_has_no_lower_case_twin() {
    let corpus = corpus();
    let mut compared = 0;
    let mut saw_empty_lower_case_n = false;

    for (kind, lower) in [("decode", false), ("lowercase", true)] {
        for row in rows(&corpus, kind) {
            let encoded = unhex(row[0]);
            let reference = row[1].as_bytes()[0];
            let matrix = SubstitutionMatrix::from_encoded(encoded);
            let indexed = if lower {
                reference.to_ascii_lowercase()
            } else {
                reference
            };
            let mine: Vec<String> = (0..CODES_PER_BASE)
                .map(|code| match matrix.base(indexed, code as u8) {
                    Ok(base) => (base as char).to_string(),
                    Err(_) => "0".to_string(),
                })
                .collect();
            assert_eq!(mine.join(","), row[2], "{}/{kind}/{}", row[0], row[1]);
            if lower && reference == b'N' {
                assert_eq!(row[2], "0,0,0,0", "n is never filled");
                saw_empty_lower_case_n = true;
            }
            compared += 1;
        }
    }
    assert_eq!(compared, 40, "decoded rows compared");
    assert!(saw_empty_lower_case_n);
}

/// `code`: what it answers and what it refuses, including the lower case reference it will not take.
#[test]
fn code_refuses_the_lower_case_reference_that_base_accepts() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "code") {
        let matrix = SubstitutionMatrix::from_encoded(unhex(row[0]));
        let reference = row[1].parse::<i32>().expect("reference") as u8;
        let read = row[2].parse::<i32>().expect("read") as u8;
        let mine = match matrix.code(reference, read) {
            Ok(code) => code.to_string(),
            Err(error) => format!(
                "ERR java.lang.IllegalArgumentException: {}",
                escape(&error.message())
            ),
        };
        assert_eq!(mine, row[3], "code({reference}, {read})");
        compared += 1;
    }
    assert_eq!(compared, 27, "code answers compared");
}

/// `base`, whose failure names the reference base rather than the code that was actually wrong.
#[test]
fn base_blames_the_reference_base_in_its_message() {
    let corpus = corpus();
    let mut compared = 0;
    let mut saw_lower_case_success = false;
    for row in rows(&corpus, "base") {
        let matrix = SubstitutionMatrix::from_encoded(unhex(row[0]));
        let reference = row[1].parse::<i32>().expect("reference") as u8;
        let code = row[2].parse::<i32>().expect("code") as u8;
        let mine = match matrix.base(reference, code) {
            Ok(base) => (base as char).to_string(),
            Err(error) => format!(
                "ERR java.lang.IllegalArgumentException: {}",
                escape(&error.message())
            ),
        };
        assert_eq!(mine, row[3], "base({reference}, {code})");
        if reference.is_ascii_lowercase() && !row[3].starts_with("ERR") {
            saw_lower_case_success = true;
        }
        compared += 1;
    }
    assert_eq!(compared, 18, "base answers compared");
    assert!(
        saw_lower_case_success,
        "a lower case reference base decodes, which is the asymmetry"
    );

    // And the message names the reference base, not the code.
    let matrix = SubstitutionMatrix::from_encoded([0x1b; BASES_SIZE]);
    let error = matrix.base(b'A', 100).expect_err("no such code");
    assert_eq!(error, MatrixError::NoSubstitutionBase(b'A'));
    assert!(error.message().contains("invalid base 'A'"));
}

/// `toString`, which is the matrix a person reads.
#[test]
fn the_display_form_is_the_reference_display_form() {
    let corpus = corpus();
    let rows = rows(&corpus, "tostring");
    assert_eq!(rows.len(), 1);
    let matrix = SubstitutionMatrix::from_encoded(unhex(rows[0][0]));
    assert_eq!(escape(&matrix.display()), rows[0][1]);
}
