//! Conformance for canonical Huffman, against `HuffmanCanoncialCodeGenerator`,
//! `CanonicalHuffmanIntegerEncoding` and `CanonicalHuffmanByteEncoding`.
//!
//! Goldens from `tools/cram-conformance/CramHuffmanDump.java` in the pinned oracle, which drives
//! the codec through the slice-blocks streams the reader and writer really use and takes the raw
//! core block back out.
//!
//! The rows that justify the suite:
//!
//! ```text
//! canon  1,2,3  1,1,1   1:0/1,2:1/1,3:2/1
//! canon  bytes:65,-128,127,0  2,2,2,2  -128:0/2,0:1/2,65:2/2,127:3/2
//! round  int  42  0  42,42,42  -  42,42,42
//! cross  1,2,3,4  2,2,2,2  4  1,2,3  2,2,2  ArrayIndexOutOfBoundsException: Index 3 out of bounds for length 3
//! ```
//!
//! Three symbols at one bit are accepted and the third is given code word `2`, because the check
//! counts set bits rather than width. Byte symbols sort signed, so `0x80` takes the first code
//! word. A one-symbol alphabet writes no bits at all. And a foreign stream runs off the end of a
//! table sized to the largest code word.

use std::io::Read;

use htsjdk_cram::bit_stream::{BitInputStream, BitOutputStream};
use htsjdk_cram::huffman::{
    parse_byte_params, parse_integer_params, serialize_byte_params, serialize_integer_params,
    CanonicalHuffman, HuffmanError, Symbol,
};

fn corpus() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/cram_huffman.txt.gz");
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

fn hex(bytes: &[u8]) -> String {
    if bytes.is_empty() {
        return "-".to_string();
    }
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

/// A dump column of numbers, where a single `-` is the empty list.
fn numbers<T: std::str::FromStr>(column: &str) -> Vec<T>
where
    T::Err: std::fmt::Debug,
{
    if column == "-" {
        return Vec::new();
    }
    column
        .split(',')
        .map(|part| part.parse().expect("number"))
        .collect()
}

/// The code words as the dump prints them: `symbol:code/length`, in the order reading walks.
fn code_words<S: Symbol>(huffman: &CanonicalHuffman<S>) -> String {
    if huffman.code_words().is_empty() {
        return "-".to_string();
    }
    huffman
        .code_words()
        .iter()
        .map(|code| format!("{}:{}/{}", code.symbol, code.code_word, code.bit_length))
        .collect::<Vec<_>>()
        .join(",")
}

/// The alphabet and its lengths, rebuilt from a `canon` row's first two columns.
enum Alphabet {
    Integer(Vec<i32>, Vec<i32>),
    Byte(Vec<i8>, Vec<i32>),
}

fn alphabet(symbols: &str, lengths: &str) -> Alphabet {
    let lengths = numbers::<i32>(lengths);
    match symbols.strip_prefix("bytes:") {
        Some(rest) => Alphabet::Byte(numbers::<i8>(rest), lengths),
        None => Alphabet::Integer(numbers::<i32>(symbols), lengths),
    }
}

/// The code words come from the lengths alone, and the alphabet's order in the file is not part of
/// the answer.
#[test]
fn the_code_words_are_rebuilt_the_way_the_reference_rebuilt_them() {
    let corpus = corpus();
    let mut compared = 0;
    let mut byte_alphabets = 0;
    for row in rows(&corpus, "canon") {
        let expected = row[2];
        let actual = match alphabet(row[0], row[1]) {
            Alphabet::Integer(symbols, lengths) => {
                code_words(&CanonicalHuffman::new(&symbols, &lengths).expect("built"))
            }
            Alphabet::Byte(symbols, lengths) => {
                byte_alphabets += 1;
                code_words(&CanonicalHuffman::new(&symbols, &lengths).expect("built"))
            }
        };
        assert_eq!(actual, expected, "canon {} {}", row[0], row[1]);
        compared += 1;
    }
    assert_eq!(compared, 21, "length tables compared");
    assert_eq!(byte_alphabets, 4, "of them byte alphabets");

    // Two alphabets differing only in order get the same code words, which is the whole point of
    // deriving them.
    let forwards = CanonicalHuffman::new(&[1, 2, 3, 4], &[2, 2, 2, 2]).expect("built");
    let backwards = CanonicalHuffman::new(&[4, 3, 2, 1], &[2, 2, 2, 2]).expect("built");
    assert_eq!(code_words(&forwards), code_words(&backwards));

    // And a byte alphabet sorts signed, so 0x80 takes the first code word rather than the last.
    let bytes = CanonicalHuffman::new(&[0x41i8, -128, 127, 0], &[2, 2, 2, 2]).expect("built");
    assert_eq!(bytes.code_words()[0].symbol, -128);
}

/// The check that refuses an impossible length table counts set bits, not width, so it fires later
/// than a check on the code word would.
#[test]
fn the_length_table_is_refused_where_the_reference_refused_it() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "canonerr") {
        let (class, message) = (row[2], row[3]);
        let error = match alphabet(row[0], row[1]) {
            Alphabet::Integer(symbols, lengths) => {
                CanonicalHuffman::new(&symbols, &lengths).expect_err("refused")
            }
            Alphabet::Byte(symbols, lengths) => {
                CanonicalHuffman::new(&symbols, &lengths).expect_err("refused")
            }
        };
        assert_eq!(error.java_exception(), class, "canonerr {}", row[0]);
        assert_eq!(error.message(), message, "canonerr {}", row[0]);
        compared += 1;
    }
    assert_eq!(compared, 2, "refused length tables compared");

    // Three symbols at one bit are accepted and the third given a code word that does not fit in
    // one bit; four are refused. Nothing between those two is a rounder statement of the rule.
    let three = CanonicalHuffman::new(&[1, 2, 3], &[1, 1, 1]).expect("built");
    assert_eq!(three.code_words()[2].code_word, 2);
    assert!(matches!(
        CanonicalHuffman::new(&[1, 2, 3, 4], &[1, 1, 1, 1]),
        Err(HuffmanError::BitLengthOutOfRange {
            bit_count: 2,
            symbol: 4
        })
    ));
}

/// The encoding parameters, ITF8 throughout for the integer flavour and raw bytes for the symbols
/// of the byte one.
#[test]
fn the_encoding_parameters_are_the_reference_bytes() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "ser") {
        let (flavour, expected, reparsed) = (row[0], row[3], row[4]);
        let lengths = numbers::<i32>(row[2]);
        let bytes = match flavour {
            "int" => serialize_integer_params(&numbers::<i32>(row[1]), &lengths),
            "byte" => serialize_byte_params(&numbers::<i8>(row[1]), &lengths),
            other => panic!("{other}"),
        };
        assert_eq!(hex(&bytes), expected, "ser {flavour} {}", row[1]);

        // And parsing them back gives parameters that serialize to the same bytes.
        let again = match flavour {
            "int" => {
                let params = parse_integer_params(&bytes).expect("parsed");
                serialize_integer_params(&params.symbols, &params.bit_lengths)
            }
            "byte" => {
                let params = parse_byte_params(&bytes).expect("parsed");
                serialize_byte_params(&params.symbols, &params.bit_lengths)
            }
            other => panic!("{other}"),
        };
        assert_eq!(hex(&again), reparsed, "ser {flavour} {} reparsed", row[1]);
        compared += 1;
    }
    assert_eq!(compared, 7, "parameter sets compared");
}

/// What lands in the core block, and what comes back out of it.
#[test]
fn a_round_trip_lands_the_bytes_the_reference_landed() {
    let corpus = corpus();
    let mut compared = 0;
    let mut empty_blocks = 0;
    for row in rows(&corpus, "round") {
        let (flavour, expected, back) = (row[0], row[4], row[5]);
        let lengths = numbers::<i32>(row[2]);
        let bytes = match flavour {
            "int" => {
                let huffman =
                    CanonicalHuffman::new(&numbers::<i32>(row[1]), &lengths).expect("built");
                let mut out = BitOutputStream::new();
                for value in numbers::<i32>(row[3]) {
                    huffman.write(&mut out, value).expect("written");
                }
                let bytes = out.into_bytes();
                let mut input = BitInputStream::new(&bytes);
                let read: Vec<String> = numbers::<i32>(row[3])
                    .iter()
                    .map(|_| huffman.read(&mut input).expect("read").to_string())
                    .collect();
                assert_eq!(read.join(","), back, "round int {}", row[1]);
                bytes
            }
            "byte" => {
                let huffman =
                    CanonicalHuffman::new(&numbers::<i8>(row[1]), &lengths).expect("built");
                let mut out = BitOutputStream::new();
                for value in numbers::<i8>(row[3]) {
                    huffman.write(&mut out, value).expect("written");
                }
                let bytes = out.into_bytes();
                let mut input = BitInputStream::new(&bytes);
                let read: Vec<String> = numbers::<i8>(row[3])
                    .iter()
                    .map(|_| huffman.read(&mut input).expect("read").to_string())
                    .collect();
                assert_eq!(read.join(","), back, "round byte {}", row[1]);
                bytes
            }
            other => panic!("{other}"),
        };
        assert_eq!(hex(&bytes), expected, "round {flavour} {}", row[1]);
        if expected == "-" {
            empty_blocks += 1;
        }
        compared += 1;
    }
    assert_eq!(compared, 14, "round trips compared");

    // A one-symbol alphabet has length zero, so it writes no bits at all: both flavours produced a
    // core block with nothing in it, and how many symbols were written is not recoverable from it.
    assert_eq!(empty_blocks, 2, "empty core blocks");
}

/// A stream written with one alphabet and read with another. Nothing marks the mismatch, so it
/// either returns the wrong symbol or runs off the end of the code-word table.
#[test]
fn a_foreign_stream_comes_out_where_the_reference_put_it() {
    let corpus = corpus();
    let mut compared = 0;
    let mut refused = 0;
    for row in rows(&corpus, "cross") {
        let (values, expected) = (numbers::<i32>(row[2]), row[5]);
        let writer =
            CanonicalHuffman::new(&numbers::<i32>(row[0]), &numbers::<i32>(row[1])).expect("built");
        let mut out = BitOutputStream::new();
        for value in values {
            writer.write(&mut out, value).expect("written");
        }
        let bytes = out.into_bytes();

        let reader =
            CanonicalHuffman::new(&numbers::<i32>(row[3]), &numbers::<i32>(row[4])).expect("built");
        let mut input = BitInputStream::new(&bytes);
        let outcome = match reader.read(&mut input) {
            Ok(symbol) => symbol.to_string(),
            Err(error) => {
                refused += 1;
                format!("{}: {}", error.java_exception(), error.message())
            }
        };
        assert_eq!(outcome, expected, "cross {} into {}", row[0], row[3]);
        compared += 1;
    }
    assert_eq!(compared, 7, "foreign streams compared");
    assert_eq!(refused, 3, "of them refused");

    // Only an empty alphabet reaches the codec's own message; everything else is an index past the
    // end of a table sized to the largest code word.
    let empty = CanonicalHuffman::<i32>::new(&[], &[]).expect("built");
    assert_eq!(
        empty.read(&mut BitInputStream::new(&[0x80])),
        Err(HuffmanError::UnableToMap)
    );
}

/// Everything the codec refuses, with the reference's own message.
#[test]
fn the_failures_are_the_reference_failures() {
    let corpus = corpus();
    let mut compared = 0;
    let mut skipped = 0;
    for row in rows(&corpus, "err") {
        let (what, detail, class, message) = (row[0], row[3], row[4], row[5]);
        let lengths = numbers::<i32>(row[2]);
        let error = match (what, alphabet(row[1], row[2])) {
            ("write-unknown", Alphabet::Integer(symbols, _)) => {
                let huffman = CanonicalHuffman::new(&symbols, &lengths).expect("built");
                huffman
                    .write(&mut BitOutputStream::new(), detail.parse().expect("symbol"))
                    .expect_err("refused")
            }
            ("write-unknown", Alphabet::Byte(symbols, _)) => {
                let huffman = CanonicalHuffman::new(&symbols, &lengths).expect("built");
                huffman
                    .write(
                        &mut BitOutputStream::new(),
                        detail.parse::<i8>().expect("symbol"),
                    )
                    .expect_err("refused")
            }
            ("read-empty", Alphabet::Integer(symbols, _)) => {
                let huffman = CanonicalHuffman::new(&symbols, &lengths).expect("built");
                huffman
                    .read(&mut BitInputStream::new(&[]))
                    .expect_err("refused")
            }
            // `read(length)` is the array codecs' method, and the reference's Huffman codecs
            // implement it only to throw. The port has no such method to call, so there is nothing
            // to compare beyond the row saying so.
            ("read-length", _) => {
                assert_eq!(message, "read(length) only applicable array codecs");
                skipped += 1;
                continue;
            }
            (other, _) => panic!("{other}"),
        };
        assert_eq!(error.java_exception(), class, "err {what} {}", row[1]);
        assert_eq!(error.message(), message, "err {what} {}", row[1]);
        compared += 1;
    }
    assert_eq!(compared, 4, "refusals compared");
    assert_eq!(skipped, 2, "and the two the port has no method for");
}
