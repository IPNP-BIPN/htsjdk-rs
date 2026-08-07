//! The port's side of a differential fuzzer over the byte-level parsers.
//!
//! Reads one hex string per line on stdin, runs each through the named parser, and prints what the
//! port did with it in the same shape `tools/fuzz/FuzzDriver.java` prints. The two are diffed line
//! by line, so a divergence is a line rather than a report.
//!
//! The parsers are the ones a hostile file reaches first, and each is a pure function of its
//! bytes: no state, no environment, nothing to seed. That is what makes a divergence a bug rather
//! than a difference of setup.

use std::io::{BufRead, Write};

use htsjdk_cram::crai::CraiEntry;
use htsjdk_cram::varint::{read_unsigned_itf8, read_unsigned_ltf8};

fn main() {
    let parser = std::env::args().nth(1).unwrap_or_else(|| "itf8".to_string());
    let stdin = std::io::stdin();
    let mut out = std::io::BufWriter::new(std::io::stdout());

    for line in stdin.lock().lines() {
        let line = line.expect("a line");
        let line = line.trim().to_string();
        if line.is_empty() {
            continue;
        }
        let bytes = unhex(&line);
        writeln!(out, "{line}\t{parser}\t{}", outcome(&parser, &bytes)).expect("a line");
    }
}

fn outcome(parser: &str, bytes: &[u8]) -> String {
    match parser {
        "itf8" => match read_unsigned_itf8(bytes) {
            Ok((value, _)) => format!("ok:{value}"),
            // The reference's RuntimeEOFException, which its stream form throws only on an empty
            // input: a short one gives a silently wrong number on both sides.
            Err(_) => "err:RuntimeEOFException".to_string(),
        },
        "ltf8" => match read_unsigned_ltf8(bytes) {
            Ok((value, _)) => format!("ok:{value}"),
            Err(_) => "err:RuntimeEOFException".to_string(),
        },
        // Java decodes invalid UTF-8 into replacement characters rather than refusing, so a line
        // that is not UTF-8 reaches the parser there and has to reach it here too.
        "crai" => {
            let text = String::from_utf8_lossy(bytes);
            match CraiEntry::parse(&text) {
                Ok(entry) => format!("ok:{}", entry.serialize().replace('\t', " ")),
                Err(error) => format!("err:{}", error.java_exception()),
            }
        }
        _ => "err:UnknownParser".to_string(),
    }
}

fn unhex(text: &str) -> Vec<u8> {
    (0..text.len() / 2)
        .map(|at| u8::from_str_radix(&text[at * 2..at * 2 + 2], 16).expect("hex"))
        .collect()
}
