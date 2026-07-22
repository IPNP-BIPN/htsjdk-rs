//! Conformance for `SAMUtils.calculateReadGroupRecordChecksum` against htsjdk 4.2.0.
//!
//! The corpus carries a header-only SAM with four read groups (one missing the leading PU tag, two
//! tying on PU, attributes written out of key order, a value with a space) and the 32-char MD5 hex
//! htsjdk produced. The port parses the SAM header and must reproduce the checksum exactly.

use std::io::Read;

use htsjdk_bam::read_group_checksum::calculate_read_group_record_checksum;
use htsjdk_bam::sam_file::read_sam;

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/rg_checksum.txt.gz");
    let f = std::fs::File::open(&p).expect("corpus");
    let mut s = String::new();
    flate2::read::GzDecoder::new(f)
        .read_to_string(&mut s)
        .expect("corpus is gzip");
    s
}

fn unescape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut chars = s.chars();
    while let Some(c) = chars.next() {
        if c != '\\' {
            out.push(c);
            continue;
        }
        match chars.next() {
            Some('n') => out.push('\n'),
            Some('t') => out.push('\t'),
            Some('\\') => out.push('\\'),
            Some(other) => {
                out.push('\\');
                out.push(other);
            }
            None => out.push('\\'),
        }
    }
    out
}

fn payload(kind: &str) -> String {
    corpus()
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .find_map(|l| {
            let mut it = l.splitn(3, '\t');
            let k = it.next()?;
            let _case = it.next()?;
            let p = it.next().unwrap_or("");
            (k == kind).then(|| unescape(p))
        })
        .unwrap_or_else(|| panic!("no {kind} row"))
}

#[test]
fn the_read_group_checksum_is_byte_identical() {
    let (header, _records) = read_sam(&payload("sam")).expect("parse sam");
    let ours = calculate_read_group_record_checksum(&header.read_groups);
    let theirs = payload("checksum");
    assert_eq!(ours, theirs, "read group checksum diverged");
}
