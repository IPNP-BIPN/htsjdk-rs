//! Conformance for `SamPairUtil.setMateInfo` against htsjdk 4.2.0.
//!
//! Each case is a 2-record SAM with the mate fields absent, and the same pair after
//! `setMateInfo(rec1, rec2, true)`. Every mate field lands in the SAM columns (flags, RNEXT, PNEXT,
//! TLEN) and the MC/MQ tags, so reproducing the output SAM validates the whole method. The cases
//! cover both ends mapped (FR and same-position), a cross-contig pair, one end unmapped, and both
//! unmapped.

use std::io::Read;

use htsjdk_bam::pair::set_mate_info;
use htsjdk_bam::sam_file::{read_sam_with, write_sam};
use htsjdk_bam::text_parse::ValidationStringency;

fn corpus() -> String {
    let p =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/set_mate_info.txt.gz");
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

fn rows(kind: &str) -> Vec<(String, String)> {
    corpus()
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .filter_map(|l| {
            let mut it = l.splitn(3, '\t');
            let k = it.next()?;
            let case = it.next()?.to_string();
            let payload = it.next().unwrap_or("");
            (k == kind).then(|| (case, unescape(payload)))
        })
        .collect()
}

#[test]
fn every_case_matches_htsjdk() {
    let inputs = rows("input");
    let outputs = rows("output");
    assert!(
        inputs.len() >= 5,
        "expected >=5 cases, got {}",
        inputs.len()
    );

    let mut failures = Vec::new();
    for ((case, input), (_, expected)) in inputs.iter().zip(outputs.iter()) {
        let (header, mut records) =
            read_sam_with(input, ValidationStringency::Lenient).expect("parse input");
        assert_eq!(records.len(), 2, "{case}: expected 2 records");

        let (a, b) = records.split_at_mut(1);
        set_mate_info(&mut a[0], &mut b[0], true);

        let ours = write_sam(&header, &records).expect("write");
        if &ours != expected {
            let at = ours
                .lines()
                .zip(expected.lines())
                .position(|(x, y)| x != y)
                .unwrap_or(0);
            failures.push(format!(
                "{case}: line {at}\n  htsjdk: {:?}\n  ours  : {:?}",
                expected.lines().nth(at),
                ours.lines().nth(at)
            ));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
