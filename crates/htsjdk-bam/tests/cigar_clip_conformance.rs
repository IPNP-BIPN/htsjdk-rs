//! Conformance for `CigarUtil.softClipEndOfRead` against htsjdk 4.2.0.
//!
//! Each case is a `(clipFrom, cigar)` pair and the cigar htsjdk produced by soft-clipping from that
//! position. The port parses the input cigar, soft-clips, and must reproduce the result cigar's
//! text. The cases exercise straddles at match/insertion/deletion boundaries and cigars that already
//! carry leading or trailing soft/hard clips, where the surgery's off-by-ones live.

use std::io::Read;

use htsjdk_bam::cigar::{soft_clip_end_of_read, Cigar};
use htsjdk_bam::text_parse::parse_cigar;

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/cigar_clip.txt.gz");
    let f = std::fs::File::open(&p).expect("corpus");
    let mut s = String::new();
    flate2::read::GzDecoder::new(f)
        .read_to_string(&mut s)
        .expect("corpus is gzip");
    s
}

#[test]
fn every_soft_clip_matches_htsjdk() {
    let mut checked = 0;
    let mut failures = Vec::new();
    for line in corpus().lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let mut it = line.splitn(3, '\t');
        let _kind = it.next().unwrap();
        let key = it.next().unwrap(); // "clipFrom:cigar"
        let expected = it.next().unwrap();

        let (clip_from, cigar_text) = key.split_once(':').unwrap();
        let clip_from: i32 = clip_from.parse().unwrap();
        let input = parse_cigar(cigar_text).expect("parse input cigar");

        let result = soft_clip_end_of_read(clip_from, &input.elements);
        let ours = Cigar::new(result).to_text();

        checked += 1;
        if ours != expected {
            failures.push(format!("{key}: htsjdk {expected}, ours {ours}"));
        }
    }
    assert!(checked >= 16, "expected at least 16 cases, got {checked}");
    assert!(
        failures.is_empty(),
        "{} of {checked} cases diverged:\n{}",
        failures.len(),
        failures.join("\n")
    );
}
