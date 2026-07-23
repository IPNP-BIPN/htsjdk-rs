//! Conformance for `build_bam_index` (the read side of `BuildBamIndex`) against htsjdk's
//! `BAMIndexer.createIndex`, byte-for-byte.
//!
//! Goldens from `tools/build-index-conformance/BuildIndexDump.java` in the pinned oracle container,
//! with the JDK deflater pinned. Each case carries the BAM htsjdk wrote and the `.bai` it built by
//! **reading** that BAM (which differs from the write-side `setCreateIndex` index: a chunk ending on
//! a BGZF block boundary is recorded as `(nextBlock, 0)` by the reading `getFilePointer`). The port
//! indexes the same BAM bytes and must reproduce the `.bai` byte-for-byte.

use std::io::Read;

fn golden_text() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/build_index.txt.gz");
    let f = std::fs::File::open(&p).expect("golden corpus");
    let mut s = String::new();
    flate2::read::GzDecoder::new(f)
        .read_to_string(&mut s)
        .expect("golden corpus is gzip");
    s
}

fn unhex(s: &str) -> Vec<u8> {
    (0..s.len() / 2)
        .map(|i| u8::from_str_radix(&s[i * 2..i * 2 + 2], 16).unwrap())
        .collect()
}

/// `kind -> name -> bytes`, preserving case order.
fn rows(kind: &str) -> Vec<(String, Vec<u8>)> {
    golden_text()
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .filter_map(|l| {
            let mut it = l.splitn(3, '\t');
            let k = it.next()?;
            let name = it.next()?;
            let hex = it.next()?;
            (k == kind).then(|| (name.to_string(), unhex(hex)))
        })
        .collect()
}

#[test]
fn build_bam_index_is_byte_identical_to_htsjdks_read_side_bai() {
    let bams = rows("bam");
    let bais = rows("bai");
    assert_eq!(bams.len(), 6, "case count");
    assert_eq!(bams.len(), bais.len(), "corpus rows have drifted");

    let mut failures = Vec::new();
    for ((name, bam), (bai_name, expected)) in bams.iter().zip(&bais) {
        assert_eq!(name, bai_name, "case lists out of step");
        let built = htsjdk_bam::build_bam_index(bam).expect("build index");
        if &built != expected {
            let at = built
                .iter()
                .zip(expected)
                .position(|(a, b)| a != b)
                .unwrap_or(built.len().min(expected.len()));
            failures.push(format!(
                "{name}: {} vs {} bytes, first differs at {at}",
                built.len(),
                expected.len()
            ));
        }
    }
    assert!(
        failures.is_empty(),
        "BAI divergences:\n{}",
        failures.join("\n")
    );
}
