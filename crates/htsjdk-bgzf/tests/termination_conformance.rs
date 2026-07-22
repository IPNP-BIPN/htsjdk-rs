//! Differential test of `check_termination` against htsjdk's
//! `BlockCompressedInputStream.checkTermination`.
//!
//! The goldens below were produced by the official htsjdk 4.2.0 jar via
//! `tools/bgzf-conformance/Term.java`: for each payload, htsjdk writes a BGZF file with
//! `BlockCompressedOutputStream` at level 5 (which `BgzfWriter` matches byte-for-byte), then reports
//! `checkTermination` for the full file, the file with its 28-byte terminator removed, and the file
//! truncated 5 bytes into its last real block. This test rebuilds the identical bytes with
//! `BgzfWriter`, applies the same truncation, and asserts the classification matches.

use std::io::Write;

use htsjdk_bgzf::{check_termination, BgzfWriter, FileTermination};

fn lcg(n: usize, seed: u64, shift: u32) -> Vec<u8> {
    let mut b = vec![0u8; n];
    let mut s = seed;
    for x in b.iter_mut() {
        s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        *x = (s >> shift) as u8;
    }
    b
}

fn runs(n: usize) -> Vec<u8> {
    (0..n).map(|i| ((i / 64) % 7) as u8).collect()
}

fn text(n: usize) -> Vec<u8> {
    let p = b"ACGTNacgtn\tSAMrecord\n";
    (0..n).map(|i| p[i % p.len()]).collect()
}

fn payload(name: &str) -> Vec<u8> {
    match name {
        "tiny" => text(10),
        "exact1" => runs(65498),
        "multi" => lcg(200_000, 12345, 58),
        "big" => text(500_000),
        other => panic!("unknown payload {other}"),
    }
}

/// The full BGZF bytes htsjdk would write for `name`, reproduced by `BgzfWriter`.
fn bgzf(name: &str) -> Vec<u8> {
    let mut w = BgzfWriter::new(Vec::new());
    w.write_all(&payload(name)).unwrap();
    w.finish().unwrap();
    w.into_inner().unwrap()
}

fn variant(full: &[u8], variant: &str) -> Vec<u8> {
    match variant {
        "full" => full.to_vec(),
        "no_terminator" => full[..full.len() - 28].to_vec(),
        "truncated" => full[..full.len() - 33].to_vec(),
        other => panic!("unknown variant {other}"),
    }
}

fn expected(name: &str) -> FileTermination {
    match name {
        "HAS_TERMINATOR_BLOCK" => FileTermination::HasTerminatorBlock,
        "HAS_HEALTHY_LAST_BLOCK" => FileTermination::HasHealthyLastBlock,
        "DEFECTIVE" => FileTermination::Defective,
        other => panic!("unknown termination {other}"),
    }
}

/// Oracle tuples from htsjdk 4.2.0 (`tools/bgzf-conformance/Term.java`): (payload, variant, result).
const GOLDENS: &[(&str, &str, &str)] = &[
    ("tiny", "full", "HAS_TERMINATOR_BLOCK"),
    ("tiny", "no_terminator", "HAS_HEALTHY_LAST_BLOCK"),
    ("tiny", "truncated", "DEFECTIVE"),
    ("exact1", "full", "HAS_TERMINATOR_BLOCK"),
    ("exact1", "no_terminator", "HAS_HEALTHY_LAST_BLOCK"),
    ("exact1", "truncated", "DEFECTIVE"),
    ("multi", "full", "HAS_TERMINATOR_BLOCK"),
    ("multi", "no_terminator", "HAS_HEALTHY_LAST_BLOCK"),
    ("multi", "truncated", "DEFECTIVE"),
    ("big", "full", "HAS_TERMINATOR_BLOCK"),
    ("big", "no_terminator", "HAS_HEALTHY_LAST_BLOCK"),
    ("big", "truncated", "DEFECTIVE"),
];

#[test]
fn check_termination_matches_htsjdk() {
    for (name, var, want) in GOLDENS {
        let full = bgzf(name);
        let bytes = variant(&full, var);
        assert_eq!(
            check_termination(&bytes),
            expected(want),
            "payload {name} variant {var}"
        );
    }
}
