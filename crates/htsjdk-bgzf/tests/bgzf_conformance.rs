//! Differential test of `BgzfWriter` against htsjdk's `BlockCompressedOutputStream`.
//!
//! Goldens were produced by the official htsjdk 4.2.0 jar from Maven Central
//! (sha256 `52c9eb1a568d8261767ebf888d6ebafa60911bd44e0c9242d413fbff1d1e2398`) on
//! OpenJDK 17.0.19, via `tools/bgzf-conformance/B.java`.
//!
//! These are whole-file comparisons: header, payload, footer, and terminator block.

use std::io::Write;

use htsjdk_bgzf::BgzfWriter;
use md5::{Digest, Md5};

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
        "empty" => Vec::new(),
        "tiny" => text(10),
        // Exactly one full uncompressed block (65498), and one byte past it.
        "exact1" => runs(65498),
        "over1" => runs(65499),
        "multi" => lcg(200_000, 12345, 58),
        // Incompressible: exercises the no-compression fallback path.
        "incompr" => lcg(200_000, 999, 56),
        "big" => text(500_000),
        other => panic!("unknown payload {other}"),
    }
}

fn md5_hex(bytes: &[u8]) -> String {
    let mut h = Md5::new();
    h.update(bytes);
    h.finalize().iter().map(|b| format!("{:02x}", b)).collect()
}

fn bgzf(input: &[u8], level: u32) -> Vec<u8> {
    let mut w = BgzfWriter::with_level(Vec::new(), level);
    w.write_all(input).unwrap();
    w.into_inner().unwrap()
}

/// (payload, level, total file length, md5 of the whole file) from htsjdk 4.2.0.
#[rustfmt::skip]
const GOLDEN: &[(&str, u32, usize, &str)] = &[
    ("empty", 0, 28, "709872fc2910431b1e8b7074bfe38c67"),
    ("tiny", 0, 69, "a7fb6e2b342c2f71132da9fb29a48677"),
    ("exact1", 0, 65557, "166c54e9e73767a7df1f0533e1565e06"),
    ("over1", 0, 65589, "0c8f695f147549e5be8c95710475aa66"),
    ("multi", 0, 200152, "43c86e01681a99d847a1f0b18f340351"),
    ("incompr", 0, 200152, "a8b934bfc0252ebb83f71292f113da04"),
    ("big", 0, 500276, "1e594cbc2baee10fc00459863608cac4"),
    ("empty", 1, 28, "709872fc2910431b1e8b7074bfe38c67"),
    ("tiny", 1, 66, "17586ed8b8255520c1b3bb4f1e269515"),
    ("exact1", 1, 960, "7854bd1c5845ad36aed5c4047260af7e"),
    ("over1", 1, 989, "81a31995f47ca1812af5513712c7509b"),
    ("multi", 1, 154254, "7698015d275f6601fcefec5a3e33a8cb"),
    ("incompr", 1, 200152, "a8b934bfc0252ebb83f71292f113da04"),
    ("big", 1, 3292, "7ac8c66c01b64338f5991fc492a788a7"),
    ("empty", 5, 28, "709872fc2910431b1e8b7074bfe38c67"),
    ("tiny", 5, 66, "17586ed8b8255520c1b3bb4f1e269515"),
    ("exact1", 5, 1991, "18d70d9727a7794d8dbd86bda170b4de"),
    ("over1", 5, 2020, "fe873fac7603b6b65af8aa84346f74a3"),
    ("multi", 5, 151541, "5642a912fdbda46318c224b590b97a0e"),
    ("incompr", 5, 200152, "a8b934bfc0252ebb83f71292f113da04"),
    ("big", 5, 1806, "61cb9535b8d533f648b4611f237882a7"),
    ("empty", 6, 28, "709872fc2910431b1e8b7074bfe38c67"),
    ("tiny", 6, 66, "17586ed8b8255520c1b3bb4f1e269515"),
    ("exact1", 6, 370, "04ff13371604707f173679248c05d681"),
    ("over1", 6, 399, "4714b07126a2e64097cc0484308cccc3"),
    ("multi", 6, 151541, "5642a912fdbda46318c224b590b97a0e"),
    ("incompr", 6, 200152, "a8b934bfc0252ebb83f71292f113da04"),
    ("big", 6, 1806, "61cb9535b8d533f648b4611f237882a7"),
    ("empty", 9, 28, "709872fc2910431b1e8b7074bfe38c67"),
    ("tiny", 9, 66, "17586ed8b8255520c1b3bb4f1e269515"),
    ("exact1", 9, 370, "7bc3b22eb175f4fe86fcf25f99c1cb4e"),
    ("over1", 9, 399, "6637aee52dd682ac47fdc1bf1f66f486"),
    ("multi", 9, 151541, "5642a912fdbda46318c224b590b97a0e"),
    ("incompr", 9, 200152, "a8b934bfc0252ebb83f71292f113da04"),
    ("big", 9, 1806, "61cb9535b8d533f648b4611f237882a7"),
];

#[test]
fn empty_input_is_just_the_terminator_block() {
    let out = bgzf(&[], 5);
    assert_eq!(out, htsjdk_bgzf::EMPTY_GZIP_BLOCK.to_vec());
}

#[test]
fn block_size_constants_match_htsjdk() {
    assert_eq!(htsjdk_bgzf::DEFAULT_UNCOMPRESSED_BLOCK_SIZE, 65498);
    assert_eq!(htsjdk_bgzf::COMPRESSED_BUFFER_SIZE, 65518);
    assert_eq!(htsjdk_bgzf::DEFAULT_COMPRESSION_LEVEL, 5);
}

#[test]
fn matches_htsjdk_block_compressed_output_stream() {
    let mut checked = 0usize;
    for &(name, level, expected_len, expected_md5) in GOLDEN {
        let out = bgzf(&payload(name), level);
        assert_eq!(
            (out.len(), md5_hex(&out).as_str()),
            (expected_len, expected_md5),
            "BGZF output diverged from htsjdk for payload `{name}` at level {level}"
        );
        checked += 1;
    }
    assert_eq!(checked, 35, "expected 35 golden vectors");
}
