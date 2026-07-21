//! Conformance of the deflate backend against the JDK's `java.util.zip.Deflater`.
//!
//! htsjdk's `BlockCompressedOutputStream` compresses BGZF blocks with
//! `new Deflater(level, true)` (`nowrap = true`). For byte-identical BAM output we must emit
//! the same compressed bytes, not merely a valid DEFLATE stream that decompresses correctly.
//!
//! Goldens below were produced by OpenJDK 17.0.19 via `tools/zlib-conformance/Z2.java`.
//! See `docs/decisions/0001-deflate-backend.md`.
//!
//! If this fails, the deflate backend changed. Do NOT relax the assertion. The default
//! `flate2` backend (miniz_oxide) emits valid but different bytes, so every round-trip test
//! would still pass while every BAM byte-comparison would fail.

use flate2::{Compress, Compression, FlushCompress};
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

fn payload(name: &str) -> Vec<u8> {
    match name {
        "lcg64k" => lcg(65536, 12345, 58),
        "rand64k" => lcg(65536, 999, 56),
        "zeros64k" => vec![0u8; 65536],
        "runs64k" => (0..65536).map(|i| ((i / 64) % 7) as u8).collect(),
        "text64k" => {
            let p = b"ACGTNacgtn\tSAMrecord\n";
            (0..65536).map(|i| p[i % p.len()]).collect()
        }
        "empty" => Vec::new(),
        "single" => vec![0x42],
        other => panic!("unknown payload {other}"),
    }
}

fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{:02x}", b)).collect()
}

fn md5_hex(bytes: &[u8]) -> String {
    let mut h = Md5::new();
    h.update(bytes);
    hex(&h.finalize())
}

fn deflate_nowrap(input: &[u8], level: u32) -> Vec<u8> {
    let mut c = Compress::new(Compression::new(level), false);
    let mut out = Vec::with_capacity(input.len().max(128) * 2);
    c.compress_vec(input, &mut out, FlushCompress::Finish)
        .expect("deflate failed");
    out
}

/// (payload, md5 of the uncompressed input) as computed by the Java harness.
const INPUT_GOLDEN: &[(&str, &str)] = &[
    ("lcg64k", "b986e8a303205123db89e21b302d03b4"),
    ("rand64k", "d2af7e0216db32a994981ef36c63adf3"),
    ("zeros64k", "fcd6bcb56c1689fcef28b57c22475bad"),
    ("runs64k", "3ebf99b511b0b02e0ae29f7080890b54"),
    ("text64k", "707189e117618885f3f54b369abc147e"),
    ("empty", "d41d8cd98f00b204e9800998ecf8427e"),
    ("single", "9d5ed678fe57bcca610140957afab571"),
];

/// (payload, level, compressed length, md5 of the compressed bytes) from OpenJDK 17.0.19.
#[rustfmt::skip]
const DEFLATE_GOLDEN: &[(&str, u32, usize, &str)] = &[
    ("lcg64k", 0, 65546, "42eb1cc7f03d2225c401931c6428a293"),
    ("lcg64k", 1, 50520, "ba6902264075ef2c7b61e9f7ee0c283b"),
    ("lcg64k", 2, 50520, "3035956dc22e00f6d38bb98ba535aef6"),
    ("lcg64k", 3, 50520, "3035956dc22e00f6d38bb98ba535aef6"),
    ("lcg64k", 4, 49603, "e8158c4d339fb75758482f9f90789aba"),
    ("lcg64k", 5, 49603, "e8158c4d339fb75758482f9f90789aba"),
    ("lcg64k", 6, 49603, "e8158c4d339fb75758482f9f90789aba"),
    ("lcg64k", 7, 49603, "e8158c4d339fb75758482f9f90789aba"),
    ("lcg64k", 8, 49603, "e8158c4d339fb75758482f9f90789aba"),
    ("lcg64k", 9, 49603, "e8158c4d339fb75758482f9f90789aba"),
    ("rand64k", 0, 65546, "b227bbe52f0f1270c2d570b323040e99"),
    ("rand64k", 1, 65556, "b902b610a30ccccffec56362334038f4"),
    ("rand64k", 2, 65556, "652ce5b89662b90bd22fa9280b3d615f"),
    ("rand64k", 3, 65556, "652ce5b89662b90bd22fa9280b3d615f"),
    ("rand64k", 4, 65556, "d9485a737bc6654410e003f84a144257"),
    ("rand64k", 5, 65556, "d9485a737bc6654410e003f84a144257"),
    ("rand64k", 6, 65556, "d9485a737bc6654410e003f84a144257"),
    ("rand64k", 7, 65556, "d9485a737bc6654410e003f84a144257"),
    ("rand64k", 8, 65556, "d9485a737bc6654410e003f84a144257"),
    ("rand64k", 9, 65556, "d9485a737bc6654410e003f84a144257"),
    ("zeros64k", 0, 65546, "bd2d4d55677c096f387e159907d014f2"),
    ("zeros64k", 1, 301, "ab1e3234b07cd95709a079198a0949fb"),
    ("zeros64k", 2, 301, "ab1e3234b07cd95709a079198a0949fb"),
    ("zeros64k", 3, 301, "ab1e3234b07cd95709a079198a0949fb"),
    ("zeros64k", 4, 78, "47e027ac2ad8db2f715ed545eb1d453e"),
    ("zeros64k", 5, 78, "47e027ac2ad8db2f715ed545eb1d453e"),
    ("zeros64k", 6, 78, "47e027ac2ad8db2f715ed545eb1d453e"),
    ("zeros64k", 7, 78, "47e027ac2ad8db2f715ed545eb1d453e"),
    ("zeros64k", 8, 78, "47e027ac2ad8db2f715ed545eb1d453e"),
    ("zeros64k", 9, 78, "47e027ac2ad8db2f715ed545eb1d453e"),
    ("runs64k", 0, 65546, "34cea58b93b8f796aa44bad7e89fa03d"),
    ("runs64k", 1, 907, "d12856c5b8009cf6f9b272500d444ea5"),
    ("runs64k", 2, 871, "d6d07888901d7ca91265d707dab6ab2e"),
    ("runs64k", 3, 902, "7ad3b16c04debd54857f66fdee492f32"),
    ("runs64k", 4, 2065, "35f0241cb2c24bbcb40ed546b2dc1604"),
    ("runs64k", 5, 1935, "cca2298a2510b0d9dea036aba2826004"),
    ("runs64k", 6, 316, "60803c6d7d54cdfb62e85f384e66c143"),
    ("runs64k", 7, 316, "60803c6d7d54cdfb62e85f384e66c143"),
    ("runs64k", 8, 315, "724266bcd791069e0b946fd26ed273e7"),
    ("runs64k", 9, 315, "724266bcd791069e0b946fd26ed273e7"),
    ("text64k", 0, 65546, "a2b3bd7c3482baace9551aa37f1d31a4"),
    ("text64k", 1, 398, "2df93af0d198b2eb6b97ed8b080f4045"),
    ("text64k", 2, 398, "2df93af0d198b2eb6b97ed8b080f4045"),
    ("text64k", 3, 398, "2df93af0d198b2eb6b97ed8b080f4045"),
    ("text64k", 4, 203, "7c644d1b5192ba39817cabb62d003aeb"),
    ("text64k", 5, 203, "7c644d1b5192ba39817cabb62d003aeb"),
    ("text64k", 6, 203, "7c644d1b5192ba39817cabb62d003aeb"),
    ("text64k", 7, 203, "7c644d1b5192ba39817cabb62d003aeb"),
    ("text64k", 8, 203, "7c644d1b5192ba39817cabb62d003aeb"),
    ("text64k", 9, 203, "7c644d1b5192ba39817cabb62d003aeb"),
    ("empty", 0, 5, "c146a7a9edbe218b6ed3bcb62ec4ad24"),
    ("empty", 1, 2, "598f4fe64aefab8f00bcbea4c9239abf"),
    ("empty", 2, 2, "598f4fe64aefab8f00bcbea4c9239abf"),
    ("empty", 3, 2, "598f4fe64aefab8f00bcbea4c9239abf"),
    ("empty", 4, 2, "598f4fe64aefab8f00bcbea4c9239abf"),
    ("empty", 5, 2, "598f4fe64aefab8f00bcbea4c9239abf"),
    ("empty", 6, 2, "598f4fe64aefab8f00bcbea4c9239abf"),
    ("empty", 7, 2, "598f4fe64aefab8f00bcbea4c9239abf"),
    ("empty", 8, 2, "598f4fe64aefab8f00bcbea4c9239abf"),
    ("empty", 9, 2, "598f4fe64aefab8f00bcbea4c9239abf"),
    ("single", 0, 6, "985a1734ba7a0b413ba6a1b99bd58fad"),
    ("single", 1, 3, "31a27c4bdb2fd8321a05396fc7a64f0c"),
    ("single", 2, 3, "31a27c4bdb2fd8321a05396fc7a64f0c"),
    ("single", 3, 3, "31a27c4bdb2fd8321a05396fc7a64f0c"),
    ("single", 4, 3, "31a27c4bdb2fd8321a05396fc7a64f0c"),
    ("single", 5, 3, "31a27c4bdb2fd8321a05396fc7a64f0c"),
    ("single", 6, 3, "31a27c4bdb2fd8321a05396fc7a64f0c"),
    ("single", 7, 3, "31a27c4bdb2fd8321a05396fc7a64f0c"),
    ("single", 8, 3, "31a27c4bdb2fd8321a05396fc7a64f0c"),
    ("single", 9, 3, "31a27c4bdb2fd8321a05396fc7a64f0c"),
];

/// The payload generators must agree byte for byte with the Java harness, otherwise any
/// agreement in the deflate test below would be comparing two different things.
#[test]
fn payload_generators_match_the_java_harness() {
    for &(name, expected) in INPUT_GOLDEN {
        assert_eq!(
            md5_hex(&payload(name)),
            expected,
            "payload generator `{name}` drifted from the Java harness"
        );
    }
}

#[test]
fn deflate_matches_jdk17_java_util_zip() {
    let mut checked = 0usize;
    for &(name, level, expected_len, expected_md5) in DEFLATE_GOLDEN {
        let out = deflate_nowrap(&payload(name), level);
        assert_eq!(
            (out.len(), md5_hex(&out).as_str()),
            (expected_len, expected_md5),
            "deflate diverged from the JDK for payload `{name}` at level {level}. Check that \
             flate2 is built with default-features = false, features = [\"zlib\"], and that \
             the vendored zlib version has not moved."
        );
        checked += 1;
    }
    assert_eq!(checked, 70, "expected 70 golden vectors");
}
