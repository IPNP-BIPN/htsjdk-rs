//! The Rust side of the zlib conformance comparison.
//!
//! Mirrors `Z2.java` exactly: the same payloads, the same levels, the same printed lines. The two
//! outputs are compared by `diff` rather than by a golden, because neither side is a reference
//! here -- the question is whether the JDK's `Deflater` in `nowrap` mode and this backend emit the
//! same bytes, which is what decision 0001 rests on.
//!
//! A line must therefore be produced identically by both programs, so the format is fixed here and
//! in `Z2.java` together, and any change to one is a change to the other.

use flate2::{Compress, Compression, FlushCompress};
use md5::{Digest, Md5};

fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

/// `Z2.lcg`: the same 64-bit LCG, taking the top byte at `shift`.
fn lcg(n: usize, seed: u64, shift: u32) -> Vec<u8> {
    let mut out = vec![0u8; n];
    let mut s = seed;
    for byte in out.iter_mut() {
        s = s
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        *byte = (s >> shift) as u8;
    }
    out
}

/// `Z2.runs`: sixty-four byte runs cycling through seven values.
fn runs(n: usize) -> Vec<u8> {
    (0..n).map(|i| ((i / 64) % 7) as u8).collect()
}

/// `Z2.text`: a repeating SAM-shaped string, which compresses well and carries a tab and a newline.
fn text(n: usize) -> Vec<u8> {
    let pattern = b"ACGTNacgtn\tSAMrecord\n";
    (0..n).map(|i| pattern[i % pattern.len()]).collect()
}

fn main() {
    let payloads: [(&str, Vec<u8>); 7] = [
        ("lcg64k", lcg(65536, 12345, 58)),
        ("rand64k", lcg(65536, 999, 56)),
        ("zeros64k", vec![0u8; 65536]),
        ("runs64k", runs(65536)),
        ("text64k", text(65536)),
        ("empty", Vec::new()),
        ("single", vec![0x42]),
    ];

    for (name, input) in &payloads {
        let mut hasher = Md5::new();
        hasher.update(input);
        println!(
            "PAYLOAD {name} len={} md5={}",
            input.len(),
            hex(&hasher.finalize())
        );
        for level in 0u32..=9 {
            // `new Deflater(level, true)`: raw deflate, no zlib header, which is what BGZF writes.
            let mut compress = Compress::new(Compression::new(level), false);
            let mut out = Vec::with_capacity(std::cmp::max(input.len() * 2, 256));
            compress
                .compress_vec(input, &mut out, FlushCompress::Finish)
                .expect("deflate never fails on an in-memory buffer");
            let mut hasher = Md5::new();
            hasher.update(&out);
            println!(
                "  (\"{name}\", {level}, {}, \"{}\"),",
                out.len(),
                hex(&hasher.finalize())
            );
        }
    }
}
