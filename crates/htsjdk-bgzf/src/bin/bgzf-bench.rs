//! The port's side of the I/O floor benchmark: the same payloads, the same levels, the same
//! `name=value` lines as `tools/benchmark/BgzfBench.java`.
//!
//! It prints the digest of every framed stream it writes, so the runner compares bytes before it
//! compares seconds. Issue #78's rule is that a path with no golden gets no benchmark row; the BGZF
//! bytes have one (the `bgzf` and `zlib` suites), and this program re-establishes it on the very
//! payload it is timing rather than trusting that.

use std::io::Write;

use htsjdk_bgzf::{decompress_all, BgzfWriter};
use md5::{Digest, Md5};

fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

fn md5(bytes: &[u8]) -> String {
    let mut hasher = Md5::new();
    hasher.update(bytes);
    hex(&hasher.finalize())
}

/// The same 64-bit LCG the conformance harnesses use.
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

fn text(n: usize) -> Vec<u8> {
    let pattern = b"ACGTNacgtn\tSAMrecord\tRG:Z:rg1\n";
    (0..n).map(|i| pattern[i % pattern.len()]).collect()
}

fn deflate(input: &[u8], level: u32) -> Vec<u8> {
    let mut writer = BgzfWriter::with_level(Vec::with_capacity(input.len()), level);
    writer
        .write_all(input)
        .expect("in-memory writes never fail");
    writer.into_inner().expect("in-memory writes never fail")
}

fn main() {
    let mut args = std::env::args().skip(1);
    let megabytes: usize = args
        .next()
        .map(|a| a.parse().expect("megabytes"))
        .unwrap_or(64);
    let reps: usize = args.next().map(|a| a.parse().expect("reps")).unwrap_or(3);
    let size = megabytes * 1024 * 1024;

    for (name, input) in [("text", text(size)), ("lcg", lcg(size, 12345, 58))] {
        println!("payload_{name}_md5={}", md5(&input));
        for level in [1u32, 5, 6, 9] {
            // One untimed run each way, so both sides pay their warm-up before the clock starts.
            let framed = deflate(&input, level);
            for run in 0..reps {
                let start = std::time::Instant::now();
                let framed = deflate(&input, level);
                let seconds = start.elapsed().as_secs_f64();
                println!(
                    "rust_deflate_{name}_level{level}_run{run}_mbps={:.2}",
                    megabytes as f64 / seconds
                );
                std::hint::black_box(framed);
            }
            println!(
                "rust_deflate_{name}_level{level}_bytes={} md5={}",
                framed.len(),
                md5(&framed)
            );

            let _ = decompress_all(&framed).expect("the port reads what it wrote");
            for run in 0..reps {
                let start = std::time::Instant::now();
                let back = decompress_all(&framed).expect("the port reads what it wrote");
                let seconds = start.elapsed().as_secs_f64();
                assert_eq!(back.len(), input.len(), "round trip lost bytes");
                println!(
                    "rust_inflate_{name}_level{level}_run{run}_mbps={:.2}",
                    megabytes as f64 / seconds
                );
            }
        }
    }
}
