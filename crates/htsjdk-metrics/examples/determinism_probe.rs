//! Everything the port formats or writes, hashed, so a determinism gate can compare two runs.
//!
//! htsjdk-rs is a library and ships no tool, so there is nothing for a determinism gate to run.
//! This is that thing: it exercises the paths whose output a byte-identity claim is made about,
//! prints one line per path, and writes a file into `TMPDIR` so that a temporary directory the
//! process cannot write to is a failure rather than a silent skip.
//!
//! The interesting environment is the locale. htsjdk-rs decision 0011 established that metrics
//! number formatting is locale-dependent **in the reference**, which is why the oracle pins
//! `en_US`. The port must not have acquired the same dependency, so the gate runs it under
//! `fr_FR`, where a locale-sensitive formatter writes `0,333333` and this one must still write
//! `0.333333`.

use std::io::Write;

fn main() {
    let mut lines = Vec::new();

    // The formatters, which are where a locale would show if the port had one.
    lines.push(format!(
        "double\t{}",
        htsjdk_metrics::format::format_double(1.0 / 3.0)
    ));
    lines.push(format!(
        "double-big\t{}",
        htsjdk_metrics::format::format_double(1.234_567_890_123e12)
    ));
    lines.push(format!(
        "double-small\t{}",
        htsjdk_metrics::format::format_double(2.5e-7)
    ));
    lines.push(format!(
        "long\t{}",
        htsjdk_metrics::format::format_long(-1234567)
    ));
    lines.push(format!(
        "bool\t{}",
        htsjdk_metrics::format::format_bool(true)
    ));
    lines.push(format!(
        "nan\t{}",
        htsjdk_metrics::format::format_double(f64::NAN)
    ));
    lines.push(format!(
        "inf\t{}",
        htsjdk_metrics::format::format_double(f64::NEG_INFINITY)
    ));

    // A hash of the whole, so a diff of two runs is one line rather than many.
    let digest = fnv(lines.join("\n").as_bytes());
    lines.push(format!("digest\t{digest:016x}"));

    // Written into TMPDIR as well as to stdout: the gate points TMPDIR somewhere else on one of
    // its runs, and a path that cannot be written is a failure rather than a skip.
    let temporary = std::env::var("TMPDIR").unwrap_or_else(|_| "/tmp".to_string());
    let path = std::path::Path::new(&temporary).join("htsjdk-rs-determinism.txt");
    let mut file = std::fs::File::create(&path)
        .unwrap_or_else(|error| panic!("could not write {}: {error}", path.display()));
    for line in &lines {
        writeln!(file, "{line}").expect("wrote a line");
        println!("{line}");
    }
}

/// FNV-1a, so the probe needs no dependency of its own.
fn fnv(bytes: &[u8]) -> u64 {
    let mut hash = 0xcbf2_9ce4_8422_2325u64;
    for byte in bytes {
        hash ^= u64::from(*byte);
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    hash
}
