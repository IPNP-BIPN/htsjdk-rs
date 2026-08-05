//! Every byte this crate emits, checked against the C zlib the workspace already links.
//!
//! This is the test that makes the port meaningful before any of GKL's own behaviour is in place.
//! A deflate stream that merely decompresses correctly proves almost nothing: miniz_oxide passes
//! that bar and produces different bytes, which is the whole subject of decision 0001. So the
//! assertion here is byte equality, on inputs chosen to reach the parts of the algorithm where two
//! implementations diverge:
//!
//! - **larger than the window**, so the slide and `slide_hash` run and the hash chains have to
//!   survive it;
//! - **more than 16384 symbols**, so a block is closed mid-stream and the trees are rebuilt;
//! - **incompressible**, so the stored block wins the three-way comparison in `flush_block`;
//! - **long runs**, where many equally long matches exist and the chain order decides which one is
//!   taken.
//!
//! Levels 1 to 9 are all compared. Level 0 is not implemented and is not tested.

use std::io::Write;

fn zlib_deflate(data: &[u8], level: u32) -> Vec<u8> {
    let mut encoder =
        flate2::write::DeflateEncoder::new(Vec::new(), flate2::Compression::new(level));
    encoder.write_all(data).unwrap();
    encoder.finish().unwrap()
}

/// A linear congruential generator, so the fixtures are the same on every machine and in every
/// language. `rand` would tie this test to a crate's version.
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u8 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 33) as u8
    }
}

fn fixtures() -> Vec<(&'static str, Vec<u8>)> {
    let mut cases: Vec<(&'static str, Vec<u8>)> = Vec::new();

    cases.push(("empty", Vec::new()));
    cases.push(("one-byte", b"A".to_vec()));
    cases.push((
        "short",
        b"the quick brown fox jumps over the lazy dog".to_vec(),
    ));

    // Four symbols, the shape of a BAM's sequence column.
    let mut lcg = Lcg(7);
    cases.push((
        "acgt-60k",
        (0..60_000)
            .map(|_| b"ACGT"[(lcg.next() % 4) as usize])
            .collect(),
    ));

    // Past the 32 KiB window, so fill_window slides and slide_hash rewrites every chain.
    let mut lcg = Lcg(13);
    cases.push((
        "acgt-200k",
        (0..200_000)
            .map(|_| b"ACGT"[(lcg.next() % 4) as usize])
            .collect(),
    ));

    // Incompressible: the stored block should win, and it is emitted by a different path.
    let mut lcg = Lcg(11);
    cases.push(("noise-60k", (0..60_000).map(|_| lcg.next()).collect()));

    // Long runs: many matches of equal length, so which one the chain yields is what is compared.
    cases.push((
        "runs-60k",
        (0..60_000usize).map(|i| b"ACGT"[(i / 300) % 4]).collect(),
    ));

    // One byte repeated: every position hashes the same, so the chain is as long as it can be and
    // max_chain is what stops the search.
    cases.push(("all-same-100k", vec![b'N'; 100_000]));

    // Text, which is what a SAM header or a read name column looks like.
    let line = b"@SRR12345.678 HWI-ST1234:56:C7ABCDEF:1:1101:1234:5678 length=151\n";
    let mut text = Vec::new();
    while text.len() < 120_000 {
        text.extend_from_slice(line);
    }
    cases.push(("readnames-120k", text));

    cases
}

#[test]
fn byte_identical_to_zlib_at_every_level() {
    let mut compared = 0usize;
    let mut failures = Vec::new();
    for (name, data) in fixtures() {
        for level in 1..=9u32 {
            let ours = gkl_deflate::deflate(&data, level as usize);
            let theirs = zlib_deflate(&data, level);
            compared += 1;
            if ours != theirs {
                let at = ours
                    .iter()
                    .zip(theirs.iter())
                    .position(|(a, b)| a != b)
                    .map(|p| p.to_string())
                    .unwrap_or_else(|| "length only".to_string());
                failures.push(format!(
                    "{name} level {level}: ours {} bytes, zlib {} bytes, first difference at {at}",
                    ours.len(),
                    theirs.len()
                ));
            }
        }
    }
    println!("gkl-deflate: {compared} (fixture, level) pairs compared against C zlib");
    assert!(
        failures.is_empty(),
        "{} of {compared} differ:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

/// Byte equality already implies this, but a round trip localises a failure: if this passes and
/// the comparison above fails, the bug is in the encoding choices rather than in the match finder.
#[test]
fn round_trips() {
    for (name, data) in fixtures() {
        for level in 1..=9usize {
            let compressed = gkl_deflate::deflate(&data, level);
            let mut out = Vec::new();
            let mut decoder = flate2::write::DeflateDecoder::new(&mut out);
            decoder.write_all(&compressed).unwrap();
            decoder.finish().unwrap();
            assert_eq!(out, data, "{name} level {level} did not round trip");
        }
    }
}
