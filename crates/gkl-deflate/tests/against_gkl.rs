//! The GKL flavour, checked against GKL itself.
//!
//! `tools/gkl-probe/emulated.txt` was produced by running Intel's own `IntelDeflater` in the
//! pinned oracle container. This test recomputes the same fixtures, compresses them with this
//! crate, and compares sha256 against that file. So the assertion is against the real library's
//! bytes, not against a reading of its source, and it fails if the port drifts *or* if the
//! recorded column is regenerated with a different GKL.
//!
//! **The fixtures are rebuilt rather than shipped.** `java.util.Random` is specified exactly, so
//! it can be reimplemented, and the file carries a sha256 per fixture: the test asserts those
//! first. If the reimplementation were wrong, the fixture assertion fails and says so, instead of
//! a deflate comparison failing for a reason that has nothing to do with deflate.
//!
//! Levels 3 to 9 only. Levels 1 and 2 are igzip inside GKL and are not implemented; the test
//! asserts that they are refused rather than silently answered.

use sha2::{Digest, Sha256};
use std::collections::HashMap;

/// `java.util.Random`, which the `java.util.Random` Javadoc specifies down to the constants. This
/// is a reimplementation of a *specification*, not a transcription of the JDK's source.
struct JavaRandom(u64);

impl JavaRandom {
    fn new(seed: u64) -> Self {
        JavaRandom((seed ^ 0x5DEECE66D) & ((1 << 48) - 1))
    }

    fn next(&mut self, bits: u32) -> i32 {
        self.0 = self.0.wrapping_mul(0x5DEECE66D).wrapping_add(0xB) & ((1 << 48) - 1);
        (self.0 >> (48 - bits)) as i32
    }

    /// The power-of-two branch, which is the only one the fixtures reach (`nextInt(4)`).
    fn next_int_pow2(&mut self, bound: i32) -> i32 {
        (((bound as i64) * (self.next(31) as i64)) >> 31) as i32
    }

    fn next_bytes(&mut self, out: &mut [u8]) {
        let mut i = 0;
        while i < out.len() {
            let mut rnd = self.next(32);
            let mut n = (out.len() - i).min(4);
            while n > 0 {
                out[i] = rnd as u8;
                i += 1;
                rnd >>= 8;
                n -= 1;
            }
        }
    }
}

fn bases(len: usize, seed: u64) -> Vec<u8> {
    let mut random = JavaRandom::new(seed);
    (0..len)
        .map(|_| b"ACGT"[random.next_int_pow2(4) as usize])
        .collect()
}

fn fixtures() -> Vec<(&'static str, Vec<u8>)> {
    let mut noise = vec![0u8; 60_000];
    JavaRandom::new(11).next_bytes(&mut noise);
    vec![
        ("acgt", bases(60_000, 7)),
        (
            "runs",
            (0..60_000usize).map(|i| b"ACGT"[(i / 300) % 4]).collect(),
        ),
        ("random", noise),
        ("acgt-2blocks", bases(200_000, 13)),
    ]
}

fn sha256(data: &[u8]) -> String {
    let mut hasher = Sha256::new();
    hasher.update(data);
    hasher
        .finalize()
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}

/// The oracle column, keyed the way the file writes it.
fn recorded() -> (HashMap<String, String>, HashMap<(String, usize), String>) {
    let text = include_str!("../../../tools/gkl-probe/emulated.txt");
    let mut inputs = HashMap::new();
    let mut outputs = HashMap::new();
    for line in text.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        match parts.as_slice() {
            ["fixture", name, _, hash] => {
                inputs.insert(name.to_string(), hash.to_string());
            }
            ["deflate", name, level, "gkl", _, hash] => {
                outputs.insert((name.to_string(), level.parse().unwrap()), hash.to_string());
            }
            _ => {}
        }
    }
    (inputs, outputs)
}

#[test]
fn the_fixtures_are_the_ones_the_oracle_compressed() {
    let (inputs, _) = recorded();
    assert!(!inputs.is_empty(), "no fixture rows in emulated.txt");
    for (name, data) in fixtures() {
        assert_eq!(
            sha256(&data),
            inputs[name],
            "fixture {name} was rebuilt wrongly, so any deflate comparison using it is meaningless"
        );
    }
}

#[test]
fn byte_identical_to_gkl_at_levels_three_to_nine() {
    let (_, outputs) = recorded();
    let mut compared = 0usize;
    let mut failures = Vec::new();
    for (name, data) in fixtures() {
        for level in 3..=9usize {
            let ours = sha256(&gkl_deflate::deflate_gkl(&data, level));
            let theirs = &outputs[&(name.to_string(), level)];
            compared += 1;
            if &ours != theirs {
                failures.push(format!("{name} level {level}: ours {ours}, GKL {theirs}"));
            }
        }
    }
    println!("gkl-deflate: {compared} (fixture, level) pairs compared against GKL's own output");
    assert!(
        failures.is_empty(),
        "{} of {compared} differ:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

/// Levels 1 and 2, which go through ISA-L rather than through a port of it.
///
/// This is the test that makes linking safe rather than hopeful. ISA-L falls back to its readable
/// C when it is built without an assembler or run on a CPU without SSE4.2, and that C emits
/// **different bytes** (decision 0034: 19749 where the assembly gives 19044). A round-trip test
/// would pass in that state and a length check would too. Comparing sha256 against the column the
/// real library produced in the pinned container is what fails instead.
///
/// Both Java levels are checked, because GKL producing identical bytes for them is a claim about
/// the level not being passed through, not an assumption to inherit.
#[cfg(feature = "isal")]
#[test]
fn levels_one_and_two_go_through_isal_and_match_gkl() {
    if !gkl_deflate::igzip_available() {
        // A skip on a host that *should* be able to prove it is a failure, not a skip.
        //
        // Without this, the two outcomes are indistinguishable from outside: a CI run that
        // compared eight fixtures and a CI run that quietly compared none both report a green
        // test. Since the comparison exists precisely because ISA-L can be silently wrong, an
        // untested green here would be the same defect one level up.
        //
        // x86-64 with SSE4.2 is the condition under which decision 0033 says the kernels run, so
        // it is the condition under which the canary must pass. If it does not, the build lost
        // its assembler.
        #[cfg(target_arch = "x86_64")]
        assert!(
            !std::arch::is_x86_feature_detected!("sse4.2"),
            "this host is x86-64 with SSE4.2, so ISA-L should reproduce GKL and does not. The \
             likeliest cause is a build without nasm, which leaves ISA-L on its readable C."
        );
        // Anywhere else the refusal is the correct behaviour and is asserted by the next test.
        println!("gkl-deflate: igzip unavailable on this host; the refusal is tested instead");
        return;
    }
    let (_, outputs) = recorded();
    let mut compared = 0usize;
    let mut failures = Vec::new();
    for (name, data) in fixtures() {
        for level in 1..=2usize {
            let ours = sha256(&gkl_deflate::deflate_gkl(&data, level));
            let theirs = &outputs[&(name.to_string(), level)];
            compared += 1;
            if &ours != theirs {
                failures.push(format!("{name} level {level}: ours {ours}, GKL {theirs}"));
            }
        }
    }
    // Asserted rather than printed, because `cargo test` captures output and a count nobody reads
    // cannot distinguish "compared eight" from "compared none".
    assert_eq!(
        compared, 8,
        "the igzip comparison ran on fewer pairs than it claims"
    );
    println!("gkl-deflate: {compared} (fixture, level) pairs compared against GKL's igzip");
    assert!(
        failures.is_empty(),
        "{} of {compared} differ. A build without nasm, or a CPU without SSE4.2, produces exactly \
         this failure:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

/// A build that cannot reproduce GKL refuses those levels rather than approximating them.
///
/// Three states reach this: the crate built without the `isal` feature, ISA-L built without an
/// assembler, and any host whose CPU has no kernels to dispatch to. All three still produce valid
/// deflate, which is what makes the refusal worth testing: the failure mode being guarded against
/// is a right-looking wrong answer, not a crash.
#[test]
fn a_build_that_cannot_reproduce_gkl_refuses_rather_than_approximates() {
    if gkl_deflate::igzip_available() {
        println!("gkl-deflate: igzip available; the byte comparison covers it");
        return;
    }
    let refused = std::panic::catch_unwind(|| gkl_deflate::deflate_gkl(b"anything", 1)).is_err();
    assert!(
        refused,
        "level 1 answered on a host that cannot reproduce GKL's igzip"
    );
    println!("gkl-deflate: igzip unavailable, and levels 1 and 2 refuse");
}

/// The other branch of GKL's own CPU check.
///
/// `Flavour::Gkl { sse42: false }` is not a hypothetical: it is what GKL does on a host without
/// SSE4.2, and it emits different bytes at the levels whose hash chains are short. The reference
/// is `tools/gkl-probe/no-sse42.txt`, which that file's header is careful to say is built from
/// Intel's source rather than measured from the library, because no host here lacks SSE4.2.
#[test]
fn the_no_sse42_branch_matches_intels_source() {
    let text = include_str!("../../../tools/gkl-probe/no-sse42.txt");
    let mut expected = HashMap::new();
    for line in text.lines().filter(|l| !l.starts_with('#')) {
        let parts: Vec<&str> = line.split('\t').collect();
        if let ["deflate", name, level, "gkl-no-sse42", _, hash] = parts.as_slice() {
            expected.insert(
                (name.to_string(), level.parse::<usize>().unwrap()),
                hash.to_string(),
            );
        }
    }
    assert_eq!(expected.len(), 28, "no-sse42.txt lost rows");

    let mut differ_from_sse42 = 0;
    for (name, data) in fixtures() {
        for level in 3..=9usize {
            let out = gkl_deflate::deflate_flavour(
                &data,
                level,
                gkl_deflate::Flavour::Gkl { sse42: false },
            );
            assert_eq!(
                sha256(&out),
                expected[&(name.to_string(), level)],
                "{name} level {level} on the no-SSE4.2 path"
            );
            if sha256(&gkl_deflate::deflate_gkl(&data, level)) != sha256(&out) {
                differ_from_sse42 += 1;
            }
        }
    }
    // The point of the branch existing at all. If this ever reaches zero, the two hashes have
    // stopped mattering and `sse42` can go.
    assert!(
        differ_from_sse42 > 0,
        "the two CPU branches produced identical bytes everywhere, which contradicts the measurement"
    );
    println!("gkl-deflate: {differ_from_sse42} of 28 rows change when SSE4.2 is absent");
}
