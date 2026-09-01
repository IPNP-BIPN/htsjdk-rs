//! The constant-memory downsampler against the reference's own decisions and counts.
//!
//! The corpus is forty records named `read0`..`read39`, rebuilt here rather than read from a
//! fixture, because the decision is a function of the name and of nothing else.
//!
//! While the suite is `golden-pending` the dump is named by `DOWNSAMPLING_DUMP` (decision 0008).

use std::path::Path;

use htsjdk_bam::downsampling::ConstantMemoryDownsampler;

#[test]
fn every_decision_and_count_matches_the_reference() {
    let golden = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/downsampling.txt.gz");
    let dump = match std::env::var("DOWNSAMPLING_DUMP") {
        Ok(path) => std::fs::read_to_string(path).expect("the dump named by DOWNSAMPLING_DUMP"),
        Err(_) if golden.exists() => {
            panic!("the golden landed: read it here instead of skipping, and drop this branch")
        }
        Err(_) => {
            println!(
                "skipped: the downsampling golden is still pending. Run the suite and point \
                 DOWNSAMPLING_DUMP at tools/conformance/pending/downsampling.DownsamplingDump.txt"
            );
            return;
        }
    };

    let (mut decisions, mut stats) = (0, 0);
    for line in dump.lines() {
        let fields: Vec<&str> = line.trim().split('\t').collect();
        match fields.as_slice() {
            ["kept", proportion, seed, index, expected] => {
                let sampler = ConstantMemoryDownsampler::new(
                    proportion.parse().expect("a proportion"),
                    seed.parse().expect("a seed"),
                );
                let ours = sampler.keep(&format!("read{index}"));
                assert_eq!(
                    ours.to_string(),
                    *expected,
                    "kept {proportion} seed={seed} read{index}"
                );
                decisions += 1;
            }
            ["stats", proportion, seed, counts] => {
                let mut sampler = ConstantMemoryDownsampler::new(
                    proportion.parse().expect("a proportion"),
                    seed.parse().expect("a seed"),
                );
                for i in 0..40 {
                    sampler.accept(&format!("read{i}"));
                }
                let ours = format!(
                    "seen={} accepted={} discarded={}",
                    sampler.seen_count(),
                    sampler.accepted_count(),
                    sampler.discarded_count()
                );
                assert_eq!(ours, *counts, "stats {proportion} seed={seed}");
                stats += 1;
            }
            _ => panic!("unrecognized dump line: {line}"),
        }
    }
    assert!(decisions > 0 && stats > 0, "both families ran");
    println!("kept={decisions} stats={stats}");
}
