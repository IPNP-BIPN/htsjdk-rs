//! Conformance for commons-math3's MersenneTwister against 3.5.
//!
//! Golden from `tools/jmath-conformance/MersenneTwisterDump.java`.
//!
//! # What this suite is for
//!
//!  * **the seeding, which is not the reference algorithm's**;
//!  * **each accessor consuming its own words, so the order of calls is part of the answer**;
//!  * **`nextInt(n)` taking two different paths for a power of two and for anything else**;
//!  * **and the permutation, of which the first k entries are the answer.**

use std::io::Read;

use jmath::mersenne_twister::{next_permutation, permute, MersenneTwister};

fn corpus() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/mersenne_twister.txt.gz");
    let file = std::fs::File::open(path).expect("the golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("the golden decompresses");
    text
}

fn row(text: &str, kind: &str, name: &str) -> String {
    let prefix = format!("{kind}\t{name}\t");
    text.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| line[prefix.len()..].to_string())
        .unwrap_or_else(|| panic!("{kind}/{name}"))
}

const SEEDS: [i32; 5] = [42, 0, 1, -1, 2147483647];

/// Every integer of every seed is the golden's.
#[test]
fn the_stream_is_the_references() {
    let text = corpus();
    for seed in SEEDS {
        let mut rng = MersenneTwister::new(seed);
        let written: Vec<String> = (0..12).map(|_| rng.next_int().to_string()).collect();
        assert_eq!(
            written.join(","),
            row(&text, "ints", &seed.to_string()),
            "{seed}"
        );
    }
    // Seeds 0 and -1 agree from the second draw on, which is what makes checking one seed unsafe.
    let zero_row = row(&text, "ints", "0");
    let minus_row = row(&text, "ints", "-1");
    let zero: Vec<&str> = zero_row.split(',').collect();
    let minus: Vec<&str> = minus_row.split(',').collect();
    assert_ne!(zero[0], minus[0]);
    assert_eq!(zero[1..], minus[1..]);
}

/// Each accessor consumes its own words.
#[test]
fn each_accessor_consumes_its_own_words() {
    let text = corpus();
    for seed in SEEDS {
        let name = seed.to_string();
        let mut rng = MersenneTwister::new(seed);
        let longs: Vec<String> = (0..12).map(|_| rng.next_long().to_string()).collect();
        assert_eq!(longs.join(","), row(&text, "longs", &name), "{seed}");

        let mut rng = MersenneTwister::new(seed);
        let booleans: Vec<String> = (0..12).map(|_| rng.next_boolean().to_string()).collect();
        assert_eq!(booleans.join(","), row(&text, "booleans", &name), "{seed}");
    }
    // And the order of calls is part of the answer: one of each, in the golden's own order. The
    // numbers are PARSED rather than formatted, because `Double.toString` is a renderer this
    // crate does not own and the claim here is about the values.
    for seed in SEEDS {
        let mixed_row = row(&text, "mixed", &seed.to_string());
        let written: Vec<&str> = mixed_row.split(',').collect();
        let mut rng = MersenneTwister::new(seed);
        assert_eq!(rng.next_int().to_string(), written[0], "{seed}");
        assert_eq!(rng.next_boolean().to_string(), written[1], "{seed}");
        assert_eq!(
            rng.next_double(),
            written[2].parse::<f64>().expect("a double"),
            "{seed}"
        );
        assert_eq!(rng.next_int_bounded(7).to_string(), written[3], "{seed}");
        assert_eq!(rng.next_long().to_string(), written[4], "{seed}");
        assert_eq!(
            rng.next_float(),
            written[5].parse::<f32>().expect("a float"),
            "{seed}"
        );
        assert_eq!(rng.next_int().to_string(), written[6], "{seed}");
    }
}

/// The two floating accessors, compared as values rather than as renderings.
#[test]
fn the_doubles_and_the_floats_are_the_references() {
    let text = corpus();
    for seed in SEEDS {
        let name = seed.to_string();
        let mut rng = MersenneTwister::new(seed);
        for written in row(&text, "doubles", &name).split(',') {
            assert_eq!(
                rng.next_double(),
                written.parse::<f64>().expect("a double"),
                "{seed}"
            );
        }
        let mut rng = MersenneTwister::new(seed);
        for written in row(&text, "floats", &name).split(',') {
            assert_eq!(
                rng.next_float(),
                written.parse::<f32>().expect("a float"),
                "{seed}"
            );
        }
    }
}

/// `nextInt(n)` takes two different paths.
#[test]
fn a_power_of_two_takes_the_other_path() {
    let text = corpus();
    for seed in SEEDS {
        for bound in [2, 7, 64, 1000] {
            let mut rng = MersenneTwister::new(seed);
            let drawn: Vec<String> = (0..12)
                .map(|_| rng.next_int_bounded(bound).to_string())
                .collect();
            assert_eq!(
                drawn.join(","),
                row(&text, "bounded", &format!("{seed}/{bound}")),
                "{seed}/{bound}"
            );
        }
    }
}

/// The permutation, and the direction a caller reads it in.
#[test]
fn the_permutation_is_a_partial_shuffle() {
    let text = corpus();
    for seed in SEEDS {
        for (n, k) in [(3, 3), (4, 4), (10, 10), (10, 3), (1, 1)] {
            let mut rng = MersenneTwister::new(seed);
            let permutation = next_permutation(&mut rng, n, k).expect("a permutation");
            let written: Vec<String> = permutation.iter().map(|i| i.to_string()).collect();
            assert_eq!(
                written.join(","),
                row(&text, "permutation", &format!("{seed}/{n}/{k}")),
                "{seed}/{n}/{k}"
            );
        }
        // `out[i] = in[perm[i]]`, which is the direction that decides where a value lands.
        let mut rng = MersenneTwister::new(seed);
        let permutation = next_permutation(&mut rng, 4, 4).expect("a permutation");
        let values = [10.0, 20.0, 30.0, 40.0];
        let permuted = permute(&values, &permutation);
        let written: Vec<f64> = row(&text, "permuted", &seed.to_string())
            .split(',')
            .map(|value| value.parse::<f64>().expect("a double"))
            .collect();
        assert_eq!(permuted, written, "{seed}");
    }
    // The two refusals the golden holds: k greater than n, and a size of zero.
    let mut rng = MersenneTwister::new(42);
    assert_eq!(next_permutation(&mut rng, 3, 5), None);
    assert!(row(&text, "error", "3/5").contains("NumberIsTooLargeException"));
    let mut rng = MersenneTwister::new(42);
    assert_eq!(next_permutation(&mut rng, 5, 0), None);
    assert!(row(&text, "error", "5/0").contains("NotStrictlyPositiveException"));
}
