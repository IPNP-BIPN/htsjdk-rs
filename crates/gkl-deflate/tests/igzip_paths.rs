//! Levels 1 and 2 on inputs chosen to reach each path of the igzip port, against the bytes ISA-L
//! 2.30's x86 kernels produced for them (amd64, SSE4.2 and AVX2, the configuration `igzip.rs`
//! documents). The four fixtures of `against_gkl.rs` reach the body's main loop and a block split;
//! these reach what they do not.

use sha2::{Digest, Sha256};

/// A small LCG, so the inputs are rebuilt rather than shipped.
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u8 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 33) as u8
    }

    fn bytes(&mut self, n: usize) -> Vec<u8> {
        (0..n).map(|_| self.next()).collect()
    }

    fn text(&mut self, n: usize) -> Vec<u8> {
        let mut out = Vec::new();
        while out.len() < n {
            let row = format!(
                "chr{}\t{}\t.\t{}\t{}\t{}\tPASS\tDP={};AF=0.{}\tGT:AD\t0/1:{},{}\n",
                self.next() % 23,
                self.next() as u32 * 997,
                b"ACGT"[(self.next() & 3) as usize] as char,
                b"ACGT"[(self.next() & 3) as usize] as char,
                self.next() % 99,
                self.next(),
                self.next() % 10,
                self.next() % 50,
                self.next() % 50,
            );
            out.extend_from_slice(row.as_bytes());
        }
        out.truncate(n);
        out
    }
}

fn inputs() -> Vec<(&'static str, Vec<u8>)> {
    let mut r = Lcg(20261113);
    let noise = r.bytes(6000);
    let unit = r.bytes(300);
    let text = r.text(15_000);
    // Over 2 KB without a match, so the kernel starts skipping, then text, then the noise again.
    let mut skip = noise[..4000].to_vec();
    skip.extend_from_slice(&text[..10_000]);
    skip.extend_from_slice(&noise[..4000]);
    let mut zeros = vec![0u8; 5000];
    zeros.extend_from_slice(&text[..3000]);
    let mut ones = vec![0xffu8; 64];
    ones.extend_from_slice(&text[..200]);
    vec![
        ("empty", Vec::new()),
        ("one", b"A".to_vec()),
        ("tiny", text[..40].to_vec()),
        ("short", text[..300].to_vec()),
        ("text", text),
        ("skip", skip),
        ("long", unit.repeat(200)),
        ("zeros", vec![0; 65_498]),
        ("zeros-then-text", zeros),
        ("ones-then-text", ones),
        ("noise", r.bytes(65_498)),
    ]
}

fn sha256(data: &[u8]) -> String {
    Sha256::digest(data)
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}

/// (input, level 1 sha256, level 2 sha256), from `igzip-oracle` in the pinned amd64 setup the
/// decision for this port describes. `text` and `long` are where the two levels differ.
const EXPECTED: [(&str, &str, &str); 11] = [
    (
        "empty",
        "f067985d352d2da6dfaef4844a66d06c5371ecbd9530a4d195ac599fef8b3427",
        "f067985d352d2da6dfaef4844a66d06c5371ecbd9530a4d195ac599fef8b3427",
    ),
    (
        "one",
        "347fd581892eff18f1849acf7d6b75e3a0f0ab3915d2ff304368ccf6202d9393",
        "347fd581892eff18f1849acf7d6b75e3a0f0ab3915d2ff304368ccf6202d9393",
    ),
    (
        "tiny",
        "e39d52788d4548bf6429228806871766d07240d8b38d523650f6d4c5f7ffaf85",
        "e39d52788d4548bf6429228806871766d07240d8b38d523650f6d4c5f7ffaf85",
    ),
    (
        "short",
        "46369129af4dd04ad43abdeccd24ea21aa02e4f09f223bcf884c8a83b664a21b",
        "46369129af4dd04ad43abdeccd24ea21aa02e4f09f223bcf884c8a83b664a21b",
    ),
    (
        "text",
        "80d2f278cba96d008f976a1a819ae3dc932565dd53e9da747fca7d900a3f2f8e",
        "f9986fee00964a0e26d7a18e04bb0102c893a88429b4ddb203d11ec9820c7c13",
    ),
    (
        "skip",
        "a1107b82b87c1faa67f93d304ff5b34aea828a17ae4b5931ae271f45c189705c",
        "c39538df9ef03b06d1f3ebfc89db8af8e06f243e7da25867f77dd67d0e89d7c0",
    ),
    (
        "long",
        "d7802e2f2c2713c244fed2a44abac52841b0c306acc7d7acd8bd8fbe5113cc73",
        "f4e07ea0691a28545fdc77260301873a00cf6928f2ef42e16f511e5e1c2c9df4",
    ),
    (
        "zeros",
        "707b06c02ef501baa1558dc0559a390272956fb3c9053fb0dc940906c4da0ab9",
        "707b06c02ef501baa1558dc0559a390272956fb3c9053fb0dc940906c4da0ab9",
    ),
    (
        "zeros-then-text",
        "f9b20f659fa2f544caef9e4eb24e6c08911e1daa95372c4fc3f62128b5ab4f6c",
        "f9b20f659fa2f544caef9e4eb24e6c08911e1daa95372c4fc3f62128b5ab4f6c",
    ),
    (
        "ones-then-text",
        "88935840bfe85547f16a4b5ff322c4e6821e5be4aaf764c7c0dc26a5b1cd17b9",
        "88935840bfe85547f16a4b5ff322c4e6821e5be4aaf764c7c0dc26a5b1cd17b9",
    ),
    (
        "noise",
        "810414400201de82f8f140e491c75b42196d0e97bf18fb5003046e098bcfdaa9",
        "810414400201de82f8f140e491c75b42196d0e97bf18fb5003046e098bcfdaa9",
    ),
];

#[test]
fn every_path_matches_isal() {
    let inputs = inputs();
    assert_eq!(inputs.len(), EXPECTED.len());
    let mut failures = Vec::new();
    for ((name, data), (expected_name, one, two)) in inputs.iter().zip(EXPECTED) {
        assert_eq!(*name, expected_name);
        for (level, expected) in [(1, one), (2, two)] {
            let ours = gkl_deflate::deflate_gkl(data, level);
            if sha256(&ours) != expected {
                failures.push(format!("{name} level {level}: {} bytes", ours.len()));
            }
        }
    }
    assert!(
        failures.is_empty(),
        "differ from ISA-L:\n{}",
        failures.join("\n")
    );
}

/// Every pinned output inflates back to its input, so a wrong expectation cannot hide behind a
/// right-looking hash.
#[test]
fn every_output_inflates_to_its_input() {
    use std::io::Read;
    for (name, data) in inputs() {
        for level in 1..=2 {
            let mut back = Vec::new();
            flate2::read::DeflateDecoder::new(&gkl_deflate::deflate_gkl(&data, level)[..])
                .read_to_end(&mut back)
                .unwrap();
            assert!(back == data, "{name} level {level}");
        }
    }
}
