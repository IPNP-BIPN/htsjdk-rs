//! Conformance for the block-copy gather (`gather_bam_files`) against htsjdk's
//! `BamFileIoUtils.gatherWithBlockCopying`, byte-for-byte.
//!
//! Goldens from `tools/bam-gather-conformance/GatherDump.java` in the pinned oracle container, with
//! the JDK deflater pinned per the oracle contract. Each case carries the input BAMs htsjdk gathered
//! (in order) and the gathered BAM it produced. The port gathers the same input bytes and must
//! reproduce the output byte-for-byte, framing included.

use std::io::Read;

fn golden_text() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/bam_gather.txt.gz");
    let f = std::fs::File::open(&p).expect("golden corpus");
    let mut s = String::new();
    flate2::read::GzDecoder::new(f)
        .read_to_string(&mut s)
        .expect("golden corpus is gzip");
    s
}

fn unhex(s: &str) -> Vec<u8> {
    (0..s.len() / 2)
        .map(|i| u8::from_str_radix(&s[i * 2..i * 2 + 2], 16).unwrap())
        .collect()
}

struct Case {
    name: String,
    inputs: Vec<Vec<u8>>,
    gathered: Vec<u8>,
}

fn cases() -> Vec<Case> {
    let text = golden_text();
    let mut order: Vec<String> = Vec::new();
    let mut inputs: std::collections::HashMap<String, Vec<Vec<u8>>> =
        std::collections::HashMap::new();
    let mut gathered = std::collections::HashMap::new();

    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let mut it = line.splitn(3, '\t');
        let kind = it.next().unwrap();
        let name = it.next().unwrap().to_string();
        let payload = it.next().unwrap_or("");
        match kind {
            "input" => {
                if !order.contains(&name) {
                    order.push(name.clone());
                }
                inputs.entry(name).or_default().push(unhex(payload));
            }
            "gathered" => {
                gathered.insert(name, unhex(payload));
            }
            other => panic!("unexpected row kind {other}"),
        }
    }

    order
        .into_iter()
        .map(|name| Case {
            inputs: inputs[&name].clone(),
            gathered: gathered[&name].clone(),
            name,
        })
        .collect()
}

#[test]
fn every_gather_case_is_byte_identical() {
    let cases = cases();
    assert_eq!(cases.len(), 5, "gather case count");
    for case in &cases {
        let refs: Vec<&[u8]> = case.inputs.iter().map(|v| v.as_slice()).collect();
        let out = htsjdk_bam::gather_bam_files(&refs).expect("gather");
        assert_eq!(
            out.len(),
            case.gathered.len(),
            "{}: length {} vs {}",
            case.name,
            out.len(),
            case.gathered.len()
        );
        let at = out.iter().zip(&case.gathered).position(|(a, b)| a != b);
        assert!(
            at.is_none(),
            "{}: first byte differs at {}",
            case.name,
            at.unwrap()
        );
    }
}
