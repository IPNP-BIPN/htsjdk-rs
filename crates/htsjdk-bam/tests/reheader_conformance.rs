//! Conformance for the block-copy reheader (`reheader_bam`) against htsjdk's
//! `BamFileIoUtils.reheaderBamFile`, byte-for-byte.
//!
//! Goldens from `tools/bam-reheader-conformance/ReheaderDump.java` in the pinned oracle container,
//! with the JDK deflater pinned per the oracle contract. Each case carries the input BAM htsjdk
//! built, the comments it added, and the reheadered BAM it produced. The port parses the input's
//! header, adds the same comments, reheaders the same input bytes, and must reproduce the output
//! byte-for-byte, framing included.

use std::io::Read;

use htsjdk_bam::reader::BamReader;
use htsjdk_bam::reheader::reheader_bam;

fn golden_text() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/bam_reheader.txt.gz");
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

/// Rows are `kind \t case \t payload`. `input`/`reheadered` are hex, `comment` repeats in order.
struct Case {
    name: String,
    input: Vec<u8>,
    comments: Vec<String>,
    reheadered: Vec<u8>,
}

fn cases() -> Vec<Case> {
    let text = golden_text();
    let mut order: Vec<String> = Vec::new();
    let mut inputs = std::collections::HashMap::new();
    let mut reheadered = std::collections::HashMap::new();
    let mut comments: std::collections::HashMap<String, Vec<String>> =
        std::collections::HashMap::new();

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
                inputs.insert(name, unhex(payload));
            }
            "reheadered" => {
                reheadered.insert(name, unhex(payload));
            }
            "comment" => comments.entry(name).or_default().push(payload.to_string()),
            other => panic!("unexpected row kind {other}"),
        }
    }

    order
        .into_iter()
        .map(|name| Case {
            input: inputs[&name].clone(),
            comments: comments.get(&name).cloned().unwrap_or_default(),
            reheadered: reheadered[&name].clone(),
            name,
        })
        .collect()
}

#[test]
fn every_reheader_case_is_byte_identical() {
    let cases = cases();
    assert_eq!(cases.len(), 6, "reheader case count");
    for case in &cases {
        // Parse the input header, add the same comments, reheader the same input bytes.
        let decoded = htsjdk_bgzf::decompress_all(&case.input).expect("input BGZF");
        let reader = BamReader::new(&decoded).expect("input BAM parses");
        let mut header = reader.header.text.clone();
        for c in &case.comments {
            header.add_comment(c);
        }

        let out = reheader_bam(&header, &case.input).expect("reheader");
        assert_eq!(
            out.len(),
            case.reheadered.len(),
            "{}: length {} vs {}",
            case.name,
            out.len(),
            case.reheadered.len()
        );
        let at = out
            .iter()
            .zip(&case.reheadered)
            .position(|(a, b)| a != b);
        assert!(
            at.is_none(),
            "{}: first byte differs at {}",
            case.name,
            at.unwrap()
        );
    }
}
