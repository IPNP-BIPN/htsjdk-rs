//! `TextualBAMIndexWriter`'s text against the reference's, over the same `.bai` files the
//! `build-index` suite measures as bytes.
//!
//! The dump carries each index twice: `bai` is the bytes, which this test parses, and `text` is
//! what the reference printed for them, which this test reproduces. A blank line is part of the
//! format and travels as `<blank>`.
//!
//! The golden is committed and re-derived by the `textual-index` suite on every run; the dump can
//! still be overridden with an environment variable while a harness change is being checked.

use std::collections::BTreeMap;
use std::io::Read;
use std::path::Path;

use htsjdk_bam::index::read_bai;
use htsjdk_bam::textual_index::render;

fn unhex(text: &str) -> Vec<u8> {
    (0..text.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&text[i..i + 2], 16).expect("a hex pair"))
        .collect()
}

#[test]
fn every_index_prints_what_the_reference_prints() {
    // The golden was produced by the pinned container on real x86-64 and is re-derived on every
    // run; `TEXTUAL_INDEX_DUMP` still overrides it, which is how a harness change is checked before CI
    // sees it.
    let dump = match std::env::var("TEXTUAL_INDEX_DUMP") {
        Ok(path) => std::fs::read_to_string(path).expect("the dump named by TEXTUAL_INDEX_DUMP"),
        Err(_) => {
            let golden =
                Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/textual_index.txt.gz");
            let file = std::fs::File::open(&golden).expect("the committed golden");
            let mut text = String::new();
            flate2::read::GzDecoder::new(file)
                .read_to_string(&mut text)
                .expect("the golden decompresses");
            text
        }
    };

    let mut indexes: BTreeMap<String, Vec<u8>> = BTreeMap::new();
    let mut expected: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for line in dump.lines() {
        let mut fields = line.splitn(3, '\t');
        match (fields.next(), fields.next(), fields.next()) {
            (Some("bai"), Some(name), Some(hex)) => {
                indexes.insert(name.to_string(), unhex(hex));
            }
            (Some("text"), Some(name), Some(text)) => {
                let text = if text == "<blank>" { "" } else { text };
                expected
                    .entry(name.to_string())
                    .or_default()
                    .push(text.to_string());
            }
            _ => panic!("unrecognized dump line: {line}"),
        }
    }

    assert!(!indexes.is_empty(), "the dump carried no index");
    for (name, bytes) in &indexes {
        let index = read_bai(bytes).unwrap_or_else(|e| panic!("{name} parses: {e:?}"));
        let ours: Vec<String> = render(&index).lines().map(str::to_string).collect();
        let theirs = expected
            .get(name)
            .unwrap_or_else(|| panic!("{name} has no text in the dump"));
        assert_eq!(&ours, theirs, "text for {name}");
    }
    println!("indexes={}", indexes.len());
}
