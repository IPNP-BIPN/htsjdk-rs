//! `TextualBAMIndexWriter`'s text against the reference's, over the same `.bai` files the
//! `build-index` suite measures as bytes.
//!
//! The dump carries each index twice: `bai` is the bytes, which this test parses, and `text` is
//! what the reference printed for them, which this test reproduces. A blank line is part of the
//! format and travels as `<blank>`.
//!
//! While the suite is `golden-pending` the dump is named by `TEXTUAL_INDEX_DUMP` (decision 0008).

use std::collections::BTreeMap;
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
    let golden = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/textual_index.txt.gz");
    let dump = match std::env::var("TEXTUAL_INDEX_DUMP") {
        Ok(path) => std::fs::read_to_string(path).expect("the dump named by TEXTUAL_INDEX_DUMP"),
        Err(_) if golden.exists() => {
            panic!("the golden landed: read it here instead of skipping, and drop this branch")
        }
        Err(_) => {
            println!(
                "skipped: the textual-index golden is still pending. Run the suite and point \
                 TEXTUAL_INDEX_DUMP at tools/conformance/pending/textual-index.TextualIndexDump.txt"
            );
            return;
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
