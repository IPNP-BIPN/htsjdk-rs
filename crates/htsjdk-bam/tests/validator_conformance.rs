//! `SamFileValidator` against the reference, over the corpus both sides read from disk.
//!
//! Each case is one SAM file broken in one way, under `tools/validator-conformance/cases`. The
//! harness prints one `error` line per `SAMValidationError.toString`, in the order the validator
//! found them, or a `clean` line for a file with nothing wrong.
//!
//! The comparison is per case and ordered: the order errors come out in is part of the validator's
//! behaviour (the header's first, then each record's in htsjdk's own sequence, then the pairs that
//! never matched), and a port that produced the same set in another order would be a different
//! program.
//!
//! The golden is committed and re-derived by the `validator` suite on every run; the dump can
//! still be overridden with an environment variable while a harness change is being checked.

use std::collections::BTreeMap;
use std::io::Read;
use std::path::{Path, PathBuf};

use htsjdk_bam::sam_file::read_sam_with;
use htsjdk_bam::text_parse::ValidationStringency;
use htsjdk_bam::validation::validate;

fn cases_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tools/validator-conformance/cases")
        .canonicalize()
        .expect("the corpus directory")
}

#[test]
fn every_case_reports_what_the_reference_reports() {
    // The golden was produced by the pinned container on real x86-64 and is re-derived on every
    // run; `VALIDATOR_DUMP` still overrides it, which is how a harness change is checked before CI
    // sees it.
    let dump = match std::env::var("VALIDATOR_DUMP") {
        Ok(path) => std::fs::read_to_string(path).expect("the dump named by VALIDATOR_DUMP"),
        Err(_) => {
            let golden = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/validator.txt.gz");
            let file = std::fs::File::open(&golden).expect("the committed golden");
            let mut text = String::new();
            flate2::read::GzDecoder::new(file)
                .read_to_string(&mut text)
                .expect("the golden decompresses");
            text
        }
    };

    let mut expected: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for line in dump.lines() {
        // Not trimmed: two of the reference's messages end with a trailing space, which is where
        // htsjdk's own string concatenation stops, and trimming would hide the difference.
        let fields: Vec<&str> = line.splitn(3, '\t').collect();
        match fields.as_slice() {
            ["clean", name] => {
                expected.entry(name.to_string()).or_default();
            }
            ["error", name, text] => {
                expected
                    .entry(name.to_string())
                    .or_default()
                    .push(text.to_string());
            }
            _ => panic!("unrecognized dump line: {line}"),
        }
    }
    assert!(!expected.is_empty(), "the dump carried no case");

    for (name, theirs) in &expected {
        let path = cases_dir().join(format!("{name}.sam"));
        let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{name}: {e}"));
        let (header, records) =
            read_sam_with(&text, ValidationStringency::Silent).expect("the corpus parses");
        let ours: Vec<String> = validate(&header, &records, None)
            .expect("no reference, so nothing to resolve")
            .iter()
            .map(ToString::to_string)
            .collect();
        assert_eq!(&ours, theirs, "errors for {name}");
    }
    println!("cases={}", expected.len());
}
