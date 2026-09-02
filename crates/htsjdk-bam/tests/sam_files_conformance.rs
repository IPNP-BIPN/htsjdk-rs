//! `SamFiles.findIndex`'s answer against the reference's, over layouts the dump carries with it.
//!
//! Each row is a directory laid out around one data file, the path the search was given, and the
//! index the reference picked. The layout travels with the answer so this test rebuilds the same
//! directory rather than keeping a second copy of the case table: a case added to the harness is
//! measured here without a line changing on this side.
//!
//! The golden is committed and re-derived by the `samfiles` suite on every run; the dump can still
//! be overridden with an environment variable while a harness change is being checked.

use std::io::Read;
use std::path::Path;

use htsjdk_bam::sam_files::find_index;

/// One row's layout, rebuilt under `root`.
fn lay_out(root: &Path, entries: &str) {
    for entry in entries.split(';') {
        let (kind, body) = entry.split_at(2);
        let body = body.to_string();
        match kind {
            "f:" => {
                let file = root.join(&body);
                std::fs::create_dir_all(file.parent().expect("a parent")).expect("the directory");
                std::fs::write(&file, [0u8]).expect("the file");
            }
            "d:" => std::fs::create_dir_all(root.join(&body)).expect("the directory"),
            "l:" => {
                let (name, target) = body.split_once("->").expect("a link target");
                let link = root.join(name);
                std::fs::create_dir_all(link.parent().expect("a parent")).expect("the directory");
                #[cfg(unix)]
                std::os::unix::fs::symlink(root.join(target), &link).expect("the symlink");
            }
            _ => panic!("unrecognized layout entry: {entry}"),
        }
    }
}

/// The answer as the dump renders it: a name under the root, both sides canonicalized.
fn rel(root: &Path, answer: Option<&Path>) -> String {
    let Some(answer) = answer else {
        return "null".to_string();
    };
    let root = root.canonicalize().expect("the root resolves");
    let answer = answer.canonicalize().expect("the answer resolves");
    answer
        .strip_prefix(&root)
        .unwrap_or_else(|_| panic!("{answer:?} is not under {root:?}"))
        .to_string_lossy()
        .into_owned()
}

#[test]
fn every_layout_finds_what_the_reference_finds() {
    // The golden was produced by the pinned container on real x86-64 and is re-derived on every
    // run; `SAMFILES_DUMP` still overrides it, which is how a harness change is checked before CI
    // sees it.
    let dump = match std::env::var("SAMFILES_DUMP") {
        Ok(path) => std::fs::read_to_string(path).expect("the dump named by SAMFILES_DUMP"),
        Err(_) => {
            let golden = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/sam_files.txt.gz");
            let file = std::fs::File::open(&golden).expect("the committed golden");
            let mut text = String::new();
            flate2::read::GzDecoder::new(file)
                .read_to_string(&mut text)
                .expect("the golden decompresses");
            text
        }
    };

    let scratch =
        std::env::temp_dir().join(format!("sam-files-conformance-{}", std::process::id()));
    let mut rows = 0;
    for line in dump.lines() {
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let mut fields = line.split('\t');
        let (Some("case"), Some(id), Some(query), Some(entries), Some(expected)) = (
            fields.next(),
            fields.next(),
            fields.next(),
            fields.next(),
            fields.next(),
        ) else {
            panic!("unrecognized dump line: {line}");
        };

        let root = scratch.join(id);
        let _ = std::fs::remove_dir_all(&root);
        std::fs::create_dir_all(&root).expect("the case's root");
        lay_out(&root, entries);

        let answer = find_index(&root.join(query));
        assert_eq!(rel(&root, answer.as_deref()), expected, "{id}");
        rows += 1;
    }

    assert!(rows > 0, "the dump carried no case");
    let _ = std::fs::remove_dir_all(&scratch);
    println!("cases={rows}");
}
