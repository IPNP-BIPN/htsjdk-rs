//! The BGZF codec against a **second, independent** decompressor.
//!
//! `rapidgzip-core` is a pure-Rust, decoder-only implementation of the rapidgzip parallel-gzip
//! algorithm, dual-licensed BSD-3-Clause or MIT. It is a **dev-dependency only**: nothing it does
//! reaches a shipped code path, and this file is the whole of its use.
//!
//! # Why it is not on the read path, and cannot be
//!
//! `BgzfReader` exists to match `BlockCompressedInputStream`'s **acceptance**, not merely to
//! decompress. Its `check_crcs` defaults to `false` because `BlockGunzipper.checkCrcs` does, and
//! `read.rs` says why: "a reader that always verified CRCs would reject files htsjdk reads
//! happily." `rapidgzip-core` verifies every member's CRC32 and uncompressed size by construction
//! and errors on a mismatch. Substituting it would change which files this library accepts, which
//! is the one property the module is there to preserve.
//!
//! It also declines seeking in decoded output and random-access indexes, which is exactly what a
//! `.bai` query needs from `vfp.rs`. So the indexed path is out on a second, independent ground.
//!
//! # What it is good for here
//!
//! Agreement. The port's writer produces BGZF and the port's reader consumes it; a bug shared by
//! both would be invisible to a round-trip test. A decompressor written by other people from the
//! same specification is an oracle for the compressed bytes in the same sense the pinned container
//! is an oracle for the numbers, and it costs one dev-dependency.
//!
//! Every golden in the workspace is compressed at each level the writer offers, then decompressed
//! three ways: by the port, by `rapidgzip-core`, and back to the original bytes.

use std::io::{Read, Write};

use htsjdk_bgzf::{decompress_all, BgzfWriter};

/// Every `.gz` golden in the workspace, which is the largest corpus of real bytes to hand.
fn corpus() -> Vec<(String, Vec<u8>)> {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("crates/")
        .parent()
        .expect("the workspace root")
        .to_path_buf();
    let mut out = Vec::new();
    let mut stack = vec![root.join("crates")];
    while let Some(directory) = stack.pop() {
        let Ok(entries) = std::fs::read_dir(&directory) else {
            continue;
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                stack.push(path);
            } else if path.extension().is_some_and(|e| e == "gz") {
                let Ok(bytes) = std::fs::read(&path) else {
                    continue;
                };
                let mut plain = Vec::new();
                if flate2::read::GzDecoder::new(&bytes[..])
                    .read_to_end(&mut plain)
                    .is_ok()
                {
                    out.push((path.display().to_string(), plain));
                }
            }
        }
    }
    out.sort_by(|a, b| a.0.cmp(&b.0));
    out
}

/// `rapidgzip-core`'s answer for a whole stream.
fn independent(bytes: &[u8]) -> std::io::Result<Vec<u8>> {
    let decoder = rapidgzip_core::Decoder::builder()
        .build()
        .map_err(|e| std::io::Error::other(e.to_string()))?;
    // `ReadAt`, not `Read`: the parallel decoder needs positioned reads, and a `Vec<u8>` provides
    // them. That requirement is the same one that makes it unusable for a non-seekable input.
    let mut reader = decoder
        .reader(bytes.to_vec())
        .map_err(|e| std::io::Error::other(e.to_string()))?;
    let mut out = Vec::new();
    reader.read_to_end(&mut out)?;
    Ok(out)
}

#[test]
fn a_second_implementation_reads_what_this_one_writes() {
    let corpus = corpus();
    assert!(!corpus.is_empty(), "no gzip goldens found to cross-check");
    let mut checked = 0;
    let mut bytes = 0usize;
    for (name, plain) in &corpus {
        // Level 5 is `BlockCompressedOutputStream`'s default; 1 and 9 bracket it.
        for level in [1u32, 5, 9] {
            let mut writer = BgzfWriter::with_level(Vec::new(), level);
            writer.write_all(plain).expect("the writer accepts bytes");
            let compressed = writer.into_inner().expect("a terminated stream");

            let ours = decompress_all(&compressed).expect("the port reads its own output");
            assert_eq!(&ours, plain, "the port's round trip lost bytes on {name}");

            let theirs = independent(&compressed)
                .unwrap_or_else(|e| panic!("rapidgzip refused {name} at level {level}: {e}"));
            assert_eq!(
                &theirs, plain,
                "the two implementations disagree on {name} at level {level}"
            );
            bytes += compressed.len();
        }
        checked += 1;
    }
    println!(
        "{checked} goldens at three levels, {bytes} compressed bytes, two implementations agreeing"
    );
}

#[test]
fn the_terminator_block_alone_decodes_to_nothing() {
    // htsjdk's empty-file marker: a BGZF member carrying an empty deflate stream. A decoder that
    // treated it as a truncated member rather than as an empty one would report an error here.
    let writer = BgzfWriter::new(Vec::new());
    let compressed = writer.into_inner().expect("a terminated stream");
    assert_eq!(decompress_all(&compressed).expect("the port"), Vec::new());
    assert_eq!(independent(&compressed).expect("rapidgzip"), Vec::new());
}
