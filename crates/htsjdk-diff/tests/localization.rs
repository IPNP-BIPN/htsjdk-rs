//! Does the harness actually find the corruption?
//!
//! The plan's Phase 0 gate is exactly this: deliberately corrupt one byte of a golden and
//! confirm the harness reports the right position. A harness that says "files differ" passes
//! a naive test and is useless on a multi-gigabyte BAM.
//!
//! These tests use real BGZF produced by the ported writer, and corrupt known positions.

use std::io::Write;

use htsjdk_bgzf::BgzfWriter;
use htsjdk_diff::{compare_bgzf, compare_bytes, Comparison};

/// Several blocks of moderately compressible content.
fn sample_bgzf() -> Vec<u8> {
    let payload: Vec<u8> = (0..250_000).map(|i| ((i / 64) % 7) as u8).collect();
    let mut w = BgzfWriter::with_level(Vec::new(), 5);
    w.write_all(&payload).unwrap();
    w.into_inner().unwrap()
}

fn divergence(c: &Comparison) -> &htsjdk_diff::Divergence {
    match c {
        Comparison::Divergent { first, .. } => first,
        other => panic!("expected a divergence, got {other}"),
    }
}

#[test]
fn identical_bgzf_files_are_byte_identical() {
    let a = sample_bgzf();
    let b = a.clone();
    assert!(compare_bgzf(&a, &b).is_byte_identical());
}

#[test]
fn a_single_corrupted_byte_is_located_exactly() {
    let a = sample_bgzf();
    // Positions are relative: the payload is highly compressible, so the file is a few
    // kilobytes even though it carries 250 kB of content.
    let positions = [0usize, 17, 100, a.len() / 2, a.len() - 1];
    for &pos in &positions {
        let mut b = a.clone();
        b[pos] ^= 0xFF;
        let c = compare_bgzf(&a, &b);
        let d = divergence(&c);
        assert_eq!(
            d.file_offset, pos as u64,
            "corruption at {pos} reported at {}",
            d.file_offset
        );
        assert_eq!(d.left, a[pos]);
        assert_eq!(d.right, b[pos]);
    }
}

#[test]
fn corruption_is_attributed_to_the_right_bgzf_block() {
    let a = sample_bgzf();

    // Find real block boundaries by walking the BSIZE fields.
    let mut boundaries = Vec::new();
    let mut off = 0usize;
    while off + 18 <= a.len() {
        let bsize = u16::from_le_bytes([a[off + 16], a[off + 17]]) as usize + 1;
        if bsize < 18 || off + bsize > a.len() {
            break;
        }
        boundaries.push(off);
        off += bsize;
    }
    assert!(
        boundaries.len() >= 4,
        "test needs several blocks, found {}",
        boundaries.len()
    );

    // Corrupt a byte inside the third block, past its header.
    let target_block = 2usize;
    let pos = boundaries[target_block] + 20;
    let mut b = a.clone();
    b[pos] ^= 0xFF;

    let c = compare_bgzf(&a, &b);
    let d = divergence(&c);
    assert_eq!(d.file_offset, pos as u64);
    assert_eq!(
        d.block_index,
        Some(target_block),
        "corruption in block {target_block} attributed to block {:?}",
        d.block_index
    );
    assert_eq!(d.block_offset, Some(boundaries[target_block] as u64));
    assert!(
        d.uncompressed_offset.is_some(),
        "the uncompressed offset is what points at the offending record"
    );
}

/// The uncompressed offset must accumulate across preceding blocks, not restart each time.
/// Getting this wrong yields a plausible small number that points at the wrong record.
#[test]
fn uncompressed_offset_accumulates_across_blocks() {
    let a = sample_bgzf();
    let mut boundaries = Vec::new();
    let mut off = 0usize;
    while off + 18 <= a.len() {
        let bsize = u16::from_le_bytes([a[off + 16], a[off + 17]]) as usize + 1;
        if bsize < 18 || off + bsize > a.len() {
            break;
        }
        boundaries.push(off);
        off += bsize;
    }

    let mut seen = Vec::new();
    for block in 1..4.min(boundaries.len()) {
        let mut b = a.clone();
        b[boundaries[block] + 20] ^= 0xFF;
        let c = compare_bgzf(&a, &b);
        seen.push(divergence(&c).uncompressed_offset.unwrap());
    }
    assert!(
        seen.windows(2).all(|w| w[0] < w[1]),
        "uncompressed offsets must increase across blocks, got {seen:?}"
    );
    // The uncompressed block size is 65498, so block 1 starts there.
    assert_eq!(seen[0], 65_498);
}

#[test]
fn truncated_file_is_reported_with_both_lengths() {
    let a = sample_bgzf();
    let mut b = a.clone();
    b.truncate(a.len() - 40);
    match compare_bgzf(&a, &b) {
        Comparison::Divergent {
            left_len,
            right_len,
            ..
        } => {
            assert_eq!(left_len as usize, a.len());
            assert_eq!(right_len as usize, a.len() - 40);
        }
        other => panic!("truncation must be a divergence, got {other}"),
    }
}

/// A one-bit change must never be reported as equal, however deep in the file.
#[test]
fn no_corruption_escapes_detection() {
    let a = sample_bgzf();
    let step = a.len() / 37;
    for pos in (1..a.len()).step_by(step.max(1)) {
        let mut b = a.clone();
        b[pos] ^= 0x01;
        assert!(
            !compare_bytes(&a, &b).is_equal(),
            "a flipped bit at {pos} was reported as equal"
        );
    }
}
