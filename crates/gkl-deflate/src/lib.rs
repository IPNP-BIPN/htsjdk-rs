//! A deflate compressor written to reproduce another implementation's bytes, not to be fast.
//!
//! ## Why this exists
//!
//! GATK writes BAM files through Intel's GKL, and GKL does not produce the bytes the JDK's
//! deflater produces. htsjdk-rs decision 0029 read the level branch out of `libgkl_compression.so`
//! and found two backends behind one `Deflater`:
//!
//! | level | backend |
//! |---|---|
//! | 1, 2 | ISA-L igzip |
//! | 3 to 9 | zlib 1.2.13 carrying Intel's `deflate_medium` patch |
//!
//! htsjdk's BGZF default is **level 5**, so the second one decides the bytes of every BAM GATK
//! writes without `--use-jdk-deflater`.
//!
//! ## What is here so far
//!
//! The zlib deflate algorithm itself: the window, the hash chains, `longest_match`, the fast and
//! slow block functions, and the whole of `trees.c`. That is the foundation both remaining pieces
//! need, and it is **verifiable on its own**: at levels 4 to 9 this crate must agree byte for byte
//! with the C zlib the workspace already links, on any input. The test does exactly that, and a
//! divergence anywhere in the match finder or the tree builder shows up as a differing byte rather
//! than as a slightly worse ratio.
//!
//! What is **not** here yet, and what H.4 still needs:
//!
//! - **`deflate_medium`**, Intel's third block function, which covers levels 4 to 6 and therefore
//!   the default;
//! - **the CRC32-based hash** Intel's fork substitutes for zlib's multiplicative one when the CPU
//!   reports SSE4.2, which changes chain order and so changes which match is found at the levels
//!   whose chains are short;
//! - **igzip**, for levels 1 and 2.
//!
//! Until those land, this crate reproduces *stock zlib*, which is the JDK's deflater and not GKL's.
//! Saying which is the point: a deflate claim that does not name its implementation is not a claim.

mod deflate;
mod trees;

/// Compress `data` as a raw deflate stream, the `nowrap` mode BGZF blocks use.
///
/// `level` is 1 to 9. Level 0 (store only) is not implemented: htsjdk never asks for it, and a
/// path nothing exercises is a path nothing checks.
pub fn deflate(data: &[u8], level: usize) -> Vec<u8> {
    assert!(
        (1..=9).contains(&level),
        "level must be 1 to 9, got {level}"
    );
    deflate::Deflater::new(data, level).finish()
}
