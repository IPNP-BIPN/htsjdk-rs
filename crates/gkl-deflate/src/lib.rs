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
//! ## What is here
//!
//! Both flavours, verified two different ways.
//!
//! [`Flavour::Jdk`] is stock zlib, which is what `java.util.zip.Deflater` reaches. It is checked
//! against the C zlib the workspace already links, byte for byte, at every level.
//!
//! [`Flavour::Gkl`] is Intel's fork: `deflate_medium` at levels 4 to 6, and a CRC-32C positional
//! hash in place of zlib's multiplicative rolling one. It is checked against the hashes GKL itself
//! produced in the pinned container, so the assertion is against the real library rather than
//! against a reading of its source.
//!
//! **GKL's output depends on the CPU, and this is where that becomes visible.** The fork selects
//! its hash on `x86_cpu_has_sse42` at load time, so a host without SSE4.2 fills the chains
//! differently and emits different bytes at every level from 3 up. `Flavour::Gkl { sse42 }` makes
//! that a parameter rather than a hidden assumption; the default of `true` is the column every
//! oracle run has been measured in, because every host so far reports SSE4.2.
//!
//! ## What is not here
//!
//! **igzip, for levels 1 and 2.** Those are the only levels GKL does not route through zlib, and
//! nothing in htsjdk asks for them by default. [`deflate_gkl`] refuses them rather than quietly
//! answering with zlib's bytes, which would be a wrong answer that looks like a right one.

mod deflate;
mod trees;

pub use deflate::Flavour;

/// Compress `data` as a raw deflate stream the way `java.util.zip.Deflater` does, which is the
/// `nowrap` mode BGZF blocks use.
///
/// `level` is 1 to 9. Level 0 (store only) is not implemented: htsjdk never asks for it, and a
/// path nothing exercises is a path nothing checks.
pub fn deflate(data: &[u8], level: usize) -> Vec<u8> {
    assert!(
        (1..=9).contains(&level),
        "level must be 1 to 9, got {level}"
    );
    deflate::Deflater::new(data, level, Flavour::Jdk).finish()
}

/// Compress `data` with an explicit flavour, for callers that need the non-default CPU branch.
pub fn deflate_flavour(data: &[u8], level: usize, flavour: Flavour) -> Vec<u8> {
    match flavour {
        Flavour::Jdk => assert!(
            (1..=9).contains(&level),
            "level must be 1 to 9, got {level}"
        ),
        Flavour::Gkl { .. } => assert!(
            (3..=9).contains(&level),
            "GKL routes levels 1 and 2 through igzip, which is not implemented; got {level}"
        ),
    }
    deflate::Deflater::new(data, level, flavour).finish()
}

/// Compress `data` the way GKL's `IntelDeflater` does, on a CPU reporting SSE4.2.
///
/// `level` is **3 to 9**. Levels 1 and 2 are igzip inside GKL, not zlib, and are not implemented;
/// answering them with zlib's bytes would be a wrong answer wearing a right one's shape.
pub fn deflate_gkl(data: &[u8], level: usize) -> Vec<u8> {
    assert!(
        (3..=9).contains(&level),
        "GKL routes levels 1 and 2 through igzip, which is not implemented; got {level}"
    );
    deflate::Deflater::new(data, level, Flavour::Gkl { sse42: true }).finish()
}
