//! Differential comparison with first-divergence localization.
//!
//! "The files differ" is useless when the file is a 40 GB BAM. What the port needs is *where*,
//! expressed in terms of the format rather than the byte stream: which BGZF block, which
//! record inside it, which field, and only then which byte.
//!
//! Two things this deliberately does **not** do:
//!
//! - It never reports "equal" for files that merely decompress to the same content. Byte
//!   equality is the claim; payload equality is a weaker, separately named result.
//! - It never applies a canonicalization rule silently. Every comparison returns which rules
//!   fired, because canonicalization is how a bit-identity claim quietly weakens.

use std::fmt;

pub mod canon;

/// Where two byte streams first differ, described at the most specific level available.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Divergence {
    /// Byte offset into the file as stored on disk.
    pub file_offset: u64,
    /// Which BGZF block, when the file is BGZF.
    pub block_index: Option<usize>,
    /// Byte offset of that block's first byte in the file.
    pub block_offset: Option<u64>,
    /// Offset into the block's *uncompressed* payload.
    pub uncompressed_offset: Option<u64>,
    /// Line number, 1-based, when the content is text.
    pub line: Option<usize>,
    /// A named position in the format, filled in by format-aware comparators.
    pub context: Option<String>,
    pub left: u8,
    pub right: u8,
}

impl fmt::Display for Divergence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "byte {}", self.file_offset)?;
        if let (Some(b), Some(off)) = (self.block_index, self.uncompressed_offset) {
            write!(f, " (BGZF block {b}, +{off} uncompressed)")?;
        }
        if let Some(l) = self.line {
            write!(f, " (line {l})")?;
        }
        if let Some(c) = &self.context {
            write!(f, " [{c}]")?;
        }
        write!(f, ": 0x{:02x} vs 0x{:02x}", self.left, self.right)
    }
}

/// The result of comparing two artefacts.
#[derive(Debug, Clone, PartialEq)]
pub enum Comparison {
    /// Byte-for-byte equal, with no canonicalization applied. The strongest result.
    ByteIdentical,
    /// Equal only after the listed canonicalization rules were applied. Weaker, and the rules
    /// are carried so a reader can see exactly how much was excused.
    IdenticalAfterCanonicalization { rules_applied: Vec<String> },
    /// Different, with the first divergence located.
    Divergent {
        first: Divergence,
        /// Total differing bytes, when both sides are the same length.
        differing_bytes: Option<u64>,
        left_len: u64,
        right_len: u64,
    },
}

impl Comparison {
    pub fn is_byte_identical(&self) -> bool {
        matches!(self, Self::ByteIdentical)
    }

    /// True for either exactness level. Use [`Self::is_byte_identical`] when the distinction
    /// matters, which is whenever a claim is being made.
    pub fn is_equal(&self) -> bool {
        !matches!(self, Self::Divergent { .. })
    }
}

impl fmt::Display for Comparison {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ByteIdentical => write!(f, "byte-identical"),
            Self::IdenticalAfterCanonicalization { rules_applied } => write!(
                f,
                "identical after canonicalization ({})",
                rules_applied.join(", ")
            ),
            Self::Divergent {
                first,
                differing_bytes,
                left_len,
                right_len,
            } => {
                write!(f, "divergent at {first}")?;
                if left_len != right_len {
                    write!(f, "; lengths {left_len} vs {right_len}")?;
                } else if let Some(n) = differing_bytes {
                    write!(f, "; {n} of {left_len} bytes differ")?;
                }
                Ok(())
            }
        }
    }
}

/// Raw byte comparison. The base case, and what every claim ultimately reduces to.
pub fn compare_bytes(left: &[u8], right: &[u8]) -> Comparison {
    if left == right {
        return Comparison::ByteIdentical;
    }
    let common = left.len().min(right.len());
    let mut first = None;
    let mut differing = 0u64;
    for i in 0..common {
        if left[i] != right[i] {
            differing += 1;
            if first.is_none() {
                first = Some(i);
            }
        }
    }
    // Falling off the end of the shorter side is itself the first divergence.
    let idx = first.unwrap_or(common);
    let (l, r) = (
        left.get(idx).copied().unwrap_or(0),
        right.get(idx).copied().unwrap_or(0),
    );
    Comparison::Divergent {
        first: Divergence {
            file_offset: idx as u64,
            block_index: None,
            block_offset: None,
            uncompressed_offset: None,
            line: None,
            context: None,
            left: l,
            right: r,
        },
        differing_bytes: (left.len() == right.len()).then_some(differing),
        left_len: left.len() as u64,
        right_len: right.len() as u64,
    }
}

/// Byte comparison that additionally reports the 1-based line of the first divergence.
pub fn compare_text(left: &[u8], right: &[u8]) -> Comparison {
    match compare_bytes(left, right) {
        Comparison::Divergent {
            mut first,
            differing_bytes,
            left_len,
            right_len,
        } => {
            let upto = (first.file_offset as usize).min(left.len());
            first.line = Some(1 + left[..upto].iter().filter(|&&b| b == b'\n').count());
            Comparison::Divergent {
                first,
                differing_bytes,
                left_len,
                right_len,
            }
        }
        other => other,
    }
}

/// Byte comparison that locates the divergence inside the BGZF block structure.
///
/// Reports both the compressed position and the offset into the decompressed payload, because
/// the two answer different questions: the first says which block to re-examine, the second
/// says which record went wrong.
pub fn compare_bgzf(left: &[u8], right: &[u8]) -> Comparison {
    let base = compare_bytes(left, right);
    let Comparison::Divergent {
        mut first,
        differing_bytes,
        left_len,
        right_len,
    } = base
    else {
        return base;
    };

    // Walk the block index of the left side up to the divergence.
    let mut off = 0usize;
    let mut index = 0usize;
    let mut uncompressed = 0u64;
    while off + 18 <= left.len() {
        let bsize = u16::from_le_bytes([left[off + 16], left[off + 17]]) as usize + 1;
        if bsize < 18 || off + bsize > left.len() {
            break;
        }
        if (first.file_offset as usize) < off + bsize {
            first.block_index = Some(index);
            first.block_offset = Some(off as u64);
            first.uncompressed_offset = Some(uncompressed);
            break;
        }
        let isize_off = off + bsize - 4;
        let block_uncompressed = u32::from_le_bytes([
            left[isize_off],
            left[isize_off + 1],
            left[isize_off + 2],
            left[isize_off + 3],
        ]) as u64;
        uncompressed += block_uncompressed;
        off += bsize;
        index += 1;
    }

    Comparison::Divergent {
        first,
        differing_bytes,
        left_len,
        right_len,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identical_inputs_are_byte_identical() {
        assert_eq!(compare_bytes(b"abc", b"abc"), Comparison::ByteIdentical);
        assert!(compare_bytes(b"abc", b"abc").is_byte_identical());
    }

    #[test]
    fn locates_the_first_differing_byte_not_just_any() {
        let c = compare_bytes(b"aXcYe", b"abcde");
        let Comparison::Divergent {
            first,
            differing_bytes,
            ..
        } = c
        else {
            panic!("expected divergence");
        };
        assert_eq!(first.file_offset, 1, "must report the first, not the last");
        assert_eq!((first.left, first.right), (b'X', b'b'));
        assert_eq!(differing_bytes, Some(2));
    }

    #[test]
    fn truncation_diverges_at_the_truncation_point() {
        let c = compare_bytes(b"abcdef", b"abc");
        let Comparison::Divergent {
            first,
            left_len,
            right_len,
            ..
        } = c
        else {
            panic!("expected divergence");
        };
        assert_eq!(first.file_offset, 3);
        assert_eq!((left_len, right_len), (6, 3));
    }

    #[test]
    fn text_comparison_reports_the_line() {
        let c = compare_text(b"one\ntwo\nthree\n", b"one\ntwo\nTHREE\n");
        let Comparison::Divergent { first, .. } = c else {
            panic!("expected divergence");
        };
        assert_eq!(first.line, Some(3));
    }

    /// An empty-vs-empty comparison must not be reported as divergent.
    #[test]
    fn empty_inputs_are_equal() {
        assert_eq!(compare_bytes(b"", b""), Comparison::ByteIdentical);
    }
}
