//! The BAI indexing bin.
//!
//! Ported from `htsjdk.samtools.GenomicIndexUtil.regionToBin` and
//! `htsjdk.samtools.SAMRecord.computeIndexingBin`.
//!
//! The bin is stored in the fixed part of every BAM record, so getting it wrong changes the
//! bytes of every single read while leaving a file that every tool still reads happily: the
//! bin is a search hint, and readers that scan linearly never consult it. It is exactly the
//! class of divergence that only byte comparison catches.

/// `GenomicIndexUtil.BIN_GENOMIC_SPAN`. The largest coordinate the 6-level scheme can address.
pub const BIN_GENOMIC_SPAN: i32 = 512 * 1024 * 1024;

/// `GenomicIndexUtil.LEVEL_STARTS`.
pub const LEVEL_STARTS: [i32; 6] = [0, 1, 9, 73, 585, 4681];

/// `GenomicIndexUtil.MAX_BINS`, `(8^6-1)/7+1`.
pub const MAX_BINS: i32 = 37450;

/// `SAMRecord.NO_ALIGNMENT_START`, which is `GenomicIndexUtil.UNSET_GENOMIC_LOCATION`.
pub const NO_ALIGNMENT_START: i32 = 0;

/// `SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX`.
pub const NO_ALIGNMENT_REFERENCE_INDEX: i32 = -1;

/// `GenomicIndexUtil.regionToBin(beg, end)`, for a 0-based half-open interval `[beg, end)`.
///
/// The `--end` on entry and the exact constants are reproduced rather than simplified. The
/// literals `((1<<15)-1)/7` and friends are the level starts computed inline, and they agree
/// with [`LEVEL_STARTS`]; they are written the same way here so the correspondence to the Java
/// is checkable line by line.
// `((1 << 3) - 1) / 7` is 7/7, which clippy flags as a suspicious self-division. It is kept
// because the five level starts are written as one visible pattern in the Java, and rewriting
// the degenerate one as `1` would break the correspondence that makes the transcription
// checkable by eye. The test `inline_level_starts_match_the_declared_table` pins the values.
#[allow(clippy::eq_op)]
pub fn region_to_bin(beg: i32, end: i32) -> i32 {
    let end = end - 1;
    if beg >> 14 == end >> 14 {
        return ((1 << 15) - 1) / 7 + (beg >> 14);
    }
    if beg >> 17 == end >> 17 {
        return ((1 << 12) - 1) / 7 + (beg >> 17);
    }
    if beg >> 20 == end >> 20 {
        return ((1 << 9) - 1) / 7 + (beg >> 20);
    }
    if beg >> 23 == end >> 23 {
        return ((1 << 6) - 1) / 7 + (beg >> 23);
    }
    if beg >> 26 == end >> 26 {
        return ((1 << 3) - 1) / 7 + (beg >> 26);
    }
    0
}

/// Why [`compute_indexing_bin`] declined to produce a bin.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BinError {
    /// `computeIndexingBin` throws `IllegalStateException` past [`BIN_GENOMIC_SPAN`].
    PositionTooHigh { start: i32, end: i32 },
}

/// `SAMRecord.computeIndexingBin()`.
///
/// `alignment_start` is 1-based as htsjdk holds it; `alignment_end` is the 1-based inclusive
/// end, or `<= 0` when it cannot be determined.
///
/// Two details carry all the risk:
///
/// - The end is **not** converted to 0-based. htsjdk subtracts one from the start only, then
///   passes the still-1-based end as the half-open exclusive end, which happens to be the
///   right number. Converting both, the symmetric-looking thing to do, shifts every bin.
/// - An undeterminable end becomes `alignment_start + 1` *after* the start was made 0-based,
///   so an unmapped-but-placed read is binned as a single base at its 0-based start.
pub fn compute_indexing_bin(alignment_start: i32, alignment_end: i32) -> Result<i32, BinError> {
    let alignment_start = alignment_start - 1; // BIN uses 0-based half-open
    let alignment_end = if alignment_end <= 0 {
        alignment_start + 1
    } else {
        alignment_end
    };

    if alignment_start > BIN_GENOMIC_SPAN || alignment_end > BIN_GENOMIC_SPAN {
        return Err(BinError::PositionTooHigh {
            start: alignment_start,
            end: alignment_end,
        });
    }

    // `& (int) BinaryCodec.MAX_USHORT` in the Java.
    Ok(region_to_bin(alignment_start, alignment_end) & 0xFFFF)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The inline literals in `regionToBin` must equal the declared level starts. If they ever
    /// disagree, one of the two transcriptions is wrong.
    #[test]
    fn inline_level_starts_match_the_declared_table() {
        assert_eq!(((1 << 15) - 1) / 7, LEVEL_STARTS[5]);
        assert_eq!(((1 << 12) - 1) / 7, LEVEL_STARTS[4]);
        assert_eq!(((1 << 9) - 1) / 7, LEVEL_STARTS[3]);
        assert_eq!(((1 << 6) - 1) / 7, LEVEL_STARTS[2]);
        assert_eq!(((1 << 3) - 1) / 7, LEVEL_STARTS[1]);
    }

    #[test]
    fn a_short_read_lands_on_the_deepest_level() {
        // 0-based [0, 100): entirely inside the first 16 kb bin.
        assert_eq!(region_to_bin(0, 100), 4681);
        assert_eq!(region_to_bin(16_384, 16_484), 4682);
    }

    #[test]
    fn spanning_a_boundary_promotes_to_a_coarser_level() {
        // Straddles the 16 kb boundary, so it cannot sit at level 5.
        let b = region_to_bin(16_300, 16_500);
        assert_eq!(
            b, 585,
            "a read crossing a 16 kb boundary belongs to level 4"
        );
    }

    /// A read that spans everything falls all the way to bin 0.
    #[test]
    fn spanning_the_whole_span_is_bin_zero() {
        assert_eq!(region_to_bin(0, BIN_GENOMIC_SPAN), 0);
    }

    /// The end is passed through 1-based-as-exclusive. Feeding a 0-based end instead would be
    /// off by one and would silently produce a different bin for reads ending on a boundary.
    #[test]
    fn a_read_ending_exactly_on_a_boundary_stays_in_the_lower_bin() {
        // 1-based start 1, 1-based inclusive end 16384: 0-based [0, 16384), still bin 4681.
        assert_eq!(compute_indexing_bin(1, 16_384).unwrap(), 4681);
        // One base further and it crosses.
        assert_eq!(compute_indexing_bin(1, 16_385).unwrap(), 585);
    }

    #[test]
    fn an_undeterminable_end_is_treated_as_one_base() {
        // start 1 (0-based 0), end unknown -> [0, 1).
        assert_eq!(compute_indexing_bin(1, 0).unwrap(), 4681);
        assert_eq!(compute_indexing_bin(1, -1).unwrap(), 4681);
    }

    #[test]
    fn positions_past_the_span_are_refused_not_wrapped() {
        assert!(compute_indexing_bin(BIN_GENOMIC_SPAN + 2, BIN_GENOMIC_SPAN + 10).is_err());
    }

    /// Every bin the scheme can produce must fit the 16-bit field it is written into.
    #[test]
    fn every_bin_fits_the_unsigned_short_field() {
        for start in (0..BIN_GENOMIC_SPAN).step_by(1_000_003) {
            for len in [1, 100, 20_000, 200_000, 2_000_000] {
                let end = (start + len).min(BIN_GENOMIC_SPAN);
                let b = region_to_bin(start, end);
                assert!(
                    (0..=MAX_BINS).contains(&b),
                    "bin {b} out of range for [{start},{end})"
                );
                assert_eq!(b & 0xFFFF, b, "bin {b} does not fit a ushort");
            }
        }
    }
}
