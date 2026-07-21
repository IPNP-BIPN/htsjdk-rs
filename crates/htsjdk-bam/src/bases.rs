//! Nibble-packed sequence bases.
//!
//! Ported from `htsjdk.samtools.SAMUtils.bytesToCompressedBases` /
//! `compressedBasesToBytes` and the `charToCompressedBase{Low,High}` tables.
//!
//! The packing itself is specified by the BAM format. What is *not* specified, and what makes
//! this a port rather than a reimplementation, is the accepted input alphabet and what happens
//! to the trailing nibble of an odd-length read.

/// The 16 nibble values, in the order htsjdk assigns them.
///
/// This is `SAMUtils`' `COMPRESSED_*_LOW` constants read out in numeric order, and it is also
/// the decode table `COMPRESSED_BASES_TO_CHARS`. The two being the same array is what makes
/// encode/decode round-trip for uppercase input.
pub const NIBBLE_TO_BASE: [u8; 16] = *b"=ACMGRSVTWYHKDBN";

/// Encoding failed because the read contained a byte htsjdk does not accept.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BadBase {
    pub base: u8,
    pub index: usize,
}

/// `SAMUtils.charToCompressedBaseLow`, as a value in `0..16`.
///
/// The accepted alphabet is `=AaCcGgTtNnMmRrSsVvWwYyHhKkDdBb` plus `.`, which htsjdk folds
/// onto `N`. Everything else raises `IllegalArgumentException`; a port that silently
/// substituted `N` would produce a valid BAM with different bases.
pub fn base_to_nibble(base: u8) -> Option<u8> {
    let n = match base {
        b'=' => 0,
        b'A' | b'a' => 1,
        b'C' | b'c' => 2,
        b'M' | b'm' => 3,
        b'G' | b'g' => 4,
        b'R' | b'r' => 5,
        b'S' | b's' => 6,
        b'V' | b'v' => 7,
        b'T' | b't' => 8,
        b'W' | b'w' => 9,
        b'Y' | b'y' => 10,
        b'H' | b'h' => 11,
        b'K' | b'k' => 12,
        b'D' | b'd' => 13,
        b'B' | b'b' => 14,
        // htsjdk folds '.' onto N, in both the low and the high table.
        b'N' | b'n' | b'.' => 15,
        _ => return None,
    };
    Some(n)
}

/// `SAMUtils.bytesToCompressedBases`.
///
/// The final nibble of an odd-length read is left as zero, because htsjdk allocates
/// `(len + 1) / 2` zeroed bytes and only ever writes the high nibble of the last one. That
/// zero is `=` when read back, so it is not a neutral padding value; it is only invisible
/// because the record's `l_seq` says to stop before it.
pub fn bytes_to_compressed_bases(read_bases: &[u8]) -> Result<Vec<u8>, BadBase> {
    let mut out = vec![0u8; read_bases.len().div_ceil(2)];
    for (i, pair) in read_bases.chunks(2).enumerate() {
        let hi = base_to_nibble(pair[0]).ok_or(BadBase {
            base: pair[0],
            index: i * 2,
        })?;
        let lo = match pair.get(1) {
            Some(&b) => base_to_nibble(b).ok_or(BadBase {
                base: b,
                index: i * 2 + 1,
            })?,
            None => 0,
        };
        out[i] = (hi << 4) | lo;
    }
    Ok(out)
}

/// `SAMUtils.compressedBasesToBytes`.
///
/// Not the inverse of [`bytes_to_compressed_bases`]: the encode table accepts lower case and
/// `.`, the decode table emits only upper case and `N`. Round-tripping a lowercase read
/// through a BAM changes it, in htsjdk exactly as here.
pub fn compressed_bases_to_bytes(length: usize, compressed: &[u8], offset: usize) -> Vec<u8> {
    let mut out = vec![0u8; length];
    for i in 0..length {
        let byte = compressed[offset + i / 2];
        let nibble = if i % 2 == 0 { byte >> 4 } else { byte & 0x0F };
        out[i] = NIBBLE_TO_BASE[nibble as usize];
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_nibble_table_matches_the_java_constants() {
        // Spot-check the values that are actually named in SAMUtils.
        assert_eq!(base_to_nibble(b'=').unwrap(), 0);
        assert_eq!(base_to_nibble(b'A').unwrap(), 1);
        assert_eq!(base_to_nibble(b'C').unwrap(), 2);
        assert_eq!(base_to_nibble(b'G').unwrap(), 4);
        assert_eq!(base_to_nibble(b'T').unwrap(), 8);
        assert_eq!(base_to_nibble(b'N').unwrap(), 15);
        // And that the decode table is the same ordering.
        for (n, &b) in NIBBLE_TO_BASE.iter().enumerate() {
            assert_eq!(base_to_nibble(b).unwrap(), n as u8);
        }
    }

    #[test]
    fn an_even_read_packs_two_bases_per_byte() {
        assert_eq!(bytes_to_compressed_bases(b"ACGT").unwrap(), [0x12, 0x48]);
    }

    /// The trailing nibble is zero, which decodes as `=`. A port that padded with `N` (0x0F)
    /// would produce a file that reads back identically and hashes differently.
    #[test]
    fn an_odd_read_leaves_the_last_nibble_zero() {
        let packed = bytes_to_compressed_bases(b"ACG").unwrap();
        assert_eq!(packed, [0x12, 0x40]);
        assert_eq!(packed[1] & 0x0F, 0, "padding must be 0 (=), not 0x0F (N)");
    }

    #[test]
    fn empty_read_packs_to_nothing() {
        assert_eq!(bytes_to_compressed_bases(b"").unwrap(), Vec::<u8>::new());
    }

    #[test]
    fn lower_case_and_dot_are_accepted_on_the_way_in() {
        assert_eq!(
            bytes_to_compressed_bases(b"acgt").unwrap(),
            bytes_to_compressed_bases(b"ACGT").unwrap()
        );
        assert_eq!(base_to_nibble(b'.').unwrap(), base_to_nibble(b'N').unwrap());
    }

    /// An unknown base is refused, not substituted. This is the difference between a port and
    /// a lenient reimplementation.
    #[test]
    fn an_unknown_base_is_refused() {
        assert_eq!(
            bytes_to_compressed_bases(b"ACXT"),
            Err(BadBase {
                base: b'X',
                index: 2
            })
        );
    }

    #[test]
    fn round_trip_holds_for_upper_case() {
        for read in [
            &b"ACGT"[..],
            b"A",
            b"ACGTN=MRSVWYHKDB",
            b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
        ] {
            let packed = bytes_to_compressed_bases(read).unwrap();
            assert_eq!(compressed_bases_to_bytes(read.len(), &packed, 0), read);
        }
    }

    /// Decoding is lossy in a specific, matching way: case and `.` do not survive.
    #[test]
    fn round_trip_upper_cases_and_folds_dot_onto_n() {
        let packed = bytes_to_compressed_bases(b"acg.").unwrap();
        assert_eq!(compressed_bases_to_bytes(4, &packed, 0), b"ACGN");
    }
}
