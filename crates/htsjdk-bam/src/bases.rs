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
    match NIBBLE[base as usize] {
        NOT_A_BASE => None,
        n => Some(n),
    }
}

const NOT_A_BASE: u8 = 0xff;

// A table rather than a match: the match compiled to a branch per base and was three quarters
// of BAM record encoding.
const NIBBLE: [u8; 256] = {
    let mut t = [NOT_A_BASE; 256];
    let codes = b"=ACMGRSVTWYHKDBN";
    let mut i = 0;
    while i < 16 {
        t[codes[i] as usize] = i as u8;
        t[codes[i].to_ascii_lowercase() as usize] = i as u8;
        i += 1;
    }
    // htsjdk folds '.' onto N, in both the low and the high table.
    t[b'.' as usize] = 15;
    t
};

/// `SAMUtils.bytesToCompressedBases`.
///
/// The final nibble of an odd-length read is left as zero, because htsjdk allocates
/// `(len + 1) / 2` zeroed bytes and only ever writes the high nibble of the last one. That
/// zero is `=` when read back, so it is not a neutral padding value; it is only invisible
/// because the record's `l_seq` says to stop before it.
pub fn bytes_to_compressed_bases(read_bases: &[u8]) -> Result<Vec<u8>, BadBase> {
    let mut out = Vec::with_capacity(read_bases.len().div_ceil(2));
    write_compressed_bases(read_bases, &mut out)?;
    Ok(out)
}

/// The same packing, appended to a buffer the caller owns.
///
/// An encoder writes these bytes into the record it is building and drops them immediately after,
/// so the `Vec` above is a per-record allocation with no reader. This is the same loop without it,
/// and the refusal is identical: the first base that is not one of the sixteen, with its index.
pub fn write_compressed_bases(read_bases: &[u8], out: &mut Vec<u8>) -> Result<(), BadBase> {
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
        out.push((hi << 4) | lo);
    }
    Ok(())
}

/// `SAMUtils.compressedBasesToBytes`.
///
/// Not the inverse of [`bytes_to_compressed_bases`]: the encode table accepts lower case and
/// `.`, the decode table emits only upper case and `N`. Round-tripping a lowercase read
/// through a BAM changes it, in htsjdk exactly as here.
pub fn compressed_bases_to_bytes(length: usize, compressed: &[u8], offset: usize) -> Vec<u8> {
    let packed = &compressed[offset..offset + length.div_ceil(2)];
    let mut out = vec![0u8; length];
    // A byte at a time, both its bases from one table: the nibble loop was a hot spot on every read.
    let mut pairs = out.chunks_exact_mut(2);
    for (pair, &byte) in (&mut pairs).zip(packed) {
        pair.copy_from_slice(&BYTE_TO_BASES[byte as usize]);
    }
    if let [last] = pairs.into_remainder() {
        *last = NIBBLE_TO_BASE[(packed[length / 2] >> 4) as usize];
    }
    out
}

/// Both bases of each packed byte, high nibble first.
const BYTE_TO_BASES: [[u8; 2]; 256] = {
    let mut table = [[0u8; 2]; 256];
    let mut i = 0;
    while i < 256 {
        table[i] = [NIBBLE_TO_BASE[i >> 4], NIBBLE_TO_BASE[i & 0x0F]];
        i += 1;
    }
    table
};

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

    /// The table gives what the nibble loop gave, for every byte, from an offset, odd or even.
    #[test]
    fn the_pair_table_is_the_nibble_loop() {
        let packed: Vec<u8> = (0..=255u8).collect();
        for length in [0, 1, 2, 301, 510] {
            let naive: Vec<u8> = (0..length)
                .map(|i| {
                    let byte = packed[1 + i / 2];
                    NIBBLE_TO_BASE[usize::from(if i % 2 == 0 { byte >> 4 } else { byte & 0x0F })]
                })
                .collect();
            assert_eq!(compressed_bases_to_bytes(length, &packed, 1), naive);
        }
    }

    /// Decoding is lossy in a specific, matching way: case and `.` do not survive.
    #[test]
    fn round_trip_upper_cases_and_folds_dot_onto_n() {
        let packed = bytes_to_compressed_bases(b"acg.").unwrap();
        assert_eq!(compressed_bases_to_bytes(4, &packed, 0), b"ACGN");
    }
}
