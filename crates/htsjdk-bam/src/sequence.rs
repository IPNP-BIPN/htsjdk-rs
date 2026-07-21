//! IUPAC base comparison.
//!
//! Ported from `htsjdk.samtools.util.SequenceUtil.basesEqual`, `isNoCall` and the static `bases`
//! table at htsjdk 4.2.0.
//!
//! The table maps each byte to a 4-bit mask of the nucleotides it stands for, and `basesEqual`
//! compares the **masks**, not the bytes. Three consequences, none of them obvious from the name:
//!
//!  1. Comparison is **case-insensitive**, because the static block copies every entry from `A`
//!     to `Z` into its lowercase counterpart. `a` equals `A`.
//!  2. `.` is given the same mask as `N`, so `.` equals `N` and equals `n`.
//!  3. Every byte *not* in the table keeps mask 0, and 0 equals 0. So two bytes that are not
//!     bases at all - `*` against `#`, say - compare **equal**, while either against a real base
//!     compares unequal. A reimplementation that compared uppercased bytes would agree on every
//!     valid input and disagree here.
//!
//! It is set equality and not set overlap: `M` is `A|C` and `A` is `A`, and those masks differ,
//! so `basesEqual(A, M)` is false. `basesEqualAmbiguous` is the overlap test, and Picard's
//! alignment metrics do not use it.

const A_MASK: u8 = 1;
const C_MASK: u8 = 2;
const G_MASK: u8 = 4;
const T_MASK: u8 = 8;

/// `SequenceUtil.BASES_ARRAY_LENGTH`.
const BASES_ARRAY_LENGTH: usize = 127;

/// The static `bases` table, built exactly as htsjdk's static block builds it.
fn bases_table() -> [u8; BASES_ARRAY_LENGTH] {
    let mut t = [0u8; BASES_ARRAY_LENGTH];
    t[b'A' as usize] = A_MASK;
    t[b'C' as usize] = C_MASK;
    t[b'G' as usize] = G_MASK;
    t[b'T' as usize] = T_MASK;
    t[b'M' as usize] = A_MASK | C_MASK;
    t[b'R' as usize] = A_MASK | G_MASK;
    t[b'W' as usize] = A_MASK | T_MASK;
    t[b'S' as usize] = C_MASK | G_MASK;
    t[b'Y' as usize] = C_MASK | T_MASK;
    t[b'K' as usize] = G_MASK | T_MASK;
    t[b'V' as usize] = A_MASK | C_MASK | G_MASK;
    t[b'H' as usize] = A_MASK | C_MASK | T_MASK;
    t[b'D' as usize] = A_MASK | G_MASK | T_MASK;
    t[b'B' as usize] = C_MASK | G_MASK | T_MASK;
    t[b'N' as usize] = A_MASK | C_MASK | G_MASK | T_MASK;
    // The lowercase copy runs over the whole A..Z range, so it also copies the zeros of the
    // letters that are not IUPAC codes. Then '.' is set, *after* the loop, so it has no
    // lowercase twin - which does not matter, since '.' has no case.
    let mut i = b'A';
    while i <= b'Z' {
        t[(i + 32) as usize] = t[i as usize];
        i += 1;
    }
    t[b'.' as usize] = A_MASK | C_MASK | G_MASK | T_MASK;
    t
}

static BASES: std::sync::LazyLock<[u8; BASES_ARRAY_LENGTH]> = std::sync::LazyLock::new(bases_table);

/// `SequenceUtil.basesEqual(lhs, rhs)`.
///
/// The bounds check is htsjdk's: a byte outside the table's range is **unequal to everything**,
/// including to another out-of-range byte. That is the one case where two unknown bytes do not
/// compare equal, and it turns on the table's length rather than on anything about bases.
pub fn bases_equal(lhs: u8, rhs: u8) -> bool {
    if lhs as usize >= BASES_ARRAY_LENGTH || rhs as usize >= BASES_ARRAY_LENGTH {
        return false;
    }
    BASES[lhs as usize] == BASES[rhs as usize]
}

/// `SequenceUtil.isNoCall(base)`.
///
/// Note this is a plain byte test and not a table lookup, so it is *not* the same predicate as
/// "has the N mask". `.` and `N` and `n` are no-calls; `B`, whose mask is also ambiguous, is not.
pub fn is_no_call(base: u8) -> bool {
    base == b'N' || base == b'n' || base == b'.'
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn comparison_ignores_case() {
        assert!(bases_equal(b'a', b'A'));
        assert!(bases_equal(b'g', b'G'));
        assert!(!bases_equal(b'a', b'g'));
    }

    #[test]
    fn dot_is_the_same_as_n() {
        assert!(bases_equal(b'.', b'N'));
        assert!(bases_equal(b'.', b'n'));
        assert!(!bases_equal(b'.', b'A'));
    }

    /// The mask test is equality, not overlap, so an ambiguity code does not match the bases it
    /// stands for.
    #[test]
    fn an_ambiguity_code_does_not_match_its_members() {
        assert!(!bases_equal(b'M', b'A'), "M is A|C, which is not A");
        assert!(!bases_equal(b'N', b'A'));
        assert!(bases_equal(b'M', b'm'));
    }

    /// Two bytes that are not bases share mask 0 and therefore compare equal. This is the case a
    /// byte-comparing reimplementation gets wrong in the direction of being *more* correct.
    #[test]
    fn two_non_bases_compare_equal_to_each_other() {
        assert!(bases_equal(b'*', b'#'));
        assert!(bases_equal(b'Z', b'Q'), "neither letter is an IUPAC code");
        assert!(!bases_equal(b'*', b'A'));
    }

    /// ...unless one of them is outside the table, where the bounds check rejects first.
    #[test]
    fn a_byte_past_the_table_is_equal_to_nothing() {
        assert!(!bases_equal(200, 200), "the same byte, still unequal");
        assert!(!bases_equal(200, b'A'));
        assert!(bases_equal(126, 125), "both inside the table, both mask 0");
    }

    /// `isNoCall` is a byte test, not a mask test, so it disagrees with the table.
    #[test]
    fn is_no_call_is_not_the_ambiguity_mask() {
        assert!(is_no_call(b'N') && is_no_call(b'n') && is_no_call(b'.'));
        assert!(!is_no_call(b'B'), "ambiguous, but not a no-call");
        assert!(
            bases_equal(b'.', b'N') && is_no_call(b'.') && is_no_call(b'N'),
            "here the two predicates happen to agree"
        );
    }
}
