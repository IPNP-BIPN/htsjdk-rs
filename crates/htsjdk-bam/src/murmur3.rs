//! The Murmur3_32 hash, as htsjdk uses it for downsampling.
//!
//! Ports `htsjdk.samtools.util.Murmur3.hashUnencodedChars` at tag 4.2.0 (a Guava-derived
//! MurmurHash3, public domain). `DownsampleSam`'s constant-memory strategy hashes each read name to a
//! 32-bit value and keeps the read when the hash falls below a probability-derived threshold, so the
//! keep/discard decision is a pure function of the read name and the seed. Reproducing that decision
//! byte-for-byte requires reproducing this hash exactly, which means matching Java's 32-bit **wrapping**
//! integer arithmetic and its **UTF-16** view of the string (`charAt` is a code unit, so a non-ASCII
//! name hashes over its surrogate units, not its scalar values).

/// `c1` in `mixK1`.
const C1: i32 = 0xcc9e2d51u32 as i32;
/// `c2` in `mixK1`.
const C2: i32 = 0x1b873593;

/// A seeded Murmur3_32 hasher. `htsjdk.samtools.util.Murmur3`.
#[derive(Debug, Clone, Copy)]
pub struct Murmur3 {
    seed: i32,
}

impl Murmur3 {
    /// `new Murmur3(seed)`.
    pub fn new(seed: i32) -> Self {
        Murmur3 { seed }
    }

    /// `hashUnencodedChars(input)`: the Murmur3_32 hash of the string's UTF-16 code units.
    pub fn hash_unencoded_chars(&self, input: &str) -> i32 {
        let units: Vec<u16> = input.encode_utf16().collect();
        let length = units.len();
        let mut h1 = self.seed;

        // Step through the code units two at a time, low unit in the low half.
        let mut i = 1;
        while i < length {
            let k1 = ((units[i - 1] as u32) | ((units[i] as u32) << 16)) as i32;
            h1 = mix_h1(h1, mix_k1(k1));
            i += 2;
        }

        // A trailing odd unit.
        if length & 1 == 1 {
            let k1 = mix_k1(units[length - 1] as i32);
            h1 ^= k1;
        }

        fmix(h1, 2 * length as i32)
    }
}

fn mix_k1(mut k1: i32) -> i32 {
    k1 = k1.wrapping_mul(C1);
    k1 = k1.rotate_left(15);
    k1 = k1.wrapping_mul(C2);
    k1
}

fn mix_h1(mut h1: i32, k1: i32) -> i32 {
    h1 ^= k1;
    h1 = h1.rotate_left(13);
    h1 = h1.wrapping_mul(5).wrapping_add(0xe6546b64u32 as i32);
    h1
}

/// Finalization mix: force all bits of a hash block to avalanche.
fn fmix(mut h1: i32, length: i32) -> i32 {
    h1 ^= length;
    h1 ^= ((h1 as u32) >> 16) as i32;
    h1 = h1.wrapping_mul(0x85ebca6bu32 as i32);
    h1 ^= ((h1 as u32) >> 13) as i32;
    h1 = h1.wrapping_mul(0xc2b2ae35u32 as i32);
    h1 ^= ((h1 as u32) >> 16) as i32;
    h1
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Values dumped from htsjdk's `Murmur3` at tag 4.2.0, covering the empty string, odd and even
    /// lengths, a realistic read name, and a non-ASCII name (which exercises the UTF-16 view).
    #[test]
    fn hashes_match_htsjdk() {
        let cases: &[(i32, &str, i32)] = &[
            (1, "", 1_364_076_727),
            (1, "A", 116_685_281),
            (1, "read0", 1_137_924_452),
            (1, "read1", -1_253_081_825),
            (1, "read12345", 587_747_829),
            (1, "HWI-ST807:461:C2P0JACXX:4:2107:8467:34718", 29_913_811),
            (1, "café", -40_327_187),
            (0, "", 0),
            (0, "A", 814_537_616),
            (0, "read0", 907_747_536),
            (42, "read1", -958_943_510),
            (-7, "read12345", 1_571_483_815),
            (-7, "café", -1_841_479_757),
        ];
        for &(seed, input, expected) in cases {
            assert_eq!(
                Murmur3::new(seed).hash_unencoded_chars(input),
                expected,
                "seed={seed} input={input:?}"
            );
        }
    }
}
