//! commons-math3's `MersenneTwister`, and the permutation a caller draws from it.
//!
//! MT19937 with commons-math3's own seeding, which is not the reference algorithm's: an int seed
//! is spread through the state with Knuth's multiplier and the array-seeding step the original
//! runs afterwards is skipped, so the stream differs from a textbook implementation seeded the
//! same way.
//!
//! Every accessor consumes its own words, which makes the ORDER of calls part of the answer: a
//! `next_long` takes two draws, a `next_boolean` takes a whole one, and a `next_double` takes a
//! 26-bit draw and a 27-bit one.
//!
//! Ported from `org.apache.commons.math3.random.MersenneTwister` and
//! `org.apache.commons.math3.random.RandomDataGenerator` (commons-math3 3.5).

const N: usize = 624;
const M: usize = 397;
const MAGIC: [u32; 2] = [0x0, 0x9908b0df];
/// The multiplier Knuth's initialisation uses, which commons-math3 spells `1812433253`.
const MAGIC_FACTOR: u32 = 1812433253;

/// The generator's state: the words and the index of the next one.
#[derive(Debug, Clone)]
pub struct MersenneTwister {
    mt: [u32; N],
    mti: usize,
}

impl MersenneTwister {
    /// `new MersenneTwister(int seed)`.
    ///
    /// The state is filled from the seed with `mt[i] = 1812433253 * (mt[i-1] ^ (mt[i-1] >> 30)) + i`
    /// and nothing else happens: commons-math3 stops here, where the reference algorithm would go
    /// on to mix an array of seeds into the state.
    pub fn new(seed: i32) -> Self {
        // The reference carries the running value in a SIGNED long and masks it back to
        // thirty-two bits each round, so the shift by thirty is arithmetic: for a negative seed
        // the first round sees `-1 >> 30 == -1` where an unsigned shift would see a large
        // positive. That is the whole of why seeds 0 and -1 differ in one value and agree in
        // every other.
        let mut running = i64::from(seed);
        let mut mt = [0u32; N];
        mt[0] = running as u32;
        for (i, word) in mt.iter_mut().enumerate().skip(1) {
            running = (i64::from(MAGIC_FACTOR)
                .wrapping_mul(running ^ (running >> 30))
                .wrapping_add(i as i64))
                & 0xffff_ffff;
            *word = running as u32;
        }
        MersenneTwister { mt, mti: N }
    }

    /// `next(int bits)`: the tempered word, kept to its top `bits`.
    fn next(&mut self, bits: u32) -> u32 {
        if self.mti >= N {
            self.twist();
        }
        let mut y = self.mt[self.mti];
        self.mti += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c5680;
        y ^= (y << 15) & 0xefc60000;
        y ^= y >> 18;
        y >> (32 - bits)
    }

    /// The generation of the next block of `N` words.
    fn twist(&mut self) {
        let mut mt = self.mt;
        for i in 0..N {
            let x = (mt[i] & 0x80000000) | (mt[(i + 1) % N] & 0x7fffffff);
            let mut next = mt[(i + M) % N] ^ (x >> 1);
            next ^= MAGIC[(x & 0x1) as usize];
            mt[i] = next;
        }
        self.mt = mt;
        self.mti = 0;
    }

    /// `nextInt()`, which is the whole tempered word read as a signed integer.
    pub fn next_int(&mut self) -> i32 {
        self.next(32) as i32
    }

    /// `nextInt(n)`, which rejects and redraws rather than taking a remainder.
    ///
    /// A power of two takes the HIGH bits instead, which is a different value from the same draw.
    pub fn next_int_bounded(&mut self, n: i32) -> i32 {
        if n <= 0 {
            return 0;
        }
        if (n & -n) == n {
            return ((i64::from(n) * i64::from(self.next(31))) >> 31) as i32;
        }
        loop {
            let bits = self.next(31) as i32;
            let value = bits % n;
            if bits - value + (n - 1) >= 0 {
                return value;
            }
        }
    }

    /// `nextDouble()`, built from two 26-bit draws and scaled by `0x1.0p-52`.
    pub fn next_double(&mut self) -> f64 {
        let high = i64::from(self.next(26)) << 26;
        let low = i64::from(self.next(26));
        ((high | low) as f64) * f64::from_bits(0x3CB0_0000_0000_0000)
    }

    /// `nextLong()`, which is two words with the first shifted left by thirty-two.
    ///
    /// The first is widened as a SIGNED int and then shifted, so a negative first word makes the
    /// high half all ones above the shift; the second is masked to its low thirty-two bits.
    pub fn next_long(&mut self) -> i64 {
        let high = (i64::from(self.next(32) as i32)) << 32;
        let low = i64::from(self.next(32) as i32) & 0xffff_ffff;
        high | low
    }

    /// `nextBoolean()`, which consumes a whole word for one bit.
    pub fn next_boolean(&mut self) -> bool {
        self.next(1) != 0
    }

    /// `nextFloat()`, which consumes a whole word for twenty-three bits of mantissa.
    pub fn next_float(&mut self) -> f32 {
        self.next(23) as f32 * f32::from_bits(0x3400_0000)
    }
}

/// `RandomDataGenerator.nextPermutation(n, k)`: a partial Fisher-Yates over `natural(n)`, of which
/// the first `k` entries are the answer.
pub fn next_permutation(rng: &mut MersenneTwister, n: usize, k: usize) -> Option<Vec<usize>> {
    if k > n || k == 0 {
        return None;
    }
    let mut index: Vec<usize> = (0..n).collect();
    // The WHOLE array is shuffled and the first k entries are then copied out, so a partial
    // permutation consumes as many draws as a full one: asking for three of ten and asking for
    // ten of ten leave the generator in the same place.
    shuffle(rng, &mut index);
    Some(index[..k].to_vec())
}

/// `MathArrays.shuffle(list, rng)`: from the end towards the front, each entry swapped with one
/// drawn from `UniformIntegerDistribution(rng, 0, i)`, which is `nextInt(i + 1)`.
fn shuffle(rng: &mut MersenneTwister, list: &mut [usize]) {
    for i in (0..list.len()).rev() {
        let target = if i == 0 {
            0
        } else {
            rng.next_int_bounded(i as i32 + 1) as usize
        };
        list.swap(i, target);
    }
}

/// `MathUtil.permute`: an INDEX map and not a destination map, `out[i] = in[perm[i]]`.
pub fn permute(values: &[f64], permutation: &[usize]) -> Vec<f64> {
    permutation.iter().map(|index| values[*index]).collect()
}
