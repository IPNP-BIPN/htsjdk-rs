//! The match-finding side of deflate: the window, the hash chains, and the two block functions.
//!
//! Ported from zlib's `deflate.c` (zlib licence, Jean-loup Gailly and Mark Adler).
//!
//! This is a **one-shot** deflater: the whole input is present from the start and the output grows
//! as needed, which is exactly how BGZF uses it (one block, `setInput`, `finish`, one `deflate`
//! call). zlib's streaming state machine, its `need_more`/`avail_out` bookkeeping and its
//! `deflate_stored` path all exist to serve callers that feed and drain in pieces, and none of
//! that changes the bytes a one-shot caller gets. What is kept is everything that does: the hash,
//! the chain order, the lazy-match rules, and the window slide.

use crate::trees::Trees;

const MIN_MATCH: usize = 3;
const MAX_MATCH: usize = 258;
const MIN_LOOKAHEAD: usize = MAX_MATCH + MIN_MATCH + 1;
const TOO_FAR: usize = 4096;
const NIL: u16 = 0;

/// The window is 15 bits and the memory level 8, which is what `deflateInit2_` is called with in
/// both htsjdk's `Deflater` and GKL's `IntelDeflater`. Nothing here reads a different value, so
/// they are constants rather than fields: a hidden second configuration is exactly the kind of
/// thing that produces bytes nobody can explain.
const W_BITS: usize = 15;
const W_SIZE: usize = 1 << W_BITS;
const W_MASK: usize = W_SIZE - 1;
const MEM_LEVEL: usize = 8;
const HASH_BITS: usize = MEM_LEVEL + 7;
const HASH_SIZE: usize = 1 << HASH_BITS;
const HASH_MASK: usize = HASH_SIZE - 1;
/// zlib writes this as `(hash_bits + MIN_MATCH - 1) / MIN_MATCH`, which is this rounding up.
const HASH_SHIFT: usize = HASH_BITS.div_ceil(MIN_MATCH);
const LIT_BUFSIZE: usize = 1 << (MEM_LEVEL + 6);
const MAX_DIST: usize = W_SIZE - MIN_LOOKAHEAD;
/// zlib's `window_size`, which is what `fill_window` measures free space against. The buffer is
/// allocated larger than this so that `longest_match` can read a few bytes past the data without
/// an out-of-bounds panic, exactly as the C reads past it without an out-of-bounds fault. Those
/// bytes never reach the output: a match is clamped to `lookahead`.
const WINDOW_SIZE: usize = 2 * W_SIZE;
const WINDOW_PADDING: usize = MAX_MATCH + 8;

/// How many of `n` leading bytes `a` and `b` share: the byte loop of zlib's `longest_match`,
/// compared a word at a time.
///
/// # Safety
/// `a` and `b` must each be valid for reads of `n` bytes.
#[inline(always)]
unsafe fn common_prefix_raw(a: *const u8, b: *const u8, n: usize) -> usize {
    let mut i = 0;
    while i + 8 <= n {
        let x = a.add(i).cast::<u64>().read_unaligned() ^ b.add(i).cast::<u64>().read_unaligned();
        if x != 0 {
            return i + (u64::from_le(x).trailing_zeros() / 8) as usize;
        }
        i += 8;
    }
    while i < n && *a.add(i) == *b.add(i) {
        i += 1;
    }
    i
}

#[cfg(test)]
#[test]
fn common_prefix_is_the_byte_loop() {
    let mut seed = 7u64;
    for _ in 0..20_000 {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let len = (seed >> 33) as usize % 300;
        let a: Vec<u8> = (0..len).map(|i| (i % 3) as u8).collect();
        let mut b = a.clone();
        if len > 0 && seed & 1 == 0 {
            b[(seed >> 40) as usize % len] ^= 1;
        }
        let expected = a.iter().zip(&b).take_while(|(x, y)| x == y).count();
        assert_eq!(
            unsafe { common_prefix_raw(a.as_ptr(), b.as_ptr(), len) },
            expected
        );
    }
}

/// zlib's `configuration_table`, verbatim. `good_length` and `max_lazy` mean different things on
/// the fast path and the slow one, which is why the table is read through [`Config`] rather than
/// inlined.
#[derive(Clone, Copy)]
struct Config {
    good_length: usize,
    max_lazy: usize,
    nice_length: usize,
    max_chain: usize,
    slow: bool,
}

const CONFIGURATION_TABLE: [Config; 10] = [
    Config {
        good_length: 0,
        max_lazy: 0,
        nice_length: 0,
        max_chain: 0,
        slow: false,
    },
    Config {
        good_length: 4,
        max_lazy: 4,
        nice_length: 8,
        max_chain: 4,
        slow: false,
    },
    Config {
        good_length: 4,
        max_lazy: 5,
        nice_length: 16,
        max_chain: 8,
        slow: false,
    },
    Config {
        good_length: 4,
        max_lazy: 6,
        nice_length: 32,
        max_chain: 32,
        slow: false,
    },
    Config {
        good_length: 4,
        max_lazy: 4,
        nice_length: 16,
        max_chain: 16,
        slow: true,
    },
    Config {
        good_length: 8,
        max_lazy: 16,
        nice_length: 32,
        max_chain: 32,
        slow: true,
    },
    Config {
        good_length: 8,
        max_lazy: 16,
        nice_length: 128,
        max_chain: 128,
        slow: true,
    },
    Config {
        good_length: 8,
        max_lazy: 32,
        nice_length: 128,
        max_chain: 256,
        slow: true,
    },
    Config {
        good_length: 32,
        max_lazy: 128,
        nice_length: 258,
        max_chain: 1024,
        slow: true,
    },
    Config {
        good_length: 32,
        max_lazy: 258,
        nice_length: 258,
        max_chain: 4096,
        slow: true,
    },
];

/// Which deflater is being reproduced.
///
/// The two differ in three ways that all change the bytes: which block function each level uses,
/// which hash feeds the chains, and, on the Intel side, whether the CPU supports SSE4.2 at all.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Flavour {
    /// Stock zlib, which is what `java.util.zip.Deflater` routes to. Levels 1 to 3 fast, 4 to 9
    /// slow, multiplicative rolling hash throughout.
    Jdk,
    /// Intel's zlib 1.2.13 fork inside `libgkl_compression.so`, the one GKL uses at levels 3 to 9.
    /// Adds `deflate_medium` at levels 4 to 6 and swaps the hash for a CRC-32C of the bytes at the
    /// position, which is not a rolling hash at all.
    ///
    /// `sse42` is the value of the fork's `x86_cpu_has_sse42`, read from CPUID at load time. **It
    /// changes the output**, because the two hashes fill the chains differently, so a GKL claim
    /// is a claim about a CPU as much as about a level. The default is `true`: every host the
    /// oracle has run on reports SSE4.2, and that is the column the goldens were measured in.
    Gkl { sse42: bool },
}

/// The hash a flavour uses to index the chains.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Hash {
    /// zlib's `UPDATE_HASH`: `h = ((h << 5) ^ c) & mask`, rolling over three bytes.
    Multiplicative,
    /// Intel's `UPDATE_HASH_CRC`: a CRC-32C of the four bytes at the position, or of three when
    /// the level is 6 or above. Not rolling: it depends only on the window, never on the previous
    /// value, which is why the priming assignments around it are skipped.
    Crc32c { three_byte: bool },
}

/// The CRC-32C the SSE4.2 `crc32` instruction computes: reflected, polynomial 0x1EDC6F41, and no
/// final inversion. `_mm_crc32_u32(0, val)` processes `val`'s four bytes least-significant first.
#[inline]
fn crc32c_u32(val: u32) -> u32 {
    // The instruction where the CPU has it: Arm's `crc32cw` is the same CRC as SSE4.2's.
    #[cfg(target_arch = "aarch64")]
    if std::arch::is_aarch64_feature_detected!("crc") {
        // SAFETY: the CPU was just found to have the instruction.
        return unsafe { crc32c_u32_arm(val) };
    }
    #[cfg(target_arch = "x86_64")]
    if std::arch::is_x86_feature_detected!("sse4.2") {
        // SAFETY: the CPU was just found to have the instruction.
        return unsafe { crc32c_u32_x86(val) };
    }
    crc32c_u32_table(val)
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "crc")]
unsafe fn crc32c_u32_arm(val: u32) -> u32 {
    std::arch::aarch64::__crc32cw(0, val)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse4.2")]
unsafe fn crc32c_u32_x86(val: u32) -> u32 {
    std::arch::x86_64::_mm_crc32_u32(0, val)
}

/// [`crc32c_u32`] without the instruction.
fn crc32c_u32_table(val: u32) -> u32 {
    // Slicing-by-4: the same value as feeding the four bytes through the bitwise CRC.
    CRC32C[3][(val & 0xFF) as usize]
        ^ CRC32C[2][((val >> 8) & 0xFF) as usize]
        ^ CRC32C[1][((val >> 16) & 0xFF) as usize]
        ^ CRC32C[0][(val >> 24) as usize]
}

/// `CRC32C[0]` is the byte table of the reflected Castagnoli polynomial; `CRC32C[k]` is a byte
/// followed by `k` zero bytes.
const CRC32C: [[u32; 256]; 4] = {
    let mut t = [[0u32; 256]; 4];
    let mut i = 0;
    while i < 256 {
        let mut crc = i as u32;
        let mut k = 0;
        while k < 8 {
            crc = if crc & 1 != 0 {
                (crc >> 1) ^ 0x82F6_3B78
            } else {
                crc >> 1
            };
            k += 1;
        }
        t[0][i] = crc;
        i += 1;
    }
    let mut i = 0;
    while i < 256 {
        let mut k = 1;
        while k < 4 {
            let prev = t[k - 1][i];
            t[k][i] = (prev >> 8) ^ t[0][(prev & 0xFF) as usize];
            k += 1;
        }
        i += 1;
    }
    t
};

#[cfg(test)]
#[test]
fn crc32c_u32_is_the_bitwise_crc() {
    let bitwise = |val: u32| {
        let mut crc = 0u32;
        for byte in val.to_le_bytes() {
            crc ^= byte as u32;
            for _ in 0..8 {
                crc = if crc & 1 != 0 {
                    (crc >> 1) ^ 0x82F6_3B78
                } else {
                    crc >> 1
                };
            }
        }
        crc
    };
    let mut v = 0x1234_5678u32;
    for _ in 0..200_000 {
        v = v.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
        assert_eq!(crc32c_u32(v), bitwise(v));
        assert_eq!(crc32c_u32(v & 0x00FF_FFFF), bitwise(v & 0x00FF_FFFF));
        assert_eq!(crc32c_u32_table(v), bitwise(v));
    }
}

/// One match, as `deflate_medium.c` carries it: where it starts in the window, where it matched,
/// how long it is, and how far its span has already been inserted into the chains.
#[derive(Clone, Copy, Default)]
struct Match {
    match_start: usize,
    match_length: usize,
    strstart: usize,
    orgstart: usize,
}

/// The whole state of a deflater, and every buffer it owns: one per thread, reset between blocks.
pub struct Deflater {
    next_in: usize,

    window: Box<[u8; WINDOW_SIZE + WINDOW_PADDING]>,
    prev: Box<[u16; W_SIZE]>,
    head: Box<[u16; HASH_SIZE]>,

    ins_h: usize,
    strstart: usize,
    lookahead: usize,
    /// Where the current block starts in the window, or negative once a slide has moved it out.
    block_start: isize,
    match_start: usize,
    match_length: usize,
    prev_length: usize,
    prev_match: usize,
    match_available: bool,
    insert: usize,

    config: Config,
    hash: Hash,
    /// True for Intel's levels 4 to 6, where `deflate_medium` replaces `deflate_slow`.
    medium: bool,
    trees: Trees,
    pub(crate) out: Vec<u8>,
}

/// A mask that keeps a window index in the allocation, so the hot reads carry no bounds check. It
/// changes nothing: every index in range is already below `WINDOW_SIZE`.
const WINDOW_INDEX: usize = WINDOW_SIZE - 1;

/// [`Deflater::insert_string`] on the deflate loops' local hash.
#[inline(always)]
fn insert_at(
    hash: Hash,
    window: &[u8; WINDOW_SIZE + WINDOW_PADDING],
    head: &mut [u16; HASH_SIZE],
    prev: &mut [u16; W_SIZE],
    ins_h: &mut usize,
    str_pos: usize,
) -> u16 {
    // SAFETY: `str_pos` is below `WINDOW_SIZE - MIN_MATCH`, `fill_window`'s invariant, and the
    // window holds `WINDOW_PADDING` bytes more.
    match hash {
        Hash::Multiplicative => {
            let c = unsafe { *window.get_unchecked(str_pos + MIN_MATCH - 1) };
            *ins_h = ((*ins_h << HASH_SHIFT) ^ c as usize) & HASH_MASK;
        }
        Hash::Crc32c { three_byte } => {
            let mut val = unsafe { window.as_ptr().add(str_pos).cast::<u32>().read_unaligned() };
            val = u32::from_le(val);
            if three_byte {
                val &= 0x00FF_FFFF;
            }
            *ins_h = crc32c_u32(val) as usize & HASH_MASK;
        }
    }
    let at = *ins_h & HASH_MASK;
    let old = head[at];
    prev[str_pos & W_MASK] = old;
    head[at] = str_pos as u16;
    old
}

/// A zeroed array on the heap, without passing through the stack on the way there.
fn zeroed<T: Copy + Default, const N: usize>() -> Box<[T; N]> {
    vec![T::default(); N]
        .into_boxed_slice()
        .try_into()
        .unwrap_or_else(|_| unreachable!("the vector has N elements"))
}

impl Deflater {
    pub fn new(level: usize, flavour: Flavour) -> Self {
        let (hash, medium) = Self::rules(level, flavour);
        Deflater {
            next_in: 0,
            window: zeroed(),
            prev: zeroed(),
            head: zeroed(),
            ins_h: 0,
            strstart: 0,
            lookahead: 0,
            block_start: 0,
            match_start: 0,
            match_length: MIN_MATCH - 1,
            prev_length: MIN_MATCH - 1,
            prev_match: 0,
            match_available: false,
            insert: 0,
            config: CONFIGURATION_TABLE[level],
            hash,
            medium,
            trees: Trees::new(LIT_BUFSIZE),
            out: Vec::new(),
        }
    }

    fn rules(level: usize, flavour: Flavour) -> (Hash, bool) {
        match flavour {
            Flavour::Jdk => (Hash::Multiplicative, false),
            Flavour::Gkl { sse42: false } => (Hash::Multiplicative, (4..=6).contains(&level)),
            Flavour::Gkl { sse42: true } => (
                Hash::Crc32c {
                    three_byte: level >= 6,
                },
                (4..=6).contains(&level),
            ),
        }
    }

    /// Back to the state [`Self::new`] answers. `prev` is left as it is, as zlib's `lm_init` leaves
    /// it: a chain only reaches an entry through a position hashed since the reset. The window is
    /// cleared because `longest_match` reads past the data, and those bytes must be a fresh
    /// deflater's zeros for a reused one to answer as a fresh one does.
    pub fn reset(&mut self) {
        self.window.fill(0);
        self.head.fill(NIL);
        self.next_in = 0;
        self.ins_h = 0;
        self.strstart = 0;
        self.lookahead = 0;
        self.block_start = 0;
        self.match_start = 0;
        self.match_length = MIN_MATCH - 1;
        self.prev_length = MIN_MATCH - 1;
        self.prev_match = 0;
        self.match_available = false;
        self.insert = 0;
        self.trees.reset();
        self.out.clear();
    }

    /// Deflates `input` into [`Self::out`], from a state that is new or was just reset.
    pub fn run(&mut self, input: &[u8]) {
        self.out.reserve(input.len() / 2 + 64);
        if self.medium {
            self.deflate_medium(input);
        } else if self.config.slow {
            self.deflate_slow(input);
        } else {
            self.deflate_fast(input);
        }
    }

    /// zlib's `UPDATE_HASH`, and only reachable on the multiplicative path. Intel's macro has the
    /// same name and the same call sites but ignores the running value entirely, which is why the
    /// hash is recomputed from `str_pos` rather than fed a byte.
    #[inline]
    fn update_hash(&mut self, c: u8) {
        self.ins_h = ((self.ins_h << HASH_SHIFT) ^ c as usize) & HASH_MASK;
    }

    /// The hash of the string starting at `str_pos`, by whichever rule this flavour uses.
    ///
    /// Intel's macro reaches the position by subtracting `MIN_MATCH - 1` from the address of the
    /// byte it was handed, which is how a rolling-hash call site ends up computing a positional
    /// hash. The four-byte load is a plain little-endian read of the window.
    #[inline]
    fn hash_at(&mut self, str_pos: usize) {
        match self.hash {
            Hash::Multiplicative => self.update_hash(self.window[str_pos + MIN_MATCH - 1]),
            Hash::Crc32c { three_byte } => {
                let mut val =
                    u32::from_le_bytes(self.window[str_pos..str_pos + 4].try_into().unwrap());
                if three_byte {
                    val &= 0x00FF_FFFF;
                }
                self.ins_h = crc32c_u32(val) as usize & HASH_MASK;
            }
        }
    }

    /// Ported from zlib's `INSERT_STRING`. Returns the previous head of this hash's chain, which
    /// is where a match search starts.
    #[inline]
    fn insert_string(&mut self, str_pos: usize) -> u16 {
        self.hash_at(str_pos);
        let at = self.ins_h & HASH_MASK;
        let head = self.head[at];
        self.prev[str_pos & W_MASK] = head;
        self.head[at] = str_pos as u16;
        head
    }

    /// Ported from zlib's `slide_hash`. Every recorded position moves down by one window, and
    /// anything that falls off the bottom becomes NIL.
    fn slide_hash(&mut self) {
        // `m >= W_SIZE ? m - W_SIZE : NIL` is a saturating subtraction, which vectorizes.
        for slot in self.head.iter_mut() {
            *slot = slot.saturating_sub(W_SIZE as u16);
        }
        for slot in self.prev.iter_mut() {
            *slot = slot.saturating_sub(W_SIZE as u16);
        }
    }

    /// Ported from zlib's `fill_window`. The slide is the part that matters for byte-identity:
    /// it happens at a fixed point, and `insert` carries across it the positions that were read
    /// but never hashed.
    fn fill_window(&mut self, input: &[u8]) {
        loop {
            let mut more = WINDOW_SIZE - self.lookahead - self.strstart;

            if self.strstart >= W_SIZE + MAX_DIST {
                self.window.copy_within(W_SIZE..WINDOW_SIZE - more, 0);
                // zlib's `match_start` is unsigned and wraps here when the last match lies below
                // the slide. A release build already wrapped, so this changes no byte; a debug
                // build panicked on overflow instead.
                self.match_start = self.match_start.wrapping_sub(W_SIZE);
                self.strstart -= W_SIZE;
                self.block_start -= W_SIZE as isize;
                if self.insert > self.strstart {
                    self.insert = self.strstart;
                }
                self.slide_hash();
                more += W_SIZE;
            }
            if self.next_in == input.len() {
                break;
            }

            let n = more.min(input.len() - self.next_in);
            let dst = self.strstart + self.lookahead;
            self.window[dst..dst + n].copy_from_slice(&input[self.next_in..self.next_in + n]);
            self.next_in += n;
            self.lookahead += n;

            // Hash the positions that were already in the window but too close to its end to be
            // hashed then. Without this the chains would have holes and the matches would differ.
            if self.lookahead + self.insert >= MIN_MATCH {
                let mut str_pos = self.strstart - self.insert;
                // Priming the rolling value, which the CRC hash does not have. Intel guards this
                // with `if (!x86_cpu_has_sse42)` for exactly that reason.
                if self.hash == Hash::Multiplicative {
                    self.ins_h = self.window[str_pos] as usize;
                    self.update_hash(self.window[str_pos + 1]);
                }
                while self.insert != 0 {
                    self.hash_at(str_pos);
                    self.prev[str_pos & W_MASK] = self.head[self.ins_h];
                    self.head[self.ins_h] = str_pos as u16;
                    str_pos += 1;
                    self.insert -= 1;
                    if self.lookahead + self.insert < MIN_MATCH {
                        break;
                    }
                }
            }

            if self.lookahead >= MIN_LOOKAHEAD || self.next_in == input.len() {
                break;
            }
        }
    }

    /// Ported from zlib's `longest_match`. The two `UNALIGNED_OK` variants in the C differ only in
    /// how fast they compare; this is the plain one, and it finds the same match.
    #[inline]
    fn longest_match(&mut self, mut cur_match: usize) -> usize {
        let mut chain_length = self.config.max_chain;
        let mut best_len = self.prev_length;
        let mut nice_match = self.config.nice_length;
        // zlib's NIL is 0, and a saturating subtraction lands on it for the same reason.
        let limit = self.strstart.saturating_sub(MAX_DIST);
        let scan = self.strstart;
        let window = self.window.as_ptr();
        let prev: &[u16; W_SIZE] = &self.prev;
        let mut match_start = None;

        // Already sitting on a good match: spend a quarter of the budget looking for a better one.
        if self.prev_length >= self.config.good_length {
            chain_length >>= 2;
        }
        // Never look past the end of the input, which is what makes the output deterministic.
        if nice_match > self.lookahead {
            nice_match = self.lookahead;
        }

        // SAFETY for every read below. `scan` is below `WINDOW_SIZE - MIN_LOOKAHEAD`, which is
        // `fill_window`'s invariant; `cur_match` and `best_len` are each below `scan` and
        // `MAX_MATCH`; so the furthest byte read, `scan + MAX_MATCH` and its 16-bit neighbour,
        // lies inside the window's `WINDOW_PADDING` bytes of slack.
        let pair = |at: usize| unsafe { window.add(at).cast::<u16>().read_unaligned() };
        // The cheap rejections zlib makes before comparing anything: the last byte of the current
        // best and the one past it, then the first two bytes, each pair one 16-bit compare as
        // `UNALIGNED_OK` makes it. A best of 0 has no last byte, so only the one past it is asked.
        let end_of = |best_len: usize| best_len.saturating_sub(1);
        let mut scan_end = pair(scan + end_of(best_len));
        let scan_start = pair(scan);
        loop {
            let m = cur_match;
            let end_matches = if best_len == 0 {
                unsafe { *window.add(m) == *window.add(scan) }
            } else {
                pair(m + best_len - 1) == scan_end
            };
            if end_matches && pair(m) == scan_start {
                // scan[2] == match[2] is implied by the hash, so the comparison starts at 3 and
                // stops at the first differing byte or at strend, eight bytes at a time.
                let len = 3 + unsafe {
                    common_prefix_raw(window.add(scan + 3), window.add(m + 3), MAX_MATCH - 3)
                };
                if len > best_len {
                    match_start = Some(cur_match);
                    best_len = len;
                    if len >= nice_match {
                        break;
                    }
                    scan_end = pair(scan + end_of(best_len));
                }
            }
            cur_match = prev[cur_match & W_MASK] as usize;
            chain_length -= 1;
            if cur_match <= limit || chain_length == 0 {
                break;
            }
        }
        if let Some(start) = match_start {
            self.match_start = start;
        }

        if best_len <= self.lookahead {
            best_len
        } else {
            self.lookahead
        }
    }

    /// Close the current block. `buf` is `None` once a slide has pushed the block's start out of
    /// the window, which is how zlib decides a stored block is no longer possible.
    fn flush_block(&mut self, last: bool) {
        let stored_len = (self.strstart as isize - self.block_start) as usize;
        let buf = if self.block_start >= 0 {
            let start = self.block_start as usize;
            Some(&self.window[start..start + stored_len])
        } else {
            None
        };
        self.trees.flush_block(&mut self.out, buf, stored_len, last);
        self.block_start = self.strstart as isize;
    }

    /// Ported from zlib's `deflate_fast`, used at levels 1 to 3: take the first match the chain
    /// offers, never look for a better one starting a byte later.
    fn deflate_fast(&mut self, input: &[u8]) {
        // In locals for the same reason as in `deflate_slow`: on data with no matches this loop is
        // one literal per turn, and every store through `self` is one the compiler has to keep.
        let hash = self.hash;
        let max_lazy = self.config.max_lazy;
        let mut strstart = self.strstart;
        let mut lookahead = self.lookahead;
        let mut ins_h = self.ins_h;
        let mut match_length = self.match_length;
        loop {
            if lookahead < MIN_LOOKAHEAD {
                (self.strstart, self.lookahead, self.ins_h) = (strstart, lookahead, ins_h);
                self.match_length = match_length;
                self.fill_window(input);
                (strstart, lookahead, ins_h) = (self.strstart, self.lookahead, self.ins_h);
                if lookahead == 0 {
                    break;
                }
            }

            let mut hash_head = NIL;
            if lookahead >= MIN_MATCH {
                hash_head = insert_at(
                    hash,
                    &self.window,
                    &mut self.head,
                    &mut self.prev,
                    &mut ins_h,
                    strstart,
                );
            }

            if hash_head != NIL && strstart - hash_head as usize <= MAX_DIST {
                self.strstart = strstart;
                self.lookahead = lookahead;
                self.match_length = match_length;
                match_length = self.longest_match(hash_head as usize);
            }

            let bflush;
            if match_length >= MIN_MATCH {
                bflush = self
                    .trees
                    .tally_dist(strstart - self.match_start, match_length - MIN_MATCH);
                lookahead -= match_length;

                // Insert every position the match covers, unless the match is long enough that
                // zlib judges the insertions not worth it.
                if match_length <= max_lazy && lookahead >= MIN_MATCH {
                    match_length -= 1;
                    loop {
                        strstart += 1;
                        insert_at(
                            hash,
                            &self.window,
                            &mut self.head,
                            &mut self.prev,
                            &mut ins_h,
                            strstart,
                        );
                        match_length -= 1;
                        if match_length == 0 {
                            break;
                        }
                    }
                    strstart += 1;
                } else {
                    strstart += match_length;
                    match_length = 0;
                    // The same priming, and Intel skips it on the same condition.
                    if hash == Hash::Multiplicative {
                        self.ins_h = self.window[strstart] as usize;
                        self.update_hash(self.window[strstart + 1]);
                        ins_h = self.ins_h;
                    }
                }
            } else {
                bflush = self.trees.tally_lit(self.window[strstart]);
                lookahead -= 1;
                strstart += 1;
            }
            if bflush {
                (self.strstart, self.lookahead, self.ins_h) = (strstart, lookahead, ins_h);
                self.match_length = match_length;
                self.flush_block(false);
            }
        }
        self.insert = self.strstart.min(MIN_MATCH - 1);
        self.flush_block(true);
    }

    /// Ported from zlib 1.2.13's `deflate_medium.c` in Intel's fork (zlib licence), the block
    /// function that fork uses at levels 4 to 6. **htsjdk's BGZF default is level 5, so this is the function that decides the
    /// bytes of every BAM GATK writes** without `--use-jdk-deflater`.
    ///
    /// It is a different idea from `deflate_slow`, not a tuning of it. `deflate_slow` looks one
    /// byte ahead and keeps whichever of the two matches is longer. This looks one *match* ahead:
    /// it finds the match at the current position, then the match at the position that one would
    /// end at, and emits the current one regardless, carrying the second forward so it is never
    /// searched twice. The lookahead is therefore about avoiding work, not about choosing better.
    ///
    /// [`Self::fizzle_matches`] runs between the two searches and can rewrite both.
    fn deflate_medium(&mut self, input: &[u8]) {
        // `current` is always assigned before it is read; the initialiser only satisfies the
        // compiler, as `memset` does in the C.
        let mut current;
        let mut next = Match::default();

        loop {
            let mut hash_head: u16 = 0;
            if self.lookahead < MIN_LOOKAHEAD {
                self.fill_window(input);
                if self.lookahead == 0 {
                    break;
                }
                next.match_length = 0;
            }

            // Every search starts from a best-so-far of 2, unlike deflate_slow where the previous
            // match's length carries in.
            self.prev_length = 2;

            if next.match_length > 0 {
                current = next;
                next.match_length = 0;
            } else {
                if self.lookahead >= MIN_MATCH {
                    hash_head = self.insert_string(self.strstart);
                }
                if hash_head != 0 && hash_head as usize == self.strstart {
                    hash_head -= 1;
                }
                current = Match {
                    match_start: 0,
                    match_length: 1,
                    strstart: self.strstart,
                    orgstart: self.strstart,
                };
                if hash_head != 0 && self.strstart - hash_head as usize <= MAX_DIST {
                    current.match_length = self.longest_match(hash_head as usize);
                    current.match_start = self.match_start;
                    if current.match_length < MIN_MATCH {
                        current.match_length = 1;
                    }
                    if current.match_start >= current.strstart {
                        current.match_length = 1;
                    }
                }
            }

            self.insert_match(current);

            // Look one match ahead, and only if there is room. The search leaves `strstart`
            // where it found it, so the emit below is unaffected.
            if self.lookahead - current.match_length > MIN_LOOKAHEAD {
                self.strstart = current.strstart + current.match_length;
                hash_head = self.insert_string(self.strstart);
                if hash_head != 0 && hash_head as usize == self.strstart {
                    hash_head -= 1;
                }
                next = Match {
                    match_start: 0,
                    match_length: 1,
                    strstart: self.strstart,
                    orgstart: self.strstart,
                };
                if hash_head != 0 && self.strstart - hash_head as usize <= MAX_DIST {
                    next.match_length = self.longest_match(hash_head as usize);
                    next.match_start = self.match_start;
                    if next.match_start >= next.strstart {
                        next.match_length = 1;
                    }
                    if next.match_length < MIN_MATCH {
                        next.match_length = 1;
                    } else {
                        self.fizzle_matches(&mut current, &mut next);
                    }
                }
                // A three-byte match from very far away is dropped, on a threshold of its own
                // that has nothing to do with `deflate_slow`'s TOO_FAR of 4096.
                if next.match_length == 3 && next.strstart - next.match_start > 12000 {
                    next.match_length = 1;
                }
                self.strstart = current.strstart;
            } else {
                next.match_length = 0;
            }

            let bflush = self.emit_match(current);
            self.strstart += current.match_length;
            if bflush {
                self.flush_block(false);
            }
        }
        self.insert = self.strstart.min(MIN_MATCH - 1);
        self.flush_block(true);
    }

    /// Ported from zlib 1.2.13's `deflate_medium.c`, function `fizzle_matches`: slide the *next* match backwards, one
    /// byte at a time, for as long as the byte before it still matches, shortening the current one
    /// to pay for it.
    ///
    /// It only commits when the current match has been shortened to nothing, so the trade it is
    /// looking for is "two matches" becoming "one longer match". Anything less is abandoned, which
    /// is why both sides are walked on copies and written back only at the end.
    ///
    /// Worth being explicit about, because it is easy to read this function as dead code: the two
    /// assignments that commit the result are `*current = c;` and `*next = n;` on the last branch,
    /// and every other path leaves both arguments alone. Reading it as a no-op costs six differing
    /// bytes in 18 kilobytes, which is what it cost here before the symbol-stream trace found it.
    fn fizzle_matches(&mut self, current: &mut Match, next: &mut Match) {
        if current.match_length <= 1
            || current.match_length > 1 + next.match_start
            || current.match_length > 1 + next.strstart
        {
            return;
        }
        // The cheap rejection: the byte each match would have to grow into, on both sides.
        let back = current.match_length - 1;
        if self.window[next.match_start - back] != self.window[next.strstart - back] {
            return;
        }
        // Overlapping matches are given up on rather than handled.
        if next.match_start + next.match_length >= current.strstart {
            return;
        }

        let mut c = *current;
        let mut n = *next;
        let limit = next.strstart.saturating_sub(MAX_DIST);
        let mut match_at = n.match_start as isize - 1;
        let mut orig_at = n.strstart as isize - 1;
        let mut changed = 0;

        while match_at >= 0
            && orig_at >= 0
            && self.window[match_at as usize] == self.window[orig_at as usize]
        {
            if c.match_length < 1
                || n.strstart <= limit
                || n.match_length >= 256
                || n.match_start == 0
            {
                break;
            }
            n.strstart -= 1;
            n.match_start -= 1;
            n.match_length += 1;
            c.match_length -= 1;
            match_at -= 1;
            orig_at -= 1;
            changed += 1;
            if match_at < 0 || orig_at < 0 {
                break;
            }
        }

        if changed == 0 {
            return;
        }
        // Committed only when the current match has been consumed entirely. `orgstart` moves with
        // it so the chain insertion that follows does not redo the span this just claimed.
        if c.match_length <= 1 && n.match_length != 2 {
            n.orgstart += 1;
            *current = c;
            *next = n;
        }
    }

    /// Ported from zlib 1.2.13's `deflate_medium.c`, function `emit_match`. A "match" shorter than MIN_MATCH is a run of
    /// literals rather than a match, which is how the function represents having found nothing.
    fn emit_match(&mut self, mut m: Match) -> bool {
        let mut flush = false;
        if m.match_length < MIN_MATCH {
            while m.match_length != 0 {
                flush |= self.trees.tally_lit(self.window[m.strstart]);
                self.lookahead -= 1;
                m.strstart += 1;
                m.match_length -= 1;
            }
            return flush;
        }
        flush |= self
            .trees
            .tally_dist(m.strstart - m.match_start, m.match_length - MIN_MATCH);
        self.lookahead -= m.match_length;
        flush
    }

    /// Ported from zlib 1.2.13's `deflate_medium.c`, function `insert_match`: hash the positions a match covers, so the
    /// chains stay complete even though the match finder skipped over them.
    ///
    /// The `strstart >= orgstart` guard is what keeps the lookahead honest. A match carried over
    /// from the previous iteration has already had part of its span inserted, and `orgstart`
    /// remembers where that stopped.
    fn insert_match(&mut self, mut m: Match) {
        if self.lookahead <= m.match_length + MIN_MATCH {
            return;
        }
        if m.match_length < MIN_MATCH {
            while m.match_length != 0 {
                m.strstart += 1;
                m.match_length -= 1;
                if m.match_length != 0 && m.strstart >= m.orgstart {
                    self.insert_string(m.strstart);
                }
            }
            return;
        }
        // Sixteen times `max_insert_length`, where `deflate_fast` uses one: this function is
        // willing to pay much more for complete chains.
        if m.match_length <= 16 * self.config.max_lazy && self.lookahead >= MIN_MATCH {
            m.match_length -= 1;
            loop {
                m.strstart += 1;
                if m.strstart >= m.orgstart {
                    self.insert_string(m.strstart);
                }
                m.match_length -= 1;
                if m.match_length == 0 {
                    break;
                }
            }
        } else {
            m.strstart += m.match_length;
            self.ins_h = self.window[m.strstart] as usize;
            if m.strstart >= 1 {
                self.insert_string(m.strstart - 1);
            }
        }
    }

    /// Ported from zlib's `deflate_slow`, used at levels 4 to 9: hold each match for one byte to
    /// see whether the next position offers a longer one.
    ///
    /// The position and the hash live in locals and go back to `self` only around the calls that
    /// read them there: a store through `head` or `prev` would otherwise make the compiler reload
    /// every field of `self` at every position.
    fn deflate_slow(&mut self, input: &[u8]) {
        let hash = self.hash;
        let max_lazy = self.config.max_lazy;
        let mut strstart = self.strstart;
        let mut lookahead = self.lookahead;
        let mut ins_h = self.ins_h;
        let mut match_length = self.match_length;
        let mut match_start = self.match_start;
        let mut match_available = self.match_available;
        let mut prev_length;
        let mut prev_match;
        loop {
            if lookahead < MIN_LOOKAHEAD {
                (self.strstart, self.lookahead, self.ins_h) = (strstart, lookahead, ins_h);
                self.match_start = match_start;
                self.fill_window(input);
                (strstart, lookahead, ins_h) = (self.strstart, self.lookahead, self.ins_h);
                match_start = self.match_start;
                if lookahead == 0 {
                    break;
                }
            }

            let mut hash_head = NIL;
            if lookahead >= MIN_MATCH {
                hash_head = insert_at(
                    hash,
                    &self.window,
                    &mut self.head,
                    &mut self.prev,
                    &mut ins_h,
                    strstart,
                );
            }

            prev_length = match_length;
            prev_match = match_start;
            match_length = MIN_MATCH - 1;

            if hash_head != NIL
                && prev_length < max_lazy
                && strstart - hash_head as usize <= MAX_DIST
            {
                self.prev_length = prev_length;
                self.strstart = strstart;
                self.lookahead = lookahead;
                self.match_start = match_start;
                match_length = self.longest_match(hash_head as usize);
                match_start = self.match_start;
                // A three-byte match from far away costs more to encode than three literals.
                if match_length <= 5
                    && match_length == MIN_MATCH
                    && strstart - match_start > TOO_FAR
                {
                    match_length = MIN_MATCH - 1;
                }
            }

            if prev_length >= MIN_MATCH && match_length <= prev_length {
                let max_insert = strstart + lookahead - MIN_MATCH;
                let bflush = self
                    .trees
                    .tally_dist(strstart - 1 - prev_match, prev_length - MIN_MATCH);

                lookahead -= prev_length - 1;
                prev_length -= 2;
                loop {
                    strstart += 1;
                    if strstart <= max_insert {
                        insert_at(
                            hash,
                            &self.window,
                            &mut self.head,
                            &mut self.prev,
                            &mut ins_h,
                            strstart,
                        );
                    }
                    prev_length -= 1;
                    if prev_length == 0 {
                        break;
                    }
                }
                match_available = false;
                match_length = MIN_MATCH - 1;
                strstart += 1;
                if bflush {
                    self.strstart = strstart;
                    self.flush_block(false);
                }
            } else if match_available {
                // The previous byte was held back and has now lost: emit it as a literal.
                let bflush = self
                    .trees
                    .tally_lit(self.window[(strstart - 1) & WINDOW_INDEX]);
                if bflush {
                    self.strstart = strstart;
                    self.flush_block(false);
                }
                strstart += 1;
                lookahead -= 1;
            } else {
                match_available = true;
                strstart += 1;
                lookahead -= 1;
            }
        }
        (self.strstart, self.lookahead, self.ins_h) = (strstart, lookahead, ins_h);
        (self.match_length, self.match_start) = (match_length, match_start);
        if match_available {
            self.trees.tally_lit(self.window[self.strstart - 1]);
            self.match_available = false;
        }
        self.insert = self.strstart.min(MIN_MATCH - 1);
        self.flush_block(true);
    }
}
