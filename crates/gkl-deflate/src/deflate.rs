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
fn crc32c_u32(val: u32) -> u32 {
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

pub struct Deflater<'a> {
    input: &'a [u8],
    next_in: usize,

    window: Vec<u8>,
    prev: Vec<u16>,
    head: Vec<u16>,

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
    out: Vec<u8>,
}

impl<'a> Deflater<'a> {
    pub fn new(input: &'a [u8], level: usize, flavour: Flavour) -> Self {
        let (hash, medium) = match flavour {
            Flavour::Jdk => (Hash::Multiplicative, false),
            Flavour::Gkl { sse42: false } => (Hash::Multiplicative, (4..=6).contains(&level)),
            Flavour::Gkl { sse42: true } => (
                Hash::Crc32c {
                    three_byte: level >= 6,
                },
                (4..=6).contains(&level),
            ),
        };
        Deflater {
            input,
            next_in: 0,
            window: vec![0; WINDOW_SIZE + WINDOW_PADDING],
            prev: vec![0; W_SIZE],
            head: vec![0; HASH_SIZE],
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
            out: Vec::with_capacity(input.len() / 2 + 64),
        }
    }

    pub fn finish(mut self) -> Vec<u8> {
        if self.medium {
            self.deflate_medium();
        } else if self.config.slow {
            self.deflate_slow();
        } else {
            self.deflate_fast();
        }
        self.out
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
                let mut val = u32::from_le_bytes([
                    self.window[str_pos],
                    self.window[str_pos + 1],
                    self.window[str_pos + 2],
                    self.window[str_pos + 3],
                ]);
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
        let head = self.head[self.ins_h];
        self.prev[str_pos & W_MASK] = head;
        self.head[self.ins_h] = str_pos as u16;
        head
    }

    /// Ported from zlib's `slide_hash`. Every recorded position moves down by one window, and
    /// anything that falls off the bottom becomes NIL.
    fn slide_hash(&mut self) {
        for slot in self.head.iter_mut() {
            let m = *slot as usize;
            *slot = if m >= W_SIZE {
                (m - W_SIZE) as u16
            } else {
                NIL
            };
        }
        for slot in self.prev.iter_mut() {
            let m = *slot as usize;
            *slot = if m >= W_SIZE {
                (m - W_SIZE) as u16
            } else {
                NIL
            };
        }
    }

    /// Ported from zlib's `fill_window`. The slide is the part that matters for byte-identity:
    /// it happens at a fixed point, and `insert` carries across it the positions that were read
    /// but never hashed.
    fn fill_window(&mut self) {
        loop {
            let mut more = WINDOW_SIZE - self.lookahead - self.strstart;

            if self.strstart >= W_SIZE + MAX_DIST {
                self.window.copy_within(W_SIZE..WINDOW_SIZE - more, 0);
                self.match_start -= W_SIZE;
                self.strstart -= W_SIZE;
                self.block_start -= W_SIZE as isize;
                if self.insert > self.strstart {
                    self.insert = self.strstart;
                }
                self.slide_hash();
                more += W_SIZE;
            }
            if self.next_in == self.input.len() {
                break;
            }

            let n = more.min(self.input.len() - self.next_in);
            let dst = self.strstart + self.lookahead;
            self.window[dst..dst + n].copy_from_slice(&self.input[self.next_in..self.next_in + n]);
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

            if self.lookahead >= MIN_LOOKAHEAD || self.next_in == self.input.len() {
                break;
            }
        }
    }

    /// Ported from zlib's `longest_match`. The two `UNALIGNED_OK` variants in the C differ only in
    /// how fast they compare; this is the plain one, and it finds the same match.
    fn longest_match(&mut self, mut cur_match: usize) -> usize {
        let mut chain_length = self.config.max_chain;
        let mut best_len = self.prev_length;
        let mut nice_match = self.config.nice_length;
        // zlib's NIL is 0, and a saturating subtraction lands on it for the same reason.
        let limit = self.strstart.saturating_sub(MAX_DIST);
        let scan = self.strstart;
        let strend = self.strstart + MAX_MATCH;

        // Already sitting on a good match: spend a quarter of the budget looking for a better one.
        if self.prev_length >= self.config.good_length {
            chain_length >>= 2;
        }
        // Never look past the end of the input, which is what makes the output deterministic.
        if nice_match > self.lookahead {
            nice_match = self.lookahead;
        }

        loop {
            let m = cur_match;
            // The three cheap rejections zlib makes before comparing anything: the byte one past
            // the current best, the byte at the best, and the first two bytes.
            if self.window[m + best_len] != self.window[scan + best_len]
                || (best_len > 0
                    && self.window[m + best_len - 1] != self.window[scan + best_len - 1])
                || self.window[m] != self.window[scan]
                || self.window[m + 1] != self.window[scan + 1]
            {
                // fall through to the chain step
            } else {
                // scan[2] == match[2] is implied by the hash, so the comparison starts at 3.
                let mut s = scan + 2;
                let mut t = m + 2;
                loop {
                    s += 1;
                    t += 1;
                    if self.window[s] != self.window[t] || s >= strend {
                        break;
                    }
                }
                let len = MAX_MATCH - (strend - s);
                if len > best_len {
                    self.match_start = cur_match;
                    best_len = len;
                    if len >= nice_match {
                        break;
                    }
                }
            }
            cur_match = self.prev[cur_match & W_MASK] as usize;
            chain_length -= 1;
            if cur_match <= limit || chain_length == 0 {
                break;
            }
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
            Some(self.window[start..start + stored_len].to_vec())
        } else {
            None
        };
        self.trees
            .flush_block(&mut self.out, buf.as_deref(), stored_len, last);
        self.block_start = self.strstart as isize;
    }

    /// Ported from zlib's `deflate_fast`, used at levels 1 to 3: take the first match the chain
    /// offers, never look for a better one starting a byte later.
    fn deflate_fast(&mut self) {
        loop {
            if self.lookahead < MIN_LOOKAHEAD {
                self.fill_window();
                if self.lookahead == 0 {
                    break;
                }
            }

            let mut hash_head = NIL;
            if self.lookahead >= MIN_MATCH {
                hash_head = self.insert_string(self.strstart);
            }

            if hash_head != NIL && self.strstart - hash_head as usize <= MAX_DIST {
                self.match_length = self.longest_match(hash_head as usize);
            }

            let bflush;
            if self.match_length >= MIN_MATCH {
                bflush = self.trees.tally_dist(
                    self.strstart - self.match_start,
                    self.match_length - MIN_MATCH,
                );
                self.lookahead -= self.match_length;

                // Insert every position the match covers, unless the match is long enough that
                // zlib judges the insertions not worth it.
                if self.match_length <= self.config.max_lazy && self.lookahead >= MIN_MATCH {
                    self.match_length -= 1;
                    loop {
                        self.strstart += 1;
                        self.insert_string(self.strstart);
                        self.match_length -= 1;
                        if self.match_length == 0 {
                            break;
                        }
                    }
                    self.strstart += 1;
                } else {
                    self.strstart += self.match_length;
                    self.match_length = 0;
                    // The same priming, and Intel skips it on the same condition.
                    if self.hash == Hash::Multiplicative {
                        self.ins_h = self.window[self.strstart] as usize;
                        self.update_hash(self.window[self.strstart + 1]);
                    }
                }
            } else {
                bflush = self.trees.tally_lit(self.window[self.strstart]);
                self.lookahead -= 1;
                self.strstart += 1;
            }
            if bflush {
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
    fn deflate_medium(&mut self) {
        // `current` is always assigned before it is read; the initialiser only satisfies the
        // compiler, as `memset` does in the C.
        let mut current;
        let mut next = Match::default();

        loop {
            let mut hash_head: u16 = 0;
            if self.lookahead < MIN_LOOKAHEAD {
                self.fill_window();
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
    fn deflate_slow(&mut self) {
        loop {
            if self.lookahead < MIN_LOOKAHEAD {
                self.fill_window();
                if self.lookahead == 0 {
                    break;
                }
            }

            let mut hash_head = NIL;
            if self.lookahead >= MIN_MATCH {
                hash_head = self.insert_string(self.strstart);
            }

            self.prev_length = self.match_length;
            self.prev_match = self.match_start;
            self.match_length = MIN_MATCH - 1;

            if hash_head != NIL
                && self.prev_length < self.config.max_lazy
                && self.strstart - hash_head as usize <= MAX_DIST
            {
                self.match_length = self.longest_match(hash_head as usize);
                // A three-byte match from far away costs more to encode than three literals.
                if self.match_length <= 5
                    && self.match_length == MIN_MATCH
                    && self.strstart - self.match_start > TOO_FAR
                {
                    self.match_length = MIN_MATCH - 1;
                }
            }

            if self.prev_length >= MIN_MATCH && self.match_length <= self.prev_length {
                let max_insert = self.strstart + self.lookahead - MIN_MATCH;
                let bflush = self.trees.tally_dist(
                    self.strstart - 1 - self.prev_match,
                    self.prev_length - MIN_MATCH,
                );

                self.lookahead -= self.prev_length - 1;
                self.prev_length -= 2;
                loop {
                    self.strstart += 1;
                    if self.strstart <= max_insert {
                        self.insert_string(self.strstart);
                    }
                    self.prev_length -= 1;
                    if self.prev_length == 0 {
                        break;
                    }
                }
                self.match_available = false;
                self.match_length = MIN_MATCH - 1;
                self.strstart += 1;
                if bflush {
                    self.flush_block(false);
                }
            } else if self.match_available {
                // The previous byte was held back and has now lost: emit it as a literal.
                let bflush = self.trees.tally_lit(self.window[self.strstart - 1]);
                if bflush {
                    self.flush_block(false);
                }
                self.strstart += 1;
                self.lookahead -= 1;
            } else {
                self.match_available = true;
                self.strstart += 1;
                self.lookahead -= 1;
            }
        }
        if self.match_available {
            self.trees.tally_lit(self.window[self.strstart - 1]);
            self.match_available = false;
        }
        self.insert = self.strstart.min(MIN_MATCH - 1);
        self.flush_block(true);
    }
}
