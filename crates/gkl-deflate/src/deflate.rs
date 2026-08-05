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
    trees: Trees,
    out: Vec<u8>,
}

impl<'a> Deflater<'a> {
    pub fn new(input: &'a [u8], level: usize) -> Self {
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
            trees: Trees::new(LIT_BUFSIZE),
            out: Vec::with_capacity(input.len() / 2 + 64),
        }
    }

    pub fn finish(mut self) -> Vec<u8> {
        if self.config.slow {
            self.deflate_slow();
        } else {
            self.deflate_fast();
        }
        self.out
    }

    #[inline]
    fn update_hash(&mut self, c: u8) {
        self.ins_h = ((self.ins_h << HASH_SHIFT) ^ c as usize) & HASH_MASK;
    }

    /// Ported from zlib's `INSERT_STRING`. Returns the previous head of this hash's chain, which
    /// is where a match search starts.
    #[inline]
    fn insert_string(&mut self, str_pos: usize) -> u16 {
        self.update_hash(self.window[str_pos + MIN_MATCH - 1]);
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
                self.ins_h = self.window[str_pos] as usize;
                self.update_hash(self.window[str_pos + 1]);
                while self.insert != 0 {
                    self.update_hash(self.window[str_pos + MIN_MATCH - 1]);
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
                    self.ins_h = self.window[self.strstart] as usize;
                    self.update_hash(self.window[self.strstart + 1]);
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
