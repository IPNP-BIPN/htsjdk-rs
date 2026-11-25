//! Inflating a block fast, without accepting anything zlib refuses.
//!
//! Inflate has one correct answer for a stream zlib accepts, so speed is free to come from anywhere
//! as long as acceptance does not move: the reader must refuse what htsjdk refuses. libdeflate does
//! not qualify, and a CRC does not rescue it. Its tables read the invalid static codes 286 and 287
//! as a length of 258, so one flipped bit can leave a stream whose payload and CRC are both right
//! and which zlib refuses, as `a_damaged_block_answers_as_zlib_alone_did` found.
//!
//! So this inflater accepts a strict subset of what zlib accepts. Everything zlib would refuse is
//! refused here, and so is everything merely unusual (a code that does not fill its codespace, a
//! stream shorter or longer than declared). A refusal here is not an answer: the caller runs zlib
//! from the start of the block, and zlib's verdict and message are the ones anyone sees. Decision
//! 0049 lists zlib's refusals against the line here that covers each.

/// Bits of a code looked up directly; longer codes go through a subtable.
const LIT_ROOT: u32 = 11;
const DIST_ROOT: u32 = 8;
/// A complete code has at most `symbols / 2` prefixes that need a subtable, each of at most
/// `2^(15 - root)` entries, on top of the root.
const LIT_TABLE: usize = (1 << LIT_ROOT) + 143 * (1 << (15 - LIT_ROOT));
const DIST_TABLE: usize = (1 << DIST_ROOT) + 15 * (1 << (15 - DIST_ROOT));

/// A table entry, laid out as libdeflate's: the low byte is every bit the entry consumes, its
/// code's and its extra bits' together, so one shift takes both; bits 8 to 12 are the code's own
/// length (the extra bits sit above it) or a subtable's size; then three flags; the value is the
/// high half, and a literal sets the top bit so that one test tells it apart.
const LITERAL: u32 = 1 << 31;
/// Anything the fast loop does not decode inline: a subtable, the end of the block, or a symbol
/// zlib calls invalid, which carries only this flag.
const EXCEPTIONAL: u32 = 1 << 15;
const SUBTABLE: u32 = 1 << 14;
const END: u32 = 1 << 13;

const fn entry(code_bits: u32, extra: u32, value: u32) -> u32 {
    (code_bits + extra) | (code_bits << 8) | (value << 16)
}

/// The value an entry yields: its base plus the extra bits above its code, from `buf` as it was
/// before the entry was consumed.
#[inline(always)]
fn value_of(e: u32, buf: u64) -> usize {
    (e >> 16) as usize + ((buf & ((1u64 << (e & 0xff)) - 1)) >> ((e >> 8) & 0x1f)) as usize
}

const LENGTH_BASE: [u16; 29] = [
    3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131,
    163, 195, 227, 258,
];
const LENGTH_EXTRA: [u8; 29] = [
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0,
];
const DIST_BASE: [u16; 30] = [
    1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537,
    2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577,
];
const DIST_EXTRA: [u8; 30] = [
    0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13,
    13,
];
const CODE_ORDER: [u8; 19] = [
    16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15,
];

fn literal_entry(symbol: usize, bits: u32) -> u32 {
    match symbol {
        0..=255 => LITERAL | entry(bits, 0, symbol as u32),
        256 => END | EXCEPTIONAL | bits,
        257..=285 => entry(
            bits,
            LENGTH_EXTRA[symbol - 257] as u32,
            LENGTH_BASE[symbol - 257] as u32,
        ),
        // Reachable only through the static code, where zlib calls them invalid.
        _ => EXCEPTIONAL | bits,
    }
}

fn distance_entry(symbol: usize, bits: u32) -> u32 {
    match symbol {
        0..=29 => entry(bits, DIST_EXTRA[symbol] as u32, DIST_BASE[symbol] as u32),
        _ => EXCEPTIONAL | bits,
    }
}

/// Builds a two-level decoding table from code lengths; `false` unless they exactly fill the
/// codespace, which is stricter than zlib's rule and also covers its over-subscribed refusals.
///
/// libdeflate's `build_decode_table`: the symbols sorted by length, so codes come in canonical
/// order; the root filled for the shortest length and doubled for each longer one; a subtable
/// per root prefix, sized for the codes that share it.
#[inline(always)]
fn build(table: &mut [u32], root: u32, lens: &[u8], make: impl Fn(usize, u32) -> u32) -> bool {
    let mut count = [0u16; 16];
    for &len in lens {
        count[len as usize] += 1;
    }
    let mut offsets = [0u16; 16];
    let mut left = 1i32;
    for len in 1..16 {
        left = (left << 1) - count[len] as i32;
        if left < 0 {
            return false;
        }
        if len < 15 {
            offsets[len + 1] = offsets[len] + count[len];
        }
    }
    if left != 0 {
        return false;
    }
    let mut sorted = [0u16; 320];
    for (symbol, &len) in lens.iter().enumerate() {
        if len != 0 {
            sorted[offsets[len as usize] as usize] = symbol as u16;
            offsets[len as usize] += 1;
        }
    }
    let mut sorted = sorted.iter().map(|&s| s as usize);

    // Codewords are bit-reversed, so the next one is found from the highest bit that is clear.
    let mut codeword = 0usize;
    let mut len = 1u32;
    while count[len as usize] == 0 {
        len += 1;
    }
    let mut remaining = count[len as usize];
    let mut end = 1usize << len;
    while len <= root {
        loop {
            table[codeword] = make(sorted.next().unwrap_or(0), len);
            if codeword == end - 1 {
                // The last codeword: the root is full once the table is doubled up to its size.
                while len < root {
                    table.copy_within(..end, end);
                    end <<= 1;
                    len += 1;
                }
                return true;
            }
            let bit = 1 << (usize::BITS - 1 - (codeword ^ (end - 1)).leading_zeros());
            codeword = (codeword & (bit - 1)) | bit;
            remaining -= 1;
            if remaining == 0 {
                break;
            }
        }
        loop {
            len += 1;
            if len <= root {
                table.copy_within(..end, end);
                end <<= 1;
            }
            remaining = count[len as usize];
            if remaining != 0 {
                break;
            }
        }
    }

    let root_mask = (1usize << root) - 1;
    let mut end = 1usize << root;
    let mut prefix = usize::MAX;
    let mut start = 0usize;
    loop {
        if codeword & root_mask != prefix {
            prefix = codeword & root_mask;
            start = end;
            let mut sub_bits = len - root;
            let mut used = remaining as u32;
            while used < 1 << sub_bits {
                sub_bits += 1;
                used = (used << 1) + count[(root + sub_bits) as usize] as u32;
            }
            end = start + (1 << sub_bits);
            if end > table.len() {
                return false;
            }
            table[prefix] =
                EXCEPTIONAL | SUBTABLE | root | (sub_bits << 8) | ((start as u32) << 16);
        }
        let value = make(sorted.next().unwrap_or(0), len - root);
        let mut at = start + (codeword >> root);
        let stride = 1 << (len - root);
        while at < end {
            table[at] = value;
            at += stride;
        }
        if codeword == (1 << len) - 1 {
            return true;
        }
        let bit = 1 << (usize::BITS - 1 - (codeword ^ ((1 << len) - 1)).leading_zeros());
        codeword = (codeword & (bit - 1)) | bit;
        remaining -= 1;
        while remaining == 0 {
            len += 1;
            remaining = count[len as usize];
        }
    }
}

/// One thread's or one reader's decoding tables.
pub(crate) struct Inflater {
    literal: Box<[u32]>,
    distance: Box<[u32]>,
    /// Whether the tables already hold the static code, which many small blocks use in a row.
    static_tables: bool,
    /// Whether to run the copy compiled for BMI2, whose shifts by a register take one micro-op
    /// instead of three; libdeflate picks the same variant at run time. Every copy gives the same
    /// answer, only faster or slower.
    #[cfg(target_arch = "x86_64")]
    bmi2: bool,
}

struct Bits<'a> {
    input: &'a [u8],
    at: usize,
    buf: u64,
    /// How many bits of `buf` are real; the ones above them, if any, are the next input bytes.
    left: u32,
}

impl Bits<'_> {
    /// At least 48 bits when eight input bytes remain, otherwise as many as remain.
    #[inline(always)]
    fn refill(&mut self) {
        if self.at + 8 <= self.input.len() {
            // SAFETY: eight bytes are in bounds, just checked.
            let word = unsafe {
                self.input
                    .as_ptr()
                    .add(self.at)
                    .cast::<u64>()
                    .read_unaligned()
            };
            self.buf |= u64::from_le(word) << self.left;
            self.at += (63 - self.left as usize) >> 3;
            self.left |= 56;
        } else {
            while self.left <= 56 && self.at < self.input.len() {
                self.buf |= (self.input[self.at] as u64) << self.left;
                self.at += 1;
                self.left += 8;
            }
        }
    }

    #[inline(always)]
    fn take(&mut self, n: u32) -> Option<u32> {
        self.refill();
        if self.left < n {
            return None;
        }
        let value = (self.buf & ((1u64 << n) - 1)) as u32;
        self.buf >>= n;
        self.left -= n;
        Some(value)
    }
}

impl Inflater {
    pub(crate) fn new() -> Self {
        Inflater {
            literal: vec![0; LIT_TABLE].into_boxed_slice(),
            distance: vec![0; DIST_TABLE].into_boxed_slice(),
            static_tables: false,
            #[cfg(target_arch = "x86_64")]
            bmi2: std::arch::is_x86_feature_detected!("bmi1")
                && std::arch::is_x86_feature_detected!("bmi2"),
        }
    }

    /// One that never takes a CPU-specific path, so that tests on any machine cover the baseline.
    #[cfg(test)]
    pub(crate) fn baseline() -> Self {
        Inflater {
            #[cfg(target_arch = "x86_64")]
            bmi2: false,
            ..Self::new()
        }
    }

    /// Whether `out` now holds all of `deflated`'s inflation, which is exactly `out.len()` bytes.
    /// On `false`, `out` holds nothing anyone may read, and the caller asks zlib.
    pub(crate) fn inflate(&mut self, deflated: &[u8], out: &mut [u8]) -> bool {
        // SAFETY: the slice is `out.len()` bytes of writable memory.
        unsafe { self.inflate_raw(deflated, out.as_mut_ptr(), out.len()) }
    }

    /// [`Self::inflate`] into a vector's spare capacity, whose size is the declared size.
    pub(crate) fn inflate_vec(&mut self, deflated: &[u8], out: &mut Vec<u8>) -> bool {
        debug_assert!(out.is_empty());
        let size = out.capacity();
        // SAFETY: the capacity is `size` bytes of writable memory, and `set_len` follows only an
        // answer that wrote every one of them.
        unsafe {
            let done = self.inflate_raw(deflated, out.as_mut_ptr(), size);
            if done {
                out.set_len(size);
            }
            done
        }
    }

    /// # Safety
    /// `out` must be valid for writes of `size` bytes. Bytes are read back only after being
    /// written, so they need not be initialized.
    unsafe fn inflate_raw(&mut self, deflated: &[u8], out: *mut u8, size: usize) -> bool {
        #[cfg(target_arch = "x86_64")]
        if self.bmi2 {
            return self.inflate_bmi2(deflated, out, size);
        }
        self.inflate_any(deflated, out, size)
    }

    /// # Safety
    /// As for [`Self::inflate_raw`], on a CPU with BMI1 and BMI2.
    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "bmi1,bmi2")]
    unsafe fn inflate_bmi2(&mut self, deflated: &[u8], out: *mut u8, size: usize) -> bool {
        self.inflate_any(deflated, out, size)
    }

    /// # Safety
    /// As for [`Self::inflate_raw`].
    #[inline(always)]
    unsafe fn inflate_any(&mut self, deflated: &[u8], out: *mut u8, size: usize) -> bool {
        let mut bits = Bits {
            input: deflated,
            at: 0,
            buf: 0,
            left: 0,
        };
        let mut pos = 0usize;
        loop {
            let Some(header) = bits.take(3) else {
                return false;
            };
            let last = header & 1 != 0;
            let fine = match header >> 1 {
                0 => stored(&mut bits, out, size, &mut pos),
                1 => {
                    self.use_static_tables();
                    self.codes(&mut bits, out, size, &mut pos)
                }
                2 => {
                    self.read_dynamic_tables(&mut bits)
                        && self.codes(&mut bits, out, size, &mut pos)
                }
                _ => false,
            };
            if !fine {
                return false;
            }
            if last {
                return pos == size;
            }
        }
    }

    #[inline(always)]
    fn use_static_tables(&mut self) {
        if self.static_tables {
            return;
        }
        let mut lens = [8u8; 288];
        lens[144..256].fill(9);
        lens[256..280].fill(7);
        build(&mut self.literal, LIT_ROOT, &lens, literal_entry);
        build(&mut self.distance, DIST_ROOT, &[5u8; 32], distance_entry);
        self.static_tables = true;
    }

    /// The header of a dynamic block. Each `false` below is one of zlib's refusals or stricter.
    #[inline(always)]
    fn read_dynamic_tables(&mut self, bits: &mut Bits<'_>) -> bool {
        self.static_tables = false;
        let Some(counts) = bits.take(14) else {
            return false;
        };
        let literals = (counts & 31) as usize + 257;
        let distances = ((counts >> 5) & 31) as usize + 1;
        let code_lengths = ((counts >> 10) & 15) as usize + 4;
        // "too many length or distance symbols".
        if literals > 286 || distances > 30 {
            return false;
        }
        let mut order = [0u8; 19];
        for &at in &CODE_ORDER[..code_lengths] {
            let Some(len) = bits.take(3) else {
                return false;
            };
            order[at as usize] = len as u8;
        }
        // "invalid code lengths set": a precode must fill its codespace.
        let mut precode = [0u32; 128];
        if !build(&mut precode, 7, &order, |symbol, bits| {
            entry(bits, 0, symbol as u32)
        }) {
            return false;
        }

        let total = literals + distances;
        let mut lens = [0u8; 320];
        let mut i = 0;
        while i < total {
            bits.refill();
            let e = precode[(bits.buf & 127) as usize];
            let used = e & 0xff;
            if used > bits.left {
                return false;
            }
            bits.buf >>= used;
            bits.left -= used;
            let symbol = e >> 16;
            if symbol < 16 {
                lens[i] = symbol as u8;
                i += 1;
                continue;
            }
            let (value, extra, base) = match symbol {
                // "invalid bit length repeat": nothing to repeat.
                16 if i == 0 => return false,
                16 => (lens[i - 1], 2, 3),
                17 => (0, 3, 3),
                _ => (0, 7, 11),
            };
            let Some(more) = bits.take(extra) else {
                return false;
            };
            let run = base + more as usize;
            // "invalid bit length repeat": past the end.
            if i + run > total {
                return false;
            }
            lens[i..i + run].fill(value);
            i += run;
        }
        // "invalid code -- missing end-of-block".
        if lens[256] == 0 {
            return false;
        }
        // "invalid literal/lengths set" and "invalid distances set".
        build(
            &mut self.literal,
            LIT_ROOT,
            &lens[..literals],
            literal_entry,
        ) && build(
            &mut self.distance,
            DIST_ROOT,
            &lens[literals..total],
            distance_entry,
        )
    }

    /// The symbols of one block, up to its end-of-block code.
    ///
    /// # Safety
    /// `out` must be valid for writes of `size` bytes, and `*pos` at most `size`, with everything
    /// below it written.
    #[inline(always)]
    unsafe fn codes(
        &self,
        bits: &mut Bits<'_>,
        out: *mut u8,
        size: usize,
        pos: &mut usize,
    ) -> bool {
        // SAFETY of the table reads: an index is masked to a root's size, or is a subtable's offset
        // plus a value masked to that subtable's size, and `build` placed each inside its table.
        let literal = self.literal.as_ptr();
        let distance = self.distance.as_ptr();
        const LIT_MASK: u64 = (1 << LIT_ROOT) - 1;
        const DIST_MASK: u64 = (1 << DIST_ROOT) - 1;
        let mut at = *pos;
        loop {
            // The fast loop: input for every refill, and room for a match plus the copy's overrun,
            // so nothing inside it needs checking but the symbols and the distance.
            //
            // libdeflate's: `left` drops whole entries, so only its low byte counts; a refill
            // leaves at least 56 bits counted and 64 real, so a match (48 bits) needs one refill,
            // and the next entry can be looked up before the refill that follows, so that the
            // copy overlaps both.
            if bits.at + 16 <= bits.input.len() && at + 258 + 40 <= size {
                let input = bits.input.as_ptr();
                let in_end = bits.input.len() - 16;
                let mut inp = bits.at;
                let mut buf = bits.buf;
                let mut left = bits.left;
                macro_rules! refill {
                    () => {
                        buf |= u64::from_le(input.add(inp).cast::<u64>().read_unaligned())
                            << (left as u8);
                        inp += 7 - ((left >> 3) & 7) as usize;
                        left |= 56;
                    };
                }
                macro_rules! consume {
                    ($e:expr) => {
                        buf >>= $e as u8;
                        left = left.wrapping_sub($e);
                    };
                }
                macro_rules! leave {
                    () => {
                        bits.at = inp;
                        bits.buf = buf;
                        bits.left = left & 0xff;
                    };
                }
                refill!();
                let mut e = *literal.add((buf & LIT_MASK) as usize);
                loop {
                    if e & LITERAL != 0 {
                        // Up to three literals of at most LIT_ROOT bits each fit in one refill.
                        consume!(e);
                        *out.add(at) = (e >> 16) as u8;
                        at += 1;
                        e = *literal.add((buf & LIT_MASK) as usize);
                        if e & LITERAL != 0 {
                            consume!(e);
                            *out.add(at) = (e >> 16) as u8;
                            at += 1;
                            e = *literal.add((buf & LIT_MASK) as usize);
                            if e & LITERAL != 0 {
                                consume!(e);
                                *out.add(at) = (e >> 16) as u8;
                                at += 1;
                                e = *literal.add((buf & LIT_MASK) as usize);
                                refill!();
                                if inp > in_end || at + 258 + 40 > size {
                                    break;
                                }
                                continue;
                            }
                        }
                        // Fewer than 48 bits may be left, and a match takes up to 48.
                        refill!();
                    }
                    if e & EXCEPTIONAL != 0 {
                        if e & SUBTABLE != 0 {
                            consume!(e);
                            let sub = (buf & ((1u64 << ((e >> 8) & 0x1f)) - 1)) as usize;
                            e = *literal.add((e >> 16) as usize + sub);
                            if e & LITERAL != 0 {
                                consume!(e);
                                *out.add(at) = (e >> 16) as u8;
                                at += 1;
                                e = *literal.add((buf & LIT_MASK) as usize);
                                refill!();
                                if inp > in_end || at + 258 + 40 > size {
                                    break;
                                }
                                continue;
                            }
                        }
                        if e & EXCEPTIONAL != 0 {
                            if e & END != 0 {
                                consume!(e);
                                leave!();
                                *pos = at;
                                return true;
                            }
                            // "invalid literal/length code".
                            return false;
                        }
                    }
                    let length = value_of(e, buf);
                    consume!(e);

                    let mut d = *distance.add((buf & DIST_MASK) as usize);
                    if d & EXCEPTIONAL != 0 {
                        if d & SUBTABLE == 0 {
                            // "invalid distance code".
                            return false;
                        }
                        consume!(d);
                        let sub = (buf & ((1u64 << ((d >> 8) & 0x1f)) - 1)) as usize;
                        d = *distance.add((d >> 16) as usize + sub);
                        if d & EXCEPTIONAL != 0 {
                            return false;
                        }
                    }
                    let dist = value_of(d, buf);
                    consume!(d);
                    // "invalid distance too far back".
                    if dist > at {
                        return false;
                    }
                    e = *literal.add((buf & LIT_MASK) as usize);
                    refill!();
                    copy_match(out.add(at), dist, length);
                    at += length;
                    if inp > in_end || at + 258 + 40 > size {
                        break;
                    }
                }
                leave!();
            }

            // The tail: every step checked against what is really left of the input and the output.
            bits.refill();
            let mut e = *literal.add((bits.buf & LIT_MASK) as usize);
            if e & (EXCEPTIONAL | SUBTABLE) == EXCEPTIONAL | SUBTABLE {
                let root = e & 0xff;
                if root > bits.left {
                    return false;
                }
                bits.buf >>= root;
                bits.left -= root;
                let sub = (bits.buf & ((1u64 << ((e >> 8) & 0x1f)) - 1)) as usize;
                e = *literal.add((e >> 16) as usize + sub);
            }
            if e & 0xff > bits.left {
                return false;
            }
            if e & LITERAL != 0 {
                if at >= size {
                    return false;
                }
                bits.buf >>= e & 0xff;
                bits.left -= e & 0xff;
                *out.add(at) = (e >> 16) as u8;
                at += 1;
                continue;
            }
            if e & EXCEPTIONAL != 0 {
                if e & END != 0 {
                    bits.buf >>= e & 0xff;
                    bits.left -= e & 0xff;
                    *pos = at;
                    return true;
                }
                return false;
            }
            let length = value_of(e, bits.buf);
            bits.buf >>= e & 0xff;
            bits.left -= e & 0xff;
            bits.refill();
            let mut d = *distance.add((bits.buf & DIST_MASK) as usize);
            if d & (EXCEPTIONAL | SUBTABLE) == EXCEPTIONAL | SUBTABLE {
                let root = d & 0xff;
                if root > bits.left {
                    return false;
                }
                bits.buf >>= root;
                bits.left -= root;
                let sub = (bits.buf & ((1u64 << ((d >> 8) & 0x1f)) - 1)) as usize;
                d = *distance.add((d >> 16) as usize + sub);
            }
            if d & EXCEPTIONAL != 0 || d & 0xff > bits.left {
                return false;
            }
            let dist = value_of(d, bits.buf);
            bits.buf >>= d & 0xff;
            bits.left -= d & 0xff;
            if dist > at || at + length > size {
                return false;
            }
            for i in at..at + length {
                *out.add(i) = *out.add(i - dist);
            }
            at += length;
        }
    }
}

/// Copies a match of `length` bytes from `dist` back, writing up to 39 bytes past its end.
///
/// # Safety
/// `dst - dist` to `dst` must be written, and `dst` valid for writes of `length + 39` bytes.
#[inline(always)]
unsafe fn copy_match(dst: *mut u8, dist: usize, length: usize) {
    let end = dst.add(length);
    let mut d = dst;
    if dist >= 8 {
        // Each word's source lies wholly behind its destination, so it is all written: five words
        // cover most matches without a branch.
        let mut s = dst.sub(dist);
        loop {
            for _ in 0..5 {
                d.cast::<u64>()
                    .write_unaligned(s.cast::<u64>().read_unaligned());
                s = s.add(8);
                d = d.add(8);
            }
            if d >= end {
                return;
            }
        }
    }
    // A period shorter than a word: the period's bytes repeated across a word, written at a stride
    // that is a whole number of periods, so every word starts in phase. Only the period is read.
    let src = dst.sub(dist);
    let (word, stride) = if dist == 1 {
        (0x0101_0101_0101_0101 * *src as u64, 8)
    } else {
        let mut word = 0u64;
        let mut j = 0;
        for i in 0..8 {
            word |= (*src.add(j) as u64) << (8 * i);
            j += 1;
            if j == dist {
                j = 0;
            }
        }
        (word, 8 - 8 % dist)
    };
    while d < end {
        d.cast::<u64>().write_unaligned(word);
        d = d.add(stride);
    }
}

/// A stored block: back to the byte boundary, the length and its complement, the bytes.
///
/// # Safety
/// As for [`Inflater::codes`].
#[inline(always)]
unsafe fn stored(bits: &mut Bits<'_>, out: *mut u8, size: usize, pos: &mut usize) -> bool {
    // Whole bytes still in the buffer go back to the input; the rest of this byte is padding.
    bits.at -= (bits.left >> 3) as usize;
    bits.buf = 0;
    bits.left = 0;
    let Some(header) = bits.input.get(bits.at..bits.at + 4) else {
        return false;
    };
    let length = u16::from_le_bytes([header[0], header[1]]);
    // "invalid stored block lengths".
    if length != !u16::from_le_bytes([header[2], header[3]]) {
        return false;
    }
    let length = length as usize;
    bits.at += 4;
    let Some(body) = bits.input.get(bits.at..bits.at + length) else {
        return false;
    };
    if *pos + length > size {
        return false;
    }
    std::ptr::copy_nonoverlapping(body.as_ptr(), out.add(*pos), length);
    *pos += length;
    bits.at += length;
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BgzfWriter, Deflater};
    use flate2::{Decompress, FlushDecompress};
    use std::io::Write;

    struct Lcg(u64);

    impl Lcg {
        fn next(&mut self) -> usize {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (self.0 >> 33) as usize
        }
    }

    /// What zlib alone makes of a block's deflate stream: the bytes, or nothing for a refusal or a
    /// size that disagrees.
    fn zlib(deflated: &[u8], size: usize) -> Option<Vec<u8>> {
        let mut out = Vec::with_capacity(size);
        Decompress::new(false)
            .decompress_vec(deflated, &mut out, FlushDecompress::Finish)
            .ok()?;
        (out.len() == size).then_some(out)
    }

    fn deflate_stream(payload: &[u8], level: u32, deflater: Deflater) -> Vec<u8> {
        let mut writer = BgzfWriter::with_deflater(Vec::new(), level, deflater).parallel(false);
        writer.write_all(payload).unwrap();
        let stream = writer.into_inner_without_terminator().unwrap();
        stream[crate::BLOCK_HEADER_LENGTH..stream.len() - crate::BLOCK_FOOTER_LENGTH].to_vec()
    }

    /// Whatever this inflater accepts, zlib accepts with the same bytes: the one property the
    /// reader rests on. The reverse need not hold, so what is asserted about it is only that
    /// ordinary streams, from either encoder, do not need zlib.
    #[test]
    fn it_accepts_nothing_zlib_refuses_and_answers_as_zlib_does() {
        let mut lcg = Lcg(19);
        let mut payloads: Vec<Vec<u8>> = vec![
            b"A".to_vec(),
            b"ACGTTGCAACGT".repeat(40),
            vec![b'A'; 60_000],
            (0..2_000).map(|i| (i % 7) as u8).collect(),
        ];
        payloads.push((0..60_000).map(|_| b"ACGTN"[lcg.next() % 5]).collect());
        payloads.push((0..60_000).map(|_| (lcg.next() % 41) as u8).collect());
        payloads.push((0..3_000).map(|_| lcg.next() as u8).collect());
        let mut inflaters = [Inflater::new(), Inflater::baseline()];
        let mut trials = 0;
        for payload in &payloads {
            for (level, deflater) in [
                (1, Deflater::Jdk),
                (5, Deflater::Jdk),
                (9, Deflater::Jdk),
                (5, Deflater::Gkl),
                (1, Deflater::Libdeflate),
                (6, Deflater::Libdeflate),
                (12, Deflater::Libdeflate),
            ] {
                if level > 9 {
                    continue;
                }
                let stream = deflate_stream(payload, level, deflater);
                for inflater in &mut inflaters {
                    let mut out = vec![0u8; payload.len()];
                    assert!(inflater.inflate(&stream, &mut out), "{level} {deflater:?}");
                    assert_eq!(out, *payload);
                }

                for _ in 0..3_000 {
                    let mut damaged = stream.clone();
                    let flips = 1 + lcg.next() % 3;
                    for _ in 0..flips {
                        let at = lcg.next() % damaged.len().max(1);
                        if let Some(byte) = damaged.get_mut(at) {
                            *byte ^= 1 << (lcg.next() % 8);
                        }
                    }
                    if lcg.next().is_multiple_of(8) {
                        damaged.truncate(lcg.next() % (damaged.len() + 1));
                    }
                    let size = if lcg.next().is_multiple_of(8) {
                        lcg.next() % (payload.len() + 2)
                    } else {
                        payload.len()
                    };
                    // Every CPU path accepts the same streams, with zlib's bytes.
                    let [fast, baseline] = &mut inflaters;
                    let mut out = vec![0u8; size];
                    let mut out_baseline = vec![0u8; size];
                    let accepted = fast.inflate(&damaged, &mut out);
                    assert_eq!(accepted, baseline.inflate(&damaged, &mut out_baseline));
                    if accepted {
                        assert_eq!(out, out_baseline);
                        assert_eq!(Some(out), zlib(&damaged, size), "level {level}");
                    }
                    trials += 1;
                }
            }
        }
        assert!(trials > 100_000);
    }
}
