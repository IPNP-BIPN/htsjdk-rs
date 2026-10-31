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
const LIT_ROOT: u32 = 10;
const DIST_ROOT: u32 = 8;
/// A complete code has at most `symbols / 2` prefixes that need a subtable, each of at most
/// `2^(15 - root)` entries, on top of the root.
const LIT_TABLE: usize = (1 << LIT_ROOT) + 143 * (1 << (15 - LIT_ROOT));
const DIST_TABLE: usize = (1 << DIST_ROOT) + 15 * (1 << (15 - DIST_ROOT));

/// A table entry: bits to consume (low five), then the kind, then extra bits or a subtable's size,
/// then the value in the high half.
const KIND_LITERAL: u32 = 0 << 5;
const KIND_BASE: u32 = 1 << 5;
const KIND_END: u32 = 2 << 5;
const KIND_INVALID: u32 = 3 << 5;
const KIND_SUBTABLE: u32 = 4 << 5;
const KIND_MASK: u32 = 7 << 5;

const fn entry(bits: u32, kind: u32, extra: u32, value: u32) -> u32 {
    bits | kind | (extra << 8) | (value << 16)
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
        0..=255 => entry(bits, KIND_LITERAL, 0, symbol as u32),
        256 => entry(bits, KIND_END, 0, 0),
        257..=285 => entry(
            bits,
            KIND_BASE,
            LENGTH_EXTRA[symbol - 257] as u32,
            LENGTH_BASE[symbol - 257] as u32,
        ),
        // Reachable only through the static code, where zlib calls them invalid.
        _ => entry(bits, KIND_INVALID, 0, 0),
    }
}

fn distance_entry(symbol: usize, bits: u32) -> u32 {
    match symbol {
        0..=29 => entry(
            bits,
            KIND_BASE,
            DIST_EXTRA[symbol] as u32,
            DIST_BASE[symbol] as u32,
        ),
        _ => entry(bits, KIND_INVALID, 0, 0),
    }
}

/// Builds a two-level decoding table from code lengths; `false` unless they exactly fill the
/// codespace, which is stricter than zlib's rule and also covers its over-subscribed refusals.
fn build(table: &mut [u32], root: u32, lens: &[u8], make: fn(usize, u32) -> u32) -> bool {
    let mut count = [0u16; 16];
    for &len in lens {
        count[len as usize] += 1;
    }
    count[0] = 0;
    let mut left = 1i32;
    for &n in &count[1..] {
        left = (left << 1) - n as i32;
        if left < 0 {
            return false;
        }
    }
    if left != 0 {
        return false;
    }
    let longest = (1..16).rev().find(|&l| count[l] != 0).unwrap_or(0) as u32;
    let sub_bits = longest.saturating_sub(root);

    let mut next = [0u16; 16];
    let mut code = 0u16;
    for len in 1..16 {
        code = (code + count[len - 1]) << 1;
        next[len] = code;
    }
    let root_size = 1usize << root;
    let mut used = root_size;
    for slot in &mut table[..root_size] {
        *slot = 0;
    }
    for (symbol, &len) in lens.iter().enumerate() {
        let len = len as u32;
        if len == 0 {
            continue;
        }
        let code = next[len as usize];
        next[len as usize] += 1;
        let reversed = (code.reverse_bits() >> (16 - len)) as usize;
        if len <= root {
            let value = make(symbol, len);
            let mut at = reversed;
            while at < root_size {
                table[at] = value;
                at += 1 << len;
            }
        } else {
            let prefix = reversed & (root_size - 1);
            if table[prefix] == 0 {
                let size = 1usize << sub_bits;
                if used + size > table.len() {
                    return false;
                }
                table[prefix] = entry(root, KIND_SUBTABLE, sub_bits, used as u32);
                used += size;
            }
            let offset = (table[prefix] >> 16) as usize;
            let value = make(symbol, len - root);
            let mut at = reversed >> root;
            while at < 1 << sub_bits {
                table[offset + at] = value;
                at += 1 << (len - root);
            }
        }
    }
    true
}

/// One thread's or one reader's decoding tables.
pub(crate) struct Inflater {
    literal: Box<[u32]>,
    distance: Box<[u32]>,
    /// Whether the tables already hold the static code, which many small blocks use in a row.
    static_tables: bool,
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
            entry(bits, KIND_LITERAL, 0, symbol as u32)
        }) {
            return false;
        }

        let total = literals + distances;
        let mut lens = [0u8; 320];
        let mut i = 0;
        while i < total {
            bits.refill();
            let e = precode[(bits.buf & 127) as usize];
            let used = e & 31;
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
    #[inline]
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
            // The fast loop: input for both refills, each of which then leaves at least 56 bits,
            // and room for a whole match plus a word's overrun, so nothing inside it needs checking
            // but the distance.
            while bits.at + 16 <= bits.input.len() && at + 258 + 8 <= size {
                bits.refill();
                let mut e = *literal.add((bits.buf & LIT_MASK) as usize);
                if e & KIND_MASK == KIND_LITERAL {
                    bits.buf >>= e & 31;
                    bits.left -= e & 31;
                    *out.add(at) = (e >> 16) as u8;
                    at += 1;
                    e = *literal.add((bits.buf & LIT_MASK) as usize);
                    if e & KIND_MASK == KIND_LITERAL {
                        bits.buf >>= e & 31;
                        bits.left -= e & 31;
                        *out.add(at) = (e >> 16) as u8;
                        at += 1;
                        continue;
                    }
                    // Fewer than 48 bits may be left, and a match takes up to 48.
                    bits.refill();
                }
                if e & KIND_MASK == KIND_SUBTABLE {
                    let root = e & 31;
                    bits.buf >>= root;
                    bits.left -= root;
                    let offset = (e >> 16) as usize;
                    let sub = (bits.buf & ((1u64 << ((e >> 8) & 31)) - 1)) as usize;
                    e = *literal.add(offset + sub);
                }
                bits.buf >>= e & 31;
                bits.left -= e & 31;
                match e & KIND_MASK {
                    KIND_LITERAL => {
                        *out.add(at) = (e >> 16) as u8;
                        at += 1;
                        continue;
                    }
                    KIND_BASE => {}
                    KIND_END => {
                        *pos = at;
                        return true;
                    }
                    // "invalid literal/length code".
                    _ => return false,
                }
                let extra = (e >> 8) & 31;
                let length = (e >> 16) as usize + (bits.buf & ((1u64 << extra) - 1)) as usize;
                bits.buf >>= extra;
                bits.left -= extra;

                let mut d = *distance.add((bits.buf & DIST_MASK) as usize);
                if d & KIND_MASK == KIND_SUBTABLE {
                    let root = d & 31;
                    bits.buf >>= root;
                    bits.left -= root;
                    let offset = (d >> 16) as usize;
                    let sub = (bits.buf & ((1u64 << ((d >> 8) & 31)) - 1)) as usize;
                    d = *distance.add(offset + sub);
                }
                // "invalid distance code".
                if d & KIND_MASK != KIND_BASE {
                    return false;
                }
                bits.buf >>= d & 31;
                bits.left -= d & 31;
                let extra = (d >> 8) & 31;
                let dist = (d >> 16) as usize + (bits.buf & ((1u64 << extra) - 1)) as usize;
                bits.buf >>= extra;
                bits.left -= extra;
                // "invalid distance too far back".
                if dist > at {
                    return false;
                }
                let dst = out.add(at);
                let src = dst.sub(dist);
                if dist >= 8 {
                    // Each word's source lies wholly behind its destination, so it is all written.
                    let mut done = 0;
                    while done < length {
                        dst.add(done)
                            .cast::<u64>()
                            .write_unaligned(src.add(done).cast::<u64>().read_unaligned());
                        done += 8;
                    }
                } else if dist == 1 {
                    std::ptr::write_bytes(dst, *src, length);
                } else {
                    for i in 0..length {
                        *dst.add(i) = *src.add(i);
                    }
                }
                at += length;
            }

            // The tail: every step checked against what is really left of the input and the output.
            bits.refill();
            let mut e = *literal.add((bits.buf & LIT_MASK) as usize);
            if e & KIND_MASK == KIND_SUBTABLE {
                let root = e & 31;
                if root > bits.left {
                    return false;
                }
                bits.buf >>= root;
                bits.left -= root;
                let sub = (bits.buf & ((1u64 << ((e >> 8) & 31)) - 1)) as usize;
                e = *literal.add((e >> 16) as usize + sub);
            }
            if e & 31 > bits.left {
                return false;
            }
            bits.buf >>= e & 31;
            bits.left -= e & 31;
            match e & KIND_MASK {
                KIND_LITERAL => {
                    if at >= size {
                        return false;
                    }
                    *out.add(at) = (e >> 16) as u8;
                    at += 1;
                }
                KIND_END => {
                    *pos = at;
                    return true;
                }
                KIND_BASE => {
                    let Some(more) = bits.take((e >> 8) & 31) else {
                        return false;
                    };
                    let length = (e >> 16) as usize + more as usize;
                    bits.refill();
                    let mut d = *distance.add((bits.buf & DIST_MASK) as usize);
                    if d & KIND_MASK == KIND_SUBTABLE {
                        let root = d & 31;
                        if root > bits.left {
                            return false;
                        }
                        bits.buf >>= root;
                        bits.left -= root;
                        let sub = (bits.buf & ((1u64 << ((d >> 8) & 31)) - 1)) as usize;
                        d = *distance.add((d >> 16) as usize + sub);
                    }
                    if d & KIND_MASK != KIND_BASE || d & 31 > bits.left {
                        return false;
                    }
                    bits.buf >>= d & 31;
                    bits.left -= d & 31;
                    let Some(more) = bits.take((d >> 8) & 31) else {
                        return false;
                    };
                    let dist = (d >> 16) as usize + more as usize;
                    if dist > at || at + length > size {
                        return false;
                    }
                    for i in at..at + length {
                        *out.add(i) = *out.add(i - dist);
                    }
                    at += length;
                }
                _ => return false,
            }
        }
    }
}

/// A stored block: back to the byte boundary, the length and its complement, the bytes.
///
/// # Safety
/// As for [`Inflater::codes`].
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
        let mut inflater = Inflater::new();
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
                let mut out = vec![0u8; payload.len()];
                assert!(inflater.inflate(&stream, &mut out), "{level} {deflater:?}");
                assert_eq!(out, *payload);

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
                    let mut out = vec![0u8; size];
                    if inflater.inflate(&damaged, &mut out) {
                        assert_eq!(Some(out), zlib(&damaged, size), "level {level}");
                    }
                    trials += 1;
                }
            }
        }
        assert!(trials > 100_000);
    }
}
