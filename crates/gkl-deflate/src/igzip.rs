//! Levels 1 and 2: ISA-L igzip as GKL ships it, in Rust, on any CPU.
//!
//! Ported from ISA-L 2.31 (BSD-3-Clause, Intel): the x86 kernels `igzip_icf_body_h1_gr_bt.asm`
//! and `igzip_icf_finish.asm`, and the C around them in `igzip.c`, `huff_codes.c`,
//! `flatten_ll.c`, `encode_df.c` and `proc_heap_base.c`.
//!
//! ## Why the assembly and not the C
//!
//! ISA-L ships two compressors. Its readable C (`igzip_icf_base.c`) is not what GKL runs: on the
//! same fixtures it makes 19749 bytes where the kernels make 19044 (decision 0034). The kernels
//! hash with `crc32` (SSE4.2), look at two positions per step, start skipping blind after 2 KB
//! without a match, and emit matches longer than 258 in a loop of their own; the C does none of
//! that. Decision 0033 measured that the SSE4.2, AVX and AVX2 builds of the kernels agree byte for
//! byte, so one scalar transcription of them is the target. Everything below `body` and `finish`
//! is shared C, translated as it is.
//!
//! ## The configuration
//!
//! What GKL passes to `isal_deflate_stateless` (decisions 0031 and 0043): the Java level as the
//! ISA-L level, `ISAL_DEF_LVL1_DEFAULT` or `ISAL_DEF_LVL2_DEFAULT` bytes of level buffer,
//! `end_of_stream = 1`, no flush, raw deflate. The level buffer's size decides where blocks split,
//! which is why [`ICF_BUF_BYTES`] is derived rather than rounded.
//!
//! ## One thing the kernels do that cannot be reproduced
//!
//! On a fresh stream the body kernel's first-byte path hashes a register that still holds the
//! address of the `isal_zstream` (`rcx`, read as `tmp2`) and uses that bucket as the one for the
//! third byte. So GKL looks up and overwrites one pointer-dependent bucket of its hash table. Here
//! that bucket is a slot outside the table holding the table's initial value, which is what the
//! real bucket holds unless the address happens to hash to one of the first two bytes' buckets. In
//! that case GKL's own output varies from run to run, so no fixed answer is GKL's.

/// The bytes a level buffer leaves for intermediate tokens after its fixed part, the same at both
/// levels: `level_buf_size - (sizeof(struct level_buf) - MAX_LVL_BUF_SIZE + hash table) - 4`.
/// Level 1 is 282624 - (87800 - 83088 + 16384) - 4, level 2 is 331776 - (87800 - 83088 + 65536)
/// - 4. Decision 0031 is what a wrong figure here costs: 62 bytes on a 200 KB fixture.
const ICF_BUF_BYTES: u64 = 261_524;

/// `sizeof(state->buffer)`, which only enters the stored-block test.
const STATE_BUFFER: u64 = 2 * 32768 + 288;
const MAX_HDR_SIZE: usize = 328;

const DIST_MASK: u64 = 0x7fff;
const LIT: u32 = 30 << 10;
const NULL_DIST_SYM: u32 = 30;

// The body kernel's constants.
const SLOP: u64 = 8;
const MIN_INBUF_PADDING: i64 = 16;
const MAX_EMIT_SIZE: i64 = 258 * 16;
const LARGE_MATCH_MIN: u64 = 264;
const SKIP_SIZE_BASE: i64 = 2 << 10;
const SKIP_BASE: i64 = 32;
const SKIP_START: i64 = 512;
const SKIP_RATE: i64 = 2;
const MAX_SKIP_LEVEL: i64 = (128 - SKIP_BASE) / SKIP_RATE + SKIP_START;

const LEN_EXTRA_BITS: [u32; 29] = [
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0,
];
const DIST_EXTRA_BITS: [u32; 30] = [
    0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13,
    13,
];
const CODE_LENGTH_ORDER: [usize; 19] = [
    16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15,
];

/// `repeated_char_header`, the fixed dynamic header ISA-L writes for a run of 0x00 or 0xFF.
const REPEATED_CHAR_HEADER: [[u32; 4]; 2] = [
    [0x0121c0ec, 0xc30c0000, 0x7d57fab0, 0x49270938],
    [0x0121c0ec, 0xc30c0000, 0x7baaff30, 0x49270938],
];
const MIN_REPEAT_LEN: usize = 4 * 1024;

/// Compress `data` exactly as `IntelDeflater` does at Java level `level`, which is 1 or 2.
pub fn deflate(data: &[u8], level: usize) -> Vec<u8> {
    assert!(
        level == 1 || level == 2,
        "igzip levels are 1 and 2, got {level}"
    );
    assert!(
        u32::try_from(data.len()).is_ok(),
        "ISA-L counts input in 32 bits"
    );
    let mut stream = Stream::new(data, level);
    if stream.run() {
        stream.mem.truncate(stream.next_out);
        stream.mem
    } else {
        stored_blocks(data)
    }
}

/// `isal_deflate_stateless`'s fallback when the compressed stream does not fit in the stored
/// size: the whole input as stored blocks of at most 65535 bytes.
fn stored_blocks(data: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(data.len() + 5 * (data.len() / 65535 + 1));
    let mut chunks = data.chunks(65535).peekable();
    if chunks.peek().is_none() {
        out.extend_from_slice(&[1, 0, 0, 0xff, 0xff]);
    }
    while let Some(chunk) = chunks.next() {
        let len = chunk.len() as u16;
        out.push(chunks.peek().is_none() as u8);
        out.extend_from_slice(&len.to_le_bytes());
        out.extend_from_slice(&(!len).to_le_bytes());
        out.extend_from_slice(chunk);
    }
    out
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum State {
    NewHdr,
    Body,
    FlushReadBuffer,
    CreateHdr,
    Hdr,
    FlushIcfBuffer,
    Type0Hdr,
    Type0Body,
    Trl,
    End,
}

/// `struct BitBuf2`, writing into `Stream::mem` at indices. Like the C it stores eight bytes at a
/// time, so `mem` keeps slack past every limit.
#[derive(Clone, Copy)]
struct BitBuf {
    bits: u64,
    count: u32,
    start: usize,
    pos: usize,
    /// `m_out_end`: the last position a store may start at, `buf + len - 8`. Signed, because a
    /// buffer shorter than eight bytes puts it before the start.
    end: i64,
}

impl BitBuf {
    fn set_buf(&mut self, at: usize, len: u64) {
        self.start = at;
        self.pos = at;
        self.end = at as i64 + len as i64 - 8;
    }

    fn is_full(&self) -> bool {
        self.pos as i64 > self.end
    }

    fn used(&self) -> u64 {
        (self.pos - self.start) as u64
    }

    fn write_unsafe(&mut self, code: u64, count: u32) {
        self.bits |= code << self.count;
        self.count += count;
    }

    fn flush_bits(&mut self, mem: &mut [u8]) {
        mem[self.pos..self.pos + 8].copy_from_slice(&self.bits.to_le_bytes());
        let bits = self.count & !7;
        self.count -= bits;
        self.pos += bits as usize / 8;
        self.bits = if bits == 64 { 0 } else { self.bits >> bits };
    }

    fn write(&mut self, mem: &mut [u8], code: u64, count: u32) {
        self.write_unsafe(code, count);
        self.flush_bits(mem);
    }

    fn flush(&mut self, mem: &mut [u8]) {
        if self.count != 0 {
            mem[self.pos..self.pos + 8].copy_from_slice(&self.bits.to_le_bytes());
            self.pos += self.count.div_ceil(8) as usize;
        }
        self.bits = 0;
        self.count = 0;
    }
}

struct Stream<'a> {
    data: &'a [u8],
    state: State,

    // The input side. Positions are absolute: `file_start` is `data`, so `next_in == total_in`.
    total_in: u64,
    avail_in: u64,
    block_next: u64,
    block_end: u64,
    has_hist: bool,
    hash_mask: u64,
    /// The hash table, plus one slot past the mask for the bucket GKL derives from a pointer.
    table: Vec<u16>,

    // The output side: `mem[..out_cap]` is the caller's buffer, then the header scratch buffer.
    mem: Vec<u8>,
    hdr_at: usize,
    next_out: usize,
    avail_out: u64,
    bb: BitBuf,
    count: u64,
    has_eob_hdr: bool,

    // The level buffer.
    icf: Vec<u32>,
    icf_next: u64,
    icf_avail_out: u64,
    ll_hist: [u32; 513],
    d_hist: [u32; 30],
    tables: [u32; 544],
    hdr_count: u32,
    hdr_extra_bits: u32,
}

/// `_mm_crc32_u32(0, data)`: CRC-32C of four bytes from zero, without the final inversion. This
/// is the SSE4.2 instruction the kernels hash with, and the one thing decision 0033 says they
/// cannot do without.
fn crc32c(data: u64) -> u64 {
    static TABLES: std::sync::OnceLock<[[u32; 256]; 4]> = std::sync::OnceLock::new();
    let t = TABLES.get_or_init(|| {
        let mut t = [[0u32; 256]; 4];
        for (i, slot) in t[0].iter_mut().enumerate() {
            *slot = (0..8).fold(i as u32, |c, _| (c >> 1) ^ (0x82f6_3b78 * (c & 1)));
        }
        for k in 1..4 {
            let (done, rest) = t.split_at_mut(k);
            for (slot, &prev) in rest[0].iter_mut().zip(&done[k - 1]) {
                *slot = (prev >> 8) ^ done[0][(prev & 0xff) as usize];
            }
        }
        t
    });
    let w = data as u32;
    (t[3][(w & 0xff) as usize]
        ^ t[2][((w >> 8) & 0xff) as usize]
        ^ t[1][((w >> 16) & 0xff) as usize]
        ^ t[0][(w >> 24) as usize]) as u64
}

/// `get_dist_icf_code` on `d`, which is the distance minus one: the symbol in bits 10..19 and the
/// extra bits from bit 19, as the ICF token wants them.
fn dist_icf_code(d: u64) -> u32 {
    let d = d as u32;
    let code = if d <= 1 {
        d
    } else {
        let nbits = 31 - d.leading_zeros() - 1;
        let extra = d & ((1 << nbits) - 1);
        (extra << 9) + (d >> nbits) + 2 * nbits
    };
    code << 10
}

impl<'a> Stream<'a> {
    fn new(data: &'a [u8], level: usize) -> Self {
        let n = data.len() as u64;
        let stored_len = if n == 0 { 5 } else { 5 * n.div_ceil(65535) + n };
        let mut hash_mask = if level == 2 { 32767 } else { 8191 };
        if hash_mask > 2 * n {
            hash_mask = (1u64 << (64 - n.leading_zeros())) - 1;
        }
        let out_cap = stored_len as usize + 512;
        let hdr_at = out_cap;
        Stream {
            data,
            state: State::NewHdr,
            total_in: 0,
            avail_in: n,
            block_next: 0,
            block_end: 0,
            has_hist: false,
            hash_mask,
            table: vec![0; hash_mask as usize + 2],
            mem: vec![0; hdr_at + MAX_HDR_SIZE + 16],
            hdr_at,
            next_out: 0,
            avail_out: stored_len,
            bb: BitBuf {
                bits: 0,
                count: 0,
                start: 0,
                pos: 0,
                end: 0,
            },
            count: 0,
            has_eob_hdr: false,
            icf: vec![0; ICF_BUF_BYTES as usize / 4 + 64],
            icf_next: 0,
            icf_avail_out: 0,
            ll_hist: [0; 513],
            d_hist: [0; 30],
            tables: [0; 544],
            hdr_count: 0,
            hdr_extra_bits: 0,
        }
    }

    fn byte(&self, i: i64) -> u32 {
        self.data[i as usize] as u32
    }

    /// An unaligned little-endian load of `N` bytes. The kernels never read past the input in a
    /// stateless call; the zero padding only keeps a mistake here from being a panic.
    fn load<const N: usize>(&self, i: i64) -> u128 {
        debug_assert!(i >= 0);
        let i = i as usize;
        let mut b = [0u8; 16];
        if let Some(s) = self.data.get(i..(i + N).min(self.data.len())) {
            b[..s.len()].copy_from_slice(s);
        }
        u128::from_le_bytes(b)
    }

    fn ld64(&self, i: i64) -> u64 {
        self.load::<8>(i) as u64
    }

    /// How many bytes from `start` match the bytes `dist` back, up to `max`, knowing the first
    /// `known` already do: the `compare_large` and `compare` macros, which all stop exactly at
    /// `max`.
    fn compare(&self, src: i64, cur: i64, known: u64, max: u64) -> u64 {
        let (src, cur) = (src as usize, cur as usize);
        let mut len = known as usize;
        let max = max as usize;
        while len + 8 <= max {
            let a = u64::from_le_bytes(self.data[src + len..src + len + 8].try_into().unwrap());
            let b = u64::from_le_bytes(self.data[cur + len..cur + len + 8].try_into().unwrap());
            if a != b {
                return (len + ((a ^ b).trailing_zeros() / 8) as usize) as u64;
            }
            len += 8;
        }
        while len < max && self.data[src + len] == self.data[cur + len] {
            len += 1;
        }
        len as u64
    }

    fn set(&mut self, bucket: u64, pos: i64) {
        self.table[bucket as usize] = pos as u16;
    }

    fn emit(&mut self, mo: &mut u64, token: u32) {
        self.icf[(*mo / 4) as usize] = token;
        *mo += 4;
    }

    /// `isal_deflate_stateless` from the point where it has checked its arguments. True when the
    /// stream ended normally; false is `STATELESS_OVERFLOW`, after which the caller stores.
    fn run(&mut self) -> bool {
        if self.avail_in >= 8 {
            let w = self.ld64(0);
            if w == 0 || w == u64::MAX {
                let repeated = self.repeated_length();
                if self.avail_in == repeated || repeated >= MIN_REPEAT_LEN as u64 {
                    self.write_constant(repeated);
                }
            }
        }
        if self.state == State::NewHdr {
            self.reset_match_history();
        }
        self.count = 0;
        let start_in = self.total_in;
        loop {
            if self.state == State::NewHdr {
                self.init_new_icf_block();
            }
            if self.state == State::Body {
                self.body();
            }
            if self.state == State::FlushReadBuffer {
                self.finish();
            }
            if self.state == State::CreateHdr {
                self.create_icf_block_hdr(start_in);
            }
            if self.state == State::Hdr {
                self.write_header();
            }
            if self.state == State::FlushIcfBuffer {
                self.flush_icf_block();
            }
            if matches!(self.state, State::Type0Hdr | State::Type0Body) {
                self.write_stored_block();
            }
            if self.state != State::NewHdr {
                break;
            }
        }
        if self.state == State::Trl {
            self.write_trailer();
        }
        self.state == State::End
    }

    fn repeated_length(&self) -> u64 {
        let c = self.data[0];
        let mut i = 0;
        while i + 8 <= self.data.len() && self.data[i..i + 8] == self.data[..8] {
            i += 8;
        }
        while i < self.data.len() && self.data[i] == c {
            i += 1;
        }
        i as u64
    }

    /// `write_constant_compressed_stateless`.
    fn write_constant(&mut self, repeated: u64) {
        let rep_bits = ((repeated - 1) / 258) * 2;
        let rep_bytes = rep_bits / 8;
        let mut rep_extra = (repeated - 1) % 258;
        let c = self.data[0];
        if self.avail_out < 16 + 8 + rep_bytes + 8 {
            return;
        }
        for (k, w) in REPEATED_CHAR_HEADER[(c & 1) as usize].iter().enumerate() {
            self.mem[self.next_out + 4 * k..self.next_out + 4 * k + 4]
                .copy_from_slice(&w.to_le_bytes());
        }
        if self.avail_in == repeated {
            self.mem[self.next_out] |= 1;
            self.has_eob_hdr = true;
            self.state = State::Trl;
        } else {
            self.state = State::NewHdr;
        }
        let zeros = self.next_out + 16;
        self.mem[zeros..zeros + rep_bytes as usize].fill(0);
        self.avail_out -= 16 + rep_bytes;
        self.next_out += 16 + rep_bytes as usize;
        self.bb.set_buf(self.next_out, self.avail_out);
        self.bb.bits = 0;
        self.bb.count = (rep_bits % 8) as u32;
        let mem = &mut self.mem;
        let bb = &mut self.bb;
        if rep_extra > 115 {
            while rep_extra > 130 && rep_extra < 230 {
                bb.write(mem, 0x3, 4);
                rep_extra -= 10;
            }
            if rep_extra >= 230 {
                bb.write(mem, 0x0f | ((rep_extra / 2 - 115) << 4), 9);
                rep_extra -= rep_extra / 2;
            }
            bb.write(mem, 0x0f | ((rep_extra - 115) << 4), 9);
        } else {
            while rep_extra >= 10 {
                bb.write(mem, 0x3, 4);
                rep_extra -= 10;
            }
            for _ in 0..rep_extra {
                bb.write(mem, 0x1, 2);
            }
        }
        bb.write(mem, 0x7, 4);
        self.total_in += repeated;
        self.avail_in -= repeated;
        self.block_end += repeated;
        let bytes = self.bb.used();
        self.next_out = self.bb.pos;
        self.avail_out -= bytes;
    }

    fn reset_match_history(&mut self) {
        self.has_hist = false;
        let fill = self.total_in as u16;
        self.table.fill(fill);
    }

    fn init_new_icf_block(&mut self) {
        self.block_next = self.block_end;
        self.icf_next = 0;
        self.icf_avail_out = ICF_BUF_BYTES;
        self.ll_hist = [0; 513];
        self.d_hist = [0; 30];
        self.state = State::Body;
    }

    /// `isal_deflate_icf_body_hash_hist_04`, the kernel both levels dispatch to. Labels are the
    /// assembly's; locals carry its register names where that helps a comparison.
    fn body(&mut self) {
        if self.avail_in == 0 {
            self.state = State::CreateHdr;
            return;
        }
        #[derive(Clone, Copy)]
        enum L {
            Loop2,
            LenDistLitHuffman,
            LenDistHuffman,
            WriteLitBits,
            CompareLoop,
            CompareLoop2,
            DoEmit,
            WriteFirstByte,
            InputEnd,
            OutputEnd,
        }
        let mask = self.hash_mask;
        let phantom = self.hash_mask + 1;
        let m_out_start = self.icf_next;
        let mut mo = self.icf_next;
        let m_out_end = self.icf_next + self.icf_avail_out - SLOP;
        let mut f_i = self.total_in as i64;
        let slop = MIN_INBUF_PADDING;
        let file_length = f_i + self.avail_in as i64 - slop;
        let mut skip_match = f_i + SKIP_SIZE_BASE;
        let mut skip_level: i64 = 0;

        let mut xdata: u128 = 0;
        let mut curr_data: u64 = 0;
        let (mut hash, mut hash2) = (0u64, 0u64);
        let (mut d1, mut d2) = (0u64, 0u64);
        let (mut len, mut len2) = (0u64, 0u64);
        let (mut lit_code, mut lit_code2) = (0u32, 0u32);

        let mut at = if file_length <= f_i {
            L::InputEnd
        } else {
            xdata = self.load::<16>(f_i);
            curr_data = self.ld64(f_i);
            hash = crc32c(curr_data) & mask;
            hash2 = crc32c(curr_data >> 8) & mask;
            if self.has_hist {
                L::Loop2
            } else {
                L::WriteFirstByte
            }
        };

        loop {
            at = match at {
                L::Loop2 => {
                    if mo > m_out_end {
                        L::OutputEnd
                    } else {
                        let p = f_i;
                        d1 = ((p - 1) as u16).wrapping_sub(self.table[hash as usize]) as u64
                            & DIST_MASK;
                        self.set(hash, p);
                        f_i = p + 1;
                        let tmp2 = curr_data;
                        curr_data >>= 16;
                        hash = crc32c(curr_data) & mask;
                        d2 = (p as u16).wrapping_sub(self.table[hash2 as usize]) as u64 & DIST_MASK;
                        self.set(hash2, p + 1);
                        hash2 = crc32c(tmp2 >> 24) & mask;
                        let first = xdata as u64;
                        curr_data = first;
                        xdata >>= 8;
                        let x1 = first ^ self.ld64(p - d1 as i64 - 1);
                        let x2 = (xdata as u64) ^ self.ld64(p - d2 as i64);
                        if x1 == 0 {
                            L::CompareLoop
                        } else if x2 == 0 {
                            L::CompareLoop2
                        } else {
                            lit_code = (curr_data & 0xff) as u32;
                            curr_data >>= 8;
                            if x1 as u32 == 0 {
                                len = (x1.trailing_zeros() / 8) as u64;
                                L::LenDistHuffman
                            } else {
                                xdata >>= 8;
                                self.ll_hist[lit_code as usize] += 1;
                                lit_code2 = (curr_data & 0xff) as u32;
                                if x2 as u32 != 0 {
                                    L::WriteLitBits
                                } else {
                                    len2 = (x2.trailing_zeros() / 8) as u64;
                                    L::LenDistLitHuffman
                                }
                            }
                        }
                    }
                }
                L::LenDistLitHuffman => {
                    // A literal at f_i - 1, then a match of len2 at f_i, d2 back.
                    self.icf[(mo / 4) as usize] = lit_code | LIT;
                    let dcode = dist_icf_code(d2);
                    let k = f_i + 1;
                    f_i += len2 as i64;
                    let token = dcode + 254 + len2 as u32;
                    if f_i > file_length {
                        self.set(hash, k);
                        self.set(hash2, k + 1);
                    } else {
                        skip_match = f_i + SKIP_SIZE_BASE;
                        skip_level -= len2 as i64;
                        xdata = self.load::<16>(f_i);
                        let next = self.ld64(f_i);
                        let hash3 = crc32c(curr_data >> 24) & mask;
                        curr_data = next;
                        self.set(hash, k);
                        hash = crc32c(curr_data);
                        self.set(hash2, k + 1);
                        hash2 = crc32c(next >> 8);
                        self.set(hash3, k + 2);
                    }
                    self.ll_hist[(len2 + 254) as usize] += 1;
                    self.icf[(mo / 4) as usize + 1] = token;
                    mo += 8;
                    self.d_hist[((token >> 10) & 0x1f) as usize] += 1;
                    if f_i > file_length {
                        L::InputEnd
                    } else {
                        hash &= mask;
                        hash2 &= mask;
                        L::Loop2
                    }
                }
                L::LenDistHuffman => {
                    // A match of len at f_i - 1, d1 back.
                    f_i -= 1;
                    let k = f_i + 2;
                    let mut dcode = dist_icf_code(d1);
                    f_i += len as i64;
                    dcode |= len as u32 + 254;
                    if f_i <= file_length {
                        skip_match = f_i + SKIP_SIZE_BASE;
                        skip_level -= len as i64;
                        xdata = self.load::<16>(f_i);
                        curr_data = self.ld64(f_i);
                    }
                    self.set(hash, k);
                    self.set(hash2, k + 1);
                    self.ll_hist[(len + 254) as usize] += 1;
                    self.emit(&mut mo, dcode);
                    self.d_hist[((dcode >> 10) & 0x1f) as usize] += 1;
                    if f_i > file_length {
                        L::InputEnd
                    } else {
                        hash = crc32c(curr_data) & mask;
                        hash2 = crc32c(curr_data >> 8) & mask;
                        L::Loop2
                    }
                }
                L::WriteLitBits => {
                    // Literals at f_i - 1 and f_i, one token.
                    curr_data = xdata as u64;
                    f_i += 1;
                    self.ll_hist[lit_code2 as usize] += 1;
                    let pair = lit_code + (lit_code2 << 10) + (31 << 10);
                    self.emit(&mut mo, pair);
                    if f_i > file_length {
                        L::InputEnd
                    } else {
                        xdata = self.load::<16>(f_i);
                        if f_i <= skip_match {
                            L::Loop2
                        } else {
                            // Two kilobytes without a match: skip ahead blind, further each time.
                            let mut rcx = (skip_level + 1).clamp(0, MAX_SKIP_LEVEL);
                            skip_match = if rcx > 63 { 0 } else { SKIP_SIZE_BASE >> rcx };
                            skip_level = rcx;
                            rcx = (rcx - SKIP_START).max(0);
                            let mut skip_count = (SKIP_RATE * rcx + SKIP_BASE) & -SKIP_BASE;
                            let room = (m_out_end + 4 - mo) >> 1;
                            let short =
                                (room as i64) < skip_count || slop + file_length - f_i < skip_count;
                            if short {
                                loop {
                                    let a = self.byte(f_i);
                                    let b = self.byte(f_i + 1);
                                    self.ll_hist[a as usize] += 1;
                                    self.ll_hist[b as usize] += 1;
                                    self.emit(&mut mo, a + (b << 10) + (31 << 10));
                                    f_i += 2;
                                    if mo > m_out_end {
                                        break L::OutputEnd;
                                    }
                                    if file_length <= f_i {
                                        break L::InputEnd;
                                    }
                                }
                            } else {
                                loop {
                                    for k in 0..4 {
                                        let a = self.byte(f_i + 2 * k);
                                        let b = self.byte(f_i + 2 * k + 1);
                                        self.ll_hist[a as usize] += 1;
                                        self.ll_hist[b as usize] += 1;
                                        self.emit(&mut mo, a + (b << 10) + (31 << 10));
                                    }
                                    f_i += 8;
                                    skip_count -= 8;
                                    if skip_count <= 0 {
                                        break;
                                    }
                                }
                                if file_length <= f_i {
                                    L::InputEnd
                                } else {
                                    // The hashes are not recomputed: the next step looks up the
                                    // buckets of the bytes before the skip, as the kernel does.
                                    curr_data = self.ld64(f_i);
                                    xdata = self.load::<16>(f_i);
                                    skip_match += f_i;
                                    L::Loop2
                                }
                            }
                        }
                    }
                }
                L::CompareLoop => {
                    let p = f_i - 1;
                    let max = (file_length - f_i + slop + 1).min(MAX_EMIT_SIZE) as u64;
                    len = self.compare(p - d1 as i64 - 1, p, 8, max);
                    if len <= 258 {
                        L::LenDistHuffman
                    } else if len >= LARGE_MATCH_MIN {
                        L::DoEmit
                    } else {
                        len = 258;
                        L::LenDistHuffman
                    }
                }
                L::CompareLoop2 => {
                    let p = f_i - 1;
                    let max = (file_length - f_i + slop).min(MAX_EMIT_SIZE) as u64;
                    len2 = self.compare(p - d2 as i64, p + 1, 8, max);
                    lit_code = (curr_data & 0xff) as u32;
                    curr_data >>= 8;
                    self.ll_hist[lit_code as usize] += 1;
                    if len2 <= 258 {
                        L::LenDistLitHuffman
                    } else if len2 >= LARGE_MATCH_MIN {
                        // do_emit2: the literal, then the long match one byte on.
                        self.emit(&mut mo, lit_code | LIT);
                        f_i += 1;
                        d1 = d2;
                        len = len2;
                        L::DoEmit
                    } else {
                        len2 = 258;
                        L::LenDistLitHuffman
                    }
                }
                L::DoEmit => {
                    // A match longer than 258 at f_i - 1, d1 back, as 258s and a remainder. From
                    // do_emit2 the two buckets below are one byte early; the kernel does the same.
                    let mut dcode = dist_icf_code(d1);
                    let long = (258 + 254) | dcode;
                    let dsym = ((dcode >> 10) & 0x1f) as usize;
                    let k = f_i + 1;
                    f_i -= 1;
                    self.set(hash, k);
                    self.set(hash2, k + 1);
                    let mut full = false;
                    loop {
                        len -= 258;
                        f_i += 258;
                        self.ll_hist[512] += 1;
                        self.d_hist[dsym] += 1;
                        self.emit(&mut mo, long);
                        if mo > m_out_end {
                            full = true;
                            break;
                        }
                        if len < LARGE_MATCH_MIN {
                            break;
                        }
                    }
                    if full {
                        L::OutputEnd
                    } else {
                        len = len.min(258);
                        f_i += len as i64;
                        dcode |= len as u32 + 254;
                        self.ll_hist[(len + 254) as usize] += 1;
                        self.d_hist[dsym] += 1;
                        self.emit(&mut mo, dcode);
                        if file_length <= f_i {
                            L::InputEnd
                        } else {
                            let t = f_i - 4;
                            for k in 0..4 {
                                let h = crc32c(self.ld64(t + k)) & mask;
                                self.set(h, t + k);
                            }
                            xdata = self.load::<16>(f_i);
                            curr_data = self.ld64(f_i);
                            hash = crc32c(curr_data) & mask;
                            hash2 = crc32c(curr_data >> 8) & mask;
                            L::Loop2
                        }
                    }
                }
                L::WriteFirstByte => {
                    if mo > m_out_end {
                        L::OutputEnd
                    } else {
                        self.has_hist = true;
                        self.set(hash, f_i);
                        hash = hash2;
                        hash2 = phantom;
                        let lit = (curr_data & 0xff) as u32;
                        self.ll_hist[lit as usize] += 1;
                        self.emit(&mut mo, lit | LIT);
                        f_i += 1;
                        xdata = self.load::<16>(f_i);
                        curr_data = self.ld64(f_i);
                        if f_i < file_length {
                            L::Loop2
                        } else {
                            L::InputEnd
                        }
                    }
                }
                L::InputEnd => {
                    self.state = State::FlushReadBuffer;
                    break;
                }
                L::OutputEnd => {
                    self.state = State::CreateHdr;
                    break;
                }
            };
        }
        self.total_in = f_i as u64;
        self.block_end = f_i as u64;
        self.avail_in = (file_length + slop - f_i) as u64;
        self.icf_next = mo;
        self.icf_avail_out -= mo - m_out_start;
    }

    /// `isal_deflate_icf_finish_hash_hist_01`: the last bytes, one position at a time.
    fn finish(&mut self) {
        let mask = self.hash_mask;
        let m_out_start = self.icf_next;
        let mut mo = self.icf_next;
        let m_out_end = self.icf_next + self.icf_avail_out - 4;
        let mut f_i = self.total_in as i64;
        let f_end_i = f_i + self.avail_in as i64 - 3;

        let mut literal_first = false;
        let mut stopped = f_i >= f_end_i;
        if !stopped && !self.has_hist {
            if mo > m_out_end {
                stopped = true;
            } else {
                let h = crc32c(self.ld64(f_i)) & mask;
                self.set(h, f_i);
                self.has_hist = true;
                literal_first = true;
            }
        }
        if !stopped {
            loop {
                let mut literal = std::mem::take(&mut literal_first);
                if !literal {
                    if mo > m_out_end {
                        break;
                    }
                    let h = crc32c(self.ld64(f_i)) & mask;
                    let seen = self.table[h as usize] as i64;
                    self.set(h, f_i);
                    let dist = ((f_i - seen) as u32 & 0xffff) as u64;
                    if dist.wrapping_sub(1) as u32 >= DIST_MASK as u32 {
                        literal = true;
                    } else {
                        let max = ((f_end_i - f_i + 3) as u64).min(258);
                        let len = self.compare(f_i - dist as i64, f_i, 0, max);
                        if len < 4 {
                            literal = true;
                        } else {
                            let code = len as u32 + 254;
                            let dcode = dist_icf_code(dist - 1) | code;
                            self.ll_hist[code as usize] += 1;
                            let k = f_i + 1;
                            f_i += len as i64;
                            if (f_i as u64) < f_end_i as u64 {
                                for k in [k, k + 1] {
                                    let h = crc32c(self.ld64(k)) & mask;
                                    self.set(h, k);
                                }
                            }
                            self.emit(&mut mo, dcode);
                            self.d_hist[((dcode >> 10) & 0x1f) as usize] += 1;
                            if f_i >= f_end_i {
                                break;
                            }
                            continue;
                        }
                    }
                }
                if literal {
                    let lit = self.byte(f_i);
                    self.ll_hist[lit as usize] += 1;
                    self.emit(&mut mo, lit | LIT);
                    f_i += 1;
                    if f_i >= f_end_i {
                        break;
                    }
                }
            }
        }
        let f_end_i = f_end_i + 3;
        while f_i < f_end_i && mo <= m_out_end {
            let lit = self.byte(f_i);
            self.ll_hist[lit as usize] += 1;
            self.emit(&mut mo, lit | LIT);
            f_i += 1;
        }
        // end_of_stream is set, so every exit is towards a header.
        self.state = State::CreateHdr;
        self.total_in = f_i as u64;
        self.block_end = f_i as u64;
        self.avail_in = (f_end_i - f_i) as u64;
        self.icf_next = mo;
        self.icf_avail_out -= mo - m_out_start;
    }

    fn create_icf_block_hdr(&mut self, start_in: u64) {
        let saved = self.bb;
        let block_in_size = self.block_end - self.block_next;
        let mut block_size = 5 * block_in_size.div_ceil(65535) + block_in_size;
        if block_size == 0 {
            block_size = 5;
        }
        block_size += (self.bb.count as u64 + 2) / 8;

        self.ll_hist[256] = 1;
        self.icf[(self.icf_next / 4) as usize] = 0x100 | (NULL_DIST_SYM << 10);
        self.icf_next += 4;

        self.has_eob_hdr = self.avail_in == 0;
        let buffer_header = self.avail_out < MAX_HDR_SIZE as u64;
        if buffer_header {
            self.bb.set_buf(self.hdr_at, MAX_HDR_SIZE as u64);
        } else {
            self.bb.set_buf(self.next_out, self.avail_out);
        }
        let bit_count = self.create_hufftables_icf();

        let block_start_offset = self.total_in - self.block_next;
        let cur_in_processed = self.total_in - start_in;
        let avail_output = self.avail_out + STATE_BUFFER - (self.total_in - self.block_end);
        if bit_count / 8 >= block_size
            && cur_in_processed >= block_start_offset
            && block_size <= avail_output
        {
            self.has_eob_hdr = false;
            self.bb = saved;
            self.state = State::Type0Hdr;
        } else if buffer_header {
            self.hdr_count = self.bb.used() as u32;
            self.hdr_extra_bits = self.bb.count;
            self.bb.flush(&mut self.mem);
            self.bb = saved;
            self.bb.bits = 0;
            self.bb.count = 0;
            self.state = State::Hdr;
        } else {
            let used = self.bb.used();
            self.next_out = self.bb.pos;
            self.avail_out -= used;
            self.state = State::FlushIcfBuffer;
        }
    }

    /// `write_header` with the buffered header, `ZSTATE_FLUSH_ICF_BUFFER` next and no toggling.
    fn write_header(&mut self) {
        let hdr_count = self.hdr_count as u64;
        let extra = self.mem[self.hdr_at + hdr_count as usize] as u64;
        self.state = State::Hdr;
        if self.bb.count != 0 {
            if self.avail_out < 8 {
                return;
            }
            self.bb.set_buf(self.next_out, self.avail_out);
            self.bb.flush(&mut self.mem);
            let used = self.bb.used();
            self.next_out = self.bb.pos;
            self.avail_out -= used;
        }
        let mut count = hdr_count - self.count;
        if count != 0 {
            count = count.min(self.avail_out);
            let from = self.hdr_at + self.count as usize;
            self.mem
                .copy_within(from..from + count as usize, self.next_out);
            self.next_out += count as usize;
            self.avail_out -= count;
            self.count += count;
            count = hdr_count - self.count;
        }
        if count == 0 && self.avail_out >= 8 {
            self.bb.set_buf(self.next_out, self.avail_out);
            self.bb.write(&mut self.mem, extra, self.hdr_extra_bits);
            self.state = State::FlushIcfBuffer;
            self.count = 0;
            let used = self.bb.used();
            self.next_out = self.bb.pos;
            self.avail_out -= used;
        }
    }

    /// `flush_icf_block` around `encode_deflate_icf`, which checks for room before every token.
    fn flush_icf_block(&mut self) {
        self.bb.set_buf(self.next_out, self.avail_out);
        let end = (self.icf_next / 4) as usize;
        let mut next = self.count as usize;
        while next < end && !self.bb.is_full() {
            let e = self.icf[next];
            let lsym = self.tables[31 + (e & 0x3ff) as usize];
            let dsym = self.tables[((e >> 10) & 0x1ff) as usize];
            self.bb.write_unsafe((lsym & 0xff_ffff) as u64, lsym >> 24);
            self.bb.write_unsafe((dsym & 0xffff) as u64, dsym >> 24);
            self.bb.write_unsafe((e >> 19) as u64, (dsym >> 16) & 0xff);
            self.bb.flush_bits(&mut self.mem);
            next += 1;
        }
        self.count = next as u64;
        let used = self.bb.used();
        self.next_out = self.bb.pos;
        self.avail_out -= used;
        if next >= end {
            self.count = 0;
            self.state = if self.avail_in == 0 {
                State::Trl
            } else {
                State::NewHdr
            };
        }
    }

    fn write_type0_header(&mut self) {
        let block_in_size = self.block_end - self.block_next;
        let (header, copy_size) = if block_in_size > 65535 {
            (0xffffu64, 65535)
        } else {
            let pending = self.avail_in + (self.total_in - self.block_next);
            if pending == block_in_size {
                self.has_eob_hdr = true;
            }
            (
                ((!(block_in_size as u32) as u64) << 16) | (block_in_size & 0xffff),
                block_in_size,
            )
        };
        let written = if self.bb.count == 0 && self.avail_out >= 5 {
            let h = (header << 8) | self.has_eob_hdr as u64;
            self.mem[self.next_out..self.next_out + 5].copy_from_slice(&h.to_le_bytes()[..5]);
            5
        } else if self.avail_out >= 8 {
            self.bb.set_buf(self.next_out, self.avail_out);
            self.bb.write_unsafe(self.has_eob_hdr as u64, 3);
            self.bb.flush(&mut self.mem);
            let used = self.bb.used();
            self.next_out = self.bb.pos;
            self.avail_out -= used;
            self.mem[self.next_out..self.next_out + 4].copy_from_slice(&header.to_le_bytes()[..4]);
            4
        } else {
            self.has_eob_hdr = false;
            return;
        };
        self.next_out += written;
        self.avail_out -= written as u64;
        self.state = State::Type0Body;
        self.count = copy_size;
    }

    fn write_stored_block(&mut self) {
        loop {
            if self.state == State::Type0Hdr {
                self.write_type0_header();
                if self.state == State::Type0Hdr {
                    break;
                }
            }
            let mut copy_size = self.count;
            let pending = self.avail_in + (self.total_in - self.block_next);
            let from = self.block_next as usize;
            if copy_size > self.avail_out || copy_size > pending {
                copy_size = self.avail_out.min(pending);
                self.count -= copy_size;
            } else {
                self.count = 0;
                self.state = State::Type0Hdr;
            }
            let n = copy_size as usize;
            self.mem[self.next_out..self.next_out + n].copy_from_slice(&self.data[from..from + n]);
            self.block_next += copy_size;
            self.next_out += n;
            self.avail_out -= copy_size;
            if self.block_next == self.block_end {
                self.state = if self.has_eob_hdr {
                    State::Trl
                } else {
                    State::NewHdr
                };
            }
            if self.state != State::Type0Hdr {
                break;
            }
        }
    }

    fn write_trailer(&mut self) {
        self.bb.set_buf(self.next_out, self.avail_out);
        if !self.has_eob_hdr {
            if self.avail_out < 8 {
                return;
            }
            self.has_eob_hdr = true;
            self.bb.write(&mut self.mem, 0x003, 10);
            if self.bb.is_full() {
                let used = self.bb.used();
                self.next_out = self.bb.pos;
                self.avail_out -= used;
                return;
            }
        }
        if self.bb.count != 0 {
            if self.avail_out < 8 {
                return;
            }
            self.bb.flush(&mut self.mem);
        }
        self.next_out = self.bb.pos;
        let used = self.bb.used();
        self.avail_out -= used;
        self.state = State::End;
    }

    /// `create_hufftables_icf`: the block's codes from its histogram, the dynamic header written
    /// to the bit buffer unless the fixed codes are no longer, and the total length in bits.
    fn create_hufftables_icf(&mut self) -> u64 {
        let saved = self.bb;
        flatten_ll(&mut self.ll_hist);
        if self.ll_hist[256] == 0 {
            self.ll_hist[256] = 1;
        }
        let mut heap = Heap::new();
        let mut bl_count = [0u32; 16];

        let mut ll = [0u32; 286];
        let size = heap.init(
            &self.ll_hist[..286]
                .iter()
                .map(|&h| h as u64)
                .collect::<Vec<_>>(),
        );
        heap.code_lens(size, &mut bl_count, &mut ll, 15);
        let max_ll = set_huff_codes(&mut ll, &bl_count);

        let mut d = [0u32; 30];
        let size = heap.init(&self.d_hist.map(|h| h as u64));
        heap.code_lens(size, &mut bl_count, &mut d, 15);
        let max_d = set_huff_codes(&mut d, &bl_count);
        for (i, code) in d.iter_mut().enumerate() {
            if *code >> 24 != 0 {
                *code |= DIST_EXTRA_BITS[i] << 16;
            }
        }
        // set_huff_codes writes lengths for the whole table; the C's set_dist_huff_codes stops at
        // the last used symbol in the same way.
        self.tables[..30].copy_from_slice(&d);
        self.tables[31..31 + 286].copy_from_slice(&ll);

        let statics = static_tables();
        let length = |code: u32| code >> 24;
        let mut combined = Vec::with_capacity(286 + 30);
        let mut compressed: u64 = 0;
        let mut fixed: u64 = 3;
        for i in 0..=max_ll {
            combined.push(length(ll[i]) as u16);
            let extra = if i > 256 { LEN_EXTRA_BITS[i - 257] } else { 0 };
            compressed += ((length(ll[i]) + extra) * self.ll_hist[i]) as u64;
            fixed += ((length(statics[31 + i]) + extra) * self.ll_hist[i]) as u64;
        }
        for i in 0..=max_d {
            combined.push(length(d[i]) as u16);
            compressed += ((length(d[i]) + DIST_EXTRA_BITS[i]) * self.d_hist[i]) as u64;
            fixed += ((length(statics[i]) + DIST_EXTRA_BITS[i]) * self.d_hist[i]) as u64;
        }
        let eob = self.has_eob_hdr;
        if fixed > compressed {
            let mut counts = [0u64; 19];
            let tokens = rl_encode(&combined, &mut counts);
            create_header(
                &mut self.bb,
                &mut self.mem,
                &mut heap,
                &tokens,
                &counts,
                (max_ll - 256) as u64,
                max_d as u64,
                eob,
            );
            compressed += 8 * self.bb.used() + self.bb.count as u64;
        }
        if fixed <= compressed {
            self.tables = *statics;
            self.bb = saved;
            self.bb.write(&mut self.mem, 0x2 | eob as u64, 3);
            compressed = fixed;
        }
        expand_hufftables_icf(&mut self.tables);
        compressed
    }
}

/// `flatten_ll`: the ICF histogram counts lengths one per length; deflate counts symbols.
fn flatten_ll(s: &mut [u32; 513]) {
    s[265] += s[266];
    s[266] = s[267] + s[268];
    s[267] = s[269] + s[270];
    s[268] = s[271] + s[272];
    s[269] = s[273] + s[274] + s[275] + s[276];
    s[270] = s[277] + s[278] + s[279] + s[280];
    s[271] = s[281] + s[282] + s[283] + s[284];
    s[272] = s[285] + s[286] + s[287] + s[288];
    let mut p = 289;
    for (range, width) in [(273..277, 8), (277..281, 16), (281..285, 32)] {
        for i in range {
            s[i] = s[p..p + width].iter().sum();
            p += width;
        }
    }
    s[284] -= s[512];
    s[285] = s[512];
}

/// The fixed codes as `static_hufftables` holds them before expansion: distance codes first, with
/// their extra-bit counts, then literal/length codes. Lengths in the top byte, codes bit-reversed.
fn static_tables() -> &'static [u32; 544] {
    static T: std::sync::OnceLock<[u32; 544]> = std::sync::OnceLock::new();
    T.get_or_init(|| {
        let mut t = [0u32; 544];
        let mut lens = [0u32; 288];
        lens[..144].fill(8);
        lens[144..256].fill(9);
        lens[256..280].fill(7);
        lens[280..].fill(8);
        let mut bl = [0u32; 16];
        for &l in &lens {
            bl[l as usize] += 1;
        }
        bl[0] = 0;
        let mut next = [0u32; 16];
        for b in 1..16 {
            next[b] = (next[b - 1] + bl[b - 1]) << 1;
        }
        for (i, &l) in lens.iter().enumerate().take(286) {
            t[31 + i] = (l << 24) | reverse(next[l as usize], l);
            next[l as usize] += 1;
        }
        for (i, slot) in t.iter_mut().enumerate().take(30) {
            *slot = (5 << 24) | (DIST_EXTRA_BITS[i] << 16) | reverse(i as u32, 5);
        }
        t
    })
}

fn reverse(code: u32, len: u32) -> u32 {
    (code as u16).reverse_bits() as u32 >> (16 - len)
}

/// `expand_hufftables_icf`: one entry per length from 11 up, so a token's length indexes the
/// table directly with its extra bits already appended.
fn expand_hufftables_icf(t: &mut [u32; 544]) {
    let ll = 31;
    let orig: [u32; 21] = t[ll + 265..ll + 286].try_into().unwrap();
    let mut p = ll + 265;
    let mut i = 0;
    for eb in 1..6 {
        for _ in 0..4 {
            let len = orig[i] >> 24;
            let code = orig[i] & 0xffff;
            i += 1;
            for j in 0..(1u32 << eb) {
                t[p] = ((code | (j << len)) & 0xff_ffff) | ((len + eb) << 24);
                p += 1;
            }
        }
    }
    t[p - 1] = orig[i];
    t[30] = 0;
}

/// `set_huff_codes` (and `set_dist_huff_codes`, which computes the same codes): canonical codes
/// from the lengths, bit-reversed. Returns the last symbol with a code.
fn set_huff_codes(table: &mut [u32], count: &[u32; 16]) -> usize {
    let mut next = [0u16; 16];
    for i in 1..16 {
        next[i] = next[i - 1].wrapping_add(count[i - 1] as u16) << 1;
    }
    let mut max = 0;
    for (i, code) in table.iter_mut().enumerate() {
        let len = *code >> 24;
        if len != 0 {
            *code = (len << 24) | reverse(next[len as usize] as u32, len);
            next[len as usize] = next[len as usize].wrapping_add(1);
            max = i;
        }
    }
    max
}

/// `struct heap_tree` and the functions of `proc_heap_base.c` and `huff_codes.c` that use it. One
/// array is the heap from index 1 up and the tree from `NODE_START` down, as in the C union; a
/// tree node is `child` in the low half of a slot and `depth` in the high half.
struct Heap {
    a: [u64; HEAP_TREE_SIZE],
    /// Sixteen counts in the C, which overruns them into the dead heap below the tree when an
    /// unlimited code is deeper than 15; sized here for the deepest tree instead.
    code_len_count: [u64; MAX_HISTHEAP_SIZE as usize + 1],
}

const MAX_HISTHEAP_SIZE: u64 = 286;
const HEAP_TREE_SIZE: usize = 3 * 286 + 1;
const NODE_START: usize = HEAP_TREE_SIZE - 1;

impl Heap {
    fn new() -> Self {
        Heap {
            a: [0; HEAP_TREE_SIZE],
            code_len_count: [0; MAX_HISTHEAP_SIZE as usize + 1],
        }
    }

    fn init(&mut self, hist: &[u64]) -> usize {
        self.a = [0; HEAP_TREE_SIZE];
        let mut size = 0;
        for (i, &h) in hist.iter().enumerate() {
            if h != 0 {
                size += 1;
                self.a[size] = (h << 16) | i as u64;
            }
        }
        if size < 2 {
            if size == 0 {
                self.a[1] = 1 << 16;
                self.a[2] = (1 << 16) | 1;
            } else if hist[0] == 0 {
                self.a[2] = 1 << 16;
            } else {
                self.a[2] = (1 << 16) | 1;
            }
            size = 2;
        }
        self.a[size + 1] = u64::MAX;
        for i in (1..=size / 2).rev() {
            self.heapify(size, i);
        }
        size
    }

    fn heapify(&mut self, size: usize, mut index: usize) {
        let h = &mut self.a;
        let mut child = 2 * index;
        while child <= size {
            if h[child] > h[child + 1] {
                child += 1;
            }
            if h[index] > h[child] {
                h.swap(index, child);
                index = child;
                child = 2 * index;
            } else {
                break;
            }
        }
    }

    fn build_huff_tree(&mut self, mut size: usize) -> usize {
        let mut node = NODE_START;
        while size > 1 {
            let h1 = self.a[1];
            self.a[1] = self.a[size];
            self.a[size] = u64::MAX;
            size -= 1;
            self.heapify(size, 1);
            let h2 = self.a[1];
            self.a[1] = ((h1 + h2) & !0xffff) | node as u64;
            self.heapify(size, 1);
            self.a[node] = h1 & 0xffff;
            self.a[node - 1] = h2 & 0xffff;
            node -= 2;
        }
        self.a[node] = self.a[1] & 0xffff;
        node
    }

    fn child(&self, i: usize) -> u64 {
        self.a[i] & 0xffff_ffff
    }

    fn depth(&self, i: usize) -> u64 {
        self.a[i] >> 32
    }

    fn set_depth(&mut self, i: usize, depth: u64) {
        self.a[i] = (self.a[i] & 0xffff_ffff) | (depth << 32);
    }

    /// `fix_code_lens` and `gen_huff_code_lens`: depths from the tree, lengths over `max` pushed
    /// back under it, then each symbol's length into the top byte of its entry in `codes`.
    fn code_lens(&mut self, size: usize, bl_count: &mut [u32; 16], codes: &mut [u32], max: usize) {
        let root = self.build_huff_tree(size);
        let mut code_len = 0usize;
        let mut j = root;
        for i in root..=NODE_START {
            let child = self.child(i);
            if child > MAX_HISTHEAP_SIZE {
                let depth = 1 + self.depth(i);
                self.set_depth(child as usize, depth);
                self.set_depth(child as usize - 1, depth);
            } else {
                self.a[j] = self.a[i];
                j += 1;
                let depth = self.depth(i) as usize;
                while code_len < depth {
                    code_len += 1;
                    self.code_len_count[code_len] = 0;
                }
                self.code_len_count[depth] += 1;
            }
        }
        let limited = code_len > max;
        let c = &mut self.code_len_count;
        while code_len > max {
            let mut i = max - 1;
            while i != 0 && c[i] == 0 {
                i -= 1;
            }
            assert!(i != 0, "a Huffman code that cannot be length-limited");
            c[i] -= 1;
            c[i + 1] += 2;
            c[code_len - 1] += 1;
            c[code_len] -= 2;
            if c[code_len] == 0 {
                code_len -= 1;
            }
        }
        bl_count[0] = 0;
        for (i, slot) in bl_count.iter_mut().enumerate().skip(1) {
            *slot = if i <= code_len { c[i] as u32 } else { 0 };
        }
        if limited {
            // Leaves in tree order, shortest lengths first.
            let mut k = 1;
            for i in root..j {
                while self.code_len_count[k] == 0 {
                    k += 1;
                }
                self.set_depth(i, k as u64);
                self.code_len_count[k] -= 1;
            }
        }
        codes.fill(0);
        for i in root..j {
            codes[self.child(i) as usize] = (self.depth(i) as u32) << 24;
        }
    }
}

/// `rl_encode` and `write_rl`: the code lengths as code-length symbols 0 to 18 with their extra
/// bits, counting each symbol.
fn rl_encode(codes: &[u16], counts: &mut [u64; 19]) -> Vec<(u8, u8)> {
    let mut out = Vec::new();
    let mut write_rl = |last: u16, mut run: u32, out: &mut Vec<(u8, u8)>| {
        if last == 0 {
            while run > 138 {
                out.push((18, 127));
                run -= 138;
                counts[18] += 1;
            }
            if run > 10 {
                out.push((18, (run - 11) as u8));
                counts[18] += 1;
            } else if run > 2 {
                out.push((17, (run - 3) as u8));
                counts[17] += 1;
            } else {
                for _ in 0..run {
                    out.push((0, 0));
                    counts[0] += 1;
                }
            }
        } else {
            out.push((last as u8, 0));
            counts[last as usize] += 1;
            run -= 1;
            if run != 0 {
                while run > 6 {
                    out.push((16, 3));
                    run -= 6;
                    counts[16] += 1;
                }
                if run <= 2 {
                    for _ in 0..run {
                        out.push((last as u8, 0));
                        counts[last as usize] += 1;
                    }
                } else {
                    out.push((16, (run - 3) as u8));
                    counts[16] += 1;
                }
            }
        }
    };
    let mut last = codes[0];
    let mut run = 1;
    for &len in &codes[1..] {
        if len == last {
            run += 1;
            continue;
        }
        write_rl(last, run, &mut out);
        last = len;
        run = 1;
    }
    write_rl(last, run, &mut out);
    out
}

/// `create_header` and `create_huffman_header`: the dynamic block header.
#[allow(clippy::too_many_arguments)]
fn create_header(
    bb: &mut BitBuf,
    mem: &mut [u8],
    heap: &mut Heap,
    tokens: &[(u8, u8)],
    counts: &[u64; 19],
    hlit: u64,
    hdist: u64,
    eob: bool,
) {
    let mut lookup = [0u32; 19];
    let mut bl_count = [0u32; 16];
    let size = heap.init(counts);
    heap.code_lens(size, &mut bl_count, &mut lookup, 7);
    set_huff_codes(&mut lookup, &bl_count);
    let length = |sym: usize| (lookup[sym] >> 24) as u64;
    let mut i = 18;
    while i > 3 && length(CODE_LENGTH_ORDER[i]) == 0 {
        i -= 1;
    }
    let hclen = i as u64 - 3;
    let data = (if eob { 5 } else { 4 })
        | (hlit << 3)
        | (hdist << 8)
        | (hclen << 13)
        | (length(CODE_LENGTH_ORDER[0]) << 17);
    bb.write(mem, data, 20);
    let mut data = 0u64;
    for i in (1..=hclen as usize + 3).rev() {
        data = (data << 3) | length(CODE_LENGTH_ORDER[i]);
    }
    bb.write(mem, data, (hclen as u32 + 3) * 3);
    for &(code, extra) in tokens {
        let entry = lookup[code as usize];
        bb.write(mem, (entry & 0xffff) as u64, entry >> 24);
        if code > 15 {
            bb.write(mem, extra as u64, [2, 3, 7][code as usize - 16]);
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::igzip_canary::{EXPECTED, INPUT};

    /// The 2 KB GKL answered in the pinned container, which used to guard the ISA-L binding.
    #[test]
    fn the_canary_is_gkls_bytes() {
        assert_eq!(super::deflate(&INPUT, 1), EXPECTED);
    }

    /// With the `isal` feature, the port against the library itself on many more inputs. Only a
    /// build whose ISA-L reproduces GKL can say anything, which is x86-64 with nasm; elsewhere
    /// the library is its own readable C and the comparison would be against the wrong answer.
    #[cfg(feature = "isal")]
    #[test]
    fn agrees_with_linked_isal() {
        if !crate::igzip_isal::usable() {
            #[cfg(target_arch = "x86_64")]
            assert!(
                !std::arch::is_x86_feature_detected!("sse4.2"),
                "x86-64 with SSE4.2 and ISA-L does not reproduce GKL: was it built without nasm?"
            );
            return;
        }
        let mut seed = 0x9e37_79b9_7f4a_7c15u64;
        let mut next = move || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        let mut compared = 0;
        for round in 0..400 {
            let len = match round % 4 {
                0 => (next() % 400) as usize,
                1 => (next() % 70_000) as usize,
                2 => 65_498,
                _ => (next() % 300_000) as usize,
            };
            // Mixtures of noise, a small alphabet, zero runs and copies of earlier spans.
            let mut data = Vec::with_capacity(len);
            while data.len() < len {
                let span = (next() % 5000) as usize + 1;
                match next() % 4 {
                    0 => data.extend((0..span).map(|_| next() as u8)),
                    1 => data.extend((0..span).map(|_| b"ACGT"[(next() & 3) as usize])),
                    2 => data.extend(std::iter::repeat_n(0u8, span)),
                    _ if !data.is_empty() => {
                        let from = (next() as usize) % data.len();
                        let end = (from + span).min(data.len());
                        data.extend_from_within(from..end);
                    }
                    _ => {}
                }
            }
            data.truncate(len);
            for level in 1..=2 {
                assert!(
                    super::deflate(&data, level) == crate::igzip_isal::deflate(&data, level),
                    "round {round}, {len} bytes, level {level}"
                );
                compared += 1;
            }
        }
        assert_eq!(compared, 800);
    }
}
