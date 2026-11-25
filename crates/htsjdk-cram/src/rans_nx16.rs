//! The CRAM 3.1 rANS Nx16 codec, block compression method 5.
//!
//! Ported from `htsjdk.samtools.cram.compression.rans.RANSNx16Decode` and `RANSNx16Encode` at
//! htsjdk 5.0.0, which is the htsjdk Picard 3.5.0 bundles.
//!
//! One leading byte of format flags says which of five transformations wrap the entropy coder:
//! order 0 or 1, four or thirty-two interleaved states, a byte stripe, run-length coding, bit
//! packing, or no entropy coding at all (`CAT`). The uncompressed size follows as a uint7 unless
//! `NOSZ` says the caller knows it, which only a stripe's sub-streams do.

use std::fmt;

pub const ORDER_FLAG_MASK: u8 = 0x01;
pub const N32_FLAG_MASK: u8 = 0x04;
pub const STRIPE_FLAG_MASK: u8 = 0x08;
pub const NOSZ_FLAG_MASK: u8 = 0x10;
pub const CAT_FLAG_MASK: u8 = 0x20;
pub const RLE_FLAG_MASK: u8 = 0x40;
pub const PACK_FLAG_MASK: u8 = 0x80;

const NUMBER_OF_SYMBOLS: usize = 256;
const TOTAL_FREQ_SHIFT: u32 = 12;
const LOWER_BOUND: u64 = 1 << 15;

/// What decoding refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum NxError {
    /// `CRAMException`, with the reference's message.
    Cram(String),
    /// The stream ended before the decoder did: an `ArrayIndexOutOfBoundsException` there.
    Truncated,
}

impl fmt::Display for NxError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            NxError::Cram(message) => write!(f, "{message}"),
            NxError::Truncated => write!(f, "the rANS Nx16 stream is truncated"),
        }
    }
}

/// A cursor over the compressed bytes.
pub(crate) struct Input<'a> {
    pub(crate) bytes: &'a [u8],
    pub(crate) pos: usize,
}

impl<'a> Input<'a> {
    pub(crate) fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, pos: 0 }
    }

    #[inline(always)]
    pub(crate) fn byte(&mut self) -> Result<u8, NxError> {
        let b = *self.bytes.get(self.pos).ok_or(NxError::Truncated)?;
        self.pos += 1;
        Ok(b)
    }

    pub(crate) fn take(&mut self, n: usize) -> Result<&'a [u8], NxError> {
        let s = self
            .bytes
            .get(self.pos..self.pos + n)
            .ok_or(NxError::Truncated)?;
        self.pos += n;
        Ok(s)
    }

    /// `CompressionUtils.readUint7`.
    pub(crate) fn uint7(&mut self) -> Result<i32, NxError> {
        let mut value: i32 = 0;
        loop {
            let c = self.byte()?;
            value = (value << 7) | i32::from(c & 0x7f);
            if c & 0x80 == 0 {
                return Ok(value);
            }
        }
    }

    fn u32_le(&mut self) -> Result<u64, NxError> {
        let b = self.take(4)?;
        Ok(u64::from(u32::from_le_bytes([b[0], b[1], b[2], b[3]])))
    }
}

/// `CompressionUtils.writeUint7`.
pub(crate) fn write_uint7(value: i32, out: &mut Vec<u8>) {
    let mut s = 0;
    let mut x = value;
    loop {
        s += 7;
        x >>= 7;
        if x <= 0 {
            break;
        }
    }
    loop {
        s -= 7;
        let more = if s > 0 { 0x80 } else { 0 };
        out.push((((value >> s) & 0x7f) | more) as u8);
        if s <= 0 {
            break;
        }
    }
}

/// `RANSNx16Decode.uncompress`.
pub fn uncompress(input: &[u8]) -> Result<Vec<u8>, NxError> {
    if input.is_empty() {
        return Ok(Vec::new());
    }
    let mut cursor = Input::new(input);
    uncompress_internal(&mut cursor, 0)
}

fn uncompress_internal(input: &mut Input<'_>, out_size: usize) -> Result<Vec<u8>, NxError> {
    if input.pos >= input.bytes.len() {
        return Ok(Vec::new());
    }
    let flags = input.byte()?;
    let mut size = if flags & NOSZ_FLAG_MASK != 0 {
        out_size
    } else {
        input.uint7()? as usize
    };

    if flags & STRIPE_FLAG_MASK != 0 {
        return decode_stripe(input, size);
    }

    let mut pack: Option<(Vec<u8>, usize, usize)> = None;
    if flags & PACK_FLAG_MASK != 0 {
        let pack_length = size;
        let symbols = input.byte()? as usize;
        if symbols <= 16 && symbols != 0 {
            let table = input.take(symbols)?.to_vec();
            size = input.uint7()? as usize;
            pack = Some((table, symbols, pack_length));
        } else {
            return Err(NxError::Cram(format!(
                "Bit Packing is not permitted when number of distinct symbols is greater than 16 or equal to 0. \
                 Number of distinct symbols: {symbols}"
            )));
        }
    }

    let mut rle: Option<(Vec<u8>, [bool; NUMBER_OF_SYMBOLS], usize, usize)> = None;
    if flags & RLE_FLAG_MASK != 0 {
        let meta_length = input.uint7()? as usize;
        let rle_output_length = size;
        size = input.uint7()? as usize;
        let (meta, symbols, meta_pos) = decode_rle_meta(input, meta_length, flags)?;
        rle = Some((meta, symbols, meta_pos, rle_output_length));
    }

    let mut out = if flags & CAT_FLAG_MASK != 0 {
        input.take(size)?.to_vec()
    } else {
        if size == 0 {
            return Err(NxError::Cram(
                "Unexpected uncompressed size of 0 in RANSNx16 stream".to_string(),
            ));
        }
        let mut out = vec![0u8; size];
        let ways = if flags & N32_FLAG_MASK != 0 { 32 } else { 4 };
        if flags & ORDER_FLAG_MASK == 0 {
            order0(input, &mut out, ways)?;
        } else {
            order1(input, &mut out, ways)?;
        }
        out
    };

    if let Some((meta, symbols, meta_pos, length)) = rle {
        out = decode_rle(&out, &symbols, &meta, meta_pos, length)?;
    }
    if let Some((table, symbols, length)) = pack {
        out = decode_pack(&out, &table, symbols, length)?;
    }
    Ok(out)
}

/// The decoding symbols and the reverse lookup of one context.
struct Table {
    start: [u32; NUMBER_OF_SYMBOLS],
    freq: [u32; NUMBER_OF_SYMBOLS],
    reverse: Vec<u8>,
}

impl Table {
    fn new() -> Self {
        Self {
            start: [0; NUMBER_OF_SYMBOLS],
            freq: [0; NUMBER_OF_SYMBOLS],
            reverse: vec![0u8; 1 << TOTAL_FREQ_SHIFT],
        }
    }

    /// `RANSDecodingSymbol.advanceSymbolStep`.
    #[inline]
    fn advance(&self, symbol: u8, r: u64, shift: u32) -> u64 {
        let mask = (1u64 << shift) - 1;
        u64::from(self.freq[symbol as usize]) * (r >> shift) + (r & mask)
            - u64::from(self.start[symbol as usize])
    }
}

/// `Utils.RANSDecodeRenormalizeNx16`.
#[inline]
fn renormalize(r: u64, input: &mut Input<'_>) -> Result<u64, NxError> {
    if r < LOWER_BOUND {
        let b = input.take(2)?;
        return Ok((r << 16) | u64::from(b[0]) | (u64::from(b[1]) << 8));
    }
    Ok(r)
}

/// `readAlphabet`.
fn read_alphabet(input: &mut Input<'_>) -> Result<[bool; NUMBER_OF_SYMBOLS], NxError> {
    let mut alphabet = [false; NUMBER_OF_SYMBOLS];
    let mut rle = 0u32;
    let mut symbol = u32::from(input.byte()?);
    let mut last = symbol;
    loop {
        alphabet[(symbol & 0xFF) as usize] = true;
        if rle != 0 {
            rle -= 1;
            symbol += 1;
        } else {
            symbol = u32::from(input.byte()?);
            if symbol == last + 1 {
                rle = u32::from(input.byte()?);
            }
        }
        last = symbol;
        if symbol == 0 {
            break;
        }
        if symbol > 0xFF {
            return Err(NxError::Truncated);
        }
    }
    Ok(alphabet)
}

/// `Utils.normaliseFrequenciesOrder0Shift`.
fn normalise_shift(freq: &mut [u32; NUMBER_OF_SYMBOLS], bits: u32) {
    let mut total: u64 = freq.iter().map(|&f| u64::from(f)).sum();
    if total == 0 || total == 1 << bits {
        return;
    }
    let mut shift = 0;
    while total < 1 << bits {
        total *= 2;
        shift += 1;
    }
    for f in freq.iter_mut() {
        if *f != 0 {
            *f <<= shift;
        }
    }
}

fn fill_table(table: &mut Table, freq: &[u32; NUMBER_OF_SYMBOLS], only: Option<&[bool; 256]>) {
    let mut cumulative = 0u32;
    for j in 0..NUMBER_OF_SYMBOLS {
        if only.is_some_and(|a| !a[j]) {
            continue;
        }
        table.start[j] = cumulative;
        table.freq[j] = freq[j];
        let end = (cumulative + freq[j]).min(table.reverse.len() as u32);
        for slot in &mut table.reverse[cumulative.min(end) as usize..end as usize] {
            *slot = j as u8;
        }
        cumulative += freq[j];
    }
}

/// One slot of a context's reverse lookup, packed so that the hot loop reads one word per symbol:
/// the symbol, its frequency, and the slot's distance from the symbol's start.
///
/// With it, `RANSDecodingSymbol.advanceSymbolStep`,
/// ```java
/// freq * (r >> scaleBits) + (r & mask) - start
/// ```
/// becomes `freq * (r >> shift) + bias`, the same value: `bias` is `(r & mask) - start` worked out
/// once per slot instead of once per symbol. That is the layout htscodecs' decoders use.
///
/// There are two packings. `u64` holds any frequency a damaged table can carry. `u32` holds the
/// tables a conforming encoder writes, where every context's frequencies add up to exactly
/// `1 << shift` ([`conforming`]): each slot then belongs to a symbol of frequency 1 to 4096, which
/// fits 12 bits as `freq - 1`, and the table is half the size in cache.
trait Slot: Copy {
    fn pack(freq: u32, bias: u32, symbol: u8) -> Self;
    fn symbol(self) -> u8;
    /// The state advanced past this slot's symbol.
    fn advance(self, x: u64, shift: u32) -> u64;
}

impl Slot for u64 {
    fn pack(freq: u32, bias: u32, symbol: u8) -> Self {
        (u64::from(freq) << 32) | (u64::from(bias) << 8) | u64::from(symbol)
    }
    #[inline(always)]
    fn symbol(self) -> u8 {
        self as u8
    }
    #[inline(always)]
    fn advance(self, x: u64, shift: u32) -> u64 {
        // Wrapping, as the reference's `long` arithmetic is, for a damaged table's frequencies.
        (self >> 32)
            .wrapping_mul(x >> shift)
            .wrapping_add((self >> 8) & 0xFF_FFFF)
    }
}

impl Slot for u32 {
    fn pack(freq: u32, bias: u32, symbol: u8) -> Self {
        (freq.wrapping_sub(1) << 20) | ((bias & 0xFFF) << 8) | u32::from(symbol)
    }
    #[inline(always)]
    fn symbol(self) -> u8 {
        self as u8
    }
    #[inline(always)]
    fn advance(self, x: u64, shift: u32) -> u64 {
        u64::from((self >> 20) + 1) * (x >> shift) + u64::from((self >> 8) & 0xFFF)
    }
}

/// Whether a context's normalised frequencies add up to exactly `1 << shift`, as every table a
/// conforming encoder writes does. Then no slot is left to `Table`'s defaults, and every
/// frequency fits the `u32` packing.
fn conforming(freq: &[u32; NUMBER_OF_SYMBOLS], shift: u32) -> bool {
    freq.iter().map(|&f| u64::from(f)).sum::<u64>() == 1 << shift
}

/// Fill one context's packed lookup the way [`fill_table`] fills a [`Table`], so that every slot
/// decodes to the symbol, frequency and start the reference would read.
///
/// A slot no symbol covers keeps `Table`'s defaults, symbol 0 with `freq[0]` and start 0, and the
/// sums wrap as the reference's `int` arithmetic does. Both only matter for a damaged table.
fn fill_packed<S: Slot>(
    table: &mut [S],
    freq: &[u32; NUMBER_OF_SYMBOLS],
    only: Option<&[bool; 256]>,
) {
    for (m, s) in table.iter_mut().enumerate() {
        *s = S::pack(freq[0], m as u32, 0);
    }
    let len = table.len() as u32;
    let mut cumulative = 0u32;
    for j in 0..NUMBER_OF_SYMBOLS {
        if only.is_some_and(|a| !a[j]) {
            continue;
        }
        let end = cumulative.wrapping_add(freq[j]).min(len);
        for m in cumulative.min(end)..end {
            table[m as usize] = S::pack(freq[j], m - cumulative, j as u8);
        }
        cumulative = cumulative.wrapping_add(freq[j]);
    }
}

/// `Utils.RANSDecodeRenormalizeNx16` without a branch, for a round that cannot run out of input.
///
/// The two bytes are read whether or not the state needs them, and the cursor moves only if it
/// does. The caller has checked that `2 * ways` bytes remain, so the read is always in bounds.
#[inline(always)]
fn renormalize_unchecked(x: u64, bytes: &[u8], pos: &mut usize) -> u64 {
    let refill = (x << 16) | u64::from(u16::from_le_bytes([bytes[*pos], bytes[*pos + 1]]));
    let need = x < LOWER_BOUND;
    *pos += 2 * usize::from(need);
    if need {
        refill
    } else {
        x
    }
}

/// `uncompressOrder0WayN`.
///
/// The reference's loop, with the lookup packed ([`Slot`]) and the states in a fixed-size array.
/// Each round of `ways` symbols first checks that the input holds the `2 * ways` bytes the round
/// could read; while it does, nothing can run out, so the round skips the per-byte checks. Rounds
/// past that point go through [`renormalize`], which refuses exactly where the reference does.
pub(crate) fn order0(input: &mut Input<'_>, out: &mut [u8], ways: usize) -> Result<(), NxError> {
    #[cfg(test)]
    if reference::ON.with(|on| on.get()) {
        return order0_reference(input, out, ways);
    }
    let alphabet = read_alphabet(input)?;
    let mut freq = [0u32; NUMBER_OF_SYMBOLS];
    for j in 0..NUMBER_OF_SYMBOLS {
        if alphabet[j] {
            freq[j] = input.uint7()? as u32;
        }
    }
    normalise_shift(&mut freq, TOTAL_FREQ_SHIFT);
    // `<4, u32>` reads "four states, slots packed in 32 bits".
    match (ways == 4, conforming(&freq, TOTAL_FREQ_SHIFT)) {
        (true, true) => order0_ways::<4, u32>(input, out, &freq, &alphabet),
        (true, false) => order0_ways::<4, u64>(input, out, &freq, &alphabet),
        (false, true) => order0_ways::<32, u32>(input, out, &freq, &alphabet),
        (false, false) => order0_ways::<32, u64>(input, out, &freq, &alphabet),
    }
}

fn order0_ways<const W: usize, S: Slot>(
    input: &mut Input<'_>,
    out: &mut [u8],
    freq: &[u32; NUMBER_OF_SYMBOLS],
    alphabet: &[bool; NUMBER_OF_SYMBOLS],
) -> Result<(), NxError> {
    let mut table = vec![S::pack(0, 0, 0); 1 << TOTAL_FREQ_SHIFT];
    fill_packed(&mut table, freq, Some(alphabet));

    let mut x = [0u64; W];
    for r in x.iter_mut() {
        *r = input.u32_le()?;
    }
    let size = out.len();
    let end = size / W * W;
    let mask = (1u64 << TOTAL_FREQ_SHIFT) - 1;
    let bytes = input.bytes;

    for chunk in out[..end].chunks_exact_mut(W) {
        if input.pos + 2 * W <= bytes.len() {
            let mut pos = input.pos;
            for r in 0..W {
                let entry = table[(x[r] & mask) as usize];
                chunk[r] = entry.symbol();
                x[r] =
                    renormalize_unchecked(entry.advance(x[r], TOTAL_FREQ_SHIFT), bytes, &mut pos);
            }
            input.pos = pos;
        } else {
            for r in 0..W {
                let entry = table[(x[r] & mask) as usize];
                chunk[r] = entry.symbol();
                x[r] = renormalize(entry.advance(x[r], TOTAL_FREQ_SHIFT), input)?;
            }
        }
    }
    for (r, slot) in out[end..].iter_mut().enumerate() {
        let entry = table[(x[r] & mask) as usize];
        x[r] = renormalize(entry.advance(x[r], TOTAL_FREQ_SHIFT), input)?;
        *slot = entry.symbol();
    }
    Ok(())
}

/// `uncompressOrder1WayN`.
///
/// Each state decodes its own stretch of the output, `interleave` bytes long, with the previous
/// symbol of that stretch as the context; the last state then decodes what is left. As in
/// [`order0`], the lookups are packed and a round that cannot run out of input skips the checks.
///
/// The packed tables of the contexts the alphabet names sit one after the other in one buffer,
/// and a context the alphabet leaves out points at a shared table of `Table::new()`'s defaults.
/// That keeps memory at the few contexts a block uses rather than 256 full tables.
///
/// A shift past 12 is left to [`order1_reference`]: the reference sizes an absent context's
/// table at 4096 slots whatever the shift, and refuses a state that indexes past it.
fn order1(input: &mut Input<'_>, out: &mut [u8], ways: usize) -> Result<(), NxError> {
    #[cfg(test)]
    if reference::ON.with(|on| on.get()) {
        return order1_reference(input, out, ways);
    }
    let start = input.pos;
    let first = input.byte()?;
    let shift = u32::from(first >> 4);
    if shift > TOTAL_FREQ_SHIFT {
        input.pos = start;
        return order1_reference(input, out, ways);
    }
    let mut rows: Vec<(usize, [u32; NUMBER_OF_SYMBOLS])> = Vec::new();
    let mut add = |context: usize, freq: &[u32; NUMBER_OF_SYMBOLS]| rows.push((context, *freq));
    if first & 1 != 0 {
        let uncompressed_length = input.uint7()? as usize;
        let compressed_length = input.uint7()? as usize;
        let compressed = input.take(compressed_length)?;
        let mut table_bytes = vec![0u8; uncompressed_length];
        // `new RANSNx16Params(~(ORDER | N32))`: order 0, four ways.
        order0(&mut Input::new(compressed), &mut table_bytes, 4)?;
        read_order1_freqs(&mut Input::new(&table_bytes), shift, &mut add)?;
    } else {
        read_order1_freqs(input, shift, &mut add)?;
    }
    // Every state starts in context 0, so the `u32` packing also needs context 0 to have a row:
    // then every symbol decoded is one the alphabet names, and the shared table is never read.
    let packs = rows.first().is_some_and(|&(c, _)| c == 0)
        && rows.iter().all(|(_, freq)| conforming(freq, shift));
    match (ways == 4, packs) {
        (true, true) => order1_ways::<4, u32>(input, out, shift, &rows),
        (true, false) => order1_ways::<4, u64>(input, out, shift, &rows),
        (false, true) => order1_ways::<32, u32>(input, out, shift, &rows),
        (false, false) => order1_ways::<32, u64>(input, out, shift, &rows),
    }
}

fn order1_ways<const W: usize, S: Slot + 'static>(
    input: &mut Input<'_>,
    out: &mut [u8],
    shift: u32,
    rows: &[(usize, [u32; NUMBER_OF_SYMBOLS])],
) -> Result<(), NxError> {
    let slots = 1usize << shift.max(1);
    // The table every absent context shares sits first, at offset 0.
    let mut tables = vec![S::pack(0, 0, 0); slots * (rows.len() + 1)];
    fill_packed(&mut tables[..slots], &[0; NUMBER_OF_SYMBOLS], None);
    let mut offset = [0usize; NUMBER_OF_SYMBOLS];
    for (k, (context, freq)) in rows.iter().enumerate() {
        let at = slots * (k + 1);
        fill_packed(&mut tables[at..at + slots], freq, None);
        offset[*context] = at;
    }

    let size = out.len();
    let interleave = size / W;
    let mut x = [0u64; W];
    for r in x.iter_mut() {
        *r = input.u32_le()?;
    }
    let mut base = [offset[0]; W];
    let mask = (1u64 << shift) - 1;
    let bytes = input.bytes;

    let mut i = 0;
    #[cfg(target_arch = "x86_64")]
    if let (Some(tables), Some(x), Some(base)) = (
        (&tables as &dyn std::any::Any).downcast_ref::<Vec<u32>>(),
        (&mut x as &mut dyn std::any::Any).downcast_mut::<[u64; 32]>(),
        (&mut base as &mut dyn std::any::Any).downcast_mut::<[usize; 32]>(),
    ) {
        if std::arch::is_x86_feature_detected!("avx2") {
            // SAFETY: the CPU has AVX2, and `tables` and `offset` are the conforming tables the
            // function's contract asks for.
            i = unsafe {
                avx2::order1_rounds(
                    tables,
                    &offset,
                    shift,
                    x,
                    base,
                    bytes,
                    &mut input.pos,
                    out,
                    interleave,
                )
            };
        }
    }

    while i < interleave {
        if input.pos + 2 * W <= bytes.len() {
            let mut pos = input.pos;
            for r in 0..W {
                let entry = tables[base[r] + (x[r] & mask) as usize];
                let symbol = entry.symbol();
                out[r * interleave + i] = symbol;
                x[r] = renormalize_unchecked(entry.advance(x[r], shift), bytes, &mut pos);
                base[r] = offset[symbol as usize];
            }
            input.pos = pos;
        } else {
            for r in 0..W {
                let entry = tables[base[r] + (x[r] & mask) as usize];
                let symbol = entry.symbol();
                out[r * interleave + i] = symbol;
                x[r] = renormalize(entry.advance(x[r], shift), input)?;
                base[r] = offset[symbol as usize];
            }
        }
        i += 1;
    }
    let last = W - 1;
    for slot in out[W * interleave..].iter_mut() {
        let entry = tables[base[last] + (x[last] & mask) as usize];
        let symbol = entry.symbol();
        *slot = symbol;
        x[last] = renormalize(entry.advance(x[last], shift), input)?;
        base[last] = offset[symbol as usize];
    }
    Ok(())
}

/// The 32-state order-1 loop on AVX2, eight states to a vector.
///
/// Ported from htscodecs 1.6 `rans_uncompress_O1_32x16_avx2` (`rANS_static32x16pr_avx2.c`),
/// written over this file's tables. A decoder only, which decision 0052 allows: its output is
/// the scalar loop's.
///
/// It is the scalar loop's arithmetic lane by lane, so the two can hand over between rounds. The
/// states fit 32 bits because the tables conform ([`conforming`]): a state below `2^32` advances
/// to below `freq * 2^(32 - shift) <= 2^32`, and a refilled one is below `2^31`. Refills take the
/// next 16-bit words in state order, as the scalar loop does: a lane that needs one gets the word
/// whose index is the number of lanes before it that need one too.
#[cfg(target_arch = "x86_64")]
mod avx2 {
    use std::arch::x86_64::*;

    /// For each mask of the eight lanes that need a refill, the word each of them takes.
    static REFILL: [[u32; 8]; 256] = {
        let mut table = [[0u32; 8]; 256];
        let mut mask = 0;
        while mask < 256 {
            let (mut lane, mut taken) = (0, 0);
            while lane < 8 {
                if mask >> lane & 1 == 1 {
                    table[mask][lane] = taken;
                    taken += 1;
                }
                lane += 1;
            }
            mask += 1;
        }
        table
    };

    /// Rounds `0..interleave` of the 32 states, for as long as a round cannot run out of input
    /// (`pos + 64 <= bytes.len()`), updating the states, their context offsets and the cursor as
    /// the scalar loop does. Returns the round it stopped before.
    ///
    /// # Safety
    ///
    /// The CPU has AVX2; `tables` are conforming `u32` slots; `offset` points every context into
    /// `tables` with room for `1 << shift` slots.
    #[target_feature(enable = "avx2")]
    #[allow(clippy::too_many_arguments)]
    pub(super) unsafe fn order1_rounds(
        tables: &[u32],
        offset: &[usize; 256],
        shift: u32,
        x: &mut [u64; 32],
        base: &mut [usize; 32],
        bytes: &[u8],
        pos: &mut usize,
        out: &mut [u8],
        interleave: usize,
    ) -> usize {
        let offsets: [i32; 256] = std::array::from_fn(|c| offset[c] as i32);
        let slot_mask = _mm256_set1_epi32(((1u32 << shift) - 1) as i32);
        let count = _mm_cvtsi32_si128(shift as i32);
        let low12 = _mm256_set1_epi32(0xFFF);
        let low8 = _mm256_set1_epi32(0xFF);
        let one = _mm256_set1_epi32(1);
        let below = _mm256_set1_epi32((1 << 15) - 1);
        let lanes = |values: [i32; 8]| _mm256_loadu_si256(values.as_ptr().cast());
        let mut s: [__m256i; 4] =
            std::array::from_fn(|v| lanes(std::array::from_fn(|k| x[8 * v + k] as u32 as i32)));
        let mut b: [__m256i; 4] =
            std::array::from_fn(|v| lanes(std::array::from_fn(|k| base[8 * v + k] as i32)));
        let mut symbols = [0u32; 32];
        let mut i = 0;
        while i < interleave && *pos + 64 <= bytes.len() {
            for v in 0..4 {
                let slot = _mm256_add_epi32(b[v], _mm256_and_si256(s[v], slot_mask));
                let entry = _mm256_i32gather_epi32::<4>(tables.as_ptr().cast(), slot);
                let freq = _mm256_add_epi32(_mm256_srli_epi32::<20>(entry), one);
                let bias = _mm256_and_si256(_mm256_srli_epi32::<8>(entry), low12);
                let symbol = _mm256_and_si256(entry, low8);
                let state = _mm256_add_epi32(
                    _mm256_mullo_epi32(freq, _mm256_srl_epi32(s[v], count)),
                    bias,
                );
                // `state < 2^15`, unsigned: the minimum with 2^15 - 1 is the state itself.
                let need = _mm256_cmpeq_epi32(_mm256_min_epu32(state, below), state);
                let mask = _mm256_movemask_ps(_mm256_castsi256_ps(need)) as usize;
                let words = _mm256_cvtepu16_epi32(_mm_loadu_si128(bytes.as_ptr().add(*pos).cast()));
                let order = _mm256_loadu_si256(REFILL[mask].as_ptr().cast());
                let refill = _mm256_or_si256(
                    _mm256_slli_epi32::<16>(state),
                    _mm256_permutevar8x32_epi32(words, order),
                );
                s[v] = _mm256_blendv_epi8(state, refill, need);
                *pos += 2 * mask.count_ones() as usize;
                b[v] = _mm256_i32gather_epi32::<4>(offsets.as_ptr(), symbol);
                _mm256_storeu_si256(symbols.as_mut_ptr().add(8 * v).cast(), symbol);
            }
            for (r, &symbol) in symbols.iter().enumerate() {
                out[r * interleave + i] = symbol as u8;
            }
            i += 1;
        }
        let mut lane = [0u32; 8];
        for v in 0..4 {
            _mm256_storeu_si256(lane.as_mut_ptr().cast(), s[v]);
            for k in 0..8 {
                x[8 * v + k] = u64::from(lane[k]);
            }
            _mm256_storeu_si256(lane.as_mut_ptr().cast(), b[v]);
            for k in 0..8 {
                base[8 * v + k] = lane[k] as usize;
            }
        }
        i
    }
}

/// `uncompressOrder0WayN` as first ported, symbol by symbol: what [`order0`] is checked against.
#[cfg(test)]
fn order0_reference(input: &mut Input<'_>, out: &mut [u8], ways: usize) -> Result<(), NxError> {
    let alphabet = read_alphabet(input)?;
    let mut freq = [0u32; NUMBER_OF_SYMBOLS];
    for j in 0..NUMBER_OF_SYMBOLS {
        if alphabet[j] {
            freq[j] = input.uint7()? as u32;
        }
    }
    normalise_shift(&mut freq, TOTAL_FREQ_SHIFT);
    let mut table = Table::new();
    fill_table(&mut table, &freq, Some(&alphabet));

    let mut rans = vec![0u64; ways];
    for r in rans.iter_mut() {
        *r = input.u32_le()?;
    }
    let size = out.len();
    let interleave = if ways == 4 { size >> 2 } else { size >> 5 };
    let mut remaining = size - interleave * ways;
    let end = size - remaining;
    let mask = (1u64 << TOTAL_FREQ_SHIFT) - 1;

    let mut i = 0;
    while i < end {
        for r in 0..ways {
            let symbol = table.reverse[(rans[r] & mask) as usize];
            out[i + r] = symbol;
            rans[r] = table.advance(symbol, rans[r], TOTAL_FREQ_SHIFT);
            rans[r] = renormalize(rans[r], input)?;
        }
        i += ways;
    }
    let mut r = 0;
    let mut at = end;
    while remaining > 0 {
        let symbol = table.reverse[(rans[r] & mask) as usize];
        rans[r] = table.advance(symbol, rans[r], TOTAL_FREQ_SHIFT);
        rans[r] = renormalize(rans[r], input)?;
        out[at] = symbol;
        at += 1;
        remaining -= 1;
        r += 1;
    }
    Ok(())
}

/// `uncompressOrder1WayN` as first ported, symbol by symbol. [`order1`] falls back to it for a
/// frequency shift past 12, which no conforming encoder writes, and is checked against it.
fn order1_reference(input: &mut Input<'_>, out: &mut [u8], ways: usize) -> Result<(), NxError> {
    let first = input.byte()?;
    let compressed_table = first & 1 != 0;
    let shift = u32::from(first >> 4);

    let tables = if compressed_table {
        let uncompressed_length = input.uint7()? as usize;
        let compressed_length = input.uint7()? as usize;
        let compressed = input.take(compressed_length)?;
        let mut table_bytes = vec![0u8; uncompressed_length];
        // `new RANSNx16Params(~(ORDER | N32))`: order 0, four ways.
        order0(&mut Input::new(compressed), &mut table_bytes, 4)?;
        read_order1_tables(&mut Input::new(&table_bytes), shift)?
    } else {
        read_order1_tables(input, shift)?
    };

    let size = out.len();
    let interleave = if ways == 4 { size >> 2 } else { size >> 5 };
    let mut rans = vec![0u64; ways];
    let mut index = vec![0usize; ways];
    let mut context = vec![0usize; ways];
    for r in 0..ways {
        rans[r] = input.u32_le()?;
        index[r] = r * interleave;
    }
    let mask = (1u64 << shift) - 1;
    let empty = Table::new();
    let table_of = |c: usize| tables[c].as_ref().unwrap_or(&empty);

    while index[0] < interleave {
        for r in 0..ways {
            let table = table_of(context[r]);
            let symbol = *table
                .reverse
                .get((rans[r] & mask) as usize)
                .ok_or(NxError::Truncated)?;
            out[index[r]] = symbol;
            rans[r] = table.advance(symbol, rans[r], shift);
            rans[r] = renormalize(rans[r], input)?;
            context[r] = symbol as usize;
        }
        for i in index.iter_mut() {
            *i += 1;
        }
    }
    let last = ways - 1;
    while index[last] < size {
        let table = table_of(context[last]);
        let symbol = *table
            .reverse
            .get((rans[last] & mask) as usize)
            .ok_or(NxError::Truncated)?;
        out[index[last]] = symbol;
        rans[last] = table.advance(symbol, rans[last], shift);
        rans[last] = renormalize(rans[last], input)?;
        context[last] = symbol as usize;
        index[last] += 1;
    }
    Ok(())
}

/// `readFrequencyTableOrder1`, for [`order1_reference`].
fn read_order1_tables(input: &mut Input<'_>, shift: u32) -> Result<Vec<Option<Table>>, NxError> {
    let mut tables: Vec<Option<Table>> = (0..NUMBER_OF_SYMBOLS).map(|_| None).collect();
    read_order1_freqs(input, shift, &mut |i, freq| {
        let mut table = Table::new();
        table.reverse = vec![0u8; 1 << shift.max(1)];
        fill_table(&mut table, freq, None);
        tables[i] = Some(table);
    })?;
    Ok(tables)
}

/// `readFrequencyTableOrder1`: the alphabet, then for each context in it the normalised
/// frequencies of the symbols that follow it, handed to `context` as they are read.
fn read_order1_freqs(
    input: &mut Input<'_>,
    shift: u32,
    context: &mut dyn FnMut(usize, &[u32; NUMBER_OF_SYMBOLS]),
) -> Result<(), NxError> {
    let alphabet = read_alphabet(input)?;
    for i in 0..NUMBER_OF_SYMBOLS {
        if !alphabet[i] {
            continue;
        }
        let mut freq = [0u32; NUMBER_OF_SYMBOLS];
        let mut run = 0u32;
        for j in 0..NUMBER_OF_SYMBOLS {
            if !alphabet[j] {
                continue;
            }
            if run > 0 {
                run -= 1;
            } else {
                freq[j] = input.uint7()? as u32;
                if freq[j] == 0 {
                    run = u32::from(input.byte()?);
                }
            }
        }
        normalise_shift(&mut freq, shift);
        context(i, &freq);
    }
    Ok(())
}

/// `decodeRLEMeta`: the symbols that carry runs, and the run lengths after them.
fn decode_rle_meta(
    input: &mut Input<'_>,
    meta_length: usize,
    flags: u8,
) -> Result<(Vec<u8>, [bool; NUMBER_OF_SYMBOLS], usize), NxError> {
    let meta = if meta_length & 1 != 0 {
        input.take((meta_length - 1) / 2)?.to_vec()
    } else {
        let compressed_length = input.uint7()? as usize;
        let compressed = input.take(compressed_length)?;
        let mut meta = vec![0u8; meta_length / 2];
        let ways = if flags & N32_FLAG_MASK != 0 { 32 } else { 4 };
        order0(&mut Input::new(compressed), &mut meta, ways)?;
        meta
    };
    let mut symbols = [false; NUMBER_OF_SYMBOLS];
    let mut pos = 0;
    let mut count = *meta.get(pos).ok_or(NxError::Truncated)? as usize;
    pos += 1;
    if count == 0 {
        count = NUMBER_OF_SYMBOLS;
    }
    for _ in 0..count {
        symbols[*meta.get(pos).ok_or(NxError::Truncated)? as usize] = true;
        pos += 1;
    }
    Ok((meta, symbols, pos))
}

/// `decodeRLE`.
///
/// The output is allocated once at its declared length and filled by position. A run that would
/// overshoot that length is refused before it is written: the reference's array refuses it too,
/// and checking first keeps a damaged run length from writing past the end on the way there. A
/// run length that `readUint7` overflows to a negative number writes the symbol not even once, as
/// the reference's `for (j = 0; j <= run; j++)` does.
///
/// Most runs are short, so a run that leaves room is written as one 16-byte store and the cursor
/// moves by the run's length: the bytes past it are overwritten by what follows.
fn decode_rle(
    input: &[u8],
    symbols: &[bool; NUMBER_OF_SYMBOLS],
    meta: &[u8],
    meta_pos: usize,
    length: usize,
) -> Result<Vec<u8>, NxError> {
    #[cfg(test)]
    if reference::ON.with(|on| on.get()) {
        return decode_rle_reference(input, symbols, meta, meta_pos, length);
    }
    const SHORT: usize = 16;
    let mut out = vec![0u8; length];
    let mut runs = Input {
        bytes: meta,
        pos: meta_pos,
    };
    let mut literals = input.iter();
    let mut at = 0;
    while at < length {
        let symbol = *literals.next().ok_or(NxError::Truncated)?;
        if symbols[symbol as usize] {
            let run = runs.uint7()?;
            let count = if run < 0 { 0 } else { run as usize + 1 };
            if count > length - at {
                return Err(NxError::Truncated);
            }
            if count <= SHORT && at + SHORT <= length {
                out[at..at + SHORT].fill(symbol);
            } else {
                out[at..at + count].fill(symbol);
            }
            at += count;
        } else {
            out[at] = symbol;
            at += 1;
        }
    }
    Ok(out)
}

/// `decodeRLE` as first ported, one byte at a time: what [`decode_rle`] is checked against.
#[cfg(test)]
fn decode_rle_reference(
    input: &[u8],
    symbols: &[bool; NUMBER_OF_SYMBOLS],
    meta: &[u8],
    meta_pos: usize,
    length: usize,
) -> Result<Vec<u8>, NxError> {
    let mut out = Vec::with_capacity(length);
    let mut runs = Input {
        bytes: meta,
        pos: meta_pos,
    };
    let mut i = 0;
    while out.len() < length {
        let symbol = *input.get(i).ok_or(NxError::Truncated)?;
        i += 1;
        if symbols[symbol as usize] {
            let run = runs.uint7()?;
            for _ in 0..=run {
                out.push(symbol);
                // Past the length the answer is already the refusal below; stopping here keeps a
                // damaged run from allocating gigabytes in the tests.
                if out.len() > length {
                    break;
                }
            }
        } else {
            out.push(symbol);
        }
    }
    // A run may overshoot the declared length, which the reference's array would refuse.
    if out.len() > length {
        return Err(NxError::Truncated);
    }
    Ok(out)
}

/// `CompressionUtils.decodePack`.
pub(crate) fn decode_pack(
    input: &[u8],
    table: &[u8],
    symbols: usize,
    length: usize,
) -> Result<Vec<u8>, NxError> {
    let mut out = vec![0u8; length];
    let (bits, per_byte, mask) = match symbols {
        0..=1 => {
            out.fill(*table.first().ok_or(NxError::Truncated)?);
            return Ok(out);
        }
        2 => (1, 8, 1u8),
        3..=4 => (2, 4, 3u8),
        _ => (4, 2, 15u8),
    };
    let mut v = 0u8;
    let mut j = 0;
    for (i, slot) in out.iter_mut().enumerate() {
        if i % per_byte == 0 {
            v = *input.get(j).ok_or(NxError::Truncated)?;
            j += 1;
        }
        *slot = *table.get((v & mask) as usize).ok_or(NxError::Truncated)?;
        v >>= bits;
    }
    Ok(out)
}

/// `decodeStripe`.
fn decode_stripe(input: &mut Input<'_>, size: usize) -> Result<Vec<u8>, NxError> {
    let streams = input.byte()? as usize;
    for _ in 0..streams {
        input.uint7()?;
    }
    let mut lengths = vec![0usize; streams];
    let mut parts = Vec::with_capacity(streams);
    for (j, length) in lengths.iter_mut().enumerate() {
        *length = size / streams + usize::from(size % streams > j);
        parts.push(uncompress_internal(input, *length)?);
    }
    let mut out = vec![0u8; size];
    for (j, part) in parts.iter().enumerate() {
        for i in 0..lengths[j] {
            out[i * streams + j] = *part.get(i).ok_or(NxError::Truncated)?;
        }
    }
    Ok(out)
}

// ---- Encoding: `RANSNx16Encode` ----

/// `RANSEncodingSymbol`, the reciprocal form the encoder divides by.
#[derive(Clone, Copy, Default)]
struct EncodingSymbol {
    x_max: u64,
    rcp_freq: u32,
    bias: u32,
    cmpl_freq: u32,
    rcp_shift: u32,
}

impl EncodingSymbol {
    fn set(start: u32, freq: u32, scale_bits: u32) -> Self {
        let mut symbol = EncodingSymbol {
            x_max: (1u64 << (31 - scale_bits)) * u64::from(freq),
            cmpl_freq: (1u32 << scale_bits).wrapping_sub(freq),
            ..Default::default()
        };
        if freq < 2 {
            symbol.rcp_freq = u32::MAX;
            symbol.rcp_shift = 0;
            symbol.bias = start + (1 << scale_bits) - 1;
        } else {
            let mut shift = 0;
            while u64::from(freq) > (1u64 << shift) {
                shift += 1;
            }
            symbol.rcp_freq = (1u64 << (shift + 31)).div_ceil(u64::from(freq)) as u32;
            symbol.rcp_shift = shift - 1;
            symbol.bias = start;
        }
        symbol.rcp_shift += 32;
        symbol
    }

    /// `putSymbolNx16`, writing backwards from `pos`.
    #[inline]
    fn put(&self, r: u64, out: &mut [u8], pos: &mut usize) -> u64 {
        let mut x = r;
        if x >= self.x_max {
            *pos -= 1;
            out[*pos] = (x >> 8) as u8;
            *pos -= 1;
            out[*pos] = x as u8;
            x >>= 16;
            if x >= self.x_max {
                *pos -= 1;
                out[*pos] = (x >> 8) as u8;
                *pos -= 1;
                out[*pos] = x as u8;
                x >>= 16;
            }
        }
        let q = (x.wrapping_mul(u64::from(self.rcp_freq))) >> self.rcp_shift;
        x.wrapping_add(u64::from(self.bias))
            .wrapping_add(q.wrapping_mul(u64::from(self.cmpl_freq)))
    }
}

/// `Utils.normaliseFrequenciesOrder0`.
fn normalise_order0(f: &mut [i32; NUMBER_OF_SYMBOLS], bits: u32) {
    let total: i32 = f.iter().fold(0i32, |a, &b| a.wrapping_add(b));
    let renorm = 1i32 << bits;
    let mut m = 0;
    let mut big = 0;
    for (symbol, &freq) in f.iter().enumerate() {
        if m < freq {
            m = freq;
            big = symbol;
        }
    }
    let tr: i64 = if total > 0 {
        ((i64::from(renorm) << 31) / i64::from(total)) + i64::from((1i32 << 30) / total)
    } else {
        0
    };
    let mut sum = 0i32;
    for freq in f.iter_mut() {
        if *freq == 0 {
            continue;
        }
        *freq = ((i64::from(*freq) * tr) >> 31) as i32;
        if *freq == 0 {
            *freq = 1;
        }
        sum += *freq;
    }
    if sum < renorm {
        f[big] += renorm - sum;
    } else if sum > renorm {
        f[big] -= sum - renorm;
    }
}

/// `Utils.normaliseFrequenciesOrder0Shift`, on signed counts.
fn normalise_shift_i32(f: &mut [i32; NUMBER_OF_SYMBOLS], bits: u32) {
    let mut total: i64 = f.iter().map(|&x| i64::from(x)).sum();
    if total == 0 || total == 1 << bits {
        return;
    }
    let mut shift = 0;
    while total < 1 << bits {
        total *= 2;
        shift += 1;
    }
    for freq in f.iter_mut() {
        if *freq != 0 {
            *freq <<= shift;
        }
    }
}

fn build_symbols(f: &[i32; NUMBER_OF_SYMBOLS]) -> [EncodingSymbol; NUMBER_OF_SYMBOLS] {
    let mut symbols = [EncodingSymbol::default(); NUMBER_OF_SYMBOLS];
    let mut cumulative = 0u32;
    for (s, &freq) in f.iter().enumerate() {
        if freq != 0 {
            symbols[s] = EncodingSymbol::set(cumulative, freq as u32, TOTAL_FREQ_SHIFT);
            cumulative += freq as u32;
        }
    }
    symbols
}

/// `writeAlphabet`.
fn write_alphabet(out: &mut Vec<u8>, f: &[i32; NUMBER_OF_SYMBOLS]) {
    let mut rle = 0usize;
    for j in 0..NUMBER_OF_SYMBOLS {
        if f[j] != 0 {
            if rle != 0 {
                rle -= 1;
            } else {
                out.push(j as u8);
                if j != 0 && f[j - 1] != 0 {
                    let mut k = j + 1;
                    while k < NUMBER_OF_SYMBOLS && f[k] != 0 {
                        k += 1;
                    }
                    rle = k - (j + 1);
                    out.push(rle as u8);
                }
            }
        }
    }
    out.push(0);
}

/// `RANSNx16Encode.compress`.
pub fn compress(input: &[u8], flags: u8) -> Vec<u8> {
    if input.is_empty() {
        return Vec::new();
    }
    let mut out = vec![flags];
    if flags & NOSZ_FLAG_MASK == 0 {
        write_uint7(input.len() as i32, &mut out);
    }
    if flags & STRIPE_FLAG_MASK != 0 {
        let sizes = crate::range::stripe_sizes(input.len());
        let chunks = crate::range::stripe_transpose(input, &sizes);
        let compressed: Vec<Vec<u8>> = chunks.iter().map(|c| compress(c, NOSZ_FLAG_MASK)).collect();
        out.push(chunks.len() as u8);
        for chunk in &compressed {
            write_uint7(chunk.len() as i32, &mut out);
        }
        for chunk in &compressed {
            out.extend_from_slice(chunk);
        }
        return out;
    }

    let mut data: Vec<u8> = input.to_vec();
    if flags & PACK_FLAG_MASK != 0 {
        let mut frequencies = [0u32; 256];
        for &b in input {
            frequencies[b as usize] += 1;
        }
        let symbols = frequencies.iter().filter(|&&f| f > 0).count();
        // Unlike the arithmetic coder's, one symbol is not packed here.
        if symbols > 1 && symbols <= 16 {
            data = crate::range::encode_pack(input, &mut out).expect("2..=16 symbols pack");
        } else {
            out[0] &= !PACK_FLAG_MASK;
        }
    }
    if flags & RLE_FLAG_MASK != 0 {
        data = encode_rle(&data, &mut out, flags);
    }
    if flags & CAT_FLAG_MASK != 0 {
        out.extend_from_slice(&data);
        return out;
    }
    let ways = if out[0] & N32_FLAG_MASK != 0 { 32 } else { 4 };
    let mut order1 = out[0] & ORDER_FLAG_MASK != 0;
    if data.len() < ways && order1 {
        out[0] &= !ORDER_FLAG_MASK;
        order1 = false;
        if data.is_empty() {
            return out;
        }
    }
    let encoded = if order1 {
        compress_order1(&data, ways)
    } else {
        compress_order0(&data, ways)
    };
    out.extend_from_slice(&encoded);
    out
}

/// `compressOrder0WayN`: the frequency table, then the states and the stream.
fn compress_order0(input: &[u8], ways: usize) -> Vec<u8> {
    let size = input.len();
    let mut bit_size = if size <= 1 {
        0
    } else {
        32 - ((size - 1) as u32).leading_zeros()
    };
    if bit_size > TOTAL_FREQ_SHIFT {
        bit_size = TOTAL_FREQ_SHIFT;
    }
    let mut f = [0i32; NUMBER_OF_SYMBOLS];
    for &b in input {
        f[b as usize] += 1;
    }
    normalise_order0(&mut f, bit_size);
    let mut table = Vec::new();
    write_alphabet(&mut table, &f);
    for &freq in f.iter() {
        if freq != 0 {
            if freq < 128 {
                table.push((freq & 0x7f) as u8);
            } else {
                table.push((128 | (freq >> 7)) as u8);
                table.push((freq & 0x7f) as u8);
            }
        }
    }
    if bit_size != TOTAL_FREQ_SHIFT {
        normalise_shift_i32(&mut f, TOTAL_FREQ_SHIFT);
    }
    let symbols = build_symbols(&f);

    let interleave = if ways == 4 { size >> 2 } else { size >> 5 };
    let mut remaining = size - interleave * ways;
    let mut rans = vec![LOWER_BOUND; ways];
    let capacity = size + size / 4 + ways * 4 + 64;
    let mut buffer = vec![0u8; capacity];
    let mut pos = capacity;

    let mut reverse = 1;
    while remaining > 0 {
        let symbol = &symbols[input[size - reverse] as usize];
        rans[remaining - 1] = symbol.put(rans[remaining - 1], &mut buffer, &mut pos);
        remaining -= 1;
        reverse += 1;
    }
    let mut i = interleave * ways;
    while i > 0 {
        for r in (0..ways).rev() {
            let symbol = &symbols[input[i - (ways - r)] as usize];
            rans[r] = symbol.put(rans[r], &mut buffer, &mut pos);
        }
        i -= ways;
    }
    for r in (0..ways).rev() {
        let state = rans[r] as u32;
        for shift in [24, 16, 8, 0] {
            pos -= 1;
            buffer[pos] = (state >> shift) as u8;
        }
    }
    table.extend_from_slice(&buffer[pos..]);
    table
}

/// `compressOrder1WayN`.
fn compress_order1(input: &[u8], ways: usize) -> Vec<u8> {
    let size = input.len();
    // Row 256 is the context totals.
    let mut f = vec![[0i32; NUMBER_OF_SYMBOLS]; NUMBER_OF_SYMBOLS + 1];
    let mut context = 0usize;
    for &b in input {
        f[NUMBER_OF_SYMBOLS][context] += 1;
        f[context][b as usize] += 1;
        context = b as usize;
    }
    f[NUMBER_OF_SYMBOLS][context] += 1;
    for n in 1..ways {
        let at = if ways == 4 {
            n * (size >> 2)
        } else {
            n * (size >> 5)
        };
        f[0][input[at] as usize] += 1;
    }
    f[NUMBER_OF_SYMBOLS][0] += ways as i32 - 1;

    // `normaliseFrequenciesOrder1`, whose bit size comes from `Math.log`.
    let totals = f[NUMBER_OF_SYMBOLS];
    for (j, &total) in totals.iter().enumerate() {
        if total == 0 {
            continue;
        }
        let ratio = jmath::math::log(f64::from(total)) / jmath::math::log(2.0);
        let mut bit_size = ratio.ceil() as i32;
        if bit_size > TOTAL_FREQ_SHIFT as i32 {
            bit_size = TOTAL_FREQ_SHIFT as i32;
        }
        if bit_size == 0 {
            bit_size = 1;
        }
        normalise_order0(&mut f[j], bit_size as u32);
    }

    let mut table = Vec::new();
    let totals = f[NUMBER_OF_SYMBOLS];
    write_alphabet(&mut table, &totals);
    for i in 0..NUMBER_OF_SYMBOLS {
        if totals[i] == 0 {
            continue;
        }
        let mut run = 0;
        for j in 0..NUMBER_OF_SYMBOLS {
            if totals[j] == 0 {
                continue;
            }
            if run > 0 {
                run -= 1;
                continue;
            }
            write_uint7(f[i][j], &mut table);
            if f[i][j] == 0 {
                for k in j + 1..NUMBER_OF_SYMBOLS {
                    if totals[k] == 0 {
                        continue;
                    }
                    if f[i][k] == 0 {
                        run += 1;
                    } else {
                        break;
                    }
                }
                table.push(run as u8);
            }
        }
    }

    // `~(ORDER | N32)`: order 0, four ways.
    let compressed_table = compress_order0(&table, 4);
    let mut out = Vec::new();
    if compressed_table.len() < table.len() {
        out.push((1 | TOTAL_FREQ_SHIFT << 4) as u8);
        write_uint7(table.len() as i32, &mut out);
        write_uint7(compressed_table.len() as i32, &mut out);
        out.extend_from_slice(&compressed_table);
    } else {
        out.push((TOTAL_FREQ_SHIFT << 4) as u8);
        out.extend_from_slice(&table);
    }

    let totals = f[NUMBER_OF_SYMBOLS];
    for (i, row) in f.iter_mut().take(NUMBER_OF_SYMBOLS).enumerate() {
        if totals[i] != 0 {
            normalise_shift_i32(row, TOTAL_FREQ_SHIFT);
        }
    }
    let symbols: Vec<[EncodingSymbol; NUMBER_OF_SYMBOLS]> = (0..NUMBER_OF_SYMBOLS)
        .map(|i| build_symbols(&f[i]))
        .collect();

    let mut rans = vec![LOWER_BOUND; ways];
    let interleave = if ways == 4 { size >> 2 } else { size >> 5 };
    let mut index = vec![0i64; ways];
    let mut symbol = vec![0u8; ways];
    for r in 0..ways {
        index[r] = ((r + 1) * interleave) as i64 - 2;
        if index[r] + 1 >= 0 && r != ways - 1 {
            symbol[r] = input[(index[r] + 1) as usize];
        }
        if r == ways - 1 {
            symbol[r] = input[size - 1];
        }
    }
    let capacity = size + size / 4 + ways * 4 + 64;
    let mut buffer = vec![0u8; capacity];
    let mut pos = capacity;
    let last = ways - 1;

    index[last] = size as i64 - 2;
    while index[last] > (ways * interleave) as i64 - 2 && index[last] >= 0 {
        let context = input[index[last] as usize];
        rans[last] =
            symbols[context as usize][symbol[last] as usize].put(rans[last], &mut buffer, &mut pos);
        symbol[last] = context;
        index[last] -= 1;
    }
    while index[0] >= 0 {
        for r in 0..ways {
            let k = ways - 1 - r;
            let context = input[index[k] as usize];
            rans[k] =
                symbols[context as usize][symbol[k] as usize].put(rans[k], &mut buffer, &mut pos);
            symbol[k] = context;
        }
        for i in index.iter_mut() {
            *i -= 1;
        }
    }
    for r in 0..ways {
        let k = ways - 1 - r;
        rans[k] = symbols[0][symbol[k] as usize].put(rans[k], &mut buffer, &mut pos);
    }
    for r in (0..ways).rev() {
        let state = rans[r] as u32;
        for shift in [24, 16, 8, 0] {
            pos -= 1;
            buffer[pos] = (state >> shift) as u8;
        }
    }
    out.extend_from_slice(&buffer[pos..]);
    out
}

/// `encodeRLE`: the literals, with the run metadata appended to `out`.
fn encode_rle(input: &[u8], out: &mut Vec<u8>, flags: u8) -> Vec<u8> {
    let mut run_counts = [0i32; NUMBER_OF_SYMBOLS];
    let mut last: i32 = -1;
    for &b in input {
        let s = i32::from(b);
        run_counts[b as usize] += if s == last { 1 } else { -1 };
        last = s;
    }
    let mut count = run_counts.iter().filter(|&&c| c > 0).count();
    if count == 0 {
        count = 1;
        run_counts[0] = 1;
    }
    let mut meta = vec![count as u8];
    for (i, &c) in run_counts.iter().enumerate() {
        if c > 0 {
            meta.push(i as u8);
        }
    }
    let mut literals = Vec::with_capacity(input.len());
    let mut i = 0;
    while i < input.len() {
        literals.push(input[i]);
        if run_counts[input[i] as usize] > 0 {
            let symbol = input[i];
            let mut run = 0;
            while i + run + 1 < input.len() && input[i + run + 1] == symbol {
                run += 1;
            }
            write_uint7(run as i32, &mut meta);
            i += run;
        }
        i += 1;
    }
    let ways = if flags & N32_FLAG_MASK != 0 { 32 } else { 4 };
    let compressed_meta = compress_order0(&meta, ways);
    write_uint7((meta.len() * 2) as i32, out);
    write_uint7(literals.len() as i32, out);
    write_uint7(compressed_meta.len() as i32, out);
    out.extend_from_slice(&compressed_meta);
    literals
}

/// Switches the decoders to the ones first ported, for the tests that compare the two.
#[cfg(test)]
mod reference {
    thread_local! {
        pub(super) static ON: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    }
}

#[cfg(test)]
mod decode_tests {
    use super::*;

    /// Each side's answer, `None` where it panicked: a size that `readUint7` overflows to a
    /// negative number asks both for more memory than exists.
    type Answer = Option<Result<Vec<u8>, NxError>>;

    fn both(stream: &[u8]) -> (Answer, Answer) {
        let run = || std::panic::catch_unwind(|| uncompress(stream)).ok();
        let hook = std::panic::take_hook();
        std::panic::set_hook(Box::new(|_| {}));
        let fast = run();
        reference::ON.with(|on| on.set(true));
        let first = run();
        reference::ON.with(|on| on.set(false));
        std::panic::set_hook(hook);
        (fast, first)
    }

    struct Lcg(u64);

    impl Lcg {
        fn next(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            self.0 >> 33
        }
    }

    /// Quality-like text: a few symbols, in runs, with the order-1 structure the 32-state
    /// decoders are fastest on.
    fn qualities(random: &mut Lcg, length: usize) -> Vec<u8> {
        let mut text = Vec::with_capacity(length);
        let mut symbol = b'F';
        while text.len() < length {
            if random.next().is_multiple_of(4) {
                symbol = b"#,:FF"[(random.next() % 5) as usize];
            }
            let run = 1 + (random.next() % 6) as usize;
            text.extend(std::iter::repeat_n(symbol, run.min(length - text.len())));
        }
        text
    }

    /// The decoders of this file against the ones first ported, on every flag combination over
    /// conforming streams, every prefix of short ones, and streams with bytes changed: the same
    /// bytes or the same refusal each time. On x86-64 with AVX2 this covers the vector loop too.
    #[test]
    fn fast_decoders_answer_as_the_first_port_does() {
        let mut random = Lcg(7);
        let mut texts = vec![qualities(&mut random, 100_000), qualities(&mut random, 777)];
        texts.push((0..5000u32).map(|i| (i * 7 % 13) as u8).collect());
        texts.push(b"ACGT".repeat(300));
        let flags = [
            0u8, 1, 4, 5, 0x40, 0x41, 0x44, 0x45, 0x80, 0x81, 0xC5, 0x08, 0x09,
        ];
        let mut compared = 0;
        for text in &texts {
            for &flag in &flags {
                let stream = compress(text, flag);
                let (fast, first) = both(&stream);
                assert_eq!(fast, Some(Ok(text.clone())), "flags {flag:#x}");
                assert_eq!(fast, first, "flags {flag:#x}");
                let cuts: Vec<usize> = if stream.len() < 600 {
                    (0..stream.len()).collect()
                } else {
                    (0..200)
                        .map(|_| (random.next() as usize) % stream.len())
                        .collect()
                };
                for cut in cuts {
                    let (fast, first) = both(&stream[..cut]);
                    assert_eq!(fast, first, "flags {flag:#x}, cut at {cut}");
                    compared += 1;
                }
                for _ in 0..300 {
                    let mut damaged = stream.clone();
                    // Past the flags and the size, which would only change what is allocated.
                    let at = 6 + (random.next() as usize) % (damaged.len() - 6).max(1);
                    if at < damaged.len() {
                        damaged[at] ^= 1 << (random.next() % 8);
                    }
                    let (fast, first) = both(&damaged);
                    assert_eq!(fast, first, "flags {flag:#x}, bit changed at {at}");
                    compared += 1;
                }
            }
        }
        assert!(compared > 10_000);
    }
}

#[cfg(test)]
mod encode_tests {
    use super::*;

    #[test]
    fn every_flag_combination_round_trips() {
        let mut text =
            b"AAAAAACCCCGGGTTTTTTTTTTNNNNACGTACGTAAAAAAAAAAAAAAAAAAAAAAAAACGT".repeat(40);
        text.extend((0..3000u32).map(|i| (i * 7 % 13) as u8));
        for flags in [
            0u8, 1, 4, 5, 0x40, 0x41, 0x80, 0x81, 0xC0, 0xC1, 0x20, 0x08, 0x09,
        ] {
            let packed = compress(&text, flags);
            assert_eq!(uncompress(&packed).unwrap(), text, "flags {flags:#x}");
        }
        for short in [&b"A"[..], b"AC", b"ACG", b"ACGTA"] {
            for flags in [0u8, 1, 5] {
                assert_eq!(uncompress(&compress(short, flags)).unwrap(), short);
            }
        }
    }
}
