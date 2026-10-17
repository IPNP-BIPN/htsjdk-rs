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

/// `uncompressOrder0WayN`.
pub(crate) fn order0(input: &mut Input<'_>, out: &mut [u8], ways: usize) -> Result<(), NxError> {
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

/// `uncompressOrder1WayN`.
fn order1(input: &mut Input<'_>, out: &mut [u8], ways: usize) -> Result<(), NxError> {
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

/// `readFrequencyTableOrder1`.
fn read_order1_tables(input: &mut Input<'_>, shift: u32) -> Result<Vec<Option<Table>>, NxError> {
    let alphabet = read_alphabet(input)?;
    let mut tables: Vec<Option<Table>> = (0..NUMBER_OF_SYMBOLS).map(|_| None).collect();
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
        let mut table = Table::new();
        table.reverse = vec![0u8; 1 << shift.max(1)];
        fill_table(&mut table, &freq, None);
        tables[i] = Some(table);
    }
    Ok(tables)
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
fn decode_rle(
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
