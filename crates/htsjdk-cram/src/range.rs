//! The CRAM 3.1 adaptive arithmetic codec, block compression method 6.
//!
//! Ported from `htsjdk.samtools.cram.compression.range` (`RangeCoder`, `ByteModel`,
//! `RangeEncode`, `RangeDecode`) at htsjdk 5.0.0.
//!
//! The format flags are rANS Nx16's with one difference: bit 2 means "compressed by bzip2" rather
//! than "32 ways". bzip2 is not ported, so a stream with that flag is refused either way.

use crate::rans_nx16::{decode_pack, write_uint7, Input, NxError};

pub const ORDER_FLAG_MASK: u8 = 0x01;
pub const EXT_FLAG_MASK: u8 = 0x04;
pub const STRIPE_FLAG_MASK: u8 = 0x08;
pub const NOSZ_FLAG_MASK: u8 = 0x10;
pub const CAT_FLAG_MASK: u8 = 0x20;
pub const RLE_FLAG_MASK: u8 = 0x40;
pub const PACK_FLAG_MASK: u8 = 0x80;

const MAX_FREQ: u32 = (1 << 16) - 17;
const STEP: u32 = 16;
const MAX_RANGE: u64 = 0xFFFF_FFFF;
const STRIPE_STREAMS: usize = 4;

/// `RangeCoder`.
pub(crate) struct RangeCoder {
    low: u64,
    range: u64,
    code: u64,
    ff_count: usize,
    carry: bool,
    cache: u8,
    pub(crate) out: Vec<u8>,
}

impl RangeCoder {
    pub(crate) fn new() -> Self {
        Self {
            low: 0,
            range: MAX_RANGE,
            code: 0,
            ff_count: 0,
            carry: false,
            cache: 0,
            out: Vec::new(),
        }
    }

    pub(crate) fn decode_start(&mut self, input: &mut Input<'_>) -> Result<(), NxError> {
        for _ in 0..5 {
            self.code = (self.code << 8).wrapping_add(u64::from(input.byte()?));
        }
        self.code &= MAX_RANGE;
        Ok(())
    }

    fn decode(
        &mut self,
        input: &mut Input<'_>,
        cumulative: u32,
        frequency: u32,
    ) -> Result<(), NxError> {
        self.code = self
            .code
            .wrapping_sub(u64::from(cumulative).wrapping_mul(self.range));
        self.range = self.range.wrapping_mul(u64::from(frequency));
        while self.range < 1 << 24 {
            self.range <<= 8;
            self.code = (self.code << 8).wrapping_add(u64::from(input.byte()?));
        }
        Ok(())
    }

    fn get_frequency(&mut self, total: u32) -> u32 {
        self.range /= u64::from(total);
        (self.code / self.range.max(1)) as u32
    }

    fn encode(&mut self, cumulative: u32, frequency: u32, total: u32) {
        let old_low = self.low;
        self.range /= u64::from(total);
        self.low = self
            .low
            .wrapping_add(u64::from(cumulative).wrapping_mul(self.range));
        self.low &= 0xFFFF_FFFF;
        self.range = self.range.wrapping_mul(u64::from(frequency));
        if self.low < old_low {
            self.carry = true;
        }
        while self.range < 1 << 24 {
            self.range <<= 8;
            self.shift_low();
        }
    }

    pub(crate) fn encode_end(&mut self) {
        for _ in 0..5 {
            self.shift_low();
        }
    }

    fn shift_low(&mut self) {
        if self.low < 0xFF00_0000 || self.carry {
            if !self.carry {
                self.out.push(self.cache);
                while self.ff_count > 0 {
                    self.out.push(0xFF);
                    self.ff_count -= 1;
                }
            } else {
                self.out.push(self.cache.wrapping_add(1));
                while self.ff_count > 0 {
                    self.out.push(0x00);
                    self.ff_count -= 1;
                }
            }
            self.cache = (self.low >> 24) as u8;
            self.carry = false;
        } else {
            self.ff_count += 1;
        }
        self.low = (self.low << 8) & 0xFFFF_FFFF;
    }
}

/// `ByteModel`: frequencies and symbols, kept roughly sorted by one swap per update.
#[derive(Clone)]
pub(crate) struct ByteModel {
    total: u32,
    max_symbol: usize,
    freq: Vec<u32>,
    symbol: Vec<u32>,
}

impl ByteModel {
    pub(crate) fn new(symbols: usize) -> Self {
        Self {
            total: symbols as u32,
            max_symbol: symbols - 1,
            freq: vec![1; symbols],
            symbol: (0..symbols as u32).collect(),
        }
    }

    fn renormalize(&mut self) {
        self.total = 0;
        for f in self.freq.iter_mut() {
            *f -= *f >> 1;
            self.total += *f;
        }
    }

    fn update(&mut self, x: usize) {
        self.freq[x] += STEP;
        self.total += STEP;
        if self.total > MAX_FREQ {
            self.renormalize();
        }
        if x > 0 && self.freq[x] > self.freq[x - 1] {
            self.freq.swap(x, x - 1);
            self.symbol.swap(x, x - 1);
        }
    }

    pub(crate) fn decode(
        &mut self,
        input: &mut Input<'_>,
        coder: &mut RangeCoder,
    ) -> Result<u32, NxError> {
        let target = coder.get_frequency(self.total);
        let mut cumulative = 0;
        let mut x = 0;
        while cumulative + self.freq[x] <= target {
            cumulative += self.freq[x];
            x += 1;
            if x > self.max_symbol {
                return Err(NxError::Truncated);
            }
        }
        coder.decode(input, cumulative, self.freq[x])?;
        let symbol = self.symbol[x];
        self.update(x);
        Ok(symbol)
    }

    pub(crate) fn encode(&mut self, coder: &mut RangeCoder, symbol: u32) {
        let mut cumulative = 0;
        let mut i = 0;
        while self.symbol[i] != symbol {
            cumulative += self.freq[i];
            i += 1;
        }
        coder.encode(cumulative, self.freq[i], self.total);
        self.update(i);
    }
}

/// `RangeDecode.uncompress`.
pub fn uncompress(input: &[u8]) -> Result<Vec<u8>, NxError> {
    let mut cursor = Input::new(input);
    uncompress_sized(&mut cursor, 0)
}

fn uncompress_sized(input: &mut Input<'_>, out_size: usize) -> Result<Vec<u8>, NxError> {
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
    let mut pack = None;
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

    let mut out = if flags & CAT_FLAG_MASK != 0 {
        input.take(size)?.to_vec()
    } else if flags & EXT_FLAG_MASK != 0 {
        return Err(NxError::Cram(
            "bzip2 inside the arithmetic codec is not ported".to_string(),
        ));
    } else {
        let mut out = vec![0u8; size];
        let order1 = flags & ORDER_FLAG_MASK != 0;
        if flags & RLE_FLAG_MASK != 0 {
            decode_rle(input, &mut out, order1)?;
        } else {
            decode_plain(input, &mut out, order1)?;
        }
        out
    };
    if let Some((table, symbols, length)) = pack {
        out = decode_pack(&out, &table, symbols, length)?;
    }
    Ok(out)
}

fn max_symbols(input: &mut Input<'_>) -> Result<usize, NxError> {
    let m = input.byte()? as usize;
    Ok(if m == 0 { 256 } else { m })
}

fn decode_plain(input: &mut Input<'_>, out: &mut [u8], order1: bool) -> Result<(), NxError> {
    let symbols = max_symbols(input)?;
    let mut models = vec![ByteModel::new(symbols); if order1 { symbols } else { 1 }];
    let mut coder = RangeCoder::new();
    coder.decode_start(input)?;
    let mut last = 0usize;
    for slot in out.iter_mut() {
        let model = models
            .get_mut(if order1 { last } else { 0 })
            .ok_or(NxError::Truncated)?;
        let symbol = model.decode(input, &mut coder)?;
        *slot = symbol as u8;
        last = symbol as usize;
    }
    Ok(())
}

fn decode_rle(input: &mut Input<'_>, out: &mut [u8], order1: bool) -> Result<(), NxError> {
    let symbols = max_symbols(input)?;
    let mut literals = vec![ByteModel::new(symbols); if order1 { symbols } else { 1 }];
    let mut runs = vec![ByteModel::new(4); 258];
    let mut coder = RangeCoder::new();
    coder.decode_start(input)?;
    let mut last = 0usize;
    let mut i = 0;
    while i < out.len() {
        let model = literals
            .get_mut(if order1 { last } else { 0 })
            .ok_or(NxError::Truncated)?;
        let symbol = model.decode(input, &mut coder)? as u8;
        out[i] = symbol;
        last = symbol as usize;
        let mut part = runs[last].decode(input, &mut coder)?;
        let mut run = part as usize;
        let mut context = 256;
        while part == 3 {
            part = runs[context].decode(input, &mut coder)?;
            context = 257;
            run += part as usize;
        }
        for j in 1..=run {
            *out.get_mut(i + j).ok_or(NxError::Truncated)? = symbol;
        }
        i += run + 1;
    }
    Ok(())
}

fn decode_stripe(input: &mut Input<'_>, size: usize) -> Result<Vec<u8>, NxError> {
    let streams = input.byte()? as usize;
    for _ in 0..streams {
        input.uint7()?;
    }
    let mut lengths = vec![0usize; streams];
    let mut parts = Vec::with_capacity(streams);
    for (j, length) in lengths.iter_mut().enumerate() {
        *length = size / streams + usize::from(size % streams > j);
        parts.push(uncompress_sized(input, *length)?);
    }
    let mut out = vec![0u8; size];
    for (j, part) in parts.iter().enumerate() {
        for i in 0..lengths[j] {
            out[i * streams + j] = *part.get(i).ok_or(NxError::Truncated)?;
        }
    }
    Ok(out)
}

/// `RangeEncode.compress`.
pub fn compress(input: &[u8], flags: u8) -> Result<Vec<u8>, NxError> {
    if input.is_empty() {
        return Ok(Vec::new());
    }
    let mut out = vec![flags];
    if flags & NOSZ_FLAG_MASK == 0 {
        write_uint7(input.len() as i32, &mut out);
    }
    if flags & STRIPE_FLAG_MASK != 0 {
        let sizes = stripe_sizes(input.len());
        let chunks = stripe_transpose(input, &sizes);
        let compressed: Vec<Vec<u8>> = chunks
            .iter()
            .map(|chunk| compress(chunk, NOSZ_FLAG_MASK))
            .collect::<Result<_, _>>()?;
        out.push(STRIPE_STREAMS as u8);
        for chunk in &compressed {
            write_uint7(chunk.len() as i32, &mut out);
        }
        for chunk in &compressed {
            out.extend_from_slice(chunk);
        }
        return Ok(out);
    }

    let mut data: std::borrow::Cow<'_, [u8]> = std::borrow::Cow::Borrowed(input);
    if flags & PACK_FLAG_MASK != 0 {
        match encode_pack(input, &mut out) {
            Some(packed) => data = std::borrow::Cow::Owned(packed),
            None => out[0] &= !PACK_FLAG_MASK,
        }
    }

    if flags & CAT_FLAG_MASK != 0 {
        out.extend_from_slice(&data);
    } else if flags & EXT_FLAG_MASK != 0 {
        return Err(NxError::Cram(
            "bzip2 inside the arithmetic codec is not ported".to_string(),
        ));
    } else {
        let order1 = flags & ORDER_FLAG_MASK != 0;
        let symbols = data.iter().copied().max().unwrap_or(0) as usize + 1;
        out.push(symbols as u8);
        let mut coder = RangeCoder::new();
        if flags & RLE_FLAG_MASK != 0 {
            encode_rle(&data, symbols, order1, &mut coder);
        } else {
            let mut models = vec![ByteModel::new(symbols); if order1 { symbols } else { 1 }];
            let mut last = 0usize;
            for &b in data.iter() {
                models[if order1 { last } else { 0 }].encode(&mut coder, u32::from(b));
                last = b as usize;
            }
        }
        coder.encode_end();
        out.extend_from_slice(&coder.out);
    }
    Ok(out)
}

fn encode_rle(data: &[u8], symbols: usize, order1: bool, coder: &mut RangeCoder) {
    let mut literals = vec![ByteModel::new(symbols); if order1 { symbols } else { 1 }];
    let mut runs = vec![ByteModel::new(4); 258];
    let mut i = 0;
    let mut last = 0usize;
    while i < data.len() {
        literals[if order1 { last } else { 0 }].encode(coder, u32::from(data[i]));
        let mut run = 1;
        while i + run < data.len() && data[i + run] == data[i] {
            run += 1;
        }
        run -= 1;
        let mut context = data[i] as usize;
        last = data[i] as usize;
        i += run + 1;
        let mut part = run.min(3);
        runs[context].encode(coder, part as u32);
        run -= part;
        context = 256;
        while part == 3 {
            part = run.min(3);
            runs[context].encode(coder, part as u32);
            context = 257;
            run -= part;
        }
    }
}

/// `CompressionUtils.encodePack`: the packed bytes, with the symbol count, the mapping table and
/// the packed length appended to `out`; `None` when there are too many symbols to pack.
pub(crate) fn encode_pack(input: &[u8], out: &mut Vec<u8>) -> Option<Vec<u8>> {
    let mut frequencies = [0u32; 256];
    for &b in input {
        frequencies[b as usize] += 1;
    }
    let mut mapping = [0u8; 256];
    let mut symbols = 0usize;
    for i in 0..256 {
        if frequencies[i] > 0 {
            mapping[i] = symbols as u8;
            symbols += 1;
        }
    }
    if symbols == 0 || symbols > 16 {
        return None;
    }
    let packed: Vec<u8> = if symbols <= 1 {
        Vec::new()
    } else {
        let (per_byte, bits) = match symbols {
            2 => (8, 1),
            3..=4 => (4, 2),
            _ => (2, 4),
        };
        let mut packed = vec![0u8; input.len().div_ceil(per_byte)];
        for (i, &b) in input.iter().enumerate() {
            let j = i / per_byte;
            packed[j] = packed[j].wrapping_add(mapping[b as usize] << ((i % per_byte) * bits));
        }
        packed
    };
    out.push(symbols as u8);
    for (i, &f) in frequencies.iter().enumerate() {
        if f > 0 {
            out.push(i as u8);
        }
    }
    write_uint7(packed.len() as i32, out);
    Some(packed)
}

/// `CompressionUtils.buildStripeUncompressedSizes`.
pub(crate) fn stripe_sizes(total: usize) -> [usize; STRIPE_STREAMS] {
    let q = total / STRIPE_STREAMS;
    let r = total % STRIPE_STREAMS;
    std::array::from_fn(|i| if i < r { q + 1 } else { q })
}

/// `CompressionUtils.stripeTranspose`.
pub(crate) fn stripe_transpose(input: &[u8], sizes: &[usize; STRIPE_STREAMS]) -> Vec<Vec<u8>> {
    sizes
        .iter()
        .enumerate()
        .map(|(i, &size)| (0..size).map(|j| input[j * STRIPE_STREAMS + i]).collect())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_flag_combination_round_trips() {
        let text = b"AAAAAACCCCGGGTTTTTTTTTTNNNNACGTACGTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGT";
        for flags in [0u8, 1, 0x40, 0x41, 0x80, 0x81, 0xC0, 0x20, 0x08] {
            let packed = compress(text, flags).unwrap();
            assert_eq!(uncompress(&packed).unwrap(), text, "flags {flags:#x}");
        }
    }
}
