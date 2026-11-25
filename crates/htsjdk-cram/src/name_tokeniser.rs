//! The CRAM 3.1 name tokeniser, block compression method 8.
//!
//! Ported from `htsjdk.samtools.cram.compression.nametokenisation` at htsjdk 5.0.0. Each read name
//! is split into tokens (strings, characters, digit runs) and coded against an earlier name: a
//! token is either new, a match of the earlier name's token at the same position, or a small
//! delta from it. Every (position, token type) pair is its own byte stream, each compressed with
//! rANS Nx16 or the arithmetic coder.

use crate::rans_nx16::{Input, NxError};
use std::rc::Rc;

pub const TOKEN_TYPE: u8 = 0x00;
pub const TOKEN_STRING: u8 = 0x01;
pub const TOKEN_CHAR: u8 = 0x02;
pub const TOKEN_DIGITS0: u8 = 0x03;
pub const TOKEN_DZLEN: u8 = 0x04;
pub const TOKEN_DUP: u8 = 0x05;
pub const TOKEN_DIFF: u8 = 0x06;
pub const TOKEN_DIGITS: u8 = 0x07;
pub const TOKEN_DELTA: u8 = 0x08;
pub const TOKEN_DELTA0: u8 = 0x09;
pub const TOKEN_MATCH: u8 = 0x0A;
pub const TOKEN_NOP: u8 = 0x0B;
pub const TOKEN_END: u8 = 0x0C;
pub const TOTAL_TOKEN_TYPES: usize = 13;
pub const TYPE_TOKEN_FLAG_MASK: u8 = 0x3F;
pub const DUP_PREVIOUS_STREAM_FLAG_MASK: u8 = 0x40;
pub const NEW_POSITION_FLAG_MASK: u8 = 0x80;

/// `NameTokenisationDecode.NAME_SEPARATOR`.
pub const NAME_SEPARATOR: u8 = 0;

/// One decoded stream and a read position of its own: a duplicated stream shares the bytes and
/// starts from the beginning, as `ByteBuffer.duplicate` does.
#[derive(Clone)]
struct Stream {
    bytes: Rc<Vec<u8>>,
    pos: usize,
}

impl Stream {
    fn byte(&mut self) -> Result<u8, NxError> {
        let b = *self.bytes.get(self.pos).ok_or(NxError::Truncated)?;
        self.pos += 1;
        Ok(b)
    }

    fn int(&mut self) -> Result<u32, NxError> {
        let b = self
            .bytes
            .get(self.pos..self.pos + 4)
            .ok_or(NxError::Truncated)?;
        self.pos += 4;
        Ok(u32::from_le_bytes([b[0], b[1], b[2], b[3]]))
    }

    #[cfg(test)]
    fn string(&mut self) -> Result<Vec<u8>, NxError> {
        let rest = self.bytes.get(self.pos..).ok_or(NxError::Truncated)?;
        let end = rest
            .iter()
            .position(|&b| b == 0)
            .ok_or(NxError::Truncated)?;
        let s = rest[..end].to_vec();
        self.pos += end + 1;
        Ok(s)
    }
}

/// `TokenStreams`.
struct TokenStreams {
    streams: Vec<[Option<Stream>; TOTAL_TOKEN_TYPES]>,
}

impl TokenStreams {
    fn read(input: &mut Input<'_>, use_arith: bool, names: usize) -> Result<Self, NxError> {
        let mut streams: Vec<[Option<Stream>; TOTAL_TOKEN_TYPES]> = Vec::new();
        let mut position: i64 = -1;
        while input.pos < input.bytes.len() {
            let flags = input.byte()?;
            let new_position = flags & NEW_POSITION_FLAG_MASK != 0;
            let is_dup = flags & DUP_PREVIOUS_STREAM_FLAG_MASK != 0;
            let token_type = flags & TYPE_TOKEN_FLAG_MASK;
            if token_type > TOKEN_END {
                return Err(NxError::Cram(format!(
                    "Invalid name tokenizer token stream type: {token_type}"
                )));
            }
            if new_position {
                position += 1;
                streams.push(Default::default());
                if token_type != TOKEN_TYPE {
                    // A type stream that would be one type then all MATCH is left out by the
                    // writer and regenerated here.
                    let mut types = vec![TOKEN_MATCH; names];
                    if let Some(first) = types.first_mut() {
                        *first = token_type;
                    }
                    streams[position as usize][TOKEN_TYPE as usize] = Some(Stream {
                        bytes: Rc::new(types),
                        pos: 0,
                    });
                }
            }
            let at = usize::try_from(position).map_err(|_| NxError::Truncated)?;
            if is_dup {
                let dup_position = input.byte()? as usize;
                let dup_type = input.byte()? as usize;
                let source = streams
                    .get(dup_position)
                    .and_then(|s| s.get(dup_type))
                    .and_then(|s| s.as_ref())
                    .ok_or(NxError::Truncated)?;
                let copy = Stream {
                    bytes: Rc::clone(&source.bytes),
                    pos: 0,
                };
                streams[at][token_type as usize] = Some(copy);
            } else {
                let length = input.uint7()? as usize;
                let compressed = input.take(length)?;
                let bytes = if use_arith {
                    crate::range::uncompress(compressed)?
                } else {
                    crate::rans_nx16::uncompress(compressed)?
                };
                streams[at][token_type as usize] = Some(Stream {
                    bytes: Rc::new(bytes),
                    pos: 0,
                });
            }
        }
        Ok(Self { streams })
    }

    fn get(&mut self, position: usize, token_type: u8) -> Result<&mut Stream, NxError> {
        self.streams
            .get_mut(position)
            .and_then(|s| s[token_type as usize].as_mut())
            .ok_or(NxError::Truncated)
    }
}

/// `NameTokenisationDecode.uncompress`: the names, each followed by `separator`.
///
/// htsjdk keeps every name it has decoded as a list of token strings, for later names to match
/// against. Here the names are written straight into the output and each keeps only where its
/// tokens end in it ([`Names`]): a later name's MATCH, DELTA or DUP reads that slice of the
/// output. The streams are read in the same order, so a damaged block is refused at the same
/// point with the same error.
pub fn uncompress(bytes: &[u8], separator: u8) -> Result<Vec<u8>, NxError> {
    #[cfg(test)]
    if reference::ON.with(|on| on.get()) {
        return uncompress_reference(bytes, separator);
    }
    let mut input = Input::new(bytes);
    let header = input.take(9)?;
    let length = u32::from_le_bytes([header[0], header[1], header[2], header[3]]) as usize;
    let names = u32::from_le_bytes([header[4], header[5], header[6], header[7]]) as usize;
    let use_arith = header[8] != 0;
    let mut streams = TokenStreams::read(&mut input, use_arith, names)?;

    let mut decoded = Names {
        out: Vec::with_capacity(length),
        bounds: Vec::new(),
        first: Vec::with_capacity(names.min(length) + 1),
    };
    for i in 0..names {
        decode_name(&mut streams, &mut decoded, i)?;
        decoded.out.push(separator);
    }
    Ok(decoded.out)
}

/// The names decoded so far, and where their tokens lie in the output.
///
/// Name `i`'s boundaries are `bounds[first[i]..first[i + 1]]`: the offset of its first byte in
/// `out`, then the offset just past each of its tokens, the END token's included. Its token at
/// position `p` (counted from 1, as the streams are) is `out[bounds[first[i] + p - 1]..bounds[first[i] + p]]`.
struct Names {
    out: Vec<u8>,
    bounds: Vec<usize>,
    first: Vec<usize>,
}

impl Names {
    /// The output range of an earlier name's token at `position`, as `decoded.get(reference)
    /// .get(position - 1)` finds it in the reference: absent for a name not decoded yet, or one
    /// with fewer tokens.
    fn token(&self, reference: usize, position: usize) -> Result<std::ops::Range<usize>, NxError> {
        let (from, to) = match (self.first.get(reference), self.first.get(reference + 1)) {
            (Some(&from), Some(&to)) => (from, to),
            _ => return Err(NxError::Truncated),
        };
        if position == 0 || position >= to - from {
            return Err(NxError::Truncated);
        }
        Ok(self.bounds[from + position - 1]..self.bounds[from + position])
    }
}

/// Write `value` in decimal after `width - digits` zeros, as `leftPadWith0(String.valueOf(v))`
/// does: the zeros go before a minus sign too.
fn push_decimal(out: &mut Vec<u8>, value: i64, width: usize) {
    let mut digits = [0u8; 20];
    let mut at = digits.len();
    let mut rest = value.unsigned_abs();
    loop {
        at -= 1;
        digits[at] = b'0' + (rest % 10) as u8;
        rest /= 10;
        if rest == 0 {
            break;
        }
    }
    let length = digits.len() - at + usize::from(value < 0);
    out.resize(out.len() + width.saturating_sub(length), b'0');
    if value < 0 {
        out.push(b'-');
    }
    out.extend_from_slice(&digits[at..]);
}

/// `NameTokenisationDecode.decodeSingleName`, appending the name to `decoded.out`.
fn decode_name(
    streams: &mut TokenStreams,
    decoded: &mut Names,
    index: usize,
) -> Result<(), NxError> {
    let reference_type = streams.get(0, TOKEN_TYPE)?.byte()?;
    let distance = streams.get(0, reference_type)?.int()?;
    let reference = (index as i64 - i64::from(distance as i32)) as usize;

    // The name about to be written starts here; `first` has one entry per name begun, so an
    // earlier name's range is complete once the next one has an entry.
    let start = decoded.out.len();
    decoded.first.push(decoded.bounds.len());
    decoded.bounds.push(start);
    let is_decoded = reference < index;

    if reference_type == TOKEN_DUP {
        if !is_decoded {
            return Err(NxError::Truncated);
        }
        // The same tokens again: the earlier name's bytes, and its boundaries moved to here.
        let (from, to) = (decoded.first[reference], decoded.first[reference + 1]);
        let (name_start, name_end) = (decoded.bounds[from], decoded.bounds[to - 1]);
        decoded.out.extend_from_within(name_start..name_end);
        for k in from + 1..to {
            let bound = decoded.bounds[k] - name_start + start;
            decoded.bounds.push(bound);
        }
        return Ok(());
    } else if reference_type != TOKEN_DIFF {
        return Err(NxError::Cram(format!(
            "Invalid nameType {reference_type}. nameType must be either TOKEN_DIFF or TOKEN_DUP"
        )));
    }

    let mut position = 1usize;
    loop {
        let token_type = streams.get(position, TOKEN_TYPE)?.byte()?;
        match token_type {
            TOKEN_CHAR => {
                let c = streams.get(position, TOKEN_CHAR)?.byte()?;
                decoded.out.push(c);
            }
            TOKEN_STRING => {
                let stream = streams.get(position, TOKEN_STRING)?;
                let rest = stream.bytes.get(stream.pos..).ok_or(NxError::Truncated)?;
                let end = rest
                    .iter()
                    .position(|&b| b == 0)
                    .ok_or(NxError::Truncated)?;
                decoded.out.extend_from_slice(&rest[..end]);
                stream.pos += end + 1;
            }
            TOKEN_DIGITS => {
                let value = streams.get(position, TOKEN_DIGITS)?.int()?;
                push_decimal(&mut decoded.out, i64::from(value), 0);
            }
            TOKEN_DIGITS0 => {
                let value = streams.get(position, TOKEN_DIGITS0)?.int()?;
                let width = streams.get(position, TOKEN_DZLEN)?.byte()? as usize;
                push_decimal(&mut decoded.out, i64::from(value), width);
            }
            TOKEN_DELTA | TOKEN_DELTA0 => {
                let prior = if is_decoded {
                    decoded.token(reference, position)?
                } else {
                    return Err(NxError::Truncated);
                };
                let value: i32 = std::str::from_utf8(&decoded.out[prior.clone()])
                    .ok()
                    .and_then(|s| s.parse().ok())
                    .ok_or_else(|| {
                        NxError::Cram(format!(
                            "The token in the prior name must be of type {}",
                            if token_type == TOKEN_DELTA {
                                "DIGITS or DELTA"
                            } else {
                                "DIGITS0 or DELTA0"
                            }
                        ))
                    })?;
                let delta = i64::from(streams.get(position, token_type)?.byte()?);
                let width = if token_type == TOKEN_DELTA0 {
                    prior.len()
                } else {
                    0
                };
                push_decimal(&mut decoded.out, i64::from(value) + delta, width);
            }
            TOKEN_MATCH => {
                let prior = if is_decoded {
                    decoded.token(reference, position)?
                } else {
                    return Err(NxError::Truncated);
                };
                decoded.out.extend_from_within(prior);
            }
            TOKEN_END | TOKEN_NOP => {}
            other => {
                return Err(NxError::Cram(format!(
                    "Invalid tokenType : {other}. tokenType must be one of the valid token types"
                )))
            }
        }
        decoded.bounds.push(decoded.out.len());
        if token_type == TOKEN_END {
            break;
        }
        position += 1;
    }
    Ok(())
}

/// [`uncompress`] as first ported, a list of token strings per name: what the tests hold the
/// decoder to.
#[cfg(test)]
fn uncompress_reference(bytes: &[u8], separator: u8) -> Result<Vec<u8>, NxError> {
    let mut input = Input::new(bytes);
    let header = input.take(9)?;
    let length = u32::from_le_bytes([header[0], header[1], header[2], header[3]]) as usize;
    let names = u32::from_le_bytes([header[4], header[5], header[6], header[7]]) as usize;
    let use_arith = header[8] != 0;
    let mut streams = TokenStreams::read(&mut input, use_arith, names)?;

    let mut decoded: Vec<Rc<Vec<Vec<u8>>>> = Vec::with_capacity(names);
    let mut out = Vec::with_capacity(length);
    for i in 0..names {
        let name = decode_name_reference(&mut streams, &mut decoded, i)?;
        out.extend_from_slice(&name);
        out.push(separator);
    }
    Ok(out)
}

#[cfg(test)]
fn decode_name_reference(
    streams: &mut TokenStreams,
    decoded: &mut Vec<Rc<Vec<Vec<u8>>>>,
    index: usize,
) -> Result<Vec<u8>, NxError> {
    let reference_type = streams.get(0, TOKEN_TYPE)?.byte()?;
    let distance = streams.get(0, reference_type)?.int()?;
    let reference = (index as i64 - i64::from(distance as i32)) as usize;

    if reference_type == TOKEN_DUP {
        let tokens = Rc::clone(decoded.get(reference).ok_or(NxError::Truncated)?);
        let name = tokens.concat();
        decoded.push(tokens);
        return Ok(name);
    } else if reference_type != TOKEN_DIFF {
        return Err(NxError::Cram(format!(
            "Invalid nameType {reference_type}. nameType must be either TOKEN_DIFF or TOKEN_DUP"
        )));
    }

    let previous = |decoded: &Vec<Rc<Vec<Vec<u8>>>>, position: usize| -> Result<Vec<u8>, NxError> {
        decoded
            .get(reference)
            .and_then(|tokens| tokens.get(position - 1))
            .cloned()
            .ok_or(NxError::Truncated)
    };

    let mut tokens: Vec<Vec<u8>> = Vec::new();
    let mut position = 1usize;
    loop {
        let token_type = streams.get(position, TOKEN_TYPE)?.byte()?;
        let token: Vec<u8> = match token_type {
            TOKEN_CHAR => vec![streams.get(position, TOKEN_CHAR)?.byte()?],
            TOKEN_STRING => streams.get(position, TOKEN_STRING)?.string()?,
            TOKEN_DIGITS => streams
                .get(position, TOKEN_DIGITS)?
                .int()?
                .to_string()
                .into_bytes(),
            TOKEN_DIGITS0 => {
                let digits = streams.get(position, TOKEN_DIGITS0)?.int()?.to_string();
                let width = streams.get(position, TOKEN_DZLEN)?.byte()? as usize;
                pad(digits.into_bytes(), width)
            }
            TOKEN_DELTA | TOKEN_DELTA0 => {
                let prior = previous(decoded, position)?;
                let value: i32 = std::str::from_utf8(&prior)
                    .ok()
                    .and_then(|s| s.parse().ok())
                    .ok_or_else(|| {
                        NxError::Cram(format!(
                            "The token in the prior name must be of type {}",
                            if token_type == TOKEN_DELTA {
                                "DIGITS or DELTA"
                            } else {
                                "DIGITS0 or DELTA0"
                            }
                        ))
                    })?;
                let delta = i64::from(streams.get(position, token_type)?.byte()?);
                let text = (i64::from(value) + delta).to_string().into_bytes();
                if token_type == TOKEN_DELTA0 {
                    pad(text, prior.len())
                } else {
                    text
                }
            }
            TOKEN_MATCH => previous(decoded, position)?,
            TOKEN_END | TOKEN_NOP => Vec::new(),
            other => {
                return Err(NxError::Cram(format!(
                    "Invalid tokenType : {other}. tokenType must be one of the valid token types"
                )))
            }
        };
        tokens.push(token);
        if token_type == TOKEN_END {
            break;
        }
        position += 1;
    }
    let name = tokens.concat();
    decoded.push(Rc::new(tokens));
    Ok(name)
}

/// `leftPadWith0`.
#[cfg(test)]
fn pad(value: Vec<u8>, width: usize) -> Vec<u8> {
    if value.len() >= width {
        return value;
    }
    let mut out = vec![b'0'; width - value.len()];
    out.extend_from_slice(&value);
    out
}

// ---- Encoding: `NameTokenisationEncode` ----

use crate::rans_nx16::{
    write_uint7, CAT_FLAG_MASK, ORDER_FLAG_MASK, PACK_FLAG_MASK, RLE_FLAG_MASK, STRIPE_FLAG_MASK,
};

/// The rANS flag sets tried per token type, htslib's tok3 level -3; the smallest result wins.
const RANS_FLAG_SETS: [&[u8]; TOTAL_TOKEN_TYPES] = [
    &[PACK_FLAG_MASK | RLE_FLAG_MASK, 0],
    &[PACK_FLAG_MASK | ORDER_FLAG_MASK, ORDER_FLAG_MASK, 0],
    &[0],
    &[PACK_FLAG_MASK | STRIPE_FLAG_MASK, 0],
    &[0],
    &[PACK_FLAG_MASK | RLE_FLAG_MASK | STRIPE_FLAG_MASK],
    &[PACK_FLAG_MASK | STRIPE_FLAG_MASK],
    &[PACK_FLAG_MASK | RLE_FLAG_MASK | STRIPE_FLAG_MASK],
    &[0],
    &[PACK_FLAG_MASK],
    &[0],
    &[0],
    &[0],
];

/// `EncodeToken`: the type, the token's own text, and what is written for it.
#[derive(Clone)]
struct EncodeToken {
    kind: u8,
    actual: Vec<u8>,
    /// The text written for a STRING or CHAR, or the number for the numeric types.
    relative: Relative,
}

#[derive(Clone)]
enum Relative {
    None,
    Text(Vec<u8>),
    Int(i32),
}

impl EncodeToken {
    fn relative_int(&self) -> i32 {
        match &self.relative {
            Relative::Int(v) => *v,
            Relative::Text(t) => std::str::from_utf8(t)
                .ok()
                .and_then(|s| s.parse().ok())
                .unwrap_or(0),
            Relative::None => 0,
        }
    }

    fn relative_text(&self) -> Vec<u8> {
        match &self.relative {
            Relative::Text(t) => t.clone(),
            Relative::Int(v) => v.to_string().into_bytes(),
            Relative::None => Vec::new(),
        }
    }
}

fn is_alphanumeric(c: u8) -> bool {
    c.is_ascii_alphanumeric()
}

fn is_all_digits(s: &[u8]) -> bool {
    !s.is_empty() && s.iter().all(u8::is_ascii_digit)
}

fn is_all_digits_leading_zero(s: &[u8]) -> bool {
    !s.is_empty() && s[0] == b'0' && s[1..].iter().all(u8::is_ascii_digit)
}

fn parse_int(s: &[u8]) -> i32 {
    std::str::from_utf8(s)
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(0)
}

/// `NameTokenisationEncode.compress`, with rANS rather than the arithmetic coder as htsjdk's
/// external compressor asks: the names in `input`, each ended by `separator`.
pub fn compress(input: &[u8], separator: u8) -> Result<Vec<u8>, NxError> {
    let names: Vec<&[u8]> = {
        let mut names = Vec::new();
        let mut start = 0;
        for (i, &b) in input.iter().enumerate() {
            if b == separator {
                names.push(&input[start..i]);
                start = i + 1;
            }
        }
        names
    };
    if names.is_empty() {
        return Err(NxError::Cram(
            "Name tokenizer input format requires a separator-delimited name list. No delimited names found in input."
                .to_string(),
        ));
    }
    let mut out = Vec::new();
    out.extend_from_slice(&(input.len() as i32).to_le_bytes());
    out.extend_from_slice(&(names.len() as i32).to_le_bytes());
    out.push(0);

    let mut index: std::collections::HashMap<&[u8], usize> = std::collections::HashMap::new();
    let mut frequencies = [0i32; 256];
    let mut tokens: Vec<Vec<EncodeToken>> = Vec::with_capacity(names.len());
    let mut max_positions = 0usize;
    for (n, name) in names.iter().enumerate() {
        let encoded = tokenise(name, n, &tokens, &mut index, &mut frequencies);
        max_positions = max_positions.max(encoded.len());
        tokens.push(encoded);
    }

    let mut registry: Vec<(Vec<u8>, usize, u8)> = Vec::new();
    for position in 0..max_positions {
        let streams = distribute(&tokens, position);
        serialize(&streams, &mut out, position, &mut registry);
    }
    Ok(out)
}

fn tokenise<'a>(
    name: &'a [u8],
    n: usize,
    previous: &[Vec<EncodeToken>],
    index: &mut std::collections::HashMap<&'a [u8], usize>,
    frequencies: &mut [i32; 256],
) -> Vec<EncodeToken> {
    if let Some(&first) = index.get(name) {
        return vec![EncodeToken {
            kind: TOKEN_DUP,
            actual: Vec::new(),
            relative: Relative::Int((n - first) as i32),
        }];
    }
    let mut encoded = vec![EncodeToken {
        kind: TOKEN_DIFF,
        actual: Vec::new(),
        relative: Relative::Int(if n == 0 { 0 } else { 1 }),
    }];
    index.insert(name, n);
    let prev = n.checked_sub(1).map(|p| &previous[p]);

    let mut pos = 0;
    let mut i = 1usize;
    while pos < name.len() {
        let start = pos;
        if is_alphanumeric(name[pos]) {
            let limit = name.len().min(start + 9);
            loop {
                pos += 1;
                if !(pos < limit && is_alphanumeric(name[pos])) {
                    break;
                }
            }
        } else {
            loop {
                pos += 1;
                if !(pos < name.len() && !is_alphanumeric(name[pos])) {
                    break;
                }
            }
        }
        let fragment = &name[start..pos];
        let mut kind = TOKEN_STRING;
        let mut relative = Relative::Text(fragment.to_vec());
        if is_all_digits_leading_zero(fragment) {
            kind = TOKEN_DIGITS0;
        } else if is_all_digits(fragment) {
            kind = TOKEN_DIGITS;
        } else if fragment.len() == 1 {
            kind = TOKEN_CHAR;
        }
        let prev_token = prev.filter(|p| p.len() > i + 1).map(|p| &p[i]);
        if let Some(prev_token) = prev_token.filter(|t| t.kind != TOKEN_END) {
            if prev_token.actual == fragment {
                kind = TOKEN_MATCH;
                relative = Relative::None;
            } else if kind == TOKEN_DIGITS
                && (prev_token.kind == TOKEN_DIGITS || prev_token.kind == TOKEN_DELTA)
            {
                let d = parse_int(fragment).wrapping_sub(parse_int(&prev_token.actual));
                frequencies[i] += 1;
                if (0..256).contains(&d) && frequencies[i] > (n / 2) as i32 {
                    encoded.push(EncodeToken {
                        kind: TOKEN_DELTA,
                        actual: fragment.to_vec(),
                        relative: Relative::Int(d),
                    });
                    i += 1;
                    continue;
                }
            } else if kind == TOKEN_DIGITS0
                && prev_token.actual.len() == fragment.len()
                && (prev_token.kind == TOKEN_DIGITS0 || prev_token.kind == TOKEN_DELTA0)
            {
                let d = parse_int(fragment).wrapping_sub(parse_int(&prev_token.actual));
                frequencies[i] += 1;
                if (0..256).contains(&d) && frequencies[i] > (n / 2) as i32 {
                    encoded.push(EncodeToken {
                        kind: TOKEN_DELTA0,
                        actual: fragment.to_vec(),
                        relative: Relative::Int(d),
                    });
                    i += 1;
                    continue;
                }
            }
        }
        encoded.push(EncodeToken {
            kind,
            actual: fragment.to_vec(),
            relative,
        });
        i += 1;
    }
    encoded.push(EncodeToken {
        kind: TOKEN_END,
        actual: Vec::new(),
        relative: Relative::None,
    });
    encoded
}

/// `distributeTokensForPosition`: one byte stream per token type, `None` where nothing was put.
fn distribute(tokens: &[Vec<EncodeToken>], position: usize) -> Vec<Option<Vec<u8>>> {
    let mut streams: Vec<Option<Vec<u8>>> = vec![None; TOTAL_TOKEN_TYPES];
    streams[TOKEN_TYPE as usize] = Some(Vec::new());
    for name in tokens {
        if position > 0 && name[0].kind == TOKEN_DUP {
            continue;
        }
        let Some(token) = name.get(position) else {
            continue;
        };
        streams[TOKEN_TYPE as usize]
            .as_mut()
            .unwrap()
            .push(token.kind);
        let mut put = |kind: u8, bytes: &[u8]| {
            streams[kind as usize]
                .get_or_insert_with(Vec::new)
                .extend_from_slice(bytes);
        };
        match token.kind {
            TOKEN_DIFF | TOKEN_DUP | TOKEN_DIGITS => {
                put(token.kind, &token.relative_int().to_le_bytes())
            }
            TOKEN_STRING => {
                let mut text: Vec<u8> = token
                    .relative_text()
                    .iter()
                    .map(|&b| if b.is_ascii() { b } else { b'?' })
                    .collect();
                text.push(0);
                put(TOKEN_STRING, &text);
            }
            TOKEN_CHAR => put(TOKEN_CHAR, &token.relative_text()[..1]),
            TOKEN_DIGITS0 => {
                put(TOKEN_DIGITS0, &token.relative_int().to_le_bytes());
                put(TOKEN_DZLEN, &[token.relative_text().len() as u8]);
            }
            TOKEN_DELTA | TOKEN_DELTA0 => put(token.kind, &[token.relative_int() as u8]),
            _ => {}
        }
    }
    streams
}

/// `tryCompress`.
fn try_compress(stream: &[u8], kind: u8) -> Vec<u8> {
    let mut best: Option<Vec<u8>> = None;
    for &flags in RANS_FLAG_SETS[kind as usize] {
        if flags & ORDER_FLAG_MASK != 0 && stream.len() < 100 {
            continue;
        }
        if flags & STRIPE_FLAG_MASK != 0 && !stream.len().is_multiple_of(4) {
            continue;
        }
        let compressed = crate::rans_nx16::compress(stream, flags);
        if best.as_ref().is_none_or(|b| b.len() > compressed.len()) {
            best = Some(compressed);
        }
    }
    match best {
        Some(b) if b.len() <= stream.len() => b,
        _ => crate::rans_nx16::compress(stream, CAT_FLAG_MASK),
    }
}

/// `serializeTokenStreams`.
fn serialize(
    streams: &[Option<Vec<u8>>],
    out: &mut Vec<u8>,
    position: usize,
    registry: &mut Vec<(Vec<u8>, usize, u8)>,
) {
    let mut omit_type = false;
    if let Some(types) = &streams[TOKEN_TYPE as usize] {
        if types.len() > 1 && types[1..].iter().all(|&t| t == TOKEN_MATCH) {
            omit_type = (1..=TOKEN_END as usize)
                .any(|t| streams[t].as_ref().is_some_and(|s| !s.is_empty()));
        }
    }
    let mut first = true;
    for kind in 0..=TOKEN_END {
        if omit_type && kind == TOKEN_TYPE {
            continue;
        }
        let Some(bytes) = streams[kind as usize].as_ref().filter(|s| !s.is_empty()) else {
            continue;
        };
        let mut header = kind;
        if first {
            header |= NEW_POSITION_FLAG_MASK;
            first = false;
        }
        let compressed = try_compress(bytes, kind);
        let duplicate = (compressed.len() > 4)
            .then(|| registry.iter().find(|(b, _, _)| *b == compressed))
            .flatten();
        if let Some((_, at, of)) = duplicate {
            out.push(header | DUP_PREVIOUS_STREAM_FLAG_MASK);
            out.push(*at as u8);
            out.push(*of);
        } else {
            out.push(header);
            write_uint7(compressed.len() as i32, out);
            out.extend_from_slice(&compressed);
            registry.push((compressed, position, kind));
        }
    }
}

/// Switches [`uncompress`] to the decoder as first ported, for the tests that compare the two.
#[cfg(test)]
mod reference {
    thread_local! {
        pub(super) static ON: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    }
}

#[cfg(test)]
mod decode_tests {
    use super::*;

    /// Each side's answer, `None` where it panicked.
    type Answer = Option<Result<Vec<u8>, NxError>>;

    fn both(stream: &[u8]) -> (Answer, Answer) {
        let run = || std::panic::catch_unwind(|| uncompress(stream, 0)).ok();
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

    /// Names of every token type the encoder writes: Illumina-like with deltas, zero-padded
    /// numbers, duplicates, signs and strings.
    fn names(random: &mut Lcg, count: usize) -> Vec<u8> {
        let mut out = Vec::new();
        let mut x = 1000u64;
        for i in 0..count {
            x += random.next() % 40;
            let name = match random.next() % 6 {
                0 | 1 => format!(
                    "HWI-ST{}:8:1101:{}:{}#0/1",
                    100 + i % 3,
                    x,
                    2000 + random.next() % 900
                ),
                2 => format!("SRR{:08}.{}", 62 + i / 50, i),
                3 => "dup-name".to_string(),
                4 => format!("r-{}_{:05}", random.next() % 7, x % 100_000),
                _ => format!("m{}/ccs", x * 1_234_567),
            };
            out.extend_from_slice(name.as_bytes());
            out.push(0);
        }
        out
    }

    /// Numbers are written as `pad(value.to_string(), width)` writes them, minus sign and all:
    /// a damaged block can make a DELTA0 of a negative STRING token.
    #[test]
    fn numbers_are_padded_as_the_first_port_pads_them() {
        for value in [
            0i64,
            7,
            -7,
            42,
            -42,
            1 << 31,
            -(1 << 31) - 255,
            i64::from(u32::MAX),
        ] {
            for width in 0..14 {
                let mut out = b"x".to_vec();
                push_decimal(&mut out, value, width);
                assert_eq!(
                    out[1..],
                    pad(value.to_string().into_bytes(), width),
                    "{value} {width}"
                );
            }
        }
    }

    /// The decoder against the first port's, on whole blocks, every prefix of a short one, cuts
    /// of long ones and blocks with one bit changed: the same names or the same refusal each time.
    #[test]
    fn the_decoder_answers_as_the_first_port_does() {
        let mut random = Lcg(5);
        let mut compared = 0;
        for count in [3usize, 40, 2000] {
            let text = names(&mut random, count);
            let stream = compress(&text, 0).expect("compresses");
            let (fast, first) = both(&stream);
            assert_eq!(fast, Some(Ok(text.clone())));
            assert_eq!(fast, first);
            let cuts: Vec<usize> = if stream.len() < 800 {
                (0..stream.len()).collect()
            } else {
                (0..300)
                    .map(|_| (random.next() as usize) % stream.len())
                    .collect()
            };
            for cut in cuts {
                let (fast, first) = both(&stream[..cut]);
                assert_eq!(fast, first, "{count} names, cut at {cut}");
                compared += 1;
            }
            for _ in 0..1500 {
                let mut damaged = stream.clone();
                // Past the name count, which would only change what is allocated.
                let at = 9 + (random.next() as usize) % (damaged.len() - 9);
                damaged[at] ^= 1 << (random.next() % 8);
                let (fast, first) = both(&damaged);
                assert_eq!(fast, first, "{count} names, bit changed at {at}");
                compared += 1;
            }
        }
        assert!(compared > 4000);
    }
}

#[cfg(test)]
mod encode_tests {
    use super::*;

    #[test]
    fn names_round_trip() {
        let mut input = Vec::new();
        for i in 0..500 {
            let name = match i % 5 {
                0 => format!("RUN1:frag{:02}", i / 3),
                1 => format!("read{:04}", i),
                2 => format!(
                    "HWI-ST{}:8:1101:{}:{}#0/1",
                    100 + i % 3,
                    1000 + i * 7,
                    2000 + i
                ),
                3 => "dup-name".to_string(),
                _ => format!("x{}", i * 1234567),
            };
            input.extend_from_slice(name.as_bytes());
            input.push(0);
        }
        let packed = compress(&input, 0).unwrap();
        assert_eq!(uncompress(&packed, 0).unwrap(), input);
    }
}
