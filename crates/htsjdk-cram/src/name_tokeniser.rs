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
pub fn uncompress(bytes: &[u8], separator: u8) -> Result<Vec<u8>, NxError> {
    let mut input = Input::new(bytes);
    let header = input.take(9)?;
    let length = u32::from_le_bytes([header[0], header[1], header[2], header[3]]) as usize;
    let names = u32::from_le_bytes([header[4], header[5], header[6], header[7]]) as usize;
    let use_arith = header[8] != 0;
    let mut streams = TokenStreams::read(&mut input, use_arith, names)?;

    let mut decoded: Vec<Rc<Vec<Vec<u8>>>> = Vec::with_capacity(names);
    let mut out = Vec::with_capacity(length);
    for i in 0..names {
        let name = decode_name(&mut streams, &mut decoded, i)?;
        out.extend_from_slice(&name);
        out.push(separator);
    }
    Ok(out)
}

fn decode_name(
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
fn pad(value: Vec<u8>, width: usize) -> Vec<u8> {
    if value.len() >= width {
        return value;
    }
    let mut out = vec![b'0'; width - value.len()];
    out.extend_from_slice(&value);
    out
}
