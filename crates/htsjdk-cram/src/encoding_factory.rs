//! Which encoding a data series type and an identifier resolve to.
//!
//! Ported from `htsjdk.samtools.cram.encoding.EncodingFactory` at htsjdk 4.2.0.
//!
//! A compression header names an encoding per data series, by identifier and parameters. This is
//! what turns that pair into something that can read bytes, and it is the last thing between a
//! file and the codecs.
//!
//! # The switch falls through
//!
//! Only the `BYTE` arm of the reference's switch ends in a `break`. An `INT` that matches nothing
//! falls into the `LONG` arm and then into the `BYTE_ARRAY` arm; a `LONG` that matches nothing
//! falls into `BYTE_ARRAY`. So an `INT` data series named with `BYTE_ARRAY_LEN` gets a byte array
//! encoding rather than the refusal the method's last line promises, and the port reproduces that
//! rather than the intent.
//!
//! # The parameters are not checked against the type
//!
//! Whatever arm is reached parses the bytes its own way, so the same parameters mean a content id
//! in one and the head of a Huffman alphabet in another. Nothing rejects the mismatch.
//!
//! # `NULL` is an identifier like any other
//!
//! It matches nothing anywhere, so it always reaches the refusal, whose message names both halves
//! of what was asked for.

use crate::encoding_map::EncodingId;
use crate::golomb;
use crate::huffman::{parse_byte_params, parse_integer_params, HuffmanParams};
use crate::varint::{read_unsigned_itf8, RuntimeEof};

/// `DataSeriesType`: the kind of value a data series holds.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DataSeriesType {
    Byte,
    Int,
    Long,
    ByteArray,
}

impl DataSeriesType {
    pub fn name(&self) -> &'static str {
        match self {
            DataSeriesType::Byte => "BYTE",
            DataSeriesType::Int => "INT",
            DataSeriesType::Long => "LONG",
            DataSeriesType::ByteArray => "BYTE_ARRAY",
        }
    }
}

/// What the factory refuses: a pair it has no arm for.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FactoryError {
    NotFound {
        value_type: DataSeriesType,
        encoding_id: EncodingId,
    },
    /// Parameters that ran out before the encoding had read all of them.
    Parameters,
}

impl FactoryError {
    pub fn message(&self) -> String {
        match self {
            FactoryError::NotFound {
                value_type,
                encoding_id,
            } => format!(
                "Encoding not found: value type={}, encoding id={}",
                value_type.name(),
                encoding_id.name()
            ),
            FactoryError::Parameters => "null".to_string(),
        }
    }

    pub fn java_exception(&self) -> &'static str {
        match self {
            FactoryError::NotFound { .. } => "IllegalArgumentException",
            FactoryError::Parameters => "RuntimeEOFException",
        }
    }
}

impl From<RuntimeEof> for FactoryError {
    fn from(_: RuntimeEof) -> Self {
        FactoryError::Parameters
    }
}

/// One encoding, built from an identifier and its parameters.
///
/// The variants are the reference's classes, and [`Encoding::java_class`] names them, because what
/// the factory returns for a given pair is the whole of what it does.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Encoding {
    ExternalByte(i32),
    ExternalInteger(i32),
    ExternalLong(i32),
    ExternalByteArray(i32),
    HuffmanByte(HuffmanParams<i8>),
    HuffmanInteger(HuffmanParams<i32>),
    Golomb { offset: i32, m: i32 },
    GolombLong { offset: i32, m: i32 },
    GolombRice { offset: i32, m: i32 },
    Beta { offset: i32, bits_per_value: i32 },
    Gamma { offset: i32 },
    Subexponential { offset: i32, k: i32 },
    ByteArrayStop { stop: u8, content_id: i32 },
    ByteArrayLen(Box<Encoding>, Box<Encoding>),
}

impl Encoding {
    /// The reference's class name, which is what a measurement of the factory records.
    pub fn java_class(&self) -> &'static str {
        match self {
            Encoding::ExternalByte(_) => "ExternalByteEncoding",
            Encoding::ExternalInteger(_) => "ExternalIntegerEncoding",
            Encoding::ExternalLong(_) => "ExternalLongEncoding",
            Encoding::ExternalByteArray(_) => "ExternalByteArrayEncoding",
            Encoding::HuffmanByte(_) => "CanonicalHuffmanByteEncoding",
            Encoding::HuffmanInteger(_) => "CanonicalHuffmanIntegerEncoding",
            Encoding::Golomb { .. } => "GolombIntegerEncoding",
            Encoding::GolombLong { .. } => "GolombLongEncoding",
            Encoding::GolombRice { .. } => "GolombRiceIntegerEncoding",
            Encoding::Beta { .. } => "BetaIntegerEncoding",
            Encoding::Gamma { .. } => "GammaIntegerEncoding",
            Encoding::Subexponential { .. } => "SubexponentialIntegerEncoding",
            Encoding::ByteArrayStop { .. } => "ByteArrayStopEncoding",
            Encoding::ByteArrayLen(_, _) => "ByteArrayLenEncoding",
        }
    }

    /// The reference's `toString`, which a nested encoding prints inside its parent's.
    pub fn describe(&self) -> String {
        match self {
            Encoding::ExternalByte(id)
            | Encoding::ExternalInteger(id)
            | Encoding::ExternalLong(id)
            | Encoding::ExternalByteArray(id) => format!("Content ID: {id}"),
            Encoding::HuffmanByte(params) => format!(
                "Symbols: {} BitLengths {}",
                join(&params.symbols),
                join(&params.bit_lengths)
            ),
            Encoding::HuffmanInteger(params) => format!(
                "Symbols: {} BitLengths {}",
                join(&params.symbols),
                join(&params.bit_lengths)
            ),
            Encoding::Golomb { offset, m }
            | Encoding::GolombLong { offset, m }
            | Encoding::GolombRice { offset, m } => format!("Offset: {offset} m: {m}"),
            Encoding::Beta {
                offset,
                bits_per_value,
            } => format!("Offset: {offset} BitsPerValue: {bits_per_value}"),
            Encoding::Gamma { offset } => format!("Offset: {offset}"),
            Encoding::Subexponential { offset, k } => format!("Offset: {offset} k: {k}"),
            Encoding::ByteArrayStop { stop, content_id } => {
                format!("Content ID: {content_id} StopByte: {stop}")
            }
            Encoding::ByteArrayLen(length, bytes) => format!(
                "LenEncoding: {} ByteEncoding: {}",
                length.describe(),
                bytes.describe()
            ),
        }
    }
}

/// The reference joins these with a semicolon, which only shows on an alphabet of more than one.
fn join<T: std::fmt::Display>(values: &[T]) -> String {
    values
        .iter()
        .map(|value| value.to_string())
        .collect::<Vec<_>>()
        .join(";")
}

/// `EncodingFactory.createCRAMEncoding`.
///
/// The arms are the reference's, fall-through included: `Int` tries its own identifiers, then
/// `Long`'s, then `ByteArray`'s; `Long` tries its own and then `ByteArray`'s; `Byte` tries only its
/// own, because that is the one arm with a `break`.
pub fn create_encoding(
    value_type: DataSeriesType,
    encoding_id: EncodingId,
    params: &[u8],
) -> Result<Encoding, FactoryError> {
    let not_found = FactoryError::NotFound {
        value_type,
        encoding_id,
    };

    if value_type == DataSeriesType::Byte {
        return match encoding_id {
            EncodingId::External => Ok(Encoding::ExternalByte(content_id(params)?)),
            EncodingId::Huffman => Ok(Encoding::HuffmanByte(parse_byte_params(params)?)),
            _ => Err(not_found),
        };
    }

    if value_type == DataSeriesType::Int {
        match encoding_id {
            EncodingId::Huffman => {
                return Ok(Encoding::HuffmanInteger(parse_integer_params(params)?))
            }
            EncodingId::External => return Ok(Encoding::ExternalInteger(content_id(params)?)),
            EncodingId::Golomb => {
                let (offset, m) = golomb::parse_params(params)?;
                return Ok(Encoding::Golomb { offset, m });
            }
            EncodingId::GolombRice => {
                let (offset, m) = golomb::parse_params(params)?;
                return Ok(Encoding::GolombRice { offset, m });
            }
            EncodingId::Beta => {
                let (offset, bits_per_value) = two_itf8(params)?;
                return Ok(Encoding::Beta {
                    offset,
                    bits_per_value,
                });
            }
            EncodingId::Gamma => {
                let (offset, _) = read_unsigned_itf8(params)?;
                return Ok(Encoding::Gamma { offset });
            }
            EncodingId::Subexponential => {
                let (offset, k) = two_itf8(params)?;
                return Ok(Encoding::Subexponential { offset, k });
            }
            // Everything else falls into the LONG arm below, because the INT arm has no break.
            _ => {}
        }
    }

    if value_type == DataSeriesType::Int || value_type == DataSeriesType::Long {
        match encoding_id {
            EncodingId::Golomb => {
                let (offset, m) = golomb::parse_params(params)?;
                return Ok(Encoding::GolombLong { offset, m });
            }
            EncodingId::External => return Ok(Encoding::ExternalLong(content_id(params)?)),
            // And on into the BYTE_ARRAY arm.
            _ => {}
        }
    }

    match encoding_id {
        EncodingId::ByteArrayLen => byte_array_len(params),
        EncodingId::ByteArrayStop => {
            let stop = *params.first().ok_or(FactoryError::Parameters)?;
            let (content_id, _) = read_unsigned_itf8(params.get(1..).unwrap_or(&[]))?;
            Ok(Encoding::ByteArrayStop { stop, content_id })
        }
        EncodingId::External => Ok(Encoding::ExternalByteArray(content_id(params)?)),
        _ => Err(not_found),
    }
}

/// `ByteArrayLenEncoding.fromSerializedEncodingParams`: each half as its identifier, the length of
/// its parameters, and then those parameters.
///
/// Both halves go back through the factory, the length as an `INT` and the bytes as a
/// `BYTE_ARRAY`, so the fall-through above applies to them too.
fn byte_array_len(params: &[u8]) -> Result<Encoding, FactoryError> {
    let mut cursor = 0usize;
    let next = |cursor: &mut usize| -> Result<i32, FactoryError> {
        let (value, used) = read_unsigned_itf8(params.get(*cursor..).unwrap_or(&[]))?;
        *cursor += used;
        Ok(value)
    };

    let length_id = next(&mut cursor)?;
    let length_size = next(&mut cursor)? as usize;
    let length_params = params
        .get(cursor..cursor + length_size)
        .ok_or(FactoryError::Parameters)?
        .to_vec();
    cursor += length_size;

    let bytes_id = next(&mut cursor)?;
    let bytes_size = next(&mut cursor)? as usize;
    let bytes_params = params
        .get(cursor..cursor + bytes_size)
        .ok_or(FactoryError::Parameters)?
        .to_vec();

    let length = create_encoding(
        DataSeriesType::Int,
        EncodingId::from_id(length_id).ok_or(FactoryError::Parameters)?,
        &length_params,
    )?;
    let bytes = create_encoding(
        DataSeriesType::ByteArray,
        EncodingId::from_id(bytes_id).ok_or(FactoryError::Parameters)?,
        &bytes_params,
    )?;
    Ok(Encoding::ByteArrayLen(Box::new(length), Box::new(bytes)))
}

fn content_id(params: &[u8]) -> Result<i32, RuntimeEof> {
    Ok(read_unsigned_itf8(params)?.0)
}

fn two_itf8(params: &[u8]) -> Result<(i32, i32), RuntimeEof> {
    let (first, used) = read_unsigned_itf8(params)?;
    let (second, _) = read_unsigned_itf8(params.get(used..).unwrap_or(&[]))?;
    Ok((first, second))
}
