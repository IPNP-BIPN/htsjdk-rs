//! Reading one value through whichever codec an encoding names.
//!
//! [`crate::encoding_factory`] turns an identifier and its parameters into an [`Encoding`]. This is
//! what then reads through it, and it is the only place where the core bit stream and the external
//! blocks are reached from the same call.
//!
//! The three shapes are the reference's: `read()` for a single value, and `read(length)` for a byte
//! array whose length came from somewhere else. A codec that has no `read(length)` refuses here
//! exactly as it does there.

use crate::core_codecs::{
    BetaIntegerCodec, CodecError, GammaIntegerCodec, SubexponentialIntegerCodec,
};
use crate::encoding_factory::Encoding;
use crate::external_codecs::{
    ByteArrayStopCodec, ExternalByteArrayCodec, ExternalByteCodec, ExternalError,
    ExternalIntegerCodec, ExternalLongCodec, SliceReadStreams, SliceWriteStreams,
};
use crate::golomb::{GolombError, GolombIntegerCodec, GolombLongCodec, GolombRiceIntegerCodec};
use crate::huffman::{CanonicalHuffman, HuffmanError};

/// What reading through an encoding refuses. Each variant carries the refusal of the codec that
/// raised it, so the class and the message are the reference's own.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ReadError {
    /// An encoding asked for a value of a kind it cannot produce, which the factory's fall-through
    /// makes reachable.
    WrongKind {
        encoding: &'static str,
    },
    Core(CodecError),
    External(ExternalError),
    Golomb(GolombError),
    Huffman(HuffmanError),
}

impl ReadError {
    pub fn message(&self) -> String {
        match self {
            ReadError::WrongKind { encoding } => {
                format!("{encoding} cannot produce a value of this kind")
            }
            ReadError::Core(error) => error.message(),
            ReadError::External(error) => error.message(),
            ReadError::Golomb(error) => error.message(),
            ReadError::Huffman(error) => error.message(),
        }
    }

    pub fn java_exception(&self) -> &'static str {
        match self {
            ReadError::WrongKind { .. } => "ClassCastException",
            ReadError::Core(error) => error.java_exception(),
            ReadError::External(error) => error.java_exception(),
            ReadError::Golomb(error) => error.java_exception(),
            ReadError::Huffman(error) => error.java_exception(),
        }
    }
}

impl From<CodecError> for ReadError {
    fn from(error: CodecError) -> Self {
        ReadError::Core(error)
    }
}

impl From<ExternalError> for ReadError {
    fn from(error: ExternalError) -> Self {
        ReadError::External(error)
    }
}

impl From<GolombError> for ReadError {
    fn from(error: GolombError) -> Self {
        ReadError::Golomb(error)
    }
}

impl From<HuffmanError> for ReadError {
    fn from(error: HuffmanError) -> Self {
        ReadError::Huffman(error)
    }
}

/// One integer, through whichever of the seven integer codecs the encoding names.
pub fn read_int(encoding: &Encoding, streams: &mut SliceReadStreams<'_>) -> Result<i32, ReadError> {
    Ok(match encoding {
        Encoding::ExternalInteger(content_id) => {
            ExternalIntegerCodec::new(*content_id).read(streams)?
        }
        Encoding::HuffmanInteger(params) => {
            CanonicalHuffman::new(&params.symbols, &params.bit_lengths)?.read(streams.core())?
        }
        Encoding::Beta {
            offset,
            bits_per_value,
        } => BetaIntegerCodec::new(*offset, *bits_per_value).read(streams.core())?,
        Encoding::Gamma { offset } => GammaIntegerCodec::new(*offset).read(streams.core())?,
        Encoding::Subexponential { offset, k } => {
            SubexponentialIntegerCodec::new(*offset, *k).read(streams.core())?
        }
        Encoding::Golomb { offset, m } => {
            GolombIntegerCodec::new(*offset, *m)?.read(streams.core())?
        }
        Encoding::GolombRice { offset, m } => {
            GolombRiceIntegerCodec::new(*offset, *m).read(streams.core())?
        }
        Encoding::GolombLong { offset, m } => {
            GolombLongCodec::new(i64::from(*offset), *m)?.read(streams.core())? as i32
        }
        other => {
            return Err(ReadError::WrongKind {
                encoding: other.java_class(),
            })
        }
    })
}

/// One long, which only two encodings produce.
pub fn read_long(
    encoding: &Encoding,
    streams: &mut SliceReadStreams<'_>,
) -> Result<i64, ReadError> {
    Ok(match encoding {
        Encoding::ExternalLong(content_id) => ExternalLongCodec::new(*content_id).read(streams)?,
        Encoding::GolombLong { offset, m } => {
            GolombLongCodec::new(i64::from(*offset), *m)?.read(streams.core())?
        }
        other => {
            return Err(ReadError::WrongKind {
                encoding: other.java_class(),
            })
        }
    })
}

/// One byte. External byte cannot fail, which is why its arm has no `?`.
pub fn read_byte(encoding: &Encoding, streams: &mut SliceReadStreams<'_>) -> Result<i8, ReadError> {
    Ok(match encoding {
        Encoding::ExternalByte(content_id) => ExternalByteCodec::new(*content_id).read(streams),
        Encoding::HuffmanByte(params) => {
            CanonicalHuffman::new(&params.symbols, &params.bit_lengths)?.read(streams.core())?
        }
        other => {
            return Err(ReadError::WrongKind {
                encoding: other.java_class(),
            })
        }
    })
}

/// A byte array. `length` is `None` where the encoding finds its own end, as byte-array-stop does,
/// and `Some` where the caller already knows it.
pub fn read_byte_array(
    encoding: &Encoding,
    streams: &mut SliceReadStreams<'_>,
    length: Option<usize>,
) -> Result<Vec<u8>, ReadError> {
    Ok(match (encoding, length) {
        (Encoding::ByteArrayStop { stop, content_id }, None) => {
            ByteArrayStopCodec::new(*stop, *content_id).read(streams)
        }
        (Encoding::ByteArrayStop { stop, content_id }, Some(length)) => {
            ByteArrayStopCodec::new(*stop, *content_id).read_with_length(streams, length)?
        }
        (Encoding::ExternalByteArray(content_id), Some(length)) => {
            ExternalByteArrayCodec::new(*content_id).read_with_length(streams, length)?
        }
        (Encoding::ExternalByteArray(content_id), None) => {
            ExternalByteArrayCodec::new(*content_id).read(streams)?
        }
        (Encoding::ByteArrayLen(length_encoding, bytes_encoding), _) => {
            let read_length = read_int(length_encoding, streams)?;
            read_byte_array(bytes_encoding, streams, Some(read_length.max(0) as usize))?
        }
        (other, _) => {
            return Err(ReadError::WrongKind {
                encoding: other.java_class(),
            })
        }
    })
}

/// One integer out, through whichever of the seven integer codecs the encoding names.
pub fn write_int(
    encoding: &Encoding,
    streams: &mut SliceWriteStreams,
    value: i32,
) -> Result<(), ReadError> {
    match encoding {
        Encoding::ExternalInteger(content_id) => {
            ExternalIntegerCodec::new(*content_id).write(streams, value)
        }
        Encoding::HuffmanInteger(params) => {
            CanonicalHuffman::new(&params.symbols, &params.bit_lengths)?
                .write(streams.core(), value)?;
        }
        Encoding::Beta {
            offset,
            bits_per_value,
        } => BetaIntegerCodec::new(*offset, *bits_per_value).write(streams.core(), value)?,
        Encoding::Gamma { offset } => {
            GammaIntegerCodec::new(*offset).write(streams.core(), value)?
        }
        Encoding::Subexponential { offset, k } => {
            SubexponentialIntegerCodec::new(*offset, *k).write(streams.core(), value)?
        }
        Encoding::Golomb { offset, m } => {
            GolombIntegerCodec::new(*offset, *m)?.write(streams.core(), value)?
        }
        Encoding::GolombRice { offset, m } => {
            GolombRiceIntegerCodec::new(*offset, *m).write(streams.core(), value)?
        }
        Encoding::GolombLong { offset, m } => {
            GolombLongCodec::new(i64::from(*offset), *m)?.write(streams.core(), i64::from(value))?
        }
        other => {
            return Err(ReadError::WrongKind {
                encoding: other.java_class(),
            })
        }
    }
    Ok(())
}

/// One byte out.
pub fn write_byte(
    encoding: &Encoding,
    streams: &mut SliceWriteStreams,
    value: i8,
) -> Result<(), ReadError> {
    match encoding {
        Encoding::ExternalByte(content_id) => {
            ExternalByteCodec::new(*content_id).write(streams, value)
        }
        Encoding::HuffmanByte(params) => {
            CanonicalHuffman::new(&params.symbols, &params.bit_lengths)?
                .write(streams.core(), value)?;
        }
        other => {
            return Err(ReadError::WrongKind {
                encoding: other.java_class(),
            })
        }
    }
    Ok(())
}

/// A byte array out. `ByteArrayLen` writes its length through one codec and its bytes through
/// another, which is the one place a single value reaches two blocks.
pub fn write_byte_array(
    encoding: &Encoding,
    streams: &mut SliceWriteStreams,
    value: &[u8],
) -> Result<(), ReadError> {
    match encoding {
        Encoding::ByteArrayStop { stop, content_id } => {
            ByteArrayStopCodec::new(*stop, *content_id).write(streams, value)
        }
        Encoding::ExternalByteArray(content_id) => {
            ExternalByteArrayCodec::new(*content_id).write(streams, value)
        }
        Encoding::ByteArrayLen(length_encoding, bytes_encoding) => {
            write_int(length_encoding, streams, value.len() as i32)?;
            write_byte_array(bytes_encoding, streams, value)?;
        }
        other => {
            return Err(ReadError::WrongKind {
                encoding: other.java_class(),
            })
        }
    }
    Ok(())
}
