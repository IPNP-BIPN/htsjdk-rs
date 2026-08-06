//! The codecs written on external blocks rather than on the core bit stream.
//!
//! Ported from `htsjdk.samtools.cram.encoding.external` and
//! `htsjdk.samtools.cram.encoding.ByteArrayLenCodec` at htsjdk 4.2.0.
//!
//! Six of the encoding map's identifiers put their data in a block of its own, named by a content
//! id, and two more are built out of the first six. A block is bytes rather than bits, so these
//! codecs carry no alignment problem and no offset. What they carry instead is a set of decisions
//! about where one value ends and the next begins.
//!
//! # Each codec names a block, and two codecs can name the same one
//!
//! That is what lets [`ByteArrayLenCodec`] put its lengths in one block and its bytes in another,
//! or in the same one, where the two interleave. Both were measured.
//!
//! # External byte cannot see the end of its block
//!
//! It returns `(byte) stream.read()`, and at the end of a `ByteArrayInputStream` that is
//! `(byte) -1`. A byte of `0xFF` that is really there reads back as `-1` too, so nothing
//! distinguishes a value from the end of the data.
//!
//! # Byte array stop trusts the data
//!
//! It appends a stop byte after each array and reads until it sees one, so an array containing
//! that byte is written whole and read back split in two, with nothing reporting it. And a block
//! that ends before a stop byte does is not an error either: the read ends on the end of the
//! stream exactly as it would on a separator, and every read after that returns an empty array.
//!
//! # Byte array len is a pair of codecs, but not any pair
//!
//! The length may be an external integer or a Huffman code on the core bit stream, so a single
//! value can straddle both kinds of block. But the bytes half is read through `read(length)`, and
//! [`ByteArrayStopCodec`] does not implement that: a byte-array-len wrapping a byte-array-stop
//! writes correctly and refuses on the way back.

use std::collections::BTreeMap;

use crate::bit_stream::{BitError, BitInputStream, BitOutputStream};
use crate::encoding_map::EncodingId;
use crate::huffman::{CanonicalHuffman, HuffmanError};
use crate::varint::{
    read_unsigned_itf8, read_unsigned_ltf8, write_unsigned_itf8, write_unsigned_ltf8,
};

/// What a codec on an external block refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ExternalError {
    /// `read(length)` on a codec that has no use for one, and `read()` past the end of a stopped
    /// array's own reach.
    NotImplemented,
    /// The byte array codec's `read()`, which has no length to work from.
    UnknownArrayLength,
    /// More bytes asked for than the block holds. The reference's message is `null`: the exception
    /// is constructed with no argument.
    EndOfStream,
    /// A Huffman length, which has its own refusals.
    Huffman(HuffmanError),
    /// The bit stream underneath a Huffman length.
    Bits(BitError),
}

impl ExternalError {
    pub fn message(&self) -> String {
        match self {
            ExternalError::NotImplemented => "Not implemented.".to_string(),
            ExternalError::UnknownArrayLength => {
                "Cannot read byte array of unknown length.".to_string()
            }
            ExternalError::EndOfStream => "null".to_string(),
            ExternalError::Huffman(error) => error.message(),
            ExternalError::Bits(error) => error.message(),
        }
    }

    pub fn java_exception(&self) -> &'static str {
        match self {
            ExternalError::NotImplemented | ExternalError::UnknownArrayLength => "RuntimeException",
            ExternalError::EndOfStream => "RuntimeEOFException",
            ExternalError::Huffman(error) => error.java_exception(),
            ExternalError::Bits(error) => error.java_exception(),
        }
    }
}

impl From<HuffmanError> for ExternalError {
    fn from(error: HuffmanError) -> Self {
        ExternalError::Huffman(error)
    }
}

impl From<BitError> for ExternalError {
    fn from(error: BitError) -> Self {
        ExternalError::Bits(error)
    }
}

/// The blocks a slice's write left behind: the core bit stream, and one byte block per content id
/// that was written to.
///
/// The reference creates a stream for every data series its compression header knows, so most of
/// the blocks a single codec produces are empty. Only the ones written to are kept here.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SliceBlockBytes {
    pub core: Vec<u8>,
    pub external: BTreeMap<i32, Vec<u8>>,
}

/// `SliceBlocksWriteStreams`: the core bit stream and the external blocks, handed out by content
/// id so that two codecs naming the same id share one.
#[derive(Debug, Default)]
pub struct SliceWriteStreams {
    core: BitOutputStream,
    external: BTreeMap<i32, Vec<u8>>,
}

impl SliceWriteStreams {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn core(&mut self) -> &mut BitOutputStream {
        &mut self.core
    }

    pub fn external(&mut self, content_id: i32) -> &mut Vec<u8> {
        self.external.entry(content_id).or_default()
    }

    /// `flushStreamsToBlocks`, which is where the core bit stream's partial byte is padded out.
    pub fn finish(self) -> SliceBlockBytes {
        SliceBlockBytes {
            core: self.core.into_bytes(),
            external: self.external,
        }
    }
}

/// A block being read, which reports the end of itself as `-1` rather than refusing.
#[derive(Debug, Clone)]
struct ByteCursor<'a> {
    bytes: &'a [u8],
    at: usize,
}

impl<'a> ByteCursor<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, at: 0 }
    }

    /// `ByteArrayInputStream.read()`.
    fn read(&mut self) -> i32 {
        match self.bytes.get(self.at) {
            Some(byte) => {
                self.at += 1;
                i32::from(*byte)
            }
            None => -1,
        }
    }

    fn remaining(&self) -> &'a [u8] {
        self.bytes.get(self.at..).unwrap_or(&[])
    }

    fn advance(&mut self, count: usize) {
        self.at += count;
    }

    /// `InputStreamUtils.readFully`, which refuses rather than returning what it has.
    fn read_fully(&mut self, length: usize) -> Result<&'a [u8], ExternalError> {
        let taken = self
            .remaining()
            .get(..length)
            .ok_or(ExternalError::EndOfStream)?;
        self.advance(length);
        Ok(taken)
    }
}

/// `SliceBlocksReadStreams`, over the blocks a write produced.
#[derive(Debug)]
pub struct SliceReadStreams<'a> {
    core: BitInputStream<'a>,
    external: BTreeMap<i32, ByteCursor<'a>>,
}

impl<'a> SliceReadStreams<'a> {
    pub fn new(blocks: &'a SliceBlockBytes) -> Self {
        Self {
            core: BitInputStream::new(&blocks.core),
            external: blocks
                .external
                .iter()
                .map(|(id, bytes)| (*id, ByteCursor::new(bytes)))
                .collect(),
        }
    }

    pub fn core(&mut self) -> &mut BitInputStream<'a> {
        &mut self.core
    }

    /// A block nothing was written to still reads, as an empty one.
    fn external(&mut self, content_id: i32) -> &mut ByteCursor<'a> {
        self.external
            .entry(content_id)
            .or_insert_with(|| ByteCursor::new(&[]))
    }
}

/// `ExternalIntegerCodec`: ITF8 straight onto the block.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ExternalIntegerCodec {
    pub content_id: i32,
}

impl ExternalIntegerCodec {
    pub fn new(content_id: i32) -> Self {
        Self { content_id }
    }

    pub fn write(&self, streams: &mut SliceWriteStreams, value: i32) {
        let (bytes, _) = write_unsigned_itf8(value);
        streams.external(self.content_id).extend_from_slice(&bytes);
    }

    pub fn read(&self, streams: &mut SliceReadStreams<'_>) -> Result<i32, ExternalError> {
        let cursor = streams.external(self.content_id);
        // Past the end the reference's stream returns -1 for every byte and the arithmetic makes a
        // number of it, which is what `read_unsigned_itf8` reproduces.
        let (value, used) =
            read_unsigned_itf8(cursor.remaining()).map_err(|_| ExternalError::EndOfStream)?;
        cursor.advance(used);
        Ok(value)
    }
}

/// `ExternalLongCodec`: LTF8, which has three more widths than ITF8.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ExternalLongCodec {
    pub content_id: i32,
}

impl ExternalLongCodec {
    pub fn new(content_id: i32) -> Self {
        Self { content_id }
    }

    pub fn write(&self, streams: &mut SliceWriteStreams, value: i64) {
        let (bytes, _) = write_unsigned_ltf8(value);
        streams.external(self.content_id).extend_from_slice(&bytes);
    }

    pub fn read(&self, streams: &mut SliceReadStreams<'_>) -> Result<i64, ExternalError> {
        let cursor = streams.external(self.content_id);
        let (value, used) =
            read_unsigned_ltf8(cursor.remaining()).map_err(|_| ExternalError::EndOfStream)?;
        cursor.advance(used);
        Ok(value)
    }
}

/// `ExternalByteCodec`: one byte per value, and no way to see the end of the block.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ExternalByteCodec {
    pub content_id: i32,
}

impl ExternalByteCodec {
    pub fn new(content_id: i32) -> Self {
        Self { content_id }
    }

    pub fn write(&self, streams: &mut SliceWriteStreams, value: i8) {
        streams.external(self.content_id).push(value as u8);
    }

    /// `(byte) inputStream.read()`. At the end that is `(byte) -1`, which a byte of `0xFF` that is
    /// really there produces too. The signature has no error in it because the reference has none.
    pub fn read(&self, streams: &mut SliceReadStreams<'_>) -> i8 {
        streams.external(self.content_id).read() as i8
    }
}

/// `ExternalByteArrayCodec`: bytes with no length of their own.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ExternalByteArrayCodec {
    pub content_id: i32,
}

impl ExternalByteArrayCodec {
    pub fn new(content_id: i32) -> Self {
        Self { content_id }
    }

    pub fn write(&self, streams: &mut SliceWriteStreams, value: &[u8]) {
        streams.external(self.content_id).extend_from_slice(value);
    }

    /// `read(length)`, which refuses rather than returning a short array.
    pub fn read_with_length(
        &self,
        streams: &mut SliceReadStreams<'_>,
        length: usize,
    ) -> Result<Vec<u8>, ExternalError> {
        Ok(streams
            .external(self.content_id)
            .read_fully(length)?
            .to_vec())
    }

    /// `read()`, which has no length to work from and says so.
    pub fn read(&self, _streams: &mut SliceReadStreams<'_>) -> Result<Vec<u8>, ExternalError> {
        Err(ExternalError::UnknownArrayLength)
    }
}

/// `ByteArrayStopCodec`: a separator the data is trusted not to contain.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ByteArrayStopCodec {
    pub stop: u8,
    pub content_id: i32,
}

impl ByteArrayStopCodec {
    pub fn new(stop: u8, content_id: i32) -> Self {
        Self { stop, content_id }
    }

    pub fn write(&self, streams: &mut SliceWriteStreams, value: &[u8]) {
        let block = streams.external(self.content_id);
        block.extend_from_slice(value);
        block.push(self.stop);
    }

    /// Reads until the stop byte or the end of the block, whichever comes first, and reports
    /// neither. An array holding the stop byte comes back split; a block that runs out comes back
    /// short, and then empty for ever after.
    pub fn read(&self, streams: &mut SliceReadStreams<'_>) -> Vec<u8> {
        let cursor = streams.external(self.content_id);
        let mut out = Vec::new();
        loop {
            let byte = cursor.read();
            if byte == -1 || byte == i32::from(self.stop) {
                return out;
            }
            out.push(byte as u8);
        }
    }

    pub fn read_with_length(
        &self,
        _streams: &mut SliceReadStreams<'_>,
        _length: usize,
    ) -> Result<Vec<u8>, ExternalError> {
        Err(ExternalError::NotImplemented)
    }
}

/// The half of a [`ByteArrayLenCodec`] that carries the length: an external integer, or a Huffman
/// code on the core bit stream.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LengthCodec {
    External(ExternalIntegerCodec),
    Huffman(CanonicalHuffman<i32>),
}

impl LengthCodec {
    fn write(&self, streams: &mut SliceWriteStreams, value: i32) -> Result<(), ExternalError> {
        match self {
            LengthCodec::External(codec) => {
                codec.write(streams, value);
                Ok(())
            }
            LengthCodec::Huffman(huffman) => {
                huffman.write(streams.core(), value)?;
                Ok(())
            }
        }
    }

    fn read(&self, streams: &mut SliceReadStreams<'_>) -> Result<i32, ExternalError> {
        match self {
            LengthCodec::External(codec) => codec.read(streams),
            LengthCodec::Huffman(huffman) => Ok(huffman.read(streams.core())?),
        }
    }
}

/// The half of a [`ByteArrayLenCodec`] that carries the bytes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BytesCodec {
    External(ExternalByteArrayCodec),
    Stop(ByteArrayStopCodec),
}

impl BytesCodec {
    fn write(&self, streams: &mut SliceWriteStreams, value: &[u8]) {
        match self {
            BytesCodec::External(codec) => codec.write(streams, value),
            BytesCodec::Stop(codec) => codec.write(streams, value),
        }
    }

    fn read_with_length(
        &self,
        streams: &mut SliceReadStreams<'_>,
        length: usize,
    ) -> Result<Vec<u8>, ExternalError> {
        match self {
            BytesCodec::External(codec) => codec.read_with_length(streams, length),
            BytesCodec::Stop(codec) => codec.read_with_length(streams, length),
        }
    }
}

/// `ByteArrayLenCodec`: a length codec and a byte codec, which need not share a block or even a
/// kind of block.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ByteArrayLenCodec {
    pub length: LengthCodec,
    pub bytes: BytesCodec,
}

impl ByteArrayLenCodec {
    pub fn new(length: LengthCodec, bytes: BytesCodec) -> Self {
        Self { length, bytes }
    }

    pub fn write(
        &self,
        streams: &mut SliceWriteStreams,
        value: &[u8],
    ) -> Result<(), ExternalError> {
        self.length.write(streams, value.len() as i32)?;
        self.bytes.write(streams, value);
        Ok(())
    }

    /// The length first, then that many bytes. A stop codec in the bytes half writes correctly and
    /// refuses here, because `read(length)` is what this calls and the stop codec has none.
    pub fn read(&self, streams: &mut SliceReadStreams<'_>) -> Result<Vec<u8>, ExternalError> {
        let length = self.length.read(streams)?;
        self.bytes.read_with_length(streams, length as usize)
    }
}

/// `ExternalEncoding.toSerializedEncodingParams`: the content id, ITF8.
pub fn serialize_external_params(content_id: i32) -> Vec<u8> {
    write_unsigned_itf8(content_id).0
}

/// `ByteArrayStopEncoding.toSerializedEncodingParams`: the stop byte raw, then the content id.
pub fn serialize_stop_params(stop: u8, content_id: i32) -> Vec<u8> {
    let mut out = vec![stop];
    out.extend_from_slice(&write_unsigned_itf8(content_id).0);
    out
}

/// `ByteArrayLenEncoding.toSerializedEncodingParams`: each half as its identifier, the length of
/// its parameters, and then those parameters.
pub fn serialize_byte_array_len_params(
    length: (EncodingId, &[u8]),
    bytes: (EncodingId, &[u8]),
) -> Vec<u8> {
    let mut out = Vec::new();
    for (id, params) in [length, bytes] {
        out.extend_from_slice(&write_unsigned_itf8(id as i32).0);
        out.extend_from_slice(&write_unsigned_itf8(params.len() as i32).0);
        out.extend_from_slice(params);
    }
    out
}
