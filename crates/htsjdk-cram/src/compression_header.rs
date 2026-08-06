//! The compression header: three length-prefixed maps inside one raw block.
//!
//! Ported from `htsjdk.samtools.cram.structure.CompressionHeader` at htsjdk 4.2.0.
//!
//! [`crate::preservation_map`], [`crate::encoding_map`] and [`crate::tag_encoding_map`] each cover
//! one of the three. This is what joins them, and it is the last structure between a container's
//! bytes and the codecs that read its slices.
//!
//! # A header read and written again is byte-identical
//!
//! Measured over both CRAM versions and four combinations of the three flags. That is the property
//! a byte-identical port needs, and the one place a dropped field shows without being looked for.
//!
//! # The version changes the block and nothing inside it
//!
//! A 3.0 block carries a four-byte CRC-32 that a 2.1 block does not, over the header and the
//! content together and written little-endian. It sits **outside** the compressed size, so the
//! same header is 178 bytes in one version and 174 in the other, identical up to those four.
//!
//! # The block is always raw
//!
//! Whatever the version, and it is read back as raw content rather than through a decompressor.
//! The content type is checked, and the refusal names what was found instead.

use crate::block::{BlockHeader, CompressionMethod, ContentType, NO_CONTENT_ID};
use crate::encoding_map::{EncodingMap, EncodingMapError};
use crate::preservation_map::{PreservationMap, PreservationMapError};
use crate::tag_encoding_map::{TagEncodingMap, TagEncodingMapError};
use crate::varint::{write_unsigned_itf8, RuntimeEof};

/// What reading a compression header refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CompressionHeaderError {
    /// A block of the wrong kind. The message names what was found, by the reference's own name
    /// for it.
    WrongBlockType {
        found: String,
    },
    /// The block ended before its content did.
    Truncated,
    Preservation(PreservationMapError),
    Encoding(EncodingMapError),
    TagEncoding(TagEncodingMapError),
}

impl CompressionHeaderError {
    pub fn message(&self) -> String {
        match self {
            CompressionHeaderError::WrongBlockType { found } => {
                format!("Compression header block expected, found: {found}")
            }
            CompressionHeaderError::Truncated => "null".to_string(),
            CompressionHeaderError::Preservation(error) => error.message(),
            CompressionHeaderError::Encoding(error) => error.message(),
            CompressionHeaderError::TagEncoding(error) => error.message(),
        }
    }

    pub fn java_exception(&self) -> &'static str {
        match self {
            CompressionHeaderError::WrongBlockType { .. } => "RuntimeIOException",
            CompressionHeaderError::Truncated => "RuntimeEOFException",
            CompressionHeaderError::Preservation(error) => error.java_exception(),
            CompressionHeaderError::Encoding(error) => error.java_exception(),
            CompressionHeaderError::TagEncoding(error) => error.java_exception(),
        }
    }
}

impl From<PreservationMapError> for CompressionHeaderError {
    fn from(error: PreservationMapError) -> Self {
        CompressionHeaderError::Preservation(error)
    }
}

impl From<EncodingMapError> for CompressionHeaderError {
    fn from(error: EncodingMapError) -> Self {
        CompressionHeaderError::Encoding(error)
    }
}

impl From<TagEncodingMapError> for CompressionHeaderError {
    fn from(error: TagEncodingMapError) -> Self {
        CompressionHeaderError::TagEncoding(error)
    }
}

impl From<RuntimeEof> for CompressionHeaderError {
    fn from(_: RuntimeEof) -> Self {
        CompressionHeaderError::Truncated
    }
}

/// The three maps, in the order the block carries them.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CompressionHeader {
    pub preservation: PreservationMap,
    pub encodings: EncodingMap,
    pub tag_encodings: TagEncodingMap,
}

impl CompressionHeader {
    /// `internalRead`: the three maps, each behind its own length prefix.
    pub fn read_content(content: &[u8]) -> Result<Self, CompressionHeaderError> {
        let (preservation, used) = PreservationMap::read_prefixed(content)?;
        let rest = content.get(used..).unwrap_or(&[]);
        let (encodings, used) = EncodingMap::read_prefixed(rest)?;
        let rest = rest.get(used..).unwrap_or(&[]);
        let (tag_encodings, _) = TagEncodingMap::read_prefixed(rest)?;
        Ok(Self {
            preservation,
            encodings,
            tag_encodings,
        })
    }

    /// `internalWrite`: the same three, in the same order.
    pub fn write_content(&self) -> Vec<u8> {
        let mut out = self.preservation.write_prefixed();
        out.extend_from_slice(&self.encodings.write_prefixed());
        out.extend_from_slice(&self.tag_encodings.write_prefixed());
        out
    }

    /// The whole block: a raw block header, the content, and the checksum a version 3 file carries.
    pub fn write_block(&self, major: u8) -> Vec<u8> {
        raw_block(ContentType::CompressionHeader, &self.write_content(), major)
    }

    /// Read the block and then its content. The content type is checked before anything is parsed,
    /// and the refusal names what was found.
    pub fn read_block(bytes: &[u8], major: u8) -> Result<Self, CompressionHeaderError> {
        let header = BlockHeader::read(bytes)?;
        match header.content_type() {
            Some(ContentType::CompressionHeader) => {}
            found => {
                return Err(CompressionHeaderError::WrongBlockType {
                    found: content_type_name(found, header.content_type),
                })
            }
        }
        let start = header.byte_length;
        let end = start + header.compressed_size.max(0) as usize;
        let content = bytes
            .get(start..end)
            .ok_or(CompressionHeaderError::Truncated)?;
        // The checksum is not verified here, because the reference does not verify it either: it
        // writes one and reads past it, and `total_length` is what accounts for its four bytes.
        let _ = major;
        Self::read_content(content)
    }
}

/// The reference's `BlockContentType` name, which is what its refusal prints.
fn content_type_name(content_type: Option<ContentType>, id: i32) -> String {
    match content_type {
        Some(ContentType::FileHeader) => "FILE_HEADER".to_string(),
        Some(ContentType::CompressionHeader) => "COMPRESSION_HEADER".to_string(),
        Some(ContentType::MappedSlice) => "MAPPED_SLICE".to_string(),
        Some(ContentType::Reserved) => "RESERVED".to_string(),
        Some(ContentType::External) => "EXTERNAL".to_string(),
        Some(ContentType::Core) => "CORE".to_string(),
        None => id.to_string(),
    }
}

/// A raw block: method, content type, three ITF8s, the content, and a CRC-32 from version 3 on.
///
/// The content id is zero and both sizes are the content's length, because nothing is compressed.
pub fn raw_block(content_type: ContentType, content: &[u8], major: u8) -> Vec<u8> {
    let mut out = Vec::with_capacity(content.len() + 16);
    out.push(CompressionMethod::Raw as u8);
    out.push(content_type as u8);
    out.extend_from_slice(&write_unsigned_itf8(NO_CONTENT_ID).0);
    out.extend_from_slice(&write_unsigned_itf8(content.len() as i32).0);
    out.extend_from_slice(&write_unsigned_itf8(content.len() as i32).0);
    out.extend_from_slice(content);
    if major >= 3 {
        out.extend_from_slice(&crc32(&out).to_le_bytes());
    }
    out
}

/// `java.util.zip.CRC32`, which is CRC-32/ISO-HDLC: the reflected polynomial `0xEDB88320`, an
/// all-ones seed and an all-ones final xor.
///
/// Written out rather than taken from a compression crate, because twenty lines of arithmetic is a
/// smaller thing to own than a dependency, and this is the only checksum CRAM uses.
pub fn crc32(bytes: &[u8]) -> u32 {
    let mut crc = 0xFFFF_FFFFu32;
    for byte in bytes {
        crc ^= u32::from(*byte);
        for _ in 0..8 {
            let carry = crc & 1;
            crc >>= 1;
            if carry != 0 {
                crc ^= 0xEDB8_8320;
            }
        }
    }
    !crc
}
