//! The CRAM block: what a container's byte size actually counts.
//!
//! Ported from `htsjdk.samtools.cram.structure.block.Block`, `BlockCompressionMethod` and
//! `BlockContentType` at htsjdk 4.2.0.
//!
//! [`crate::container`] gives a size and a block count; this is what that size counts. Each block
//! is a five-field header, then the compressed content, then a CRC32 from version 3 up.
//!
//! # The CRC covers the header and the content together
//!
//! `Block.read` wraps a `CRC32InputStream` around the **whole** read, so the checksum is over the
//! five header fields *and* the content, not the content alone. A block therefore cannot be
//! verified without re-reading its own header, and the four checksum bytes sit **outside** the
//! `compressedSize` the header declares. It is absent below version 3, exactly as the container
//! header's is, so the two version-dependent lengths compound.
//!
//! # A content id is only legal on an external block
//!
//! `Block`'s constructor throws `Cannot set a Content ID for non-external blocks` when the id is
//! set on anything else, so the field is present in every block and meaningful in one kind.
//!
//! # The first two blocks of every CRAM have fixed methods
//!
//! Not a per-file choice: `createGZIPFileHeaderBlock` and `createRawCompressionHeaderBlock` fix
//! them. Measured on five files, every one begins with a **GZIP** `FILE_HEADER` block in the first
//! container and a **RAW** `COMPRESSION_HEADER` block in the second.
//!
//! # An ordinary file uses rANS
//!
//! Measured, the methods present in a four-read CRAM are RAW, GZIP and **rANS**, over 29 blocks.
//! That is what makes rANS 4x8 required rather than optional, and it is the evidence behind
//! decision 0038's scoping: the codecs a file can actually contain are the ones that must be
//! ported.

use crate::varint::{read_unsigned_itf8, RuntimeEof};

/// `BlockCompressionMethod`, by id because that is what the byte carries.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CompressionMethod {
    Raw = 0,
    Gzip = 1,
    Bzip2 = 2,
    Lzma = 3,
    Rans = 4,
    Range = 5,
}

impl CompressionMethod {
    /// `BlockCompressionMethod.byId`.
    pub fn from_id(id: i32) -> Option<Self> {
        Some(match id {
            0 => CompressionMethod::Raw,
            1 => CompressionMethod::Gzip,
            2 => CompressionMethod::Bzip2,
            3 => CompressionMethod::Lzma,
            4 => CompressionMethod::Rans,
            5 => CompressionMethod::Range,
            _ => return None,
        })
    }
}

/// `BlockContentType`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContentType {
    FileHeader = 0,
    CompressionHeader = 1,
    MappedSlice = 2,
    Reserved = 3,
    External = 4,
    Core = 5,
}

impl ContentType {
    pub fn from_id(id: i32) -> Option<Self> {
        Some(match id {
            0 => ContentType::FileHeader,
            1 => ContentType::CompressionHeader,
            2 => ContentType::MappedSlice,
            3 => ContentType::Reserved,
            4 => ContentType::External,
            5 => ContentType::Core,
            _ => return None,
        })
    }
}

/// `Block.NO_CONTENT_ID`.
pub const NO_CONTENT_ID: i32 = 0;

/// One block's header, with the content left where it is.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BlockHeader {
    pub method: i32,
    pub content_type: i32,
    pub content_id: i32,
    pub compressed_size: i32,
    pub uncompressed_size: i32,
    /// The header's own byte length, which is where the content starts.
    pub byte_length: usize,
}

impl BlockHeader {
    /// The five fields at the head of a block: two raw bytes then three ITF8s.
    pub fn read(bytes: &[u8]) -> Result<Self, RuntimeEof> {
        let method = i32::from(*bytes.first().ok_or(RuntimeEof)?);
        let content_type = i32::from(*bytes.get(1).ok_or(RuntimeEof)?);
        let mut at = 2usize;
        let mut itf8 = || -> Result<i32, RuntimeEof> {
            let (value, consumed) = read_unsigned_itf8(&bytes[at.min(bytes.len())..])?;
            at += consumed;
            Ok(value)
        };
        let content_id = itf8()?;
        let compressed_size = itf8()?;
        let uncompressed_size = itf8()?;
        Ok(Self {
            method,
            content_type,
            content_id,
            compressed_size,
            uncompressed_size,
            byte_length: at,
        })
    }

    /// How many bytes this block occupies in the file, header, content and checksum together.
    ///
    /// The checksum is **outside** `compressed_size`, which is the part a port most easily gets
    /// wrong: `compressedSize` describes the content and nothing else.
    pub fn total_length(&self, major: u8) -> usize {
        self.byte_length + self.compressed_size.max(0) as usize + if major >= 3 { 4 } else { 0 }
    }

    pub fn method(&self) -> Option<CompressionMethod> {
        CompressionMethod::from_id(self.method)
    }

    pub fn content_type(&self) -> Option<ContentType> {
        ContentType::from_id(self.content_type)
    }
}

/// One block, located: its header, and where its content sits in the file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LocatedBlock {
    pub header: BlockHeader,
    /// The offset of the block's first header byte.
    pub offset: usize,
    /// The checksum, or zero below version 3 where there is none.
    pub checksum: i32,
}

impl LocatedBlock {
    /// The compressed content, borrowed from the file.
    pub fn content<'a>(&self, cram: &'a [u8]) -> &'a [u8] {
        let start = self.offset + self.header.byte_length;
        let end = start + self.header.compressed_size.max(0) as usize;
        &cram[start.min(cram.len())..end.min(cram.len())]
    }
}

/// Walk the blocks of one container, given where its blocks begin and how many there are.
pub fn blocks_of_container(
    cram: &[u8],
    blocks_start: usize,
    blocks_byte_size: i32,
    block_count: i32,
    major: u8,
) -> Result<Vec<LocatedBlock>, RuntimeEof> {
    let end = blocks_start + blocks_byte_size.max(0) as usize;
    let mut out = Vec::with_capacity(block_count.max(0) as usize);
    let mut at = blocks_start;
    for _ in 0..block_count.max(0) {
        if at >= end || at >= cram.len() {
            break;
        }
        let header = BlockHeader::read(&cram[at..])?;
        let checksum = if major >= 3 {
            let start = at + header.byte_length + header.compressed_size.max(0) as usize;
            read_int32(cram, start)
        } else {
            0
        };
        let total = header.total_length(major);
        out.push(LocatedBlock {
            header,
            offset: at,
            checksum,
        });
        at += total;
    }
    Ok(out)
}

/// `CramInt.readInt32` at an offset, little-endian, `-1` for bytes past the end.
fn read_int32(bytes: &[u8], at: usize) -> i32 {
    let byte = |offset: usize| -> i32 {
        match bytes.get(at + offset) {
            Some(value) => i32::from(*value),
            None => -1,
        }
    };
    byte(0) | (byte(1) << 8) | (byte(2) << 16) | (byte(3) << 24)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_ids_are_the_bytes_the_format_carries() {
        assert_eq!(CompressionMethod::from_id(0), Some(CompressionMethod::Raw));
        assert_eq!(CompressionMethod::from_id(4), Some(CompressionMethod::Rans));
        assert_eq!(CompressionMethod::from_id(6), None);
        assert_eq!(ContentType::from_id(4), Some(ContentType::External));
        assert_eq!(ContentType::from_id(6), None);
    }

    /// The checksum sits outside the declared content size, which is the arithmetic a port most
    /// easily gets wrong.
    #[test]
    fn the_checksum_is_outside_the_compressed_size() {
        let header = BlockHeader {
            method: 0,
            content_type: 4,
            content_id: 1,
            compressed_size: 10,
            uncompressed_size: 10,
            byte_length: 5,
        };
        assert_eq!(header.total_length(3), 5 + 10 + 4);
        assert_eq!(header.total_length(2), 5 + 10);
    }

    #[test]
    fn a_block_header_is_two_bytes_then_three_itf8s() {
        // method 4, type 4, id 1, compressed 36, uncompressed 4.
        let header = BlockHeader::read(&[0x04, 0x04, 0x01, 0x24, 0x04]).expect("parses");
        assert_eq!(header.method, 4);
        assert_eq!(header.content_type, 4);
        assert_eq!(header.content_id, 1);
        assert_eq!(header.compressed_size, 36);
        assert_eq!(header.uncompressed_size, 4);
        assert_eq!(header.byte_length, 5);
    }
}
