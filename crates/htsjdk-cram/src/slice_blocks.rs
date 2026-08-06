//! The core block and the external blocks of one slice, in the order they are written.
//!
//! Ported from `htsjdk.samtools.cram.structure.SliceBlocks` at htsjdk 4.2.0.
//!
//! A slice's data is one core block of bits and any number of external blocks of bytes, each named
//! by a content id. [`crate::external_codecs::SliceBlockBytes`] is what the codecs fill; this is
//! how it reaches a file and comes back.
//!
//! # The order is by content id, not by insertion
//!
//! The externals live in a `TreeMap`, so writing is the core block first and then ascending content
//! id, whatever order they were added in. Measured: added 3, 2, 1 and written 1, 2, 3; added 300,
//! 2, 128 and written 2, 128, 300.
//!
//! # The reader takes a count, not an order
//!
//! It reads that many blocks and sorts them by type as it goes, so a stream whose core block comes
//! last is read exactly as one whose core block comes first. The check that there was a core block
//! at all happens after every block has been read.
//!
//! # Everything else is fatal
//!
//! A duplicate content id, a block that is neither core nor external, and a stream with no core
//! block. The first message names the id and the type of both the new block and the one already
//! there, and mentions the compression header encoding map, which is not what it is talking about.

use std::collections::BTreeMap;

use crate::block::{BlockHeader, ContentType};
use crate::compression_header::raw_block_with_id;
use crate::external_codecs::SliceBlockBytes;
use crate::varint::RuntimeEof;

/// What reading or building a slice's blocks refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SliceBlocksError {
    /// Two external blocks with the same content id.
    DuplicateContentId { content_id: i32 },
    /// A stream that carried no core block, checked after all of them were read.
    NoCoreBlock,
    /// A block that is neither core nor external.
    NotASliceBlock { name: String },
    /// The stream ended inside a block.
    Truncated,
}

impl SliceBlocksError {
    pub fn message(&self) -> String {
        match self {
            // The reference's message says "compression header encoding map", which is not what
            // this is. Kept as it is written.
            SliceBlocksError::DuplicateContentId { content_id } => format!(
                "Attempt to add a duplicate block (id {content_id} of type EXTERNAL) to \
                 compression header encoding map. Existing block is of type EXTERNAL."
            ),
            SliceBlocksError::NoCoreBlock => {
                "A core block is required in a CRAM stream but none was found.".to_string()
            }
            SliceBlocksError::NotASliceBlock { name } => {
                format!("Not a slice block, content type id {name}")
            }
            SliceBlocksError::Truncated => "null".to_string(),
        }
    }

    pub fn java_exception(&self) -> &'static str {
        match self {
            SliceBlocksError::DuplicateContentId { .. } | SliceBlocksError::NoCoreBlock => {
                "CRAMException"
            }
            SliceBlocksError::NotASliceBlock { .. } => "RuntimeException",
            SliceBlocksError::Truncated => "RuntimeEOFException",
        }
    }
}

impl From<RuntimeEof> for SliceBlocksError {
    fn from(_: RuntimeEof) -> Self {
        SliceBlocksError::Truncated
    }
}

/// `writeBlocks`: the core block, then the externals in ascending content id.
///
/// Everything is written raw here, as the measurement wrote it. What compressor a real writer
/// chooses for an external block is the compression header's business and not this one's.
pub fn write_blocks(blocks: &SliceBlockBytes, major: u8) -> Vec<u8> {
    let mut out = raw_block_with_id(ContentType::Core, 0, &blocks.core, major);
    for (content_id, content) in &blocks.external {
        out.extend_from_slice(&raw_block_with_id(
            ContentType::External,
            *content_id,
            content,
            major,
        ));
    }
    out
}

/// The `SliceBlocks(cramVersion, numberOfBlocks, inputStream)` constructor.
///
/// Reads exactly `count` blocks, sorting them by type as it goes, and only then asks whether one
/// of them was the core block.
pub fn read_blocks(
    bytes: &[u8],
    count: i32,
    major: u8,
) -> Result<SliceBlockBytes, SliceBlocksError> {
    let mut core: Option<Vec<u8>> = None;
    let mut external: BTreeMap<i32, Vec<u8>> = BTreeMap::new();
    let mut at = 0usize;

    for _ in 0..count.max(0) {
        let header = BlockHeader::read(bytes.get(at..).unwrap_or(&[]))?;
        let start = at + header.byte_length;
        let end = start + header.compressed_size.max(0) as usize;
        let content = bytes
            .get(start..end)
            .ok_or(SliceBlocksError::Truncated)?
            .to_vec();

        match header.content_type() {
            Some(ContentType::Core) => core = Some(content),
            Some(ContentType::External) => {
                if external.insert(header.content_id, content).is_some() {
                    return Err(SliceBlocksError::DuplicateContentId {
                        content_id: header.content_id,
                    });
                }
            }
            found => {
                return Err(SliceBlocksError::NotASliceBlock {
                    name: content_type_name(found, header.content_type),
                })
            }
        }

        at += header.total_length(major);
    }

    Ok(SliceBlockBytes {
        core: core.ok_or(SliceBlocksError::NoCoreBlock)?,
        external,
    })
}

/// The `SliceBlocks(coreBlock, externalBlocks)` constructor: a core block and a list of externals
/// in whatever order the caller has them.
///
/// The order is lost here, which is the point: it is not part of the file. A repeated content id is
/// refused rather than overwritten.
pub fn from_blocks(
    core: Vec<u8>,
    externals: Vec<(i32, Vec<u8>)>,
) -> Result<SliceBlockBytes, SliceBlocksError> {
    let mut external: BTreeMap<i32, Vec<u8>> = BTreeMap::new();
    for (content_id, content) in externals {
        if external.insert(content_id, content).is_some() {
            return Err(SliceBlocksError::DuplicateContentId { content_id });
        }
    }
    Ok(SliceBlockBytes { core, external })
}

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
