//! A whole CRAM file, from its first byte to the blocks its records are read from.
//!
//! Every piece has a module and a suite of its own: [`crate::container`] for the file definition
//! and the container headers, [`crate::compression_header`] for the three maps,
//! [`crate::slice_header`], [`crate::slice_blocks`], [`crate::record_read`]. This is the walk that
//! puts them in order, and it is the only thing that shows whether they fit.
//!
//! # The first container is not like the others
//!
//! It holds the SAM header, in a `FILE_HEADER` block rather than a compression header, and a
//! reader that treats it as an ordinary container refuses it with the compression header's own
//! message. That is the first thing composing the pieces shows.
//!
//! # A container's blocks are compressed and the compression header's is not
//!
//! The compression header block is raw by definition; a slice's blocks carry whichever compressor
//! the writer chose, so reading records means undoing several different ones in the same container.
//!
//! # The EOF container is a container
//!
//! It parses as one whose record count is zero, and a reader that stops on a byte pattern rather
//! than on that count is reading the wrong thing.

use crate::block::{BlockHeader, CompressionMethod, ContentType};
use crate::compression_header::{CompressionHeader, CompressionHeaderError};
use crate::container::{ContainerError, ContainerHeader, FileDefinition, FILE_DEFINITION_LENGTH};
use crate::external_codecs::SliceBlockBytes;
use crate::rans;
use crate::slice_header::SliceHeader;
use crate::varint::RuntimeEof;

/// What walking a file refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FileError {
    Container(ContainerError),
    CompressionHeader(CompressionHeaderError),
    /// A block whose compressor this crate does not undo. bzip2 and LZMA are legal in a CRAM and
    /// are not ported: a file using either is refused rather than read wrong.
    UnsupportedCompression {
        method: i32,
    },
    /// A block that did not decompress.
    Corrupt,
    Truncated,
}

impl FileError {
    pub fn message(&self) -> String {
        match self {
            FileError::Container(error) => error.message(),
            FileError::CompressionHeader(error) => error.message(),
            FileError::UnsupportedCompression { method } => {
                format!("no decompressor for block compression method {method}")
            }
            FileError::Corrupt => "a block did not decompress".to_string(),
            FileError::Truncated => "null".to_string(),
        }
    }
}

impl From<ContainerError> for FileError {
    fn from(error: ContainerError) -> Self {
        FileError::Container(error)
    }
}

impl From<CompressionHeaderError> for FileError {
    fn from(error: CompressionHeaderError) -> Self {
        FileError::CompressionHeader(error)
    }
}

impl From<RuntimeEof> for FileError {
    fn from(_: RuntimeEof) -> Self {
        FileError::Truncated
    }
}

/// One block of a slice, after its compression is undone.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DecompressedBlock {
    pub content_type: ContentType,
    pub method: CompressionMethod,
    pub content_id: i32,
    pub content: Vec<u8>,
}

/// One slice of a container: its header, and its blocks decompressed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SliceWalk {
    pub header: SliceHeader,
    pub blocks: Vec<DecompressedBlock>,
}

impl SliceWalk {
    /// The blocks as the codecs want them: the core one and the externals by content id.
    pub fn block_bytes(&self) -> SliceBlockBytes {
        let mut bytes = SliceBlockBytes::default();
        for block in &self.blocks {
            match block.content_type {
                ContentType::Core => bytes.core = block.content.clone(),
                _ => {
                    bytes
                        .external
                        .insert(block.content_id, block.content.clone());
                }
            }
        }
        bytes
    }
}

/// One container: its header, its compression header, and its slices.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ContainerWalk {
    pub header: ContainerHeader,
    /// The byte offset the container starts at, which is what a CRAI entry points to.
    pub offset: usize,
    pub compression_header: Option<CompressionHeader>,
    pub slices: Vec<SliceWalk>,
}

/// A whole file, walked.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FileWalk {
    pub definition: FileDefinition,
    /// The SAM header container's text, exactly as the block held it.
    pub sam_header: Vec<u8>,
    pub containers: Vec<ContainerWalk>,
}

/// Walk a file: the definition, the SAM header container, then every container until the EOF one.
pub fn read_file(cram: &[u8]) -> Result<FileWalk, FileError> {
    let definition = FileDefinition::read(cram)?;
    let major = definition.major;
    let mut at = FILE_DEFINITION_LENGTH;

    // The SAM header container, whose one block is a FILE_HEADER rather than a compression header.
    let header_container = ContainerHeader::read(&cram[at..], major)?;
    let blocks_start = at + header_container.byte_length;
    let sam_header = first_block_content(cram, blocks_start, major)?;
    at = blocks_start + header_container.blocks_byte_size.max(0) as usize;

    let mut containers = Vec::new();
    while at < cram.len() {
        let offset = at;
        let header = ContainerHeader::read(&cram[at..], major)?;
        let blocks_start = at + header.byte_length;
        let end = blocks_start + header.blocks_byte_size.max(0) as usize;

        if header.is_eof() {
            containers.push(ContainerWalk {
                header,
                offset,
                compression_header: None,
                slices: Vec::new(),
            });
            break;
        }

        // The compression header block, then one slice header block and its blocks per slice.
        let compression_header = CompressionHeader::read_block(&cram[blocks_start..end], major)?;
        let mut cursor = blocks_start + block_length(cram, blocks_start, major)?;

        let mut slices = Vec::new();
        while cursor < end {
            let slice_header_block = BlockHeader::read(&cram[cursor..])?;
            let content_start = cursor + slice_header_block.byte_length;
            let content_end = content_start + slice_header_block.compressed_size.max(0) as usize;
            let slice_header = SliceHeader::read(
                cram.get(content_start..content_end)
                    .ok_or(FileError::Truncated)?,
            )?;
            cursor += slice_header_block.total_length(major);

            let mut blocks = Vec::new();
            // The header's block count is the core block plus the externals, and it does not
            // count the slice header block itself.
            for _ in 0..slice_header.block_count.max(0) {
                let header = BlockHeader::read(cram.get(cursor..).ok_or(FileError::Truncated)?)?;
                let start = cursor + header.byte_length;
                let end = start + header.compressed_size.max(0) as usize;
                let raw = cram.get(start..end).ok_or(FileError::Truncated)?;
                let method = header.method().ok_or(FileError::UnsupportedCompression {
                    method: header.method,
                })?;
                blocks.push(DecompressedBlock {
                    content_type: header.content_type().unwrap_or(ContentType::External),
                    method,
                    content_id: header.content_id,
                    content: decompress(method, raw, header.uncompressed_size)?,
                });
                cursor += header.total_length(major);
            }

            slices.push(SliceWalk {
                header: slice_header,
                blocks,
            });
        }

        containers.push(ContainerWalk {
            header,
            offset,
            compression_header: Some(compression_header),
            slices,
        });
        at = end;
    }

    Ok(FileWalk {
        definition,
        sam_header,
        containers,
    })
}

/// How many bytes the block at `at` occupies, header, content and checksum together.
fn block_length(cram: &[u8], at: usize, major: u8) -> Result<usize, FileError> {
    let header = BlockHeader::read(cram.get(at..).ok_or(FileError::Truncated)?)?;
    Ok(header.total_length(major))
}

/// The content of the block at `at`, decompressed.
fn first_block_content(cram: &[u8], at: usize, major: u8) -> Result<Vec<u8>, FileError> {
    let header = BlockHeader::read(cram.get(at..).ok_or(FileError::Truncated)?)?;
    let start = at + header.byte_length;
    let end = start + header.compressed_size.max(0) as usize;
    let raw = cram.get(start..end).ok_or(FileError::Truncated)?;
    let method = header.method().ok_or(FileError::UnsupportedCompression {
        method: header.method,
    })?;
    let _ = major;
    decompress(method, raw, header.uncompressed_size)
}

/// Undo a block's compression.
///
/// Raw, GZIP and rANS are the three a htsjdk-written CRAM uses. bzip2 and LZMA are legal and are
/// not ported: a file using either is refused here rather than read wrong.
pub fn decompress(
    method: CompressionMethod,
    content: &[u8],
    uncompressed_size: i32,
) -> Result<Vec<u8>, FileError> {
    match method {
        CompressionMethod::Raw => Ok(content.to_vec()),
        CompressionMethod::Gzip => {
            use std::io::Read;
            let mut out = Vec::with_capacity(uncompressed_size.max(0) as usize);
            flate2::read::GzDecoder::new(content)
                .read_to_end(&mut out)
                .map_err(|_| FileError::Corrupt)?;
            Ok(out)
        }
        CompressionMethod::Rans => rans::uncompress(content).map_err(|_| FileError::Corrupt),
        other => Err(FileError::UnsupportedCompression {
            method: other as i32,
        }),
    }
}
