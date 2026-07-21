//! BGZF codec, byte-identical to htsjdk's `BlockCompressedOutputStream`.
//!
//! Ported from htsjdk 4.2.0, commit `4cc010022ac038fb30f26e6f9717fabff3e808c1`:
//!
//! - `src/main/java/htsjdk/samtools/util/BlockCompressedStreamConstants.java`
//! - `src/main/java/htsjdk/samtools/util/BlockCompressedOutputStream.java`
//!   (`deflateBlock`, `writeGzipBlock`, `flush`, `close`)
//!
//! The deflate backend must be zlib, never miniz_oxide. See
//! `docs/decisions/0001-deflate-backend.md`.

pub mod read;
pub mod vfp;
mod write;

pub use read::{decompress_all, BgzfError, BgzfReader, DecompressedBlock};
pub use write::BgzfWriter;

// Constants transcribed from BlockCompressedStreamConstants.
pub const BLOCK_HEADER_LENGTH: usize = 18;
pub const BLOCK_FOOTER_LENGTH: usize = 8;
pub const MAX_COMPRESSED_BLOCK_SIZE: usize = 64 * 1024;
pub const GZIP_OVERHEAD: usize = BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH + 2;
pub const NO_COMPRESSION_OVERHEAD: usize = 10;

/// Uncompressed bytes per block: 65498, which is *not* 64 KiB. htsjdk reserves the gzip and
/// no-compression overhead so that even an incompressible block still fits inside 64 KiB once
/// framed.
pub const DEFAULT_UNCOMPRESSED_BLOCK_SIZE: usize =
    64 * 1024 - (GZIP_OVERHEAD + NO_COMPRESSION_OVERHEAD);

/// Size of htsjdk's `compressedBuffer`, 65518, and therefore the threshold that decides
/// whether the no-compression fallback kicks in.
pub const COMPRESSED_BUFFER_SIZE: usize = MAX_COMPRESSED_BLOCK_SIZE - BLOCK_HEADER_LENGTH;

/// htsjdk `Defaults.COMPRESSION_LEVEL`, which is 5, not zlib's own default of 6.
pub const DEFAULT_COMPRESSION_LEVEL: u32 = 5;

/// The 28-byte empty block htsjdk appends as the BGZF end-of-file marker.
pub const EMPTY_GZIP_BLOCK: [u8; 28] = [
    31, 139, 8, 4, // ID1, ID2, CM=deflate, FLG=FEXTRA
    0, 0, 0, 0, // MTIME
    0, 255, // XFL, OS=unknown
    6, 0, // XLEN
    66, 67, // BGZF_ID1 'B', BGZF_ID2 'C'
    2, 0, // BGZF_LEN
    27, 0, // BSIZE - 1
    3, 0, // the empty deflate stream
    0, 0, 0, 0, // CRC32
    0, 0, 0, 0, // ISIZE
];

// Gzip/BGZF header field values, named as in BlockCompressedStreamConstants so the read and
// write paths can assert against the same constants rather than repeating magic numbers.
pub const GZIP_ID1: u8 = 31;
pub const GZIP_ID2: u8 = 139;
pub const GZIP_CM_DEFLATE: u8 = 8;
pub const GZIP_FLG: u8 = 4;
pub const GZIP_XFL: u8 = 0;
pub const GZIP_OS_UNKNOWN: u8 = 255;
pub const GZIP_XLEN: u16 = 6;
pub const BGZF_ID1: u8 = 66;
pub const BGZF_ID2: u8 = 67;
pub const BGZF_LEN: u16 = 2;
