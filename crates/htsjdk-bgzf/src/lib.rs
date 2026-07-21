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

use std::io::{self, Write};

use flate2::{Compress, Compression, Crc, FlushCompress, Status};

// Constants transcribed from BlockCompressedStreamConstants.
pub const BLOCK_HEADER_LENGTH: usize = 18;
pub const BLOCK_FOOTER_LENGTH: usize = 8;
pub const MAX_COMPRESSED_BLOCK_SIZE: usize = 64 * 1024;
pub const GZIP_OVERHEAD: usize = BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH + 2;
pub const NO_COMPRESSION_OVERHEAD: usize = 10;

/// 65498. Note this is *not* 64 KiB: htsjdk reserves the gzip and no-compression overhead so
/// that even an incompressible block still fits inside 64 KiB once framed.
pub const DEFAULT_UNCOMPRESSED_BLOCK_SIZE: usize =
    64 * 1024 - (GZIP_OVERHEAD + NO_COMPRESSION_OVERHEAD);

/// 65518. The size of htsjdk's `compressedBuffer`, and therefore the threshold that decides
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

/// Writes a BGZF stream whose bytes match `BlockCompressedOutputStream`.
pub struct BgzfWriter<W: Write> {
    inner: W,
    buffer: Vec<u8>,
    level: u32,
    finished: bool,
}

impl<W: Write> BgzfWriter<W> {
    pub fn new(inner: W) -> Self {
        Self::with_level(inner, DEFAULT_COMPRESSION_LEVEL)
    }

    pub fn with_level(inner: W, level: u32) -> Self {
        assert!(level <= 9, "compression level must be 0..=9, got {level}");
        Self {
            inner,
            buffer: Vec::with_capacity(DEFAULT_UNCOMPRESSED_BLOCK_SIZE),
            level,
            finished: false,
        }
    }

    /// Ports `deflateBlock`. Returns the total framed block size, or 0 when nothing is buffered.
    ///
    /// The fallback condition is subtle and must mirror Java exactly. htsjdk deflates into a
    /// fixed `compressedBuffer` of [`COMPRESSED_BUFFER_SIZE`] bytes and then tests
    /// `deflater.finished()`. When the output *exactly* fills that buffer, zlib returns with
    /// `avail_out == 0` and cannot signal end-of-stream, so `finished()` is false and the
    /// no-compression path is taken even though the data did technically fit.
    ///
    /// Testing `compressed.len() > COMPRESSED_BUFFER_SIZE` instead is wrong at precisely that
    /// boundary, and incompressible payloads land on it: 65498 bytes of random data deflate to
    /// exactly 65518 at level 1. So the condition is the stream status, never the length.
    fn deflate_block(&mut self) -> io::Result<usize> {
        if self.buffer.is_empty() {
            return Ok(0);
        }

        // Capacity is the bound, matching Java's fixed-size output array.
        let mut compressed = Vec::with_capacity(COMPRESSED_BUFFER_SIZE);
        let mut c = Compress::new(Compression::new(self.level), false);
        let status = c
            .compress_vec(&self.buffer, &mut compressed, FlushCompress::Finish)
            .map_err(io::Error::other)?;

        if status != Status::StreamEnd {
            // Matches htsjdk's `noCompressionDeflater`, explicitly the plain JDK deflater at
            // NO_COMPRESSION, which predictably yields input + 10 bytes.
            compressed.clear();
            compressed.reserve(COMPRESSED_BUFFER_SIZE);
            let mut nc = Compress::new(Compression::none(), false);
            let nc_status = nc
                .compress_vec(&self.buffer, &mut compressed, FlushCompress::Finish)
                .map_err(io::Error::other)?;
            // htsjdk throws IllegalStateException("unpossible") here. NO_COMPRESSION yields
            // input + 10 bytes, and the uncompressed block size is chosen so that always fits.
            debug_assert_eq!(nc_status, Status::StreamEnd);
        }

        let mut crc = Crc::new();
        crc.update(&self.buffer);

        let total = self.write_gzip_block(&compressed, self.buffer.len(), crc.sum())?;
        self.buffer.clear();
        Ok(total)
    }

    /// Ports `writeGzipBlock`. All multi-byte fields are little-endian.
    fn write_gzip_block(
        &mut self,
        compressed: &[u8],
        uncompressed_size: usize,
        crc: u32,
    ) -> io::Result<usize> {
        let total = compressed.len() + BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH;

        let mut header = [0u8; BLOCK_HEADER_LENGTH];
        header[0] = 31; // GZIP_ID1
        header[1] = 139; // GZIP_ID2
        header[2] = 8; // CM = deflate
        header[3] = 4; // FLG = FEXTRA
                       // header[4..8] MTIME stays zero, which is what makes the output reproducible
        header[8] = 0; // XFL
        header[9] = 255; // OS = unknown
        header[10..12].copy_from_slice(&6u16.to_le_bytes()); // XLEN
        header[12] = 66; // 'B'
        header[13] = 67; // 'C'
        header[14..16].copy_from_slice(&2u16.to_le_bytes()); // BGZF_LEN
                                                             // "I don't know why we store block size - 1, but that is what the spec says" (htsjdk)
        header[16..18].copy_from_slice(&((total - 1) as u16).to_le_bytes());

        self.inner.write_all(&header)?;
        self.inner.write_all(compressed)?;
        self.inner.write_all(&crc.to_le_bytes())?;
        self.inner
            .write_all(&(uncompressed_size as u32).to_le_bytes())?;
        Ok(total)
    }

    /// Flushes any buffered data and appends the BGZF terminator block, as
    /// `close(writeTerminatorBlock = true)` does.
    pub fn finish(&mut self) -> io::Result<()> {
        if self.finished {
            return Ok(());
        }
        self.deflate_block()?;
        self.inner.write_all(&EMPTY_GZIP_BLOCK)?;
        self.inner.flush()?;
        self.finished = true;
        Ok(())
    }

    pub fn into_inner(mut self) -> io::Result<W> {
        self.finish()?;
        Ok(self.inner)
    }
}

impl<W: Write> Write for BgzfWriter<W> {
    fn write(&mut self, mut data: &[u8]) -> io::Result<usize> {
        let n = data.len();
        while !data.is_empty() {
            let free = DEFAULT_UNCOMPRESSED_BLOCK_SIZE - self.buffer.len();
            let take = free.min(data.len());
            self.buffer.extend_from_slice(&data[..take]);
            data = &data[take..];
            if self.buffer.len() == DEFAULT_UNCOMPRESSED_BLOCK_SIZE {
                self.deflate_block()?;
            }
        }
        Ok(n)
    }

    /// Ports htsjdk's `flush()`, which is **not** a no-op on the byte stream: it emits the
    /// buffered bytes as a block, creating a block boundary at the flush point. A caller that
    /// flushes mid-stream therefore gets a different (still valid) block layout, and the port
    /// must reproduce that rather than silently deferring.
    fn flush(&mut self) -> io::Result<()> {
        while !self.buffer.is_empty() {
            self.deflate_block()?;
        }
        self.inner.flush()
    }
}
