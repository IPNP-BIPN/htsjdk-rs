//! BGZF write path, byte-identical to htsjdk's `BlockCompressedOutputStream`.
//!
//! Ported from htsjdk 4.2.0
//! `src/main/java/htsjdk/samtools/util/BlockCompressedOutputStream.java`
//! (`deflateBlock`, `writeGzipBlock`, `flush`, `close`).

use std::io::{self, Write};

use flate2::{Compress, Compression, Crc, FlushCompress, Status};

use crate::{
    vfp, BGZF_ID1, BGZF_ID2, BGZF_LEN, BLOCK_FOOTER_LENGTH, BLOCK_HEADER_LENGTH,
    COMPRESSED_BUFFER_SIZE, DEFAULT_COMPRESSION_LEVEL, DEFAULT_UNCOMPRESSED_BLOCK_SIZE,
    EMPTY_GZIP_BLOCK, GZIP_CM_DEFLATE, GZIP_FLG, GZIP_ID1, GZIP_ID2, GZIP_OS_UNKNOWN, GZIP_XFL,
    GZIP_XLEN,
};

/// Writes a BGZF stream whose bytes match `BlockCompressedOutputStream`.
pub struct BgzfWriter<W: Write> {
    inner: W,
    buffer: Vec<u8>,
    level: u32,
    finished: bool,
    /// Byte offset of the block currently being filled, `mBlockAddress` in htsjdk.
    block_address: u64,
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
            block_address: 0,
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
        self.block_address += total as u64;
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
        header[0] = GZIP_ID1;
        header[1] = GZIP_ID2;
        header[2] = GZIP_CM_DEFLATE;
        header[3] = GZIP_FLG;
        // header[4..8] MTIME stays zero, which is what makes the output reproducible
        header[8] = GZIP_XFL;
        header[9] = GZIP_OS_UNKNOWN;
        header[10..12].copy_from_slice(&GZIP_XLEN.to_le_bytes());
        header[12] = BGZF_ID1;
        header[13] = BGZF_ID2;
        header[14..16].copy_from_slice(&BGZF_LEN.to_le_bytes());
        // "I don't know why we store block size - 1, but that is what the spec says" (htsjdk)
        header[16..18].copy_from_slice(&((total - 1) as u16).to_le_bytes());

        self.inner.write_all(&header)?;
        self.inner.write_all(compressed)?;
        self.inner.write_all(&crc.to_le_bytes())?;
        self.inner
            .write_all(&(uncompressed_size as u32).to_le_bytes())?;
        Ok(total)
    }

    /// `BlockCompressedOutputStream.getFilePointer()`: the virtual file pointer at the current
    /// write position.
    ///
    /// The upper 48 bits are the compressed byte offset of the block being filled; the lower 16
    /// are the offset into that block's *uncompressed* payload. Taken before and after writing a
    /// record, the pair is exactly the chunk the BAM index stores.
    pub fn file_pointer(&self) -> u64 {
        // The buffer never exceeds DEFAULT_UNCOMPRESSED_BLOCK_SIZE (65498), which fits the
        // 16-bit offset field, so this cannot fail in practice.
        vfp::make_file_pointer(self.block_address, self.buffer.len() as u32)
            .expect("block offset is bounded by the uncompressed block size")
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

    /// Flushes any buffered bytes as a block, then returns the inner writer **without** appending
    /// the EOF terminator block. This is what `BAMFileWriter.writeHeader(OutputStream)` needs: it
    /// wraps the stream in a `BlockCompressedOutputStream`, writes the header, and calls `flush()`
    /// only, never `close()`, so the header is a complete block boundary with no terminator (the
    /// terminator comes later, from whatever is appended after it).
    pub fn into_inner_without_terminator(mut self) -> io::Result<W> {
        self.flush()?;
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{decompress_all, EMPTY_GZIP_BLOCK};
    use std::io::Write;

    #[test]
    fn into_inner_without_terminator_omits_the_eof_block_but_keeps_the_data() {
        let payload = b"BAM\x01some header bytes";

        let mut w = BgzfWriter::new(Vec::new());
        w.write_all(payload).unwrap();
        let without = w.into_inner_without_terminator().unwrap();

        let mut w2 = BgzfWriter::new(Vec::new());
        w2.write_all(payload).unwrap();
        let with = w2.into_inner().unwrap();

        // The terminated stream is exactly the un-terminated one followed by the EOF block.
        assert_eq!(with.len(), without.len() + EMPTY_GZIP_BLOCK.len());
        assert_eq!(&with[..without.len()], &without[..]);
        assert!(!without.ends_with(&EMPTY_GZIP_BLOCK));
        // Both decompress to the same payload.
        assert_eq!(decompress_all(&without).unwrap(), payload);
    }
}
