//! BGZF write path, byte-identical to htsjdk's `BlockCompressedOutputStream`.
//!
//! # Which deflater, which is not a detail of the framing
//!
//! `BlockCompressedOutputStream` deflates through a **static** `defaultDeflaterFactory`, and GATK
//! replaces it: everything it writes without `--use-jdk-deflater` is Intel's GKL and not the JDK's
//! zlib. The framing is identical either way and the bytes are not, so a writer that only knew one
//! of them could reproduce htsjdk's own output and none of GATK's. [`Deflater`] is that choice,
//! and [`Deflater::Jdk`] is the default because it is htsjdk's.
//!
//! Ported from htsjdk 4.2.0
//! `src/main/java/htsjdk/samtools/util/BlockCompressedOutputStream.java`
//! (`deflateBlock`, `writeGzipBlock`, `flush`, `close`), with [`Deflater::Libdeflate`] from htsjdk
//! 5.0.0's `LibdeflateDeflater` and jlibdeflate 0.1.0's native `deflateCompress`.

use std::io::{self, Write};

use flate2::{Compress, Compression, Crc, FlushCompress, Status};

use crate::{
    vfp, BGZF_ID1, BGZF_ID2, BGZF_LEN, BLOCK_FOOTER_LENGTH, BLOCK_HEADER_LENGTH,
    COMPRESSED_BUFFER_SIZE, DEFAULT_COMPRESSION_LEVEL, DEFAULT_UNCOMPRESSED_BLOCK_SIZE,
    EMPTY_GZIP_BLOCK, GZIP_CM_DEFLATE, GZIP_FLG, GZIP_ID1, GZIP_ID2, GZIP_OS_UNKNOWN, GZIP_XFL,
    GZIP_XLEN,
};

/// Which implementation deflates a block, which decides its bytes and nothing else.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Deflater {
    /// `java.util.zip.Deflater`, which is htsjdk's own default and what `--use-jdk-deflater` asks
    /// for.
    #[default]
    Jdk,
    /// Intel's GKL, which GATK installs as the default factory: zlib 1.2.13 carrying Intel's
    /// `deflate_medium` patch above level two, and ISA-L igzip at one and two.
    ///
    /// Levels 1 and 2 need the `gkl-igzip` feature, which links ISA-L and therefore wants an
    /// assembler at build time. Without it they refuse rather than answer with zlib's bytes, which
    /// is `gkl-deflate`'s own behaviour and the reason the feature is off by default: BGZF's level
    /// is 5, so nothing here pays a C toolchain for a pair of levels it never asks for.
    Gkl,
    /// libdeflate v1.25, which htsjdk 5.0.0's own `DeflaterFactory` returns whenever its native
    /// library loads (`samjdk.use_libdeflate`, true by default). So it is what Picard 3.5.0 writes
    /// with `USE_JDK_DEFLATER=true`, despite the name: that argument only stops Picard installing
    /// GKL, and the factory left in place is htsjdk's.
    Libdeflate,
}

/// A libdeflate compressor at one level: `LibdeflateCompressor` in jlibdeflate, which
/// `LibdeflateDeflater` keeps for the life of the stream.
struct LibdeflateState(std::ptr::NonNull<libdeflate_sys::libdeflate_compressor>);

// The compressor is plain heap state that libdeflate never shares; it may move between threads,
// it is just never used from two at once, which `&mut self` already rules out.
unsafe impl Send for LibdeflateState {}

impl LibdeflateState {
    fn new(level: u32) -> io::Result<Self> {
        // SAFETY: any level 0..=12 is valid, and `with_deflater` has already bounded it to 0..=9.
        let raw = unsafe { libdeflate_sys::libdeflate_alloc_compressor(level as std::ffi::c_int) };
        std::ptr::NonNull::new(raw)
            .map(LibdeflateState)
            .ok_or_else(|| io::Error::other("libdeflate could not allocate a compressor"))
    }

    /// `libdeflate_deflate_compress` into exactly `avail` bytes: the length written, or `None` when
    /// it does not fit, which is jlibdeflate's -1 and leaves `LibdeflateDeflater.finished()` false.
    fn compress(&mut self, input: &[u8], out: &mut Vec<u8>, avail: usize) -> Option<usize> {
        out.clear();
        out.resize(avail, 0);
        // SAFETY: both pointers are valid for the lengths passed, and the compressor is ours.
        let written = unsafe {
            libdeflate_sys::libdeflate_deflate_compress(
                self.0.as_ptr(),
                input.as_ptr().cast(),
                input.len(),
                out.as_mut_ptr().cast(),
                avail,
            )
        };
        out.truncate(written);
        (written != 0).then_some(written)
    }
}

impl Drop for LibdeflateState {
    fn drop(&mut self) {
        // SAFETY: allocated by `libdeflate_alloc_compressor` and freed exactly once, here.
        unsafe { libdeflate_sys::libdeflate_free_compressor(self.0.as_ptr()) }
    }
}

/// Writes a BGZF stream whose bytes match `BlockCompressedOutputStream`.
pub struct BgzfWriter<W: Write> {
    inner: W,
    buffer: Vec<u8>,
    level: u32,
    deflater: Deflater,
    finished: bool,
    /// Byte offset of the block currently being filled, `mBlockAddress` in htsjdk.
    block_address: u64,
    /// The deflated block, reused between blocks: `compressedBuffer` in htsjdk, allocated once in
    /// the constructor rather than once per block.
    compressed: Vec<u8>,
    /// The zlib deflater, reused and reset between blocks, which is htsjdk's `deflater` field.
    /// Built on first use because the GKL path never touches it.
    deflater_state: Option<Compress>,
    /// The libdeflate compressor, likewise built on first use.
    libdeflate_state: Option<LibdeflateState>,
}

impl<W: Write> BgzfWriter<W> {
    pub fn new(inner: W) -> Self {
        Self::with_level(inner, DEFAULT_COMPRESSION_LEVEL)
    }

    pub fn with_level(inner: W, level: u32) -> Self {
        Self::with_deflater(inner, level, Deflater::Jdk)
    }

    /// The same writer with the deflater named, which is the only thing that differs between
    /// htsjdk's output and GATK's.
    pub fn with_deflater(inner: W, level: u32, deflater: Deflater) -> Self {
        assert!(level <= 9, "compression level must be 0..=9, got {level}");
        Self {
            inner,
            buffer: Vec::with_capacity(DEFAULT_UNCOMPRESSED_BLOCK_SIZE),
            compressed: Vec::with_capacity(COMPRESSED_BUFFER_SIZE),
            deflater_state: None,
            libdeflate_state: None,
            level,
            deflater,
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

        // Capacity is the bound, matching Java's fixed-size output array. The vector itself is
        // reused between blocks: `BlockCompressedOutputStream` allocates `compressedBuffer` once
        // in its constructor, and a stream is a long run of blocks.
        let mut compressed = std::mem::take(&mut self.compressed);
        compressed.clear();
        compressed.reserve(COMPRESSED_BUFFER_SIZE);
        let fits = match self.deflater {
            Deflater::Jdk => {
                // One deflater for the life of the stream, reset between blocks, which is what
                // htsjdk's `deflater` field is. Constructing one per block allocates zlib's window
                // and its hash tables every time -- a quarter of a megabyte at level 5.
                let compressor = self
                    .deflater_state
                    .get_or_insert_with(|| Compress::new(Compression::new(self.level), false));
                compressor.reset();
                let status = compressor
                    .compress_vec(&self.buffer, &mut compressed, FlushCompress::Finish)
                    .map_err(io::Error::other)?;
                status == Status::StreamEnd
            }
            Deflater::Gkl => {
                // Levels 1 and 2 are GKL's igzip pair, and `deflate_gkl` PANICS on them when ISA-L
                // is absent or unusable rather than answering with zlib's bytes. That refusal is
                // the right one and a panic is the wrong way to deliver it: a writer is handed a
                // level by a caller, and a caller can be told.
                if self.level <= 2 && !gkl_deflate::igzip_available() {
                    return Err(io::Error::other(format!(
                        "GKL routes level {} through ISA-L igzip, which this build cannot reach: \
                         rebuild with the `gkl-igzip` feature, on a host where it is usable",
                        self.level
                    )));
                }
                compressed = gkl_deflate::deflate_gkl(&self.buffer, self.level as usize);
                // htsjdk asks the deflater whether it finished, having given it an output array of
                // exactly `COMPRESSED_BUFFER_SIZE` bytes. A deflater handed a vector never runs
                // out, so the same question is asked of the length: output that would have filled
                // that array is output the JDK path could not have declared finished, EQUAL length
                // included, which is the boundary the comment below is about.
                compressed.len() < COMPRESSED_BUFFER_SIZE
            }
            Deflater::Libdeflate => {
                // One call deflates the whole block, into an output array of exactly
                // `COMPRESSED_BUFFER_SIZE` bytes. Unlike zlib's, libdeflate's answer at that
                // boundary is not ambiguous: output that fills the array exactly was written whole.
                let state = match self.libdeflate_state.as_mut() {
                    Some(state) => state,
                    None => self
                        .libdeflate_state
                        .insert(LibdeflateState::new(self.level)?),
                };
                state
                    .compress(&self.buffer, &mut compressed, COMPRESSED_BUFFER_SIZE)
                    .is_some()
            }
        };

        if !fits {
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
        self.compressed = compressed;
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

    /// GKL's igzip levels refuse rather than panic when ISA-L is out of reach.
    ///
    /// The refusal itself is `gkl-deflate`'s, and it is the right answer: levels 1 and 2 are the
    /// only pair GKL does not route through zlib, so answering them with zlib's bytes would be a
    /// wrong answer that looks like a right one. What is asserted here is that a caller is TOLD,
    /// rather than the process ending.
    #[test]
    fn a_gkl_level_that_needs_igzip_refuses_when_it_is_out_of_reach() {
        if gkl_deflate::igzip_available() {
            // A build with ISA-L usable answers those levels, which is the other half of the
            // claim: the refusal is about reach and not about the level.
            for level in [1, 2] {
                let mut w = BgzfWriter::with_deflater(Vec::new(), level, Deflater::Gkl);
                w.write_all(b"chr1\t100\t.\tA\tC\n").unwrap();
                assert_eq!(
                    decompress_all(&w.into_inner().unwrap()).unwrap(),
                    b"chr1\t100\t.\tA\tC\n"
                );
            }
            return;
        }
        for level in [1, 2] {
            let mut w = BgzfWriter::with_deflater(Vec::new(), level, Deflater::Gkl);
            w.write_all(b"chr1\t100\t.\tA\tC\n").unwrap();
            let error = w
                .into_inner()
                .expect_err("igzip is out of reach on this build");
            assert!(error.to_string().contains("igzip"), "{error}");
        }
        // Every other level is unaffected, since GKL routes them through zlib.
        let mut w = BgzfWriter::with_deflater(Vec::new(), 5, Deflater::Gkl);
        w.write_all(b"chr1\t100\t.\tA\tC\n").unwrap();
        assert!(w.into_inner().is_ok());
    }

    /// The deflater changes the block's bytes and nothing else about the file.
    ///
    /// The GKL path's own agreement with the reference is measured where a GKL oracle exists,
    /// which is `gkl-deflate`'s suites and `gatk-rs`'s goldens; what is asserted here is that this
    /// writer frames what that crate produced, and that the framing is otherwise untouched.
    #[test]
    fn the_gkl_deflater_frames_gkl_bytes_and_changes_nothing_else() {
        // Text rather than random bytes: the two implementations agree on incompressible input.
        let payload: Vec<u8> = std::iter::repeat_n(&b"chr1\t100\t.\tA\tC\t.\t.\t.\n"[..], 40)
            .flatten()
            .copied()
            .collect();

        let mut jdk = BgzfWriter::new(Vec::new());
        jdk.write_all(&payload).unwrap();
        let jdk = jdk.into_inner().unwrap();

        let mut gkl =
            BgzfWriter::with_deflater(Vec::new(), DEFAULT_COMPRESSION_LEVEL, Deflater::Gkl);
        gkl.write_all(&payload).unwrap();
        let gkl = gkl.into_inner().unwrap();

        assert_ne!(jdk, gkl, "the two deflaters do not agree on this payload");
        assert_eq!(decompress_all(&jdk).unwrap(), payload);
        assert_eq!(decompress_all(&gkl).unwrap(), payload);
        // The block's payload is what `gkl-deflate` produced, framed: the header, the CRC and the
        // uncompressed length are the format's and are not the deflater's to change.
        let deflated = gkl_deflate::deflate_gkl(&payload, DEFAULT_COMPRESSION_LEVEL as usize);
        let framed = BLOCK_HEADER_LENGTH + deflated.len() + BLOCK_FOOTER_LENGTH;
        assert_eq!(
            &gkl[BLOCK_HEADER_LENGTH..BLOCK_HEADER_LENGTH + deflated.len()],
            &deflated[..]
        );
        assert_eq!(&gkl[framed..], &EMPTY_GZIP_BLOCK[..]);
        // Both files end with the terminator, and only the middle differs.
        assert!(jdk.ends_with(&EMPTY_GZIP_BLOCK));
    }

    /// libdeflate writes its own bytes, which decompress to the input at every level BGZF accepts,
    /// for text and for incompressible data alike. Whether they are htsjdk 5.0.0's is measured
    /// against Picard 3.5.0 by picard-rs, whose `USE_JDK_DEFLATER=true` reaches this deflater.
    #[test]
    fn the_libdeflate_deflater_round_trips_at_every_level() {
        let text: Vec<u8> = std::iter::repeat_n(&b"chr1\t100\t.\tA\tC\t.\t.\t.\n"[..], 4000)
            .flatten()
            .copied()
            .collect();
        // A linear congruential generator: incompressible enough to reach libdeflate's stored
        // blocks, and the same bytes on every run.
        let mut state = 0x2545_f491_u32;
        let noise: Vec<u8> = (0..200_000)
            .map(|_| {
                state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                (state >> 24) as u8
            })
            .collect();
        for payload in [&text, &noise] {
            for level in 0..=9 {
                let mut writer = BgzfWriter::with_deflater(Vec::new(), level, Deflater::Libdeflate);
                writer.write_all(payload).unwrap();
                let written = writer.into_inner().unwrap();
                assert_eq!(&decompress_all(&written).unwrap(), payload, "level {level}");
                assert!(written.ends_with(&EMPTY_GZIP_BLOCK));
            }
        }

        let mut jdk = BgzfWriter::new(Vec::new());
        jdk.write_all(&text).unwrap();
        let mut libdeflate =
            BgzfWriter::with_deflater(Vec::new(), DEFAULT_COMPRESSION_LEVEL, Deflater::Libdeflate);
        libdeflate.write_all(&text).unwrap();
        assert_ne!(
            jdk.into_inner().unwrap(),
            libdeflate.into_inner().unwrap(),
            "libdeflate and zlib do not agree on this payload"
        );
    }
}
