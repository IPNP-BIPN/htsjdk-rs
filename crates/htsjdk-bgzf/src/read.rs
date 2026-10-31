//! BGZF read path, matching htsjdk's `BlockCompressedInputStream`.
//!
//! Ported from htsjdk 4.2.0:
//! `src/main/java/htsjdk/samtools/util/BlockCompressedInputStream.java`
//! (`processNextBlock`, `inflateBlock`, `isValidBlockHeader`) and
//! `src/main/java/htsjdk/samtools/util/BlockGunzipper.java` (`unzipBlock`).
//!
//! Fidelity here is about *acceptance*, not bytes: the port must accept exactly the files
//! htsjdk accepts and reject exactly the ones it rejects. In particular CRC checking is
//! **off by default**, matching `BlockGunzipper.checkCrcs = false`. A reader that always
//! verified CRCs would reject files htsjdk reads happily.

use std::io::{self, Read};

use flate2::{Crc, Decompress, FlushDecompress};

use crate::inflate::Inflater;
use crate::{
    vfp, BGZF_ID1, BGZF_ID2, BLOCK_FOOTER_LENGTH, BLOCK_HEADER_LENGTH, GZIP_CM_DEFLATE, GZIP_FLG,
    GZIP_ID1, GZIP_ID2, GZIP_XLEN, MAX_COMPRESSED_BLOCK_SIZE,
};

/// Why a BGZF stream was rejected. Messages mirror htsjdk's wording.
#[derive(Debug)]
pub enum BgzfError {
    InvalidGzipHeader,
    BlockSizeDisagreement { declared: usize, actual: usize },
    UnexpectedBlockLength(usize),
    IncorrectHeaderSize,
    PrematureEnd,
    CrcMismatch { expected: u32, actual: u32 },
    Inflate(String),
}

impl std::fmt::Display for BgzfError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidGzipHeader => write!(f, "Invalid GZIP header"),
            Self::BlockSizeDisagreement { declared, actual } => {
                write!(
                    f,
                    "GZIP blocksize disagreement: declared {declared}, actual {actual}"
                )
            }
            Self::UnexpectedBlockLength(n) => write!(f, "Unexpected compressed block length: {n}"),
            Self::IncorrectHeaderSize => write!(f, "Incorrect header size for file"),
            Self::PrematureEnd => write!(f, "Premature end of file"),
            Self::CrcMismatch { expected, actual } => {
                write!(f, "CRC mismatch: expected {expected:#x}, got {actual:#x}")
            }
            Self::Inflate(m) => write!(f, "inflate failed: {m}"),
        }
    }
}

impl std::error::Error for BgzfError {}

impl From<BgzfError> for io::Error {
    fn from(e: BgzfError) -> Self {
        io::Error::new(io::ErrorKind::InvalidData, e)
    }
}

/// Ports `isValidBlockHeader`. Note it checks that the FEXTRA bit is *set*, not that FLG
/// equals 4, and that XLEN is exactly 6.
pub fn is_valid_block_header(h: &[u8]) -> bool {
    h.len() >= BLOCK_HEADER_LENGTH
        && h[0] == GZIP_ID1
        && h[1] == GZIP_ID2
        && (h[3] & GZIP_FLG) != 0
        && h[10] == GZIP_XLEN as u8
        && h[12] == BGZF_ID1
        && h[13] == BGZF_ID2
}

/// One decompressed block plus where it came from.
#[derive(Debug, Clone)]
pub struct DecompressedBlock {
    /// Byte offset of this block's first byte in the compressed stream.
    pub block_address: u64,
    /// Framed size of the compressed block.
    pub block_compressed_size: usize,
    pub data: Vec<u8>,
}

/// Reads a BGZF stream the way htsjdk does.
pub struct BgzfReader<R: Read> {
    inner: R,
    current: Vec<u8>,
    offset: usize,
    stream_offset: u64,
    block_address: u64,
    check_crcs: bool,
    eof: bool,
    /// The compressed block, reused between blocks.
    ///
    /// A BGZF stream is a long run of blocks of at most 64 KiB, and allocating -- and ZEROING --
    /// one buffer per block is a megabyte of `memset` for every sixteen blocks. htsjdk allocates
    /// `compressedBuffer` once in `BlockCompressedInputStream`'s constructor for the same reason.
    scratch: Vec<u8>,
    /// The inflater, reused between blocks. `BlockGunzipper` holds one `Inflater` for the life of
    /// the stream and calls `reset()`; constructing one per block allocates zlib's window each
    /// time.
    inflater: Option<Inflater>,
    /// zlib's, built when a block first needs an answer [`Inflater`] does not give, which is any
    /// damaged block and anything unusual.
    zlib: Option<Decompress>,
}

impl<R: Read> BgzfReader<R> {
    pub fn new(inner: R) -> Self {
        Self {
            inner,
            current: Vec::new(),
            offset: 0,
            stream_offset: 0,
            block_address: 0,
            check_crcs: false, // matches BlockGunzipper's default
            eof: false,
            scratch: Vec::new(),
            inflater: None,
            zlib: None,
        }
    }

    /// Ports `setCheckCrcs`. Off by default, as in htsjdk.
    pub fn set_check_crcs(&mut self, check: bool) {
        self.check_crcs = check;
    }

    /// Ports `getFilePointer`: the virtual pointer of the next byte to be returned.
    ///
    /// htsjdk advances to the next block address when the current block is exhausted, so that
    /// a pointer never dangles past the end of a block.
    pub fn virtual_pos(&self) -> u64 {
        if self.offset >= self.current.len() && !self.current.is_empty() {
            vfp::make_file_pointer(self.stream_offset, 0).unwrap_or(0)
        } else {
            vfp::make_file_pointer(self.block_address, self.offset as u32).unwrap_or(0)
        }
    }

    /// Reads and decompresses the next block, or returns `None` at end of stream.
    ///
    /// Ports `processNextBlock` plus `BlockGunzipper.unzipBlock`.
    pub fn next_block(&mut self) -> Result<Option<DecompressedBlock>, BgzfError> {
        let block_address = self.stream_offset;
        let mut header = [0u8; BLOCK_HEADER_LENGTH];

        let got = read_up_to(&mut self.inner, &mut header).map_err(|_| BgzfError::PrematureEnd)?;
        if got == 0 {
            // htsjdk tolerates a stream that simply ends, with no terminator block.
            return Ok(None);
        }
        self.stream_offset += got as u64;
        if got != BLOCK_HEADER_LENGTH {
            return Err(BgzfError::IncorrectHeaderSize);
        }
        if header[0] != GZIP_ID1
            || header[1] != GZIP_ID2
            || header[2] != GZIP_CM_DEFLATE
            || header[3] != GZIP_FLG
        {
            return Err(BgzfError::InvalidGzipHeader);
        }
        if u16::from_le_bytes([header[10], header[11]]) != GZIP_XLEN {
            return Err(BgzfError::InvalidGzipHeader);
        }

        let block_length = u16::from_le_bytes([header[16], header[17]]) as usize + 1;
        if !(BLOCK_HEADER_LENGTH..=MAX_COMPRESSED_BLOCK_SIZE).contains(&block_length) {
            return Err(BgzfError::UnexpectedBlockLength(block_length));
        }

        let remaining = block_length - BLOCK_HEADER_LENGTH;
        let mut rest = std::mem::take(&mut self.scratch);
        // Kept at its largest, so it is zeroed once and not once per block.
        if rest.len() < MAX_COMPRESSED_BLOCK_SIZE {
            rest.resize(MAX_COMPRESSED_BLOCK_SIZE, 0);
        }
        let got = read_up_to(&mut self.inner, &mut rest[..remaining]).map_err(|error| {
            self.scratch = std::mem::take(&mut rest);
            let _ = error;
            BgzfError::PrematureEnd
        })?;
        self.stream_offset += got as u64;
        if got != remaining {
            self.scratch = rest;
            return Err(BgzfError::PrematureEnd);
        }

        // Footer: CRC32 then ISIZE, both little-endian.
        let footer = &rest[remaining - BLOCK_FOOTER_LENGTH..remaining];
        let expected_crc = u32::from_le_bytes([footer[0], footer[1], footer[2], footer[3]]);
        let uncompressed_len =
            u32::from_le_bytes([footer[4], footer[5], footer[6], footer[7]]) as usize;

        let deflated = &rest[..remaining - BLOCK_FOOTER_LENGTH];
        let mut data = Vec::with_capacity(uncompressed_len);
        let inflated = self
            .inflater
            .get_or_insert_with(Inflater::new)
            .inflate_vec(deflated, &mut data);
        let outcome = if inflated {
            Ok(())
        } else {
            let zlib = self.zlib.get_or_insert_with(|| Decompress::new(false));
            zlib.reset(false);
            zlib.decompress_vec(deflated, &mut data, FlushDecompress::Finish)
                .map(drop)
                .map_err(|e| BgzfError::Inflate(e.to_string()))
        };
        // The buffer goes back before any refusal is returned, so a stream that fails mid-way does
        // not leave the reader without one.
        self.scratch = rest;
        outcome?;

        if data.len() != uncompressed_len {
            return Err(BgzfError::BlockSizeDisagreement {
                declared: uncompressed_len,
                actual: data.len(),
            });
        }

        if self.check_crcs {
            let mut crc = Crc::new();
            crc.update(&data);
            if crc.sum() != expected_crc {
                return Err(BgzfError::CrcMismatch {
                    expected: expected_crc,
                    actual: crc.sum(),
                });
            }
        }

        Ok(Some(DecompressedBlock {
            block_address,
            block_compressed_size: block_length,
            data,
        }))
    }

    fn fill(&mut self) -> io::Result<()> {
        while self.offset >= self.current.len() {
            match self.next_block()? {
                None => {
                    self.eof = true;
                    return Ok(());
                }
                Some(b) => {
                    // A zero-length block is the terminator; keep going in case more follows,
                    // which is what htsjdk does for concatenated BGZF files.
                    self.block_address = b.block_address;
                    self.current = b.data;
                    self.offset = 0;
                    if self.current.is_empty() && self.stream_offset == 0 {
                        self.eof = true;
                        return Ok(());
                    }
                    if self.current.is_empty() {
                        continue;
                    }
                }
            }
        }
        Ok(())
    }
}

impl<R: Read> Read for BgzfReader<R> {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        if out.is_empty() {
            return Ok(0);
        }
        self.fill()?;
        if self.eof && self.offset >= self.current.len() {
            return Ok(0);
        }
        let n = (self.current.len() - self.offset).min(out.len());
        out[..n].copy_from_slice(&self.current[self.offset..self.offset + n]);
        self.offset += n;
        Ok(n)
    }
}

/// Reads until the buffer is full or the stream ends, returning how many bytes landed.
fn read_up_to<R: Read>(r: &mut R, buf: &mut [u8]) -> io::Result<usize> {
    let mut total = 0;
    while total < buf.len() {
        match r.read(&mut buf[total..]) {
            Ok(0) => break,
            Ok(n) => total += n,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
            Err(e) => return Err(e),
        }
    }
    Ok(total)
}

/// Decompresses an entire BGZF stream into one buffer.
///
/// A convenience over [`BgzfReader::next_block`] for callers that want the payload rather than
/// the block structure. It deliberately does not hide the block structure from anything that
/// needs it: the differential harness reports positions in terms of blocks, and losing that
/// would make a divergence in a large BAM unlocatable.
///
/// The blocks are inflated on rayon's workers, each into its own slice of the output, after a
/// serial walk over the headers that applies [`BgzfReader::next_block`]'s checks in its order. The
/// error returned is the one the serial loop would have stopped at: the first failing block's, or
/// the header walk's if every block before it inflated.
pub fn decompress_all(bytes: &[u8]) -> Result<Vec<u8>, BgzfError> {
    use rayon::prelude::*;

    let mut blocks = Vec::new();
    let mut at = 0usize;
    let walk_error = loop {
        match block_frame(&bytes[at..]) {
            Ok(None) => break None,
            Ok(Some(frame)) => {
                blocks.push((at, frame));
                at += frame.0;
            }
            Err(error) => break Some(error),
        }
    };
    // A declared size past BGZF's 64 KiB is a damaged block; the serial reader sizes its buffer
    // from it one block at a time, so it is left to refuse that block in its own way.
    if blocks.len() < 2 || blocks.iter().any(|(_, (_, size))| *size > 1 << 16) {
        return decompress_serial(bytes);
    }

    let total: usize = blocks.iter().map(|(_, (_, size))| size).sum();
    let mut out = vec![0u8; total];
    let mut rest = out.as_mut_slice();
    let mut jobs = Vec::with_capacity(blocks.len());
    for &(start, (length, size)) in &blocks {
        let (slot, tail) = rest.split_at_mut(size);
        rest = tail;
        let deflated = &bytes[start + BLOCK_HEADER_LENGTH..start + length - BLOCK_FOOTER_LENGTH];
        jobs.push((deflated, slot));
    }
    let first_failure = jobs
        .into_par_iter()
        .with_min_len(4)
        .map(|(deflated, slot)| inflate_on_worker(deflated, slot))
        .collect::<Vec<_>>()
        .into_iter()
        .find_map(Result::err);
    match first_failure.or(walk_error) {
        Some(error) => Err(error),
        None => Ok(out),
    }
}

/// One block on a rayon worker, each worker with its own tables. An error is the one the serial
/// loop would have stopped at, because only zlib ever gives one.
fn inflate_on_worker(deflated: &[u8], slot: &mut [u8]) -> Result<(), BgzfError> {
    thread_local! {
        static INFLATER: std::cell::RefCell<Option<Inflater>> = const { std::cell::RefCell::new(None) };
    }
    let inflated = INFLATER.with_borrow_mut(|state| {
        state
            .get_or_insert_with(Inflater::new)
            .inflate(deflated, slot)
    });
    if inflated {
        return Ok(());
    }
    let mut zlib = Decompress::new(false);
    zlib.decompress(deflated, slot, FlushDecompress::Finish)
        .map_err(|e| BgzfError::Inflate(e.to_string()))?;
    let actual = zlib.total_out() as usize;
    if actual != slot.len() {
        return Err(BgzfError::BlockSizeDisagreement {
            declared: slot.len(),
            actual,
        });
    }
    Ok(())
}

/// One block's framing at the start of `bytes`: its compressed length and declared inflated size,
/// with the refusals [`BgzfReader::next_block`] makes before it inflates anything. `None` is the
/// end of the stream.
fn block_frame(bytes: &[u8]) -> Result<Option<(usize, usize)>, BgzfError> {
    if bytes.is_empty() {
        return Ok(None);
    }
    if bytes.len() < BLOCK_HEADER_LENGTH {
        return Err(BgzfError::IncorrectHeaderSize);
    }
    let header = &bytes[..BLOCK_HEADER_LENGTH];
    if header[0] != GZIP_ID1
        || header[1] != GZIP_ID2
        || header[2] != GZIP_CM_DEFLATE
        || header[3] != GZIP_FLG
        || u16::from_le_bytes([header[10], header[11]]) != GZIP_XLEN
    {
        return Err(BgzfError::InvalidGzipHeader);
    }
    let length = u16::from_le_bytes([header[16], header[17]]) as usize + 1;
    if !(BLOCK_HEADER_LENGTH..=MAX_COMPRESSED_BLOCK_SIZE).contains(&length) {
        return Err(BgzfError::UnexpectedBlockLength(length));
    }
    if bytes.len() < length {
        return Err(BgzfError::PrematureEnd);
    }
    let footer = &bytes[length - BLOCK_FOOTER_LENGTH..length];
    let size = u32::from_le_bytes([footer[4], footer[5], footer[6], footer[7]]) as usize;
    Ok(Some((length, size)))
}

fn decompress_serial(bytes: &[u8]) -> Result<Vec<u8>, BgzfError> {
    let mut reader = BgzfReader::new(bytes);
    let mut out = Vec::new();
    while let Some(block) = reader.next_block()? {
        out.extend_from_slice(&block.data);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn stream(len: usize) -> Vec<u8> {
        let payload: Vec<u8> = (0..len).map(|i| (i * 7 % 251) as u8).collect();
        let mut w = crate::BgzfWriter::new(Vec::new());
        w.write_all(&payload).unwrap();
        w.into_inner().unwrap()
    }

    /// What the reader answered before libdeflate was in front of zlib: one block through zlib alone.
    fn zlib_answer(block: &[u8]) -> Result<Vec<u8>, String> {
        let end = block.len() - BLOCK_FOOTER_LENGTH;
        let size = u32::from_le_bytes(block[end + 4..end + 8].try_into().unwrap()) as usize;
        let mut data = Vec::with_capacity(size);
        Decompress::new(false)
            .decompress_vec(
                &block[BLOCK_HEADER_LENGTH..end],
                &mut data,
                FlushDecompress::Finish,
            )
            .map_err(|e| BgzfError::Inflate(e.to_string()).to_string())?;
        if data.len() != size {
            return Err(BgzfError::BlockSizeDisagreement {
                declared: size,
                actual: data.len(),
            }
            .to_string());
        }
        Ok(data)
    }

    /// Every way of damaging a block answers as zlib alone did: bytes for bytes, message for
    /// message. A flipped CRC is the case that matters most, since htsjdk does not check one: the
    /// reader must still hand over the block.
    #[test]
    fn a_damaged_block_answers_as_zlib_alone_did() {
        let mut seed = 3u64;
        let mut next = move || {
            seed = seed
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (seed >> 33) as usize
        };
        for payload in [
            b"ACGTTGCAACGT".repeat(40),
            (0..600).map(|_| next() as u8).collect::<Vec<u8>>(),
            (0..600).map(|_| b"ACGTN"[next() % 5]).collect::<Vec<u8>>(),
        ] {
            let mut writer = crate::BgzfWriter::new(Vec::new());
            writer.write_all(&payload).unwrap();
            let stream = writer.into_inner().unwrap();
            let (length, _) = block_frame(&stream).unwrap().unwrap();
            for byte in BLOCK_HEADER_LENGTH..length {
                for bit in 0..8 {
                    let mut block = stream[..length].to_vec();
                    block[byte] ^= 1 << bit;
                    // A declared size this damage made absurd would size a buffer in either path.
                    let footer = length - BLOCK_FOOTER_LENGTH;
                    let size =
                        u32::from_le_bytes(block[footer + 4..footer + 8].try_into().unwrap());
                    if size > 1 << 20 {
                        continue;
                    }
                    let got = BgzfReader::new(&block[..])
                        .next_block()
                        .map(|answer| answer.unwrap().data)
                        .map_err(|e| e.to_string());
                    assert_eq!(got, zlib_answer(&block), "byte {byte} bit {bit}");
                }
            }
        }
    }

    /// The parallel path against the serial loop it replaces, on whole and damaged streams.
    #[test]
    fn parallel_inflate_answers_as_the_serial_loop_does() {
        let good = stream(1 << 20);
        assert_eq!(
            decompress_all(&good).unwrap(),
            decompress_serial(&good).unwrap()
        );

        // A corrupted deflate body in the third block, a truncated tail, and both at once: the
        // serial loop stops at the earliest, and so must the parallel one.
        let (first, _) = block_frame(&good).unwrap().unwrap();
        let (second, _) = block_frame(&good[first..]).unwrap().unwrap();
        let mut corrupt = good.clone();
        corrupt[first + second + BLOCK_HEADER_LENGTH + 3] ^= 0xff;
        let truncated = &good[..good.len() - 40];
        let mut both = corrupt.clone();
        both.truncate(both.len() - 40);
        for damaged in [&corrupt[..], truncated, &both[..]] {
            let parallel = decompress_all(damaged).unwrap_err().to_string();
            let serial = decompress_serial(damaged).unwrap_err().to_string();
            assert_eq!(parallel, serial);
        }
    }
}
