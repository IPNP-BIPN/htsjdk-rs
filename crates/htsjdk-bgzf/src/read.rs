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
        if block_length < BLOCK_HEADER_LENGTH || block_length > MAX_COMPRESSED_BLOCK_SIZE {
            return Err(BgzfError::UnexpectedBlockLength(block_length));
        }

        let remaining = block_length - BLOCK_HEADER_LENGTH;
        let mut rest = vec![0u8; remaining];
        let got = read_up_to(&mut self.inner, &mut rest).map_err(|_| BgzfError::PrematureEnd)?;
        self.stream_offset += got as u64;
        if got != remaining {
            return Err(BgzfError::PrematureEnd);
        }

        // Footer: CRC32 then ISIZE, both little-endian.
        let footer = &rest[remaining - BLOCK_FOOTER_LENGTH..];
        let expected_crc = u32::from_le_bytes([footer[0], footer[1], footer[2], footer[3]]);
        let uncompressed_len =
            u32::from_le_bytes([footer[4], footer[5], footer[6], footer[7]]) as usize;

        let deflated = &rest[..remaining - BLOCK_FOOTER_LENGTH];
        let mut data = Vec::with_capacity(uncompressed_len);
        let mut d = Decompress::new(false);
        d.decompress_vec(deflated, &mut data, FlushDecompress::Finish)
            .map_err(|e| BgzfError::Inflate(e.to_string()))?;

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
