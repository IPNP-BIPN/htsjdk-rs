//! `BlockCompressedInputStream.checkTermination`: inspect a BGZF file's tail for the empty-block
//! terminator, a healthy (non-terminator) last block, or corruption.
//!
//! Ported from htsjdk 4.2.0 `BlockCompressedInputStream.checkTermination(SeekableByteChannel)`. The
//! Java takes a seekable channel and reads only the tail; this port takes the whole file bytes,
//! which is what the one caller (`CheckTerminatorBlock`) has in hand, and reads the same tail.

use crate::{
    BGZF_ID1, BGZF_ID2, BGZF_LEN, EMPTY_GZIP_BLOCK, GZIP_CM_DEFLATE, GZIP_FLG, GZIP_ID1, GZIP_ID2,
    GZIP_OS_UNKNOWN, GZIP_XFL, GZIP_XLEN, MAX_COMPRESSED_BLOCK_SIZE,
};

/// `BlockCompressedInputStream.FileTermination`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FileTermination {
    /// The file ends with the standard empty-gzip-block terminator.
    HasTerminatorBlock,
    /// The last block is a valid BGZF block but not the terminator (an older, un-terminated file).
    HasHealthyLastBlock,
    /// The file is shorter than a terminator, or its last block is truncated / corrupt.
    Defective,
}

/// `BlockCompressedStreamConstants.GZIP_BLOCK_PREAMBLE`: the fixed 16-byte head of every BGZF block.
const GZIP_BLOCK_PREAMBLE: [u8; 16] = [
    GZIP_ID1,
    GZIP_ID2,
    GZIP_CM_DEFLATE,
    GZIP_FLG,
    0,
    0,
    0,
    0, // modification time
    GZIP_XFL,
    GZIP_OS_UNKNOWN,
    GZIP_XLEN as u8,
    0, // little-endian short
    BGZF_ID1,
    BGZF_ID2,
    BGZF_LEN as u8,
    0, // little-endian short
];

/// `checkTermination`: classify the tail of a BGZF file.
pub fn check_termination(data: &[u8]) -> FileTermination {
    let file_size = data.len();
    if file_size < EMPTY_GZIP_BLOCK.len() {
        return FileTermination::Defective;
    }

    // The end of the file is the empty gzip block used to terminate a bgzipped file.
    if data[file_size - EMPTY_GZIP_BLOCK.len()..] == EMPTY_GZIP_BLOCK {
        return FileTermination::HasTerminatorBlock;
    }

    // Otherwise, scan the tail backwards for the last block preamble and check its declared size.
    let bufsize = file_size.min(MAX_COMPRESSED_BLOCK_SIZE);
    let buf = &data[file_size - bufsize..];
    let mut i = buf.len() - EMPTY_GZIP_BLOCK.len();
    loop {
        if buf[i..i + GZIP_BLOCK_PREAMBLE.len()] == GZIP_BLOCK_PREAMBLE {
            // The 2-byte little-endian BSIZE (total block size minus one) follows the preamble.
            let off = i + GZIP_BLOCK_PREAMBLE.len();
            let total_block_size_minus_one = u16::from_le_bytes([buf[off], buf[off + 1]]) as usize;
            return if buf.len() - i == total_block_size_minus_one + 1 {
                FileTermination::HasHealthyLastBlock
            } else {
                FileTermination::Defective
            };
        }
        if i == 0 {
            break;
        }
        i -= 1;
    }
    FileTermination::Defective
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::BgzfWriter;

    /// A well-formed BGZF file (some blocks plus the empty terminator) reports the terminator block.
    #[test]
    fn a_terminated_file_has_a_terminator_block() {
        let mut w = BgzfWriter::new(Vec::new());
        std::io::Write::write_all(&mut w, &vec![7u8; 5000]).unwrap();
        w.finish().unwrap();
        let bytes = w.into_inner().unwrap();
        assert_eq!(
            check_termination(&bytes),
            FileTermination::HasTerminatorBlock
        );
    }

    /// The same file with its 28-byte terminator removed reports a healthy last block.
    #[test]
    fn a_file_without_the_terminator_has_a_healthy_last_block() {
        let mut w = BgzfWriter::new(Vec::new());
        std::io::Write::write_all(&mut w, &vec![7u8; 5000]).unwrap();
        w.finish().unwrap();
        let bytes = w.into_inner().unwrap();
        let trimmed = &bytes[..bytes.len() - EMPTY_GZIP_BLOCK.len()];
        assert_eq!(
            check_termination(trimmed),
            FileTermination::HasHealthyLastBlock
        );
    }

    /// A file whose last block is truncated mid-way is defective.
    #[test]
    fn a_truncated_last_block_is_defective() {
        let mut w = BgzfWriter::new(Vec::new());
        std::io::Write::write_all(&mut w, &vec![7u8; 5000]).unwrap();
        w.finish().unwrap();
        let bytes = w.into_inner().unwrap();
        // Drop the terminator and a few bytes of the real last block so its size no longer matches.
        let trimmed = &bytes[..bytes.len() - EMPTY_GZIP_BLOCK.len() - 5];
        assert_eq!(check_termination(trimmed), FileTermination::Defective);
    }

    /// A file shorter than a terminator block is defective.
    #[test]
    fn a_tiny_file_is_defective() {
        assert_eq!(check_termination(&[0u8; 10]), FileTermination::Defective);
    }
}
