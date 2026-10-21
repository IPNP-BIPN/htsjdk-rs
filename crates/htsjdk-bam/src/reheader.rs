//! `BamFileIoUtils.reheaderBamFile` / `blockCopyBamFile`: replace a BAM's header while copying the
//! record data blocks verbatim, so the output is byte-identical to htsjdk's block-copy reheader
//! rather than to a full re-encode.
//!
//! Ported from `htsjdk.samtools.BamFileIoUtils` (with `SAMUtils.findVirtualOffsetOfFirstRecordInBam`
//! and `BAMFileReader.findVirtualOffsetOfFirstRecord`) at tag 4.2.0, for the `skipHeader=true,
//! skipTerminator=false` shape `reheaderBamFile` uses. The output is:
//!
//! 1. [`write_bam_header_block`](crate::writer::write_bam_header_block) of the new header: fresh
//!    BGZF block(s), flushed, no terminator (`BAMFileWriter.writeHeader(OutputStream)`).
//! 2. The tail of the input's first-record block re-compressed as one new block. htsjdk seeks to the
//!    virtual offset of the first record, reads the `available()` bytes remaining in that
//!    decompressed block, and flushes them through a new `BlockCompressedOutputStream`.
//! 3. A raw byte copy of the input from the next block boundary to end of file, including the
//!    terminator (`skipTerminator=false`). These compressed blocks are copied without decoding.
//!
//! Byte-identity holds against Picard run with `USE_JDK_DEFLATER=true` (java.util.zip), matching the
//! rest of htsjdk-rs's BGZF writes; the default GKL/igzip deflater is a separate surface.

use std::io::{self, BufReader, Cursor, Read, Seek, SeekFrom, Write};

use htsjdk_bgzf::{
    check_termination, vfp, BgzfError, BgzfReader, BgzfWriter, Deflater, FileTermination,
    DEFAULT_COMPRESSION_LEVEL, EMPTY_GZIP_BLOCK, MAX_COMPRESSED_BLOCK_SIZE,
};

use crate::header::SamHeader;
use crate::writer::{write_bam_header_block_with, BAM_MAGIC};

/// Why a block-copy reheader could not run.
#[derive(Debug)]
pub enum ReheaderError {
    /// The input's magic was not `BAM\1`: it is not a BAM (block-copy does not support SAM).
    NotABam,
    /// A BGZF decode error while reading the input header or locating the first record.
    Bgzf(BgzfError),
    /// The input ended before a full header could be read.
    Truncated,
    /// The virtual offset of the first record did not land on a decoded block boundary.
    BlockNotFound,
    /// The input has no valid GZIP block at the end of the file (`checkTermination == DEFECTIVE`).
    Defective,
    /// Writing the output failed.
    Write(io::Error),
}

impl From<BgzfError> for ReheaderError {
    fn from(e: BgzfError) -> Self {
        ReheaderError::Bgzf(e)
    }
}

impl From<io::Error> for ReheaderError {
    fn from(_: io::Error) -> Self {
        // The only reader here is over an in-memory slice, so an I/O error is a short read.
        ReheaderError::Truncated
    }
}

/// Reads a little-endian `i32` from the decompressed BGZF stream.
fn read_i32<R: Read>(r: &mut R) -> Result<i32, ReheaderError> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b)?;
    Ok(i32::from_le_bytes(b))
}

/// Consumes exactly the BAM header binary content (`BAMFileReader.readHeader`: magic, length-prefixed
/// text, then the binary dictionary) from `reader`, so that `reader.virtual_pos()` afterwards is the
/// virtual offset of the first record, as `mCompressedInputStream.getFilePointer()` is in htsjdk.
fn consume_header<R: Read>(reader: &mut BgzfReader<R>) -> Result<(), ReheaderError> {
    let mut magic = [0u8; 4];
    reader.read_exact(&mut magic)?;
    if magic != BAM_MAGIC {
        return Err(ReheaderError::NotABam);
    }
    let l_text = read_i32(reader)?;
    io::copy(&mut reader.by_ref().take(l_text as u64), &mut io::sink())?;
    let n_ref = read_i32(reader)?;
    for _ in 0..n_ref {
        let l_name = read_i32(reader)?;
        io::copy(&mut reader.by_ref().take(l_name as u64), &mut io::sink())?;
        let _l_ref = read_i32(reader)?;
    }
    Ok(())
}

/// `BamFileIoUtils.blockCopyBamFile(input, out, skipHeader, skipTerminator)`: copy a BAM's
/// compressed blocks into a new byte buffer without decoding the records.
///
/// - `skip_header`: drop the header. htsjdk seeks to the first record's virtual offset, re-compresses
///   the bytes remaining in that decoded block as one new block, then raw-copies from the next block.
///   With `skip_header=false` the whole file is raw-copied from offset 0, header included.
/// - `skip_terminator`: drop the trailing EOF block if the input has one, so blocks from a following
///   file (or a freshly written terminator) can be appended.
pub fn block_copy(
    input_bam: &[u8],
    skip_header: bool,
    skip_terminator: bool,
) -> Result<Vec<u8>, ReheaderError> {
    block_copy_with(
        input_bam,
        skip_header,
        skip_terminator,
        DEFAULT_COMPRESSION_LEVEL,
        Deflater::Jdk,
    )
}

/// [`block_copy`] with the compression of the one re-compressed block named: htsjdk writes it
/// through a `BlockCompressedOutputStream`, which takes the static default level and deflater.
pub fn block_copy_with(
    input_bam: &[u8],
    skip_header: bool,
    skip_terminator: bool,
    level: u32,
    deflater: Deflater,
) -> Result<Vec<u8>, ReheaderError> {
    let mut out = Vec::new();
    block_copy_to(
        &mut Cursor::new(input_bam),
        &mut out,
        skip_header,
        skip_terminator,
        level,
        deflater,
    )?;
    Ok(out)
}

/// [`block_copy_with`] from a seekable input to `out`, as it is read: only the block the first
/// record starts in is held, and the blocks after it are copied straight through, so a file of any
/// size is copied in constant memory.
pub fn block_copy_to<R: Read + Seek, W: Write>(
    input: &mut R,
    out: &mut W,
    skip_header: bool,
    skip_terminator: bool,
    level: u32,
    deflater: Deflater,
) -> Result<(), ReheaderError> {
    // checkTermination: reject a file whose tail is not a valid GZIP block, and learn whether an
    // EOF terminator is present (so skip_terminator can drop exactly it). It reads no more than the
    // last block's worth of the file.
    let length = input.seek(SeekFrom::End(0))?;
    let tail_length = length.min(MAX_COMPRESSED_BLOCK_SIZE as u64);
    input.seek(SeekFrom::Start(length - tail_length))?;
    let mut tail = vec![0; tail_length as usize];
    input.read_exact(&mut tail)?;
    let term = check_termination(&tail);
    if term == FileTermination::Defective {
        return Err(ReheaderError::Defective);
    }

    let copy_start = if skip_header {
        // Locate the first record: consume the header, then read the virtual file pointer.
        input.seek(SeekFrom::Start(0))?;
        let mut finder = BgzfReader::new(BufReader::new(&mut *input));
        consume_header(&mut finder)?;
        let vptr = finder.virtual_pos();
        drop(finder);
        let first_block_addr = vfp::block_address(vptr);
        let first_offset = vfp::block_offset(vptr) as usize;

        // Find the decoded block the first record starts in, to grab its tail and its successor's
        // compressed address. A fresh reader is used because the finder consumed part of the stream.
        input.seek(SeekFrom::Start(0))?;
        let mut scan = BgzfReader::new(BufReader::new(&mut *input));
        let mut tail: Vec<u8> = Vec::new();
        let mut next_block_addr: Option<u64> = None;
        while let Some(block) = scan.next_block()? {
            if block.block_address == first_block_addr {
                tail = block.data[first_offset..].to_vec();
                next_block_addr = Some(block.block_address + block.block_compressed_size as u64);
                break;
            }
        }
        drop(scan);
        let next_block_addr = next_block_addr.ok_or(ReheaderError::BlockNotFound)?;

        // Re-compress the first-record tail as one flushed block (empty tail writes nothing, as
        // htsjdk's flush() of a zero-byte BlockCompressedOutputStream writes no block).
        if !tail.is_empty() {
            let mut w = BgzfWriter::with_deflater(Vec::new(), level, deflater);
            io::Write::write_all(&mut w, &tail)?;
            out.write_all(&w.into_inner_without_terminator()?)
                .map_err(ReheaderError::Write)?;
        }
        next_block_addr
    } else {
        0
    };

    let skip_last = if term == FileTermination::HasTerminatorBlock && skip_terminator {
        EMPTY_GZIP_BLOCK.len() as u64
    } else {
        0
    };
    // htsjdk copies `length - skipLast - pos` bytes. For an input with no records the first-record
    // pointer lands at EOF (`pos == length`), so with the terminator dropped this count goes
    // negative and `transferByStream` copies nothing; guard the count to match.
    let end = length - skip_last;
    if copy_start < end {
        input.seek(SeekFrom::Start(copy_start))?;
        let copied = io::copy(&mut input.by_ref().take(end - copy_start), out)
            .map_err(ReheaderError::Write)?;
        if copied != end - copy_start {
            return Err(ReheaderError::Truncated);
        }
    }
    Ok(())
}

/// `reheaderBamFile(newHeader, input, skipHeader=true, skipTerminator=false)`: the whole reheadered
/// BAM as bytes. `input_bam` is the raw BAM file (BGZF-framed, terminator included).
pub fn reheader_bam(new_header: &SamHeader, input_bam: &[u8]) -> Result<Vec<u8>, ReheaderError> {
    reheader_bam_with(
        new_header,
        input_bam,
        DEFAULT_COMPRESSION_LEVEL,
        Deflater::Jdk,
    )
}

/// [`reheader_bam`] with the compression of the blocks it writes named, which are the header and
/// the first record block's tail; every other block is copied as it was.
pub fn reheader_bam_with(
    new_header: &SamHeader,
    input_bam: &[u8],
    level: u32,
    deflater: Deflater,
) -> Result<Vec<u8>, ReheaderError> {
    // The new header as flushed BGZF block(s), no terminator, then the record blocks block-copied
    // with the header dropped and the terminator kept.
    let mut out = Vec::with_capacity(input_bam.len());
    reheader_bam_to(
        new_header,
        &mut Cursor::new(input_bam),
        &mut out,
        level,
        deflater,
    )?;
    Ok(out)
}

/// [`reheader_bam_with`] from a seekable input to `out`, as it is read, which is how htsjdk's
/// `reheaderBamFile` goes from one path to another: the new header's blocks, then the records'
/// blocks copied through.
pub fn reheader_bam_to<R: Read + Seek, W: Write>(
    new_header: &SamHeader,
    input: &mut R,
    out: &mut W,
    level: u32,
    deflater: Deflater,
) -> Result<(), ReheaderError> {
    out.write_all(&write_bam_header_block_with(new_header, level, deflater)?)
        .map_err(ReheaderError::Write)?;
    block_copy_to(input, out, true, false, level, deflater)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reader::BamReader;
    use crate::sam_file::read_sam_with;
    use crate::text_parse::ValidationStringency;
    use crate::writer::write_bam_header_block;
    use crate::writer::BamWriter;

    fn build_bam(sam: &str) -> Vec<u8> {
        let (header, records) = read_sam_with(sam, ValidationStringency::Lenient).unwrap();
        let mut w = BamWriter::new(Vec::new(), &header).unwrap();
        for r in &records {
            w.write(r).unwrap();
        }
        w.finish().unwrap()
    }

    const SAM: &str = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n\
        a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n\
        b\t0\tchr1\t20\t60\t4M\t*\t0\t0\tTTTT\tIIII\n";

    #[test]
    fn reheader_preserves_records_and_applies_the_new_header() {
        let input = build_bam(SAM);
        let (mut header, _) = read_sam_with(SAM, ValidationStringency::Lenient).unwrap();
        header.add_comment("a comment");
        header.add_comment("second comment");

        let out = reheader_bam(&header, &input).unwrap();

        // The output is a valid BAM: its header carries the new comments, and its records are
        // exactly the input's, unchanged by the block copy.
        let decoded = decompress_for_test(&out);
        let reader = BamReader::new(&decoded).unwrap();
        assert!(reader.header.raw_text.contains("@CO\ta comment"));
        assert!(reader.header.raw_text.contains("@CO\tsecond comment"));
        let out_recs: Vec<_> = reader.map(|r| r.unwrap().read_name).collect();
        assert_eq!(out_recs, vec!["a".to_string(), "b".to_string()]);
    }

    #[test]
    fn a_bgzf_stream_without_the_bam_magic_is_rejected() {
        // Valid BGZF, but the decompressed content is not a BAM: the magic check must reject it.
        let mut w = BgzfWriter::new(Vec::new());
        io::Write::write_all(&mut w, b"XXXXnot bam content").unwrap();
        let framed = w.into_inner().unwrap();
        let (header, _) = read_sam_with(SAM, ValidationStringency::Lenient).unwrap();
        assert!(matches!(
            reheader_bam(&header, &framed),
            Err(ReheaderError::NotABam)
        ));
    }

    #[test]
    fn raw_garbage_is_rejected() {
        let (header, _) = read_sam_with(SAM, ValidationStringency::Lenient).unwrap();
        assert!(reheader_bam(&header, b"not a bam at all").is_err());
    }

    /// BamReader consumes already-decompressed BAM bytes, so decode the framed output for the test.
    fn decompress_for_test(framed: &[u8]) -> Vec<u8> {
        htsjdk_bgzf::decompress_all(framed).unwrap()
    }

    #[test]
    fn write_bam_header_block_has_no_terminator() {
        let (header, _) = read_sam_with(SAM, ValidationStringency::Lenient).unwrap();
        let block = write_bam_header_block(&header).unwrap();
        assert!(!block.ends_with(&htsjdk_bgzf::EMPTY_GZIP_BLOCK));
        // It decodes to exactly the header binary content, starting with the magic.
        let decoded = htsjdk_bgzf::decompress_all(&block).unwrap();
        assert_eq!(&decoded[..4], &BAM_MAGIC);
    }
}
