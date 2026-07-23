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

use std::io::{self, Read};

use htsjdk_bgzf::{vfp, BgzfError, BgzfReader, BgzfWriter};

use crate::header::SamHeader;
use crate::writer::{write_bam_header_block, BAM_MAGIC};

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

/// `reheaderBamFile(newHeader, input, skipHeader=true, skipTerminator=false)`: the whole reheadered
/// BAM as bytes. `input_bam` is the raw BAM file (BGZF-framed, terminator included).
pub fn reheader_bam(new_header: &SamHeader, input_bam: &[u8]) -> Result<Vec<u8>, ReheaderError> {
    // 1. The new header as flushed BGZF block(s), no terminator.
    let mut out = write_bam_header_block(new_header)?;

    // 2. Locate the first record: consume the input header, then read the virtual file pointer.
    let mut finder = BgzfReader::new(input_bam);
    consume_header(&mut finder)?;
    let vptr = finder.virtual_pos();
    let first_block_addr = vfp::block_address(vptr);
    let first_offset = vfp::block_offset(vptr) as usize;

    // Find the decoded block the first record starts in, to grab its tail and its successor's
    // compressed address. A fresh reader is used because the finder consumed part of the stream.
    let mut scan = BgzfReader::new(input_bam);
    let mut tail: Vec<u8> = Vec::new();
    let mut next_block_addr: Option<u64> = None;
    while let Some(block) = scan.next_block()? {
        if block.block_address == first_block_addr {
            tail = block.data[first_offset..].to_vec();
            next_block_addr = Some(block.block_address + block.block_compressed_size as u64);
            break;
        }
    }
    let next_block_addr = next_block_addr.ok_or(ReheaderError::BlockNotFound)?;

    // 3. Re-compress the first-record tail as one flushed block (empty tail writes nothing, as
    // htsjdk's flush() of a zero-byte BlockCompressedOutputStream writes no block).
    if !tail.is_empty() {
        let mut w = BgzfWriter::new(Vec::new());
        io::Write::write_all(&mut w, &tail)?;
        out.extend_from_slice(&w.into_inner_without_terminator()?);
    }

    // 4. Raw-copy the remaining compressed blocks, including the terminator (skipTerminator=false).
    out.extend_from_slice(&input_bam[next_block_addr as usize..]);
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reader::BamReader;
    use crate::sam_file::read_sam_with;
    use crate::text_parse::ValidationStringency;
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
