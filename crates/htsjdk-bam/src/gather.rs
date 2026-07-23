//! `BamFileIoUtils.gatherWithBlockCopying`: concatenate BAM files by copying their compressed
//! blocks, keeping the first file's header and dropping the rest.
//!
//! Ported from `htsjdk.samtools.BamFileIoUtils.gatherWithBlockCopying` at tag 4.2.0. This is the
//! fast path `GatherBamFiles` takes when every input is a BAM file (`determineBlockCopyingStatus`);
//! the output is:
//!
//! - the first input block-copied whole, header included, terminator dropped
//!   (`blockCopyBamFile(skipHeader=false, skipTerminator=true)`);
//! - each subsequent input block-copied with its header and terminator dropped
//!   (`blockCopyBamFile(skipHeader=true, skipTerminator=true)`), so only its record blocks are kept;
//! - a single EOF terminator block.
//!
//! Byte-identity holds against Picard with `USE_JDK_DEFLATER=true`, transitively through
//! [`block_copy`](crate::reheader::block_copy).

use htsjdk_bgzf::EMPTY_GZIP_BLOCK;

use crate::reheader::{block_copy, ReheaderError};

/// `gatherWithBlockCopying(inputs, output)`: the concatenated BAM as bytes. Each element of `bams`
/// is a whole raw BAM file (BGZF-framed, terminator included). The first file's header is kept and
/// the others' are dropped, as htsjdk does; callers that need the headers reconciled must do that
/// upstream (Picard requires the inputs already share a header for this fast path to be correct).
pub fn gather_bam_files(bams: &[&[u8]]) -> Result<Vec<u8>, ReheaderError> {
    let mut out = Vec::new();
    for (i, bam) in bams.iter().enumerate() {
        let skip_header = i != 0;
        out.extend_from_slice(&block_copy(bam, skip_header, true)?);
    }
    // And lastly add the terminator block.
    out.extend_from_slice(&EMPTY_GZIP_BLOCK);
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reader::BamReader;
    use crate::sam_file::read_sam_with;
    use crate::text_parse::ValidationStringency;
    use crate::writer::BamWriter;

    const HEADER: &str = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n";

    fn build_bam(reads: &str) -> Vec<u8> {
        let sam = format!("{HEADER}{reads}");
        let (header, records) = read_sam_with(&sam, ValidationStringency::Lenient).unwrap();
        let mut w = BamWriter::new(Vec::new(), &header).unwrap();
        for r in &records {
            w.write(r).unwrap();
        }
        w.finish().unwrap()
    }

    #[test]
    fn gathering_concatenates_records_under_the_first_header() {
        let a = build_bam("a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n");
        let b = build_bam("b\t0\tchr1\t20\t60\t4M\t*\t0\t0\tTTTT\tIIII\n");
        let c = build_bam("c\t0\tchr1\t30\t60\t4M\t*\t0\t0\tGGGG\tIIII\n");

        let gathered = gather_bam_files(&[&a, &b, &c]).unwrap();

        let decoded = htsjdk_bgzf::decompress_all(&gathered).unwrap();
        let reader = BamReader::new(&decoded).unwrap();
        assert!(reader.header.raw_text.contains("@SQ\tSN:chr1\tLN:1000"));
        let names: Vec<_> = reader.map(|r| r.unwrap().read_name).collect();
        assert_eq!(
            names,
            vec!["a".to_string(), "b".to_string(), "c".to_string()]
        );
    }

    #[test]
    fn a_single_input_round_trips() {
        let a = build_bam("a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n");
        let gathered = gather_bam_files(&[&a]).unwrap();
        let decoded = htsjdk_bgzf::decompress_all(&gathered).unwrap();
        let reader = BamReader::new(&decoded).unwrap();
        let names: Vec<_> = reader.map(|r| r.unwrap().read_name).collect();
        assert_eq!(names, vec!["a".to_string()]);
    }

    #[test]
    fn an_input_with_no_records_contributes_nothing_but_its_header() {
        let empty = build_bam("");
        let a = build_bam("a\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n");
        // Empty file first: its header is kept, and the second file's records follow.
        let gathered = gather_bam_files(&[&empty, &a]).unwrap();
        let decoded = htsjdk_bgzf::decompress_all(&gathered).unwrap();
        let reader = BamReader::new(&decoded).unwrap();
        let names: Vec<_> = reader.map(|r| r.unwrap().read_name).collect();
        assert_eq!(names, vec!["a".to_string()]);
    }
}
