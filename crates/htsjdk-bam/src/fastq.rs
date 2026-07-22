//! FASTQ record formatting.
//!
//! Ports `htsjdk.samtools.fastq.FastqRecord`, `FastqEncoder` (`encode`/`write`/`asFastqRecord`),
//! `FastqConstants`, and the `BasicFastqWriter` newline behaviour, plus the two `SAMRecord`
//! accessors the conversion reads (`getReadString`, `getBaseQualityString` via
//! `SAMUtils.phredToFastq`), at tag 4.2.0. This is the output primitive for the read-to-FASTQ
//! transforms (`SamToFastq`) and the round trip back (`FastqToSam`).
//!
//! The record is four lines:
//!
//! ```text
//! @<name>
//! <sequence>
//! +<quality header>
//! <qualities>
//! ```
//!
//! `FastqEncoder.encode` produces exactly that with **no trailing newline**; `BasicFastqWriter`
//! adds one after each record, so a file is the records joined and terminated by newlines. A null
//! field renders as the empty string in the encoder. Two `SAMRecord` sentinels survive into the
//! conversion: empty read bases and empty qualities each render as `*` (`getReadString` /
//! `getBaseQualityString` return `"*"`), not as an empty line.

use crate::record::BamRecord;
use crate::tag::{Tag, TagValue};

const READ_PAIRED: u16 = 0x1;
const FIRST_OF_PAIR: u16 = 0x40;
const SECOND_OF_PAIR: u16 = 0x80;

/// `SAMUtils.MAX_PHRED_SCORE`: the largest score `phredToFastq` will encode.
const MAX_PHRED_SCORE: u8 = 93;

/// `SAMRecord.NULL_SEQUENCE_STRING` and `NULL_QUALS_STRING`.
const NULL_STRING: &str = "*";

/// `htsjdk.samtools.fastq.FastqRecord`, carrying the four fields the encoder writes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastqRecord {
    pub read_name: Option<String>,
    pub read_string: Option<String>,
    pub quality_header: Option<String>,
    pub quality_string: Option<String>,
}

impl FastqRecord {
    /// `FastqRecord.toFastQString()` / `FastqEncoder.encode`: the four lines, no trailing newline,
    /// with a null field written as the empty string.
    pub fn to_fastq_string(&self) -> String {
        format!(
            "@{}\n{}\n+{}\n{}",
            self.read_name.as_deref().unwrap_or(""),
            self.read_string.as_deref().unwrap_or(""),
            self.quality_header.as_deref().unwrap_or(""),
            self.quality_string.as_deref().unwrap_or(""),
        )
    }
}

/// `BasicFastqWriter.write`: the encoded record followed by a newline.
pub fn write_record(record: &FastqRecord) -> String {
    let mut s = record.to_fastq_string();
    s.push('\n');
    s
}

/// `SAMUtils.phredToFastq(byte[])`: each score offset by 33 into printable ASCII.
///
/// Panics on a score outside `0..=93`, matching htsjdk's `IllegalArgumentException`.
pub fn phred_to_fastq(quals: &[u8]) -> String {
    quals
        .iter()
        .map(|&q| {
            assert!(q <= MAX_PHRED_SCORE, "Cannot encode phred score: {q}");
            (33 + q) as char
        })
        .collect()
}

/// `SAMRecord.getReadString()`: the read bases as a string, or `*` when there are none.
fn read_string(rec: &BamRecord) -> String {
    if rec.read_bases.is_empty() {
        NULL_STRING.to_string()
    } else {
        String::from_utf8_lossy(&rec.read_bases).into_owned()
    }
}

/// `SAMRecord.getBaseQualityString()`: the FASTQ-encoded qualities, or `*` when there are none.
fn base_quality_string(rec: &BamRecord) -> String {
    if rec.base_qualities.is_empty() {
        NULL_STRING.to_string()
    } else {
        phred_to_fastq(&rec.base_qualities)
    }
}

/// `FastqEncoder.asFastqRecord(SAMRecord)`.
///
/// A paired read that is first or second of its pair gets `/1` or `/2` appended to its name. The
/// quality header comes from the `CO` (comment) tag, which is usually absent.
pub fn as_fastq_record(rec: &BamRecord) -> FastqRecord {
    let mut read_name = rec.read_name.clone();
    if rec.flags & READ_PAIRED != 0 {
        if rec.flags & FIRST_OF_PAIR != 0 {
            read_name.push_str("/1");
        } else if rec.flags & SECOND_OF_PAIR != 0 {
            read_name.push_str("/2");
        }
    }

    let quality_header = match rec.tags.get(Tag::new(b"CO")) {
        Some(TagValue::Str(s)) => Some(s.clone()),
        _ => None,
    };

    FastqRecord {
        read_name: Some(read_name),
        read_string: Some(read_string(rec)),
        quality_header,
        quality_string: Some(base_quality_string(rec)),
    }
}

/// `FastqEncoder.encode(SAMRecord)`: convert then format.
pub fn encode(rec: &BamRecord) -> String {
    as_fastq_record(rec).to_fastq_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rec(name: &str, flags: u16, bases: &[u8], quals: &[u8]) -> BamRecord {
        BamRecord {
            read_name: name.to_string(),
            flags,
            read_bases: bases.to_vec(),
            base_qualities: quals.to_vec(),
            ..Default::default()
        }
    }

    #[test]
    fn phred_encodes_by_adding_33() {
        // 0 -> '!', 30 -> '?', 40 -> 'I'.
        assert_eq!(phred_to_fastq(&[0, 30, 40]), "!?I");
    }

    #[test]
    #[should_panic(expected = "Cannot encode phred score")]
    fn phred_rejects_a_score_above_93() {
        phred_to_fastq(&[94]);
    }

    /// An unpaired record, verified against htsjdk's `FastqEncoder.encode(SAMRecord)`.
    #[test]
    fn an_unpaired_record_encodes_without_a_pair_suffix() {
        let r = rec("read1", 0, b"ACGT", &[40, 40, 30, 20]);
        assert_eq!(encode(&r), "@read1\nACGT\n+\nII?5");
    }

    /// A first-of-pair record gets `/1`; a second gets `/2`.
    #[test]
    fn a_paired_record_gets_its_pair_suffix() {
        let r1 = rec("read1", READ_PAIRED | FIRST_OF_PAIR, b"AC", &[40, 40]);
        assert_eq!(encode(&r1), "@read1/1\nAC\n+\nII");
        let r2 = rec("read1", READ_PAIRED | SECOND_OF_PAIR, b"GT", &[10, 10]);
        // qual 10 -> '+', so the qualities line is "++" under the empty "+" header line.
        assert_eq!(encode(&r2), "@read1/2\nGT\n+\n++");
    }

    #[test]
    fn empty_bases_and_quals_render_as_star() {
        let r = rec("empty", 0, b"", &[]);
        assert_eq!(encode(&r), "@empty\n*\n+\n*");
    }

    #[test]
    fn the_comment_tag_becomes_the_quality_header() {
        let mut r = rec("read1", 0, b"A", &[40]);
        r.tags
            .insert(Tag::new(b"CO"), TagValue::Str("note".to_string()));
        assert_eq!(encode(&r), "@read1\nA\n+note\nI");
    }

    #[test]
    fn the_writer_terminates_the_record_with_a_newline() {
        let r = as_fastq_record(&rec("r", 0, b"A", &[40]));
        assert_eq!(write_record(&r), "@r\nA\n+\nI\n");
    }
}
