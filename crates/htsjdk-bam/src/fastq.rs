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
const READ_UNMAPPED: u16 = 0x4;
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

/// `SAMUtils.fastqToPhred(String)`: each printable FASTQ character back to a binary score.
///
/// Panics on a character outside `33..=126`, matching htsjdk's `IllegalArgumentException`. The
/// inverse of [`phred_to_fastq`].
pub fn fastq_to_phred(fastq: &str) -> Vec<u8> {
    fastq
        .bytes()
        .map(|c| {
            assert!(
                (33..=126).contains(&c),
                "Invalid fastq character: {}",
                c as char
            );
            c - 33
        })
        .collect()
}

/// `SequenceUtil.getSamReadNameFromFastqHeader`: the SAM read name for a FASTQ header.
///
/// The header is truncated at the first space, then any trailing `/1` or `/2` pair suffixes are
/// stripped (in a loop, to trap the pathological `/1/1`), because a `/1` left on an unpaired read
/// causes trouble downstream in tools like MergeBamAlignment.
pub fn get_sam_read_name_from_fastq_header(fastq_header: &str) -> String {
    let mut name = match fastq_header.find(' ') {
        Some(idx) => &fastq_header[..idx],
        None => fastq_header,
    }
    .to_string();
    while name.ends_with("/1") || name.ends_with("/2") {
        name.truncate(name.len() - 2);
    }
    name
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

impl FastqRecord {
    /// `FastqRecord.getReadBases()`: the sequence as bytes, or empty when there is none.
    pub fn get_read_bases(&self) -> Vec<u8> {
        match &self.read_string {
            None => Vec::new(),
            Some(s) => s.as_bytes().to_vec(),
        }
    }

    /// `FastqRecord.getBaseQualities()`: the FASTQ qualities decoded to binary scores.
    pub fn get_base_qualities(&self) -> Vec<u8> {
        match &self.quality_string {
            None => Vec::new(),
            Some(s) => fastq_to_phred(s),
        }
    }
}

/// `FastqEncoder.asSAMRecord(FastqRecord, header)`, without the custom hook.
///
/// Builds an **unmapped** record: the read name is the FASTQ header cleaned by
/// [`get_sam_read_name_from_fastq_header`], and the bases and qualities are decoded from the
/// record. This is the conversion `FastqToSam` performs on every input read.
pub fn as_sam_record(record: &FastqRecord) -> BamRecord {
    BamRecord {
        read_name: get_sam_read_name_from_fastq_header(record.read_name.as_deref().unwrap_or("")),
        flags: READ_UNMAPPED,
        read_bases: record.get_read_bases(),
        base_qualities: record.get_base_qualities(),
        ..Default::default()
    }
}

/// `StringUtil.isBlank`: null or all-whitespace. Here the null case is the caller's `None`.
fn is_blank(line: &str) -> bool {
    line.trim().is_empty()
}

/// `htsjdk.samtools.fastq.FastqReader`, as a parser over already-read text.
///
/// Ports `readNextRecord`, `checkLine`, and `readLineConditionallySkippingBlanks`. Blank-line
/// skipping is off by default in htsjdk (`new FastqReader(file)` passes `skipBlankLines = false`),
/// so it is a parameter here. A record is four lines: an `@` header, the sequence, a `+` header,
/// and the qualities, with the sequence and quality lines required to be the same length. The
/// record's name and quality header are the header lines with their leading `@`/`+` removed.
pub struct FastqReader<'a> {
    lines: std::str::Lines<'a>,
    skip_blank_lines: bool,
    line: usize,
}

impl<'a> FastqReader<'a> {
    pub fn new(text: &'a str, skip_blank_lines: bool) -> Self {
        FastqReader {
            lines: text.lines(),
            skip_blank_lines,
            line: 0,
        }
    }

    /// `readLineConditionallySkippingBlanks`: the next line, skipping blanks only when configured to.
    fn next_line(&mut self) -> Option<&'a str> {
        loop {
            let line = self.lines.next()?;
            if !(self.skip_blank_lines && is_blank(line)) {
                return Some(line);
            }
        }
    }

    /// `checkLine`: a line that is missing (`None`) or blank is an error.
    fn check_line(&self, line: Option<&'a str>, kind: &str) -> Result<&'a str, FastqError> {
        match line {
            None => Err(FastqError(format!("File is too short - missing {kind}"))),
            Some(l) if is_blank(l) => Err(FastqError(format!("Missing {kind}"))),
            Some(l) => Ok(l),
        }
    }

    /// `readNextRecord`: the next record, `None` at end of input.
    pub fn next_record(&mut self) -> Result<Option<FastqRecord>, FastqError> {
        let seq_header = match self.next_line() {
            None => return Ok(None),
            Some(h) => h,
        };
        if is_blank(seq_header) {
            return Err(FastqError("Missing sequence header".to_string()));
        }
        if !seq_header.starts_with('@') {
            return Err(FastqError(format!(
                "Sequence header must start with @: {seq_header}"
            )));
        }

        let seq_line = self.next_line();
        let seq_line = self.check_line(seq_line, "SequenceLine")?;

        let qual_header = self.next_line();
        let qual_header = self.check_line(qual_header, "QualityHeader")?;
        if !qual_header.starts_with('+') {
            return Err(FastqError(format!(
                "Quality header must start with +: {qual_header}"
            )));
        }

        let qual_line = self.next_line();
        let qual_line = self.check_line(qual_line, "QualityLine")?;

        if seq_line.len() != qual_line.len() {
            return Err(FastqError(
                "Sequence and quality line must be the same length".to_string(),
            ));
        }

        self.line += 4;
        Ok(Some(FastqRecord {
            read_name: Some(seq_header[1..].to_string()),
            read_string: Some(seq_line.to_string()),
            quality_header: Some(qual_header[1..].to_string()),
            quality_string: Some(qual_line.to_string()),
        }))
    }
}

/// An error while parsing FASTQ, carrying htsjdk's message (without its line-number suffix).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastqError(pub String);

impl std::fmt::Display for FastqError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

/// Parse an entire FASTQ text into its records.
pub fn parse_fastq(text: &str, skip_blank_lines: bool) -> Result<Vec<FastqRecord>, FastqError> {
    let mut reader = FastqReader::new(text, skip_blank_lines);
    let mut out = Vec::new();
    while let Some(rec) = reader.next_record()? {
        out.push(rec);
    }
    Ok(out)
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

    #[test]
    fn a_two_record_file_parses_into_two_records() {
        let text = "@r1 desc\nACGT\n+\nIIII\n@r2\nTT\n+r2\n##\n";
        let recs = parse_fastq(text, false).unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].read_name.as_deref(), Some("r1 desc"));
        assert_eq!(recs[0].read_string.as_deref(), Some("ACGT"));
        assert_eq!(recs[0].quality_header.as_deref(), Some(""));
        assert_eq!(recs[0].quality_string.as_deref(), Some("IIII"));
        // The quality header keeps whatever followed the '+'.
        assert_eq!(recs[1].quality_header.as_deref(), Some("r2"));
    }

    #[test]
    fn a_header_without_the_at_sign_is_rejected() {
        let err = parse_fastq("r1\nACGT\n+\nIIII\n", false).unwrap_err();
        assert!(err.0.contains("must start with @"), "{}", err.0);
    }

    #[test]
    fn mismatched_sequence_and_quality_lengths_are_rejected() {
        let err = parse_fastq("@r\nACGT\n+\nII\n", false).unwrap_err();
        assert!(err.0.contains("same length"), "{}", err.0);
    }

    #[test]
    fn a_truncated_record_is_too_short() {
        let err = parse_fastq("@r\nACGT\n", false).unwrap_err();
        assert!(err.0.contains("too short"), "{}", err.0);
    }

    /// With blank-line skipping off (htsjdk's default), a blank line where the sequence is expected
    /// is a missing-line error, not a skipped line.
    #[test]
    fn a_blank_line_is_significant_by_default() {
        let err = parse_fastq("@r\n\n+\n\n", false).unwrap_err();
        assert!(err.0.contains("Missing SequenceLine"), "{}", err.0);
        // With skipping on, the blanks are consumed and the record is too short instead.
        let err2 = parse_fastq("@r\n\n+\n\n", true).unwrap_err();
        assert!(err2.0.contains("too short"), "{}", err2.0);
    }

    #[test]
    fn fastq_to_phred_is_the_inverse_of_phred_to_fastq() {
        assert_eq!(fastq_to_phred("!?I"), vec![0, 30, 40]);
        assert_eq!(phred_to_fastq(&fastq_to_phred("ABCabc")), "ABCabc");
    }

    /// Verified against htsjdk's SequenceUtil.getSamReadNameFromFastqHeader.
    #[test]
    fn the_sam_read_name_drops_the_comment_and_pair_suffix() {
        assert_eq!(get_sam_read_name_from_fastq_header("read1/1"), "read1");
        assert_eq!(get_sam_read_name_from_fastq_header("read1/2"), "read1");
        assert_eq!(
            get_sam_read_name_from_fastq_header("read1 comment"),
            "read1"
        );
        assert_eq!(get_sam_read_name_from_fastq_header("read1"), "read1");
        // The loop strips a pathological doubled suffix.
        assert_eq!(get_sam_read_name_from_fastq_header("r/1/2"), "r");
        // A space is truncated before the suffix is considered, so a suffix after a space is gone
        // with the comment.
        assert_eq!(get_sam_read_name_from_fastq_header("read1 x/1"), "read1");
    }

    #[test]
    fn as_sam_record_builds_an_unmapped_read() {
        let fq = FastqRecord {
            read_name: Some("read1/1".to_string()),
            read_string: Some("ACGT".to_string()),
            quality_header: Some(String::new()),
            quality_string: Some("II?5".to_string()),
        };
        let sam = as_sam_record(&fq);
        assert_eq!(sam.read_name, "read1"); // suffix stripped
        assert_eq!(sam.flags, READ_UNMAPPED);
        assert_eq!(sam.read_bases, b"ACGT");
        assert_eq!(sam.base_qualities, vec![40, 40, 30, 20]);
    }

    /// The full SAM -> FASTQ -> SAM round trip preserves bases and qualities, and the pair suffix
    /// that encoding added is removed again by the name cleanup.
    #[test]
    fn sam_to_fastq_to_sam_round_trips_bases_and_quals() {
        let original = rec(
            "read1",
            READ_PAIRED | FIRST_OF_PAIR,
            b"ACGTN",
            &[40, 30, 20, 10, 0],
        );
        let fastq_text = write_record(&as_fastq_record(&original));
        let parsed = &parse_fastq(&fastq_text, false).unwrap()[0];
        let back = as_sam_record(parsed);
        assert_eq!(back.read_name, "read1");
        assert_eq!(back.read_bases, original.read_bases);
        assert_eq!(back.base_qualities, original.base_qualities);
        assert_eq!(back.flags, READ_UNMAPPED);
    }

    /// The encoder and reader are inverse on the fields the reader populates.
    #[test]
    fn encode_then_parse_round_trips() {
        let original = FastqRecord {
            read_name: Some("read1/1".to_string()),
            read_string: Some("ACGTN".to_string()),
            quality_header: Some(String::new()),
            quality_string: Some("IIII#".to_string()),
        };
        let text = write_record(&original);
        let parsed = parse_fastq(&text, false).unwrap();
        assert_eq!(parsed, vec![original]);
    }
}
