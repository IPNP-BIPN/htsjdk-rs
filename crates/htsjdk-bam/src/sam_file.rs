//! Whole SAM text files: header plus records.
//!
//! Ported from `htsjdk.samtools.SAMTextWriter.writeHeader` and `SAMTextReader`.
//!
//! The file layer is four lines of Java and one property worth checking. `writeHeader` writes
//! the header text **verbatim**, and `SAMTextHeaderCodec.encode` already terminates every line
//! with a newline. So the writer adds nothing between the header and the first record, and a
//! writer that inserted a separator would produce a file with a blank line that most readers
//! tolerate and no byte comparison accepts.

use crate::header::SamHeader;
use crate::record::BamRecord;
use crate::text::write_alignment;
use crate::text_parse::{parse_line_with, ParseError, ValidationStringency};

/// Writes a complete SAM text file.
///
/// `reference_name` resolves a record's reference index against the header, because a record
/// carries indices and the text form carries names.
pub fn write_sam(header: &SamHeader, records: &[BamRecord]) -> Option<String> {
    let mut out = header.encode();
    let name_of = |index: i32| -> &str {
        if index < 0 {
            "*"
        } else {
            header
                .sequences
                .get(index as usize)
                .map(|s| s.name.as_str())
                .unwrap_or("*")
        }
    };
    for (done, rec) in records.iter().enumerate() {
        let line = write_alignment(
            rec,
            name_of(rec.reference_index),
            name_of(rec.mate_reference_index),
        )?;
        // Room for the rest from the first line, and a quarter more. A page of it nothing writes
        // costs nothing, where growing by doubling would hold the old text and twice the new.
        if done == 0 {
            out.reserve((line.len() + 1) * records.len() / 4 * 5);
        }
        out.push_str(&line);
        out.push('\n');
    }
    Some(out)
}

/// Reads a complete SAM text file.
///
/// Header lines are those starting with `@`. htsjdk stops treating lines as header at the
/// first non-`@` line rather than filtering throughout, so an `@`-prefixed line *after* a
/// record is a record, not a header line. Reproduced: the split is positional.
pub fn read_sam(text: &str) -> Result<(SamHeader, Vec<BamRecord>), ParseError> {
    read_sam_with(text, ValidationStringency::default())
}

/// [`read_sam`] at an explicit stringency.
///
/// Needed because htsjdk's writer emits records its own default-stringency reader rejects; see
/// [`ValidationStringency`].
pub fn read_sam_with(
    text: &str,
    stringency: ValidationStringency,
) -> Result<(SamHeader, Vec<BamRecord>), ParseError> {
    let mut lines = text.lines().peekable();
    let mut header_text = String::new();
    while let Some(line) = lines.next_if(|line| line.starts_with('@')) {
        header_text.push_str(line);
        header_text.push('\n');
    }

    let header = crate::reader::parse_header_text(&header_text);
    let index_of = |name: &str| -> Option<i32> {
        header
            .sequences
            .iter()
            .position(|s| s.name == name)
            .map(|i| i as i32)
    };
    // One record per line, known before the first: a growing vector would hold each size and twice
    // that as it moved.
    let mut records = Vec::with_capacity(text.bytes().filter(|byte| *byte == b'\n').count());
    for line in lines.filter(|line| !line.is_empty()) {
        let mut record = parse_line_with(line, index_of, stringency)?;
        // These are held until the caller is done with all of them. One line parsed on its own need
        // not pay for this.
        record.tags.compact();
        records.push(record);
    }
    Ok((header, records))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cigar::{Cigar, CigarElement, Op};
    use crate::header::SequenceRecord;

    fn sample() -> (SamHeader, Vec<BamRecord>) {
        let mut h = SamHeader::new();
        h.set_sort_order("coordinate");
        h.sequences.push(SequenceRecord::new("chr1", 250_000_000));
        h.sequences.push(SequenceRecord::new("chr2", 200_000_000));
        let records = (0..20)
            .map(|i| BamRecord {
                read_name: format!("read{i}"),
                flags: 0,
                reference_index: i % 2,
                alignment_start: 100 + i * 7,
                mapping_quality: 60,
                cigar: Cigar::new(vec![CigarElement {
                    length: 4,
                    op: Op::M,
                }]),
                mate_reference_index: -1,
                mate_alignment_start: 0,
                inferred_insert_size: 0,
                read_bases: b"ACGT".to_vec(),
                base_qualities: vec![30, 31, 32, 33],
                tags: Default::default(),
            })
            .collect();
        (h, records)
    }

    /// No separator between the header and the first record: the header's own trailing newline
    /// is the only one.
    #[test]
    fn the_first_record_follows_the_header_with_no_blank_line() {
        let (h, records) = sample();
        let text = write_sam(&h, &records).unwrap();
        let header_len = h.encode().len();
        assert_eq!(&text[..header_len], h.encode());
        assert!(
            text[header_len..].starts_with("read0\t"),
            "a blank line crept in: {:?}",
            &text[header_len..header_len + 20]
        );
    }

    #[test]
    fn a_file_round_trips() {
        let (h, records) = sample();
        let text = write_sam(&h, &records).unwrap();
        let (back_h, back_r) = read_sam(&text).unwrap();
        assert_eq!(back_h, h);
        assert_eq!(back_r, records);
        assert_eq!(write_sam(&back_h, &back_r).unwrap(), text);
    }

    #[test]
    fn a_records_reference_name_comes_from_the_header() {
        let (h, records) = sample();
        let text = write_sam(&h, &records).unwrap();
        let first_record = text.lines().find(|l| !l.starts_with('@')).unwrap();
        assert_eq!(first_record.split('\t').nth(2), Some("chr1"));
        let second = text.lines().filter(|l| !l.starts_with('@')).nth(1).unwrap();
        assert_eq!(second.split('\t').nth(2), Some("chr2"));
    }

    /// The header/body split is positional, not a filter. An `@`-prefixed read name after the
    /// records begin is a record.
    #[test]
    fn an_at_prefixed_line_after_the_records_is_a_record() {
        let text = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n\
                    r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t?@AB\n\
                    @weird\t0\tchr1\t2\t60\t4M\t*\t0\t0\tACGT\t?@AB\n";
        let (h, records) = read_sam(text).unwrap();
        assert_eq!(h.sequences.len(), 1);
        assert_eq!(records.len(), 2, "the @-prefixed line is a record");
        assert_eq!(records[1].read_name, "@weird");
    }

    #[test]
    fn a_header_only_file_has_no_records() {
        let (h, _) = sample();
        let text = write_sam(&h, &[]).unwrap();
        let (back, records) = read_sam(&text).unwrap();
        assert_eq!(back, h);
        assert!(records.is_empty());
    }
}
