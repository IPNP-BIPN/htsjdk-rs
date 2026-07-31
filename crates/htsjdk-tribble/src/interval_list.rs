//! `IntervalListCodec`, ported from `htsjdk.tribble.IntervalList.IntervalListCodec` (htsjdk 4.2.0).
//!
//! This is what GATK's `-L regions.interval_list` runs through. The same file has a second parser
//! in htsjdk, [`htsjdk_bam::interval::IntervalList::parse_body`], and the two do **not** agree:
//! the reader accepts what this codec refuses, so the answer to "is this interval list valid"
//! depends on which door it came through. GATK's `-L` uses this door.
//!
//! # The field count is exact, and it is counted after trailing empties are dropped
//!
//! ```java
//! final String[] fields = line.split("\t");
//! if (fields.length != 5) { throw new TribbleException(...); }
//! ```
//!
//! `String.split` with no limit drops trailing empty fields, so `"chr1\t1\t10\t+\t"` has four
//! fields, not five, and is refused for the count rather than for the empty name. An empty field
//! in the *middle* survives, so `"chr1\t1\t10\t\tname"` reaches the strand check with `""`.
//!
//! # `.` is a legal strand everywhere else and an error here
//!
//! ```java
//! Strand strand = Strand.decode(fields[STRAND_POS]);
//! if (strand == Strand.NONE) throw new IllegalArgumentException("Invalid strand field: " + ...);
//! ```
//!
//! `Strand.decode` answers `NONE` for anything that is not exactly one `+` or `-`, and `NONE` is
//! then rejected, so the codec has no way to express an unstranded interval. `BEDCodec` in the
//! same package takes `.` and keeps it.
//!
//! # An unknown contig is a dropped line, an over-long one is a dead file
//!
//! ```java
//! if (sequence == null) { log.warn(...); return null; }
//! if (sequenceLength > 0 && sequenceLength < end) { throw new IllegalArgumentException(...); }
//! ```
//!
//! Two ways to disagree with the dictionary, two different outcomes: a contig it does not hold
//! costs one line, silently as far as the return value is concerned, while an interval past the
//! end of a contig it does hold fails the whole read. And the length check is switched off by a
//! declared length of 0, which is a legal `@SQ` and turns the contig into an unbounded one.
//!
//! # `start == end + 1` is an interval, `start == end + 2` is a refusal
//!
//! ```java
//! if (start > end + 1) throw new IllegalArgumentException(... "I'm afraid I cannot let you do that.");
//! ```
//!
//! The empty interval is deliberate: `IntervalList` writes zero-length intervals that way. One
//! past it is refused with a message quoting *2001*, which is part of the observable output and
//! so part of the port.

use htsjdk_bam::header::SamHeader;

/// `Strand`, as the codec uses it: only the two directed values survive the check.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Positive,
    Negative,
}

/// One decoded record, which is `htsjdk.samtools.util.Interval`.
///
/// The name is kept exactly as the column read: `.` is the string `.`, not an absent name. The
/// *reader* in `htsjdk-bam` turns `.` into `None`, and that difference between the two parsers is
/// the reason this type does not reuse it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IntervalRecord {
    pub contig: String,
    pub start: i32,
    pub end: i32,
    pub strand: Strand,
    pub name: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IntervalListError {
    /// `decode` with no dictionary, which is what a codec that never read a header holds.
    NoDictionary,
    /// A field count other than five, counted after `String.split` dropped trailing empties.
    FieldCount { count: usize, line: String },
    /// `Integer.parseInt` on the start or end column.
    NumberFormat { input: String },
    /// A start below 1.
    CoordinateBelowOne { start: i32 },
    /// A start more than one past the end.
    StartPastEnd { start: i32, end: i32 },
    /// A strand column that is not exactly `+` or `-`.
    InvalidStrand { field: String },
    /// An end past a contig whose declared length is greater than zero.
    PastContigEnd { end: i32, length: i32 },
}

impl IntervalListError {
    pub fn class(&self) -> &'static str {
        match self {
            IntervalListError::NoDictionary | IntervalListError::FieldCount { .. } => {
                "htsjdk.tribble.TribbleException"
            }
            IntervalListError::NumberFormat { .. } => "java.lang.NumberFormatException",
            _ => "java.lang.IllegalArgumentException",
        }
    }

    pub fn message(&self) -> String {
        match self {
            IntervalListError::NoDictionary => {
                "IntervalList dictionary cannot be null when decoding a record".to_string()
            }
            IntervalListError::FieldCount { count, line } => {
                format!("Invalid interval record contains {count} fields: {line}")
            }
            IntervalListError::NumberFormat { input } => format!("For input string: \"{input}\""),
            IntervalListError::CoordinateBelowOne { start } => format!(
                "Coordinate less than 1: start value of {start} is less than 1 and thus illegal"
            ),
            IntervalListError::StartPastEnd { start, end } => format!(
                "Start value of {start} is greater than end + 1 for end of value: {end}. \
                 I'm afraid I cannot let you do that."
            ),
            IntervalListError::InvalidStrand { field } => format!("Invalid strand field: {field}"),
            IntervalListError::PastContigEnd { end, length } => format!(
                "interval with end: {end} extends beyond end of sequence with length: {length}"
            ),
        }
    }
}

/// `IntervalListCodec.canDecode`: the two extensions, matched case-sensitively on the whole path.
///
/// Unlike `BEDCodec.canDecode` this does not strip a block-compressed extension first, so
/// `.interval_list.bgz` is not a Feature file even though htsjdk reads bgzip elsewhere.
pub fn can_decode(path: &str) -> bool {
    path.ends_with(".interval_list") || path.ends_with(".interval_list.gz")
}

/// `String.trim().isEmpty()`. Java's `trim` cuts every char `<= ' '`, which is not Rust's
/// `str::trim` (Unicode whitespace, and no control chars above the space).
fn java_trim_is_empty(line: &str) -> bool {
    !line.chars().any(|c| c > ' ')
}

/// `String.split("\t")` with no limit: **every** trailing empty field is dropped, so `"\t"` splits
/// into nothing at all rather than into two empty fields. The one exception is the empty input,
/// which Java answers with a single empty field.
fn split_dropping_trailing_empties(line: &str) -> Vec<&str> {
    if line.is_empty() {
        return vec![""];
    }
    let mut fields: Vec<&str> = line.split('\t').collect();
    while fields.last() == Some(&"") {
        fields.pop();
    }
    fields
}

fn parse_int(text: &str) -> Result<i32, IntervalListError> {
    text.parse::<i32>()
        .map_err(|_| IntervalListError::NumberFormat {
            input: text.to_string(),
        })
}

/// `IntervalListCodec.decode`, with the dictionary the codec would have taken from the header it
/// read.
///
/// `Ok(None)` is a line the reference drops: a header line, a blank one, or an interval on a
/// contig the dictionary does not hold. `dictionary` is `None` for a codec that never read a
/// header, which refuses every record.
pub fn decode(
    line: &str,
    dictionary: Option<&SamHeader>,
) -> Result<Option<IntervalRecord>, IntervalListError> {
    if line.starts_with('@') {
        return Ok(None);
    }
    if java_trim_is_empty(line) {
        return Ok(None);
    }
    let Some(header) = dictionary else {
        return Err(IntervalListError::NoDictionary);
    };

    let fields = split_dropping_trailing_empties(line);
    if fields.len() != 5 {
        return Err(IntervalListError::FieldCount {
            count: fields.len(),
            line: line.to_string(),
        });
    }

    // `lastSeq` in the reference is an interning cache: it replaces an equal string with the
    // previous one so the records share it. Nothing observable depends on it, so it is not here.
    let contig = fields[0];
    let start = parse_int(fields[1])?;
    let end = parse_int(fields[2])?;
    if start < 1 {
        return Err(IntervalListError::CoordinateBelowOne { start });
    }
    // `end + 1` in Java wraps at Integer.MAX_VALUE rather than trapping, so the port wraps too:
    // an end of 2147483647 makes the comparison `start > -2147483648`, which every start passes.
    if start > end.wrapping_add(1) {
        return Err(IntervalListError::StartPastEnd { start, end });
    }

    let strand = match fields[3] {
        "+" => Strand::Positive,
        "-" => Strand::Negative,
        field => {
            return Err(IntervalListError::InvalidStrand {
                field: field.to_string(),
            })
        }
    };

    let Some(sequence) = header.sequences.iter().find(|s| s.name == contig) else {
        // The reference logs a warning and returns null. The line is dropped, and the file loads.
        return Ok(None);
    };
    if sequence.length > 0 && sequence.length < end {
        return Err(IntervalListError::PastContigEnd {
            end,
            length: sequence.length,
        });
    }

    Ok(Some(IntervalRecord {
        contig: contig.to_string(),
        start,
        end,
        strand,
        name: fields[4].to_string(),
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use htsjdk_bam::header::SequenceRecord;

    fn header() -> SamHeader {
        let mut header = SamHeader::default();
        header.sequences.push(SequenceRecord::new("chr1", 200));
        header
    }

    #[test]
    fn a_trailing_tab_costs_a_field_rather_than_emptying_the_name() {
        let error = decode("chr1\t1\t10\t+\t", Some(&header())).unwrap_err();
        assert_eq!(
            error,
            IntervalListError::FieldCount {
                count: 4,
                line: "chr1\t1\t10\t+\t".into()
            }
        );
    }

    #[test]
    fn an_empty_field_in_the_middle_survives_the_count_and_fails_the_strand() {
        let error = decode("chr1\t1\t10\t\tname", Some(&header())).unwrap_err();
        assert_eq!(
            error,
            IntervalListError::InvalidStrand {
                field: String::new()
            }
        );
    }

    #[test]
    fn the_empty_interval_is_legal_and_the_next_one_is_not() {
        assert!(decode("chr1\t11\t10\t+\tn", Some(&header()))
            .unwrap()
            .is_some());
        assert_eq!(
            decode("chr1\t12\t10\t+\tn", Some(&header())).unwrap_err(),
            IntervalListError::StartPastEnd { start: 12, end: 10 }
        );
    }

    #[test]
    fn an_unknown_contig_is_dropped_and_an_over_long_one_is_refused() {
        assert_eq!(decode("chrX\t1\t10\t+\tn", Some(&header())).unwrap(), None);
        assert_eq!(
            decode("chr1\t1\t201\t+\tn", Some(&header())).unwrap_err(),
            IntervalListError::PastContigEnd {
                end: 201,
                length: 200
            }
        );
    }

    #[test]
    fn a_blank_line_is_dropped_before_the_dictionary_is_needed() {
        assert_eq!(decode("   ", None).unwrap(), None);
        assert_eq!(decode("@HD\tVN:1.6", None).unwrap(), None);
        assert_eq!(
            decode("chr1\t1\t10\t+\tn", None).unwrap_err(),
            IntervalListError::NoDictionary
        );
    }

    #[test]
    fn can_decode_matches_on_the_extension_and_does_not_strip_bgz() {
        assert!(can_decode("a.interval_list"));
        assert!(can_decode("a.interval_list.gz"));
        assert!(!can_decode("a.interval_list.bgz"));
        assert!(!can_decode("a.INTERVAL_LIST"));
    }
}
