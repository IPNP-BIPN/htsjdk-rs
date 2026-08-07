//! The CRAM index: six numbers a line, and the arithmetic that queries them.
//!
//! Ported from `htsjdk.samtools.cram.CRAIEntry` and `CRAIIndex` at htsjdk 4.2.0.
//!
//! A `.crai` is not a structure, it is a sorted text file: one line per slice, six tab-separated
//! integers, gzipped. So what is worth pinning is not a layout but four decisions.
//!
//! # Unmapped-unplaced sorts last, and its alignment start is not consulted
//!
//! Everything else sorts by reference, then start, then container offset, then slice offset. An
//! unmapped entry skips the start entirely, so two unmapped entries at starts 900 and 100 sort by
//! their container offsets and come out in that order.
//!
//! # An unmapped entry never intersects, not even with itself
//!
//! Stated in the reference as a special case rather than falling out of the arithmetic.
//!
//! # The overlap test is a midpoint comparison
//!
//! `|a0 + b0 - a1 - b1| < span0 + span1`, which is not the same expression as
//! `a0 < b1 && a1 < b0` and does not agree with it on a zero span: two identical entries of span
//! zero do **not** intersect, and neither does a zero-span entry with one that contains it.
//!
//! # A query with a start or a span below one matches the whole sequence
//!
//! So 0 and -1 are not out-of-range values here, they are a wildcard, and the overlap test is never
//! reached for such a query.

use std::cmp::Ordering;

/// `ReferenceContext.UNMAPPED_UNPLACED_ID`.
pub const UNMAPPED_UNPLACED_ID: i32 = -1;
/// `ReferenceContext.MULTIPLE_REFERENCE_ID`, which cannot be indexed directly.
pub const MULTIPLE_REFERENCE_ID: i32 = -2;
/// `CRAI_INDEX_COLUMNS`.
pub const COLUMNS: usize = 6;

/// What building or parsing an entry refuses. All three are `CRAMException`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CraiError {
    /// A multi-reference slice, which has to be indexed by the references it holds.
    MultiReference,
    /// A line with the wrong number of columns. The count is of what was found.
    WrongColumnCount { found: usize },
    /// A column that is not a number. The reference wraps the JDK's exception, so its own message
    /// is the wrapped one's `toString`.
    NotANumber { text: String },
}

impl CraiError {
    pub fn message(&self) -> String {
        match self {
            // The two spaces after the full stop are the reference's.
            CraiError::MultiReference => "Cannot directly index a multiref slice.  Index by its \
                                          constituent references instead."
                .to_string(),
            CraiError::WrongColumnCount { found } => {
                format!("Malformed CRAI index entry: expecting {COLUMNS} columns but got {found}")
            }
            CraiError::NotANumber { text } => {
                format!("java.lang.NumberFormatException: For input string: \"{text}\"")
            }
        }
    }

    pub fn java_exception(&self) -> &'static str {
        "CRAMException"
    }
}

/// One line of the index: where a slice is, and what it covers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CraiEntry {
    pub sequence_id: i32,
    pub alignment_start: i32,
    pub alignment_span: i32,
    pub container_start_byte_offset: i64,
    pub slice_byte_offset_from_compression_header_start: i32,
    pub slice_byte_size: i32,
}

impl CraiEntry {
    /// The constructor, which refuses a multi-reference slice and checks nothing else. Negative
    /// starts, spans and offsets all go through: the corpus carries a line of them.
    pub fn new(
        sequence_id: i32,
        alignment_start: i32,
        alignment_span: i32,
        container_start_byte_offset: i64,
        slice_byte_offset_from_compression_header_start: i32,
        slice_byte_size: i32,
    ) -> Result<Self, CraiError> {
        if sequence_id == MULTIPLE_REFERENCE_ID {
            return Err(CraiError::MultiReference);
        }
        Ok(Self {
            sequence_id,
            alignment_start,
            alignment_span,
            container_start_byte_offset,
            slice_byte_offset_from_compression_header_start,
            slice_byte_size,
        })
    }

    /// `new CRAIEntry(line)`.
    ///
    /// The column count is checked before anything is parsed, and the multi-reference check the
    /// other constructor makes is **not** made here: a line naming -2 parses.
    pub fn parse(line: &str) -> Result<Self, CraiError> {
        let columns: Vec<&str> = line.split('\t').collect();
        if columns.len() != COLUMNS {
            return Err(CraiError::WrongColumnCount {
                found: columns.len(),
            });
        }
        let number = |text: &str| -> Result<i64, CraiError> {
            text.parse::<i64>().map_err(|_| CraiError::NotANumber {
                text: text.to_string(),
            })
        };
        Ok(Self {
            sequence_id: number(columns[0])? as i32,
            alignment_start: number(columns[1])? as i32,
            alignment_span: number(columns[2])? as i32,
            container_start_byte_offset: number(columns[3])?,
            slice_byte_offset_from_compression_header_start: number(columns[4])? as i32,
            slice_byte_size: number(columns[5])? as i32,
        })
    }

    /// `serializeToString`: six numbers, tab separated.
    pub fn serialize(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            self.sequence_id,
            self.alignment_start,
            self.alignment_span,
            self.container_start_byte_offset,
            self.slice_byte_offset_from_compression_header_start,
            self.slice_byte_size
        )
    }

    /// `intersect`, the midpoint comparison.
    pub fn intersects(&self, other: &CraiEntry) -> bool {
        if self.sequence_id != other.sequence_id {
            return false;
        }
        // Stated as a special case in the reference, before any arithmetic.
        if self.sequence_id == UNMAPPED_UNPLACED_ID {
            return false;
        }

        let a0 = i64::from(self.alignment_start);
        let a1 = i64::from(other.alignment_start);
        let b0 = a0 + i64::from(self.alignment_span);
        let b1 = a1 + i64::from(other.alignment_span);

        (a0 + b0 - a1 - b1).abs() < i64::from(self.alignment_span) + i64::from(other.alignment_span)
    }
}

impl Ord for CraiEntry {
    /// `compareTo`: unmapped last, then start, then container offset, then slice offset. The start
    /// is skipped entirely for an unmapped entry.
    fn cmp(&self, other: &Self) -> Ordering {
        if self.sequence_id != other.sequence_id {
            if self.sequence_id == UNMAPPED_UNPLACED_ID {
                return Ordering::Greater;
            }
            if other.sequence_id == UNMAPPED_UNPLACED_ID {
                return Ordering::Less;
            }
            return self.sequence_id.cmp(&other.sequence_id);
        }

        if self.sequence_id != UNMAPPED_UNPLACED_ID && self.alignment_start != other.alignment_start
        {
            return self.alignment_start.cmp(&other.alignment_start);
        }

        if self.container_start_byte_offset != other.container_start_byte_offset {
            return self
                .container_start_byte_offset
                .cmp(&other.container_start_byte_offset);
        }

        self.slice_byte_offset_from_compression_header_start
            .cmp(&other.slice_byte_offset_from_compression_header_start)
    }
}

impl PartialOrd for CraiEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// `CRAIIndex.writeIndex`: the entries sorted, one per line, each followed by a newline.
///
/// The sort happens here rather than on the way in, so an index holds its entries in whatever order
/// they were added and only the file is ordered.
pub fn write_index(entries: &[CraiEntry]) -> Vec<u8> {
    let mut sorted = entries.to_vec();
    sorted.sort();
    let mut out = Vec::new();
    for entry in sorted {
        out.extend_from_slice(entry.serialize().as_bytes());
        out.push(b'\n');
    }
    out
}

/// `CRAIIndex.find`.
///
/// A start or a span below one matches the whole sequence, so the overlap test is never reached for
/// such a query. The result is sorted, whatever order the list was in.
pub fn find(entries: &[CraiEntry], sequence_id: i32, start: i32, span: i32) -> Vec<CraiEntry> {
    let match_entire_sequence = start < 1 || span < 1;
    // The reference builds a query entry with a dummy value of 1 in the three offset fields, which
    // the overlap test does not look at.
    let query = CraiEntry {
        sequence_id,
        alignment_start: start,
        alignment_span: span,
        container_start_byte_offset: 1,
        slice_byte_offset_from_compression_header_start: 1,
        slice_byte_size: 1,
    };

    let mut found: Vec<CraiEntry> = entries
        .iter()
        .filter(|entry| entry.sequence_id == sequence_id)
        .filter(|entry| match_entire_sequence || entry.intersects(&query))
        .copied()
        .collect();
    found.sort();
    found
}

/// `CRAIIndex.getLeftmost`: the first after sorting, or nothing at all.
pub fn leftmost(entries: &[CraiEntry]) -> Option<CraiEntry> {
    let mut sorted = entries.to_vec();
    sorted.sort();
    sorted.first().copied()
}
