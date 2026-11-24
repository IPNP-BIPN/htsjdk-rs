//! `htsjdk.samtools.BAMIteratorFilter` and its multiple-intervals implementation.
//!
//! The chunks a query reads are approximate: a BGZF block holds records that do not overlap the
//! region, so every record that comes back is asked again, one at a time, whether it belongs. This
//! is that question, and it is the third and last piece of a BAI query.
//!
//! It is also where `picard-rs`'s hand-written overlap loop went wrong. The rule that catches a
//! reader out:
//!
//! ```java
//! if (record.getReadUnmappedFlag() && record.getAlignmentStart() != SAMRecord.NO_ALIGNMENT_START) {
//!     alignmentEnd = record.getAlignmentStart();   // Unmapped read with coordinate of mate.
//! }
//! ```
//!
//! An unmapped read placed at its mate's coordinate is **in** the query: it spans the single base
//! it starts on. A filter that dropped every unmapped read would return five fewer records on this
//! programme's own corpus, which is exactly what happened.
//!
//! The state machine is the other half. The filter walks each reference's **sorted** intervals
//! alongside that reference's records, advancing its cursor when an interval is behind the record
//! and stopping the whole iteration when every reference's cursor has run off the end. So `STOP_ITERATION` is not "this record fails", it
//! is "no later record can pass either", and a port that answered a boolean per record would read
//! the whole file.

use crate::cigar::Cigar;
use crate::query::QueryInterval;
use crate::record::BamRecord;

const READ_UNMAPPED: u16 = 0x4;
/// `SAMRecord.NO_ALIGNMENT_START`.
const NO_ALIGNMENT_START: i32 = 0;

/// `BAMIteratorFilter.IntervalComparison`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntervalComparison {
    Before,
    After,
    Overlapping,
    Contained,
}

/// `BAMIteratorFilter.FilteringIteratorState`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FilteringIteratorState {
    MatchesFilter,
    StopIteration,
    ContinueIteration,
}

/// `SAMRecord.getAlignmentEnd()` for the purpose of this comparison.
///
/// A mapped read ends where its cigar says. An unmapped read that carries a coordinate ends where
/// it starts, and an unplaced one has no end at all, which the caller reaches only through the
/// reference-index comparison above it.
fn alignment_end_for_filter(record: &BamRecord) -> i32 {
    if record.flags & READ_UNMAPPED != 0 && record.alignment_start != NO_ALIGNMENT_START {
        return record.alignment_start;
    }
    record.alignment_start + reference_span(&record.cigar) - 1
}

fn reference_span(cigar: &Cigar) -> i32 {
    cigar
        .elements
        .iter()
        .filter(|element| element.op.consumes_reference_bases())
        .map(|element| element.length as i32)
        .sum()
}

/// `BAMQueryMultipleIntervalsIteratorFilter.compareIntervalToRecord`.
pub fn compare_interval_to_record(
    interval: &QueryInterval,
    record: &BamRecord,
) -> IntervalComparison {
    let interval_end = if interval.end <= 0 {
        i32::MAX
    } else {
        interval.end
    };
    let alignment_end = alignment_end_for_filter(record);

    if interval.reference_index < record.reference_index {
        IntervalComparison::Before
    } else if interval.reference_index > record.reference_index {
        IntervalComparison::After
    } else if interval_end < record.alignment_start {
        IntervalComparison::Before
    } else if alignment_end < interval.start {
        IntervalComparison::After
    } else if record.alignment_start >= interval.start && alignment_end <= interval_end {
        // CoordMath.encloses(interval.start, intervalEnd, alignmentStart, alignmentEnd)
        IntervalComparison::Contained
    } else {
        IntervalComparison::Overlapping
    }
}

/// `BAMQueryMultipleIntervalsIteratorFilter`: stateful, because each reference's interval cursor
/// only moves forward.
///
/// Diverges from htsjdk 5.0.0 on purpose, as htsjdk 6.0.0 does (samtools/htsjdk#1877): 5.0.0
/// walked one cursor over all the intervals, so in a file whose references are not in dictionary
/// order (`samtools cat` of per-contig files) it stopped early and dropped every reference that
/// came later in the file than one sorting after it. Records are now compared with the intervals
/// of their own reference only, as htslib's multi-region iterator does.
pub struct MultipleIntervalsFilter {
    intervals: Vec<QueryInterval>,
    contained: bool,
    /// The intervals' distinct references, ascending.
    references: Vec<i32>,
    /// For each of `references`, the index of its first interval not yet behind its records, or
    /// `None` once they all are.
    next_interval: Vec<Option<usize>>,
    references_with_intervals_ahead: usize,
}

impl MultipleIntervalsFilter {
    /// `intervals` sorted by reference and then start, as [`crate::query::optimize_intervals`]
    /// sorts them; `contained` is the caller's `queryContained` rather than `queryOverlapping`.
    pub fn new(intervals: Vec<QueryInterval>, contained: bool) -> Self {
        let first_intervals: Vec<usize> = (0..intervals.len())
            .filter(|&i| i == 0 || intervals[i].reference_index != intervals[i - 1].reference_index)
            .collect();
        MultipleIntervalsFilter {
            references: first_intervals
                .iter()
                .map(|&i| intervals[i].reference_index)
                .collect(),
            references_with_intervals_ahead: first_intervals.len(),
            next_interval: first_intervals.into_iter().map(Some).collect(),
            intervals,
            contained,
        }
    }

    /// `compareToFilter(record)`.
    pub fn compare_to_filter(&mut self, record: &BamRecord) -> FilteringIteratorState {
        let reference = record.reference_index;
        if let Ok(slot) = self.references.binary_search(&reference) {
            if let Some(first) = self.next_interval[slot] {
                for i in first..self.intervals.len() {
                    if self.intervals[i].reference_index != reference {
                        break;
                    }
                    let comparison = compare_interval_to_record(&self.intervals[i], record);
                    // An interval behind the record is behind every later record of its reference.
                    if comparison == IntervalComparison::Before {
                        continue;
                    }
                    self.next_interval[slot] = Some(i);
                    let matches = comparison == IntervalComparison::Contained
                        || (comparison == IntervalComparison::Overlapping && !self.contained);
                    return if matches {
                        FilteringIteratorState::MatchesFilter
                    } else {
                        FilteringIteratorState::ContinueIteration
                    };
                }
                // The record is past its reference's last interval.
                self.next_interval[slot] = None;
                self.references_with_intervals_ahead -= 1;
            }
        }
        // Nothing later can match once every reference's intervals are behind its records.
        if self.references_with_intervals_ahead == 0 {
            FilteringIteratorState::StopIteration
        } else {
            FilteringIteratorState::ContinueIteration
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cigar::{CigarElement, Op};

    fn mapped(reference: i32, start: i32, length: u32) -> BamRecord {
        BamRecord {
            reference_index: reference,
            alignment_start: start,
            cigar: Cigar::new(vec![CigarElement { length, op: Op::M }]),
            ..Default::default()
        }
    }

    fn placed_unmapped(reference: i32, start: i32) -> BamRecord {
        BamRecord {
            reference_index: reference,
            alignment_start: start,
            flags: READ_UNMAPPED,
            ..Default::default()
        }
    }

    #[test]
    fn a_placed_unmapped_read_spans_one_base_and_is_in_the_query() {
        let interval = QueryInterval::new(0, 100, 200);
        let record = placed_unmapped(0, 150);
        assert_eq!(
            compare_interval_to_record(&interval, &record),
            IntervalComparison::Contained
        );
    }

    #[test]
    fn contained_and_overlapping_are_different_answers() {
        let interval = QueryInterval::new(0, 100, 200);
        assert_eq!(
            compare_interval_to_record(&interval, &mapped(0, 150, 10)),
            IntervalComparison::Contained
        );
        assert_eq!(
            compare_interval_to_record(&interval, &mapped(0, 195, 20)),
            IntervalComparison::Overlapping
        );
    }

    #[test]
    fn the_cursor_only_moves_forward_and_then_stops() {
        let mut filter = MultipleIntervalsFilter::new(
            vec![
                QueryInterval::new(0, 100, 200),
                QueryInterval::new(0, 400, 500),
            ],
            false,
        );
        // A record past both intervals leaves the cursor at the end, and the answer is STOP.
        assert_eq!(
            filter.compare_to_filter(&mapped(0, 900, 10)),
            FilteringIteratorState::StopIteration
        );
        // And it stays stopped: the cursor does not go back for an earlier record.
        assert_eq!(
            filter.compare_to_filter(&mapped(0, 150, 10)),
            FilteringIteratorState::StopIteration
        );
    }

    #[test]
    fn query_contained_rejects_what_query_overlapping_returns() {
        let overlapping = mapped(0, 195, 20);
        let mut lenient =
            MultipleIntervalsFilter::new(vec![QueryInterval::new(0, 100, 200)], false);
        let mut strict = MultipleIntervalsFilter::new(vec![QueryInterval::new(0, 100, 200)], true);
        assert_eq!(
            lenient.compare_to_filter(&overlapping),
            FilteringIteratorState::MatchesFilter
        );
        assert_eq!(
            strict.compare_to_filter(&overlapping),
            FilteringIteratorState::ContinueIteration
        );
    }

    #[test]
    fn references_out_of_dictionary_order_do_not_end_the_query() {
        // samtools/htsjdk#1877: the file holds reference 1's records before reference 0's, as
        // `samtools cat` of per-contig files writes them. 5.0.0 stopped at the third record.
        let mut filter = MultipleIntervalsFilter::new(
            vec![QueryInterval::new(0, 1, 100), QueryInterval::new(1, 1, 100)],
            false,
        );
        let states: Vec<_> = [
            mapped(1, 10, 5),
            mapped(1, 150, 5),
            mapped(0, 10, 5),
            mapped(0, 150, 5),
        ]
        .iter()
        .map(|record| filter.compare_to_filter(record))
        .collect();
        assert_eq!(
            states,
            vec![
                FilteringIteratorState::MatchesFilter,
                FilteringIteratorState::ContinueIteration,
                FilteringIteratorState::MatchesFilter,
                FilteringIteratorState::StopIteration,
            ]
        );
    }

    #[test]
    fn a_reference_without_intervals_is_skipped_until_every_interval_is_behind() {
        let mut filter = MultipleIntervalsFilter::new(vec![QueryInterval::new(0, 1, 100)], false);
        assert_eq!(
            filter.compare_to_filter(&mapped(1, 10, 5)),
            FilteringIteratorState::ContinueIteration
        );
        assert_eq!(
            filter.compare_to_filter(&mapped(0, 150, 5)),
            FilteringIteratorState::StopIteration
        );
    }
}
