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
//! The state machine is the other half. The filter walks a **sorted** interval list alongside the
//! records, advancing its cursor when an interval is behind the record and stopping the whole
//! iteration when the cursor runs off the end. So `STOP_ITERATION` is not "this record fails", it
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

/// `BAMQueryMultipleIntervalsIteratorFilter`: stateful, because its interval cursor only moves
/// forward.
pub struct MultipleIntervalsFilter {
    intervals: Vec<QueryInterval>,
    contained: bool,
    interval_index: usize,
}

impl MultipleIntervalsFilter {
    /// `contained` is the caller's `queryContained` rather than `queryOverlapping`.
    pub fn new(intervals: Vec<QueryInterval>, contained: bool) -> Self {
        MultipleIntervalsFilter {
            intervals,
            contained,
            interval_index: 0,
        }
    }

    /// `compareToFilter(record)`.
    pub fn compare_to_filter(&mut self, record: &BamRecord) -> FilteringIteratorState {
        while self.interval_index < self.intervals.len() {
            match compare_interval_to_record(&self.intervals[self.interval_index], record) {
                // The interval is behind the record: it can never match again, so drop it.
                IntervalComparison::Before => self.interval_index += 1,
                // The interval is ahead: this record misses, later ones may not.
                IntervalComparison::After => return FilteringIteratorState::ContinueIteration,
                IntervalComparison::Contained => return FilteringIteratorState::MatchesFilter,
                IntervalComparison::Overlapping => {
                    return if self.contained {
                        FilteringIteratorState::ContinueIteration
                    } else {
                        FilteringIteratorState::MatchesFilter
                    }
                }
            }
        }
        // Past the last interval: nothing later can match either.
        FilteringIteratorState::StopIteration
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
}
