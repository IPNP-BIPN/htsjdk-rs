//! `htsjdk.samtools.QueryInterval` and the chunk arithmetic a query closes on.
//!
//! A BAI query is three steps: turn the caller's intervals into a canonical set, turn each into the
//! index bins that can hold it, and turn those bins' chunks into a list of ranges to read. This
//! module is the first and the third; `bin` is the second, and reading the index itself is
//! `index`.
//!
//! `gatk-rs` names `QueryInterval` and `Chunk` in its read source and reaches neither from here,
//! because neither was here. `picard-rs`'s `MergeSamFiles` port wrote its own overlap loop for the
//! same reason and was five records out on its first run, because `queryOverlapping` returns
//! placed-but-unmapped reads and a hand-written filter did not.
//!
//! Two conventions decide most of the behaviour and neither is obvious:
//!
//! * **`end == 0` means "to the end of the reference"**, so a comparison has to special-case it
//!   rather than treat it as an empty interval, and merging two intervals where either is open
//!   produces an open one.
//! * **Adjacency is not overlap.** `optimizeIntervals` merges an interval that *abuts* the next
//!   one, and `optimizeChunkList` coalesces chunks in adjacent BGZF blocks, both to keep the read
//!   sequential; a port that only merged overlaps would answer the same records through more
//!   seeks, and would disagree with the reference's chunk list byte for byte.

use htsjdk_bgzf::vfp::{block_address, block_offset};

use crate::index::Chunk;

/// `htsjdk.samtools.QueryInterval`: 1-based, inclusive, with `end == 0` meaning open-ended.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct QueryInterval {
    pub reference_index: i32,
    pub start: i32,
    pub end: i32,
}

impl QueryInterval {
    /// # Panics
    ///
    /// On a negative reference index, which is htsjdk's `IllegalArgumentException`.
    pub fn new(reference_index: i32, start: i32, end: i32) -> Self {
        assert!(
            reference_index >= 0,
            "Invalid reference index {reference_index}"
        );
        QueryInterval {
            reference_index,
            start,
            end,
        }
    }

    /// `compareTo`, which answers a **difference** rather than a sign, except where the open-ended
    /// end forces it to answer 1 or -1: an open interval sorts after every closed one on the same
    /// start.
    pub fn compare_to(&self, other: &QueryInterval) -> i32 {
        let comp = self.reference_index - other.reference_index;
        if comp != 0 {
            return comp;
        }
        let comp = self.start - other.start;
        if comp != 0 {
            return comp;
        }
        if self.end == other.end {
            0
        } else if self.end == 0 {
            1
        } else if other.end == 0 {
            -1
        } else {
            self.end - other.end
        }
    }

    /// `endsAtStartOf`: this interval stops exactly one base before the other starts.
    pub fn ends_at_start_of(&self, other: &QueryInterval) -> bool {
        self.reference_index == other.reference_index && self.end + 1 == other.start
    }

    /// `overlaps`, with `end == 0` read as `Integer.MAX_VALUE`.
    ///
    /// The test is `CoordMath.overlaps`, and it is not `start <= otherEnd && otherStart <= end`.
    /// The two agree on every interval whose start is at or before its end, and disagree on the
    /// ones where it is not: `CoordMath.overlaps(500, 400, 500, 400)` is **true**, because its
    /// third clause asks whether the second interval is *enclosed* by the first, and an inverted
    /// interval encloses itself. htsjdk does not refuse such an interval, so the port reproduces
    /// the answer rather than the intent.
    pub fn overlaps(&self, other: &QueryInterval) -> bool {
        if self.reference_index != other.reference_index {
            return false;
        }
        let this_end = if self.end == 0 { i32::MAX } else { self.end };
        let other_end = if other.end == 0 { i32::MAX } else { other.end };
        coord_math_overlaps(self.start, this_end, other.start, other_end)
    }
}

/// `CoordMath.overlaps(start, end, start2, end2)`, transcribed rather than simplified.
///
/// Its third clause, `encloses(start2, end2, start, end)`, is what makes an inverted interval
/// overlap itself, and dropping it would be the obvious simplification.
fn coord_math_overlaps(start: i32, end: i32, start2: i32, end2: i32) -> bool {
    (start2 >= start && start2 <= end)
        || (end2 >= start && end2 <= end)
        || encloses(start2, end2, start, end)
}

/// `CoordMath.encloses(outerStart, outerEnd, innerStart, innerEnd)`.
fn encloses(outer_start: i32, outer_end: i32, inner_start: i32, inner_end: i32) -> bool {
    inner_start >= outer_start && inner_end <= outer_end
}

/// `QueryInterval.optimizeIntervals`: sort, then merge everything that overlaps **or abuts**.
///
/// Merging an abutting pair is what makes the result a canonical set rather than a sorted one, and
/// it is why `assertIntervalsOptimized` rejects an abutting pair as "not optimized" rather than as
/// merely redundant.
pub fn optimize_intervals(intervals: &[QueryInterval]) -> Vec<QueryInterval> {
    if intervals.is_empty() {
        return Vec::new();
    }
    let mut sorted = intervals.to_vec();
    sorted.sort_by(|a, b| a.compare_to(b).cmp(&0));

    let mut unique = Vec::new();
    let mut previous = sorted[0];
    for next in sorted.into_iter().skip(1) {
        if previous.ends_at_start_of(&next) || previous.overlaps(&next) {
            // Either end being open makes the merged interval open, whatever the other says.
            let new_end = if previous.end == 0 || next.end == 0 {
                0
            } else {
                previous.end.max(next.end)
            };
            previous = QueryInterval::new(previous.reference_index, previous.start, new_end);
        } else {
            unique.push(previous);
            previous = next;
        }
    }
    unique.push(previous);
    unique
}

/// Why `assertIntervalsOptimized` would refuse a list, or `None` when it would accept it.
///
/// The message is the reference's, because a caller comparing refusals compares text.
pub fn intervals_optimized_error(intervals: &[QueryInterval]) -> Option<String> {
    for pair in intervals.windows(2) {
        let (previous, this) = (&pair[0], &pair[1]);
        if previous.compare_to(this) >= 0 {
            return Some(format!(
                "List of intervals is not sorted: {} >= {}",
                display(previous),
                display(this)
            ));
        }
        if previous.overlaps(this) {
            return Some(format!(
                "List of intervals is not optimized: {} intersects {}",
                display(previous),
                display(this)
            ));
        }
        if previous.ends_at_start_of(this) {
            return Some(format!(
                "List of intervals is not optimized: {} abuts {}",
                display(previous),
                display(this)
            ));
        }
    }
    None
}

/// `QueryInterval.toString`: `referenceIndex:start-end`.
pub fn display(interval: &QueryInterval) -> String {
    format!(
        "{}:{}-{}",
        interval.reference_index, interval.start, interval.end
    )
}

/// `Chunk.compareTo`: by start, then by end, as a **sign** (`Long.signum`).
pub fn compare_chunks(a: &Chunk, b: &Chunk) -> std::cmp::Ordering {
    (a.start, a.end).cmp(&(b.start, b.end))
}

/// `Chunk.overlaps`, which compares BGZF block addresses and offsets rather than the raw pointers.
pub fn chunks_overlap(a: &Chunk, b: &Chunk) -> bool {
    let comparison = compare_chunks(a, b);
    if comparison == std::cmp::Ordering::Equal {
        return true;
    }
    let (left, right) = if comparison == std::cmp::Ordering::Less {
        (a, b)
    } else {
        (b, a)
    };
    let left_block = block_address(left.end);
    let right_block = block_address(right.start);
    if left_block > right_block {
        true
    } else if left_block == right_block {
        block_offset(left.end) > block_offset(right.start)
    } else {
        false
    }
}

/// `Chunk.isAdjacentTo`: one chunk's end is exactly the other's start, in both directions.
pub fn chunks_are_adjacent(a: &Chunk, b: &Chunk) -> bool {
    (block_address(a.end) == block_address(b.start) && block_offset(a.end) == block_offset(b.start))
        || (block_address(a.start) == block_address(b.end)
            && block_offset(a.start) == block_offset(b.end))
}

/// `Chunk.optimizeChunkList(chunks, minimumOffset)`.
///
/// Three things happen here, and the order matters: chunks ending at or before `minimum_offset` are
/// dropped (the linear-index optimization, which is what makes a query skip the file's beginning),
/// the rest are sorted, and neighbours that overlap **or sit in adjacent blocks** are coalesced so
/// the read is sequential. A chunk that survives the drop keeps its own start even when that start
/// is below `minimum_offset`.
pub fn optimize_chunk_list(chunks: &[Chunk], minimum_offset: u64) -> Vec<Chunk> {
    let mut sorted = chunks.to_vec();
    sorted.sort_by(compare_chunks);

    let mut result: Vec<Chunk> = Vec::new();
    for chunk in sorted {
        if chunk.end <= minimum_offset {
            continue;
        }
        match result.last_mut() {
            None => result.push(chunk),
            Some(last) => {
                if !chunks_overlap(last, &chunk) && !chunks_are_adjacent(last, &chunk) {
                    result.push(chunk);
                } else if chunk.end > last.end {
                    last.end = chunk.end;
                }
            }
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    fn vfp(block: u64, offset: u16) -> u64 {
        (block << 16) | offset as u64
    }

    #[test]
    fn an_open_ended_interval_sorts_after_a_closed_one() {
        let open = QueryInterval::new(0, 100, 0);
        let closed = QueryInterval::new(0, 100, 200);
        assert_eq!(open.compare_to(&closed), 1);
        assert_eq!(closed.compare_to(&open), -1);
    }

    #[test]
    fn abutting_intervals_are_merged_not_kept() {
        let intervals = [
            QueryInterval::new(0, 100, 199),
            QueryInterval::new(0, 200, 300),
        ];
        let optimized = optimize_intervals(&intervals);
        assert_eq!(optimized, vec![QueryInterval::new(0, 100, 300)]);
        // And the assertion rejects the unmerged pair with its own words.
        let message = intervals_optimized_error(&intervals).expect("abutting is not optimized");
        assert!(message.contains("abuts"), "{message}");
    }

    #[test]
    fn merging_with_an_open_interval_gives_an_open_one() {
        let intervals = [
            QueryInterval::new(0, 100, 500),
            QueryInterval::new(0, 200, 0),
        ];
        assert_eq!(
            optimize_intervals(&intervals),
            vec![QueryInterval::new(0, 100, 0)]
        );
    }

    #[test]
    fn an_inverted_interval_overlaps_itself() {
        // CoordMath.overlaps closes on `encloses`, and (500, 400) encloses (500, 400).
        let inverted = QueryInterval::new(2, 500, 400);
        assert!(inverted.overlaps(&inverted));
    }

    #[test]
    fn intervals_on_different_references_never_merge() {
        let intervals = [
            QueryInterval::new(0, 100, 200),
            QueryInterval::new(1, 100, 200),
        ];
        assert_eq!(optimize_intervals(&intervals).len(), 2);
    }

    #[test]
    fn chunks_in_adjacent_blocks_are_coalesced() {
        // The first ends exactly where the second begins, in block terms.
        let first = Chunk {
            start: vfp(10, 0),
            end: vfp(20, 40),
        };
        let second = Chunk {
            start: vfp(20, 40),
            end: vfp(30, 0),
        };
        assert!(chunks_are_adjacent(&first, &second));
        let optimized = optimize_chunk_list(&[first, second], 0);
        assert_eq!(optimized.len(), 1);
        assert_eq!(optimized[0].end, vfp(30, 0));
    }

    #[test]
    fn the_linear_index_drops_chunks_that_end_before_the_minimum() {
        let early = Chunk {
            start: vfp(1, 0),
            end: vfp(2, 0),
        };
        let late = Chunk {
            start: vfp(50, 0),
            end: vfp(60, 0),
        };
        let optimized = optimize_chunk_list(&[early, late], vfp(10, 0));
        assert_eq!(optimized, vec![late], "the early chunk is skipped entirely");
    }
}
