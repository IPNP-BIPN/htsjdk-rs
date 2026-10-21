//! Interval overlap queries.
//!
//! Ported from `htsjdk.samtools.util.OverlapDetector` at htsjdk 4.2.0 — **behaviourally**, not
//! structurally. htsjdk backs the detector with a 1273-line augmented red-black `IntervalTree`.
//! Decision 0020 establishes, by measurement in the oracle, that the *order* of the overlap set
//! this returns does not reach the output of the consumer it was ported for
//! (`CollectRnaSeqMetrics` folds the set commutatively). So what has to be reproduced is the
//! **set of overlapping objects**, and the tree's traversal order is an unobservable
//! implementation detail. This port therefore uses a plainly correct structure and matches
//! htsjdk on the set, not on the walk.
//!
//! The one genuinely surprising rule, kept exactly, is the buffer sign. `addLhs` stores an
//! interval as `[start + lhsBuffer, end - lhsBuffer]` and `getOverlaps` queries
//! `[start + rhsBuffer, end - rhsBuffer]`. A **positive** buffer therefore *shrinks* the
//! interval and a negative buffer *grows* it, which is the opposite of what "buffer" suggests.
//! Picard's `create` passes `(0, 0)`, so the buffers are usually inert, but they are reproduced
//! because a caller that passes a non-zero buffer depends on the sign.
//!
//! An interval whose adjusted `start > end` is dropped on add and yields no overlaps on query,
//! matching the `if (start <= end)` guard on both sides.
//!
//! A query does not scan the contig. On first use each contig's intervals are sorted by start,
//! and a query visits only those starting within the longest stored interval's length of it,
//! then hands the hits back in insertion order, which is what the plain scan returned and what a
//! caller iterating the result sees.

use std::collections::HashMap;
use std::sync::OnceLock;

/// One stored interval and the object it carries.
struct Entry<T> {
    start: i32,
    end: i32,
    object: T,
}

/// One contig's intervals in insertion order, and the start-sorted view queries search.
struct Contig<T> {
    entries: Vec<Entry<T>>,
    index: OnceLock<StartIndex>,
}

/// `order[i]` is the entry with the `i`-th smallest start, `starts[i]` that start, and `span` the
/// largest `end - start` on the contig, which bounds how far left of a query a hit can start.
struct StartIndex {
    order: Vec<u32>,
    starts: Vec<i32>,
    span: i64,
}

impl<T> Contig<T> {
    fn index(&self) -> &StartIndex {
        self.index.get_or_init(|| {
            let mut order: Vec<u32> = (0..self.entries.len() as u32).collect();
            order.sort_by_key(|&i| self.entries[i as usize].start);
            let starts = order
                .iter()
                .map(|&i| self.entries[i as usize].start)
                .collect();
            let span = self
                .entries
                .iter()
                .map(|e| e.end as i64 - e.start as i64)
                .max()
                .unwrap_or(0);
            StartIndex {
                order,
                starts,
                span,
            }
        })
    }
}

/// `OverlapDetector<T>`.
///
/// Generic over the carried object. The object is returned by reference from
/// [`get_overlaps`](OverlapDetector::get_overlaps); a caller wanting set semantics over
/// non-`Copy` objects works with the references, exactly as htsjdk works with the `Set<T>` of
/// object identities.
pub struct OverlapDetector<T> {
    lhs_buffer: i32,
    rhs_buffer: i32,
    by_contig: HashMap<String, Contig<T>>,
}

impl<T> OverlapDetector<T> {
    pub fn new(lhs_buffer: i32, rhs_buffer: i32) -> Self {
        OverlapDetector {
            lhs_buffer,
            rhs_buffer,
            by_contig: HashMap::new(),
        }
    }

    /// `OverlapDetector.create(intervals)`: buffers of zero.
    pub fn create() -> Self {
        Self::new(0, 0)
    }

    /// `addLhs(object, interval)`.
    ///
    /// The adjusted interval is `[start + lhsBuffer, end - lhsBuffer]`, and an interval with no
    /// overlappable bases after the adjustment (`start > end`) is not stored.
    pub fn add(&mut self, contig: &str, start: i32, end: i32, object: T) {
        let start = start + self.lhs_buffer;
        let end = end - self.lhs_buffer;
        if start <= end {
            let contig = self
                .by_contig
                .entry(contig.to_string())
                .or_insert_with(|| Contig {
                    entries: Vec::new(),
                    index: OnceLock::new(),
                });
            contig.entries.push(Entry { start, end, object });
            contig.index = OnceLock::new();
        }
    }

    /// `getOverlaps(locatable)`: every stored object whose adjusted interval overlaps the query.
    ///
    /// Overlap is the closed-interval test `stored.start <= query.end && stored.end >= query.start`,
    /// the same one `IntervalTree.overlappers` applies. The returned order is insertion order
    /// here and traversal order in htsjdk; decision 0020 records that neither is observable.
    pub fn get_overlaps(&self, contig: &str, start: i32, end: i32) -> Vec<&T> {
        let start = start + self.rhs_buffer;
        let end = end - self.rhs_buffer;
        if start > end {
            return Vec::new();
        }
        let Some(contig) = self.by_contig.get(contig) else {
            return Vec::new();
        };
        // shortcut: one interval as long as the contig makes every query visit every start to
        // its left, which is the plain scan again; an interval tree would bound that too.
        let index = contig.index();
        let from = index
            .starts
            .partition_point(|&s| (s as i64) < start as i64 - index.span);
        let to = index.starts.partition_point(|&s| s <= end);
        let mut hits: Vec<u32> = index.order[from..to.max(from)]
            .iter()
            .copied()
            .filter(|&i| contig.entries[i as usize].end >= start)
            .collect();
        hits.sort_unstable();
        hits.into_iter()
            .map(|i| &contig.entries[i as usize].object)
            .collect()
    }

    /// `overlapsAny(locatable)`.
    pub fn overlaps_any(&self, contig: &str, start: i32, end: i32) -> bool {
        !self.get_overlaps(contig, start, end).is_empty()
    }

    /// `getAll()`: every stored object, in no particular order.
    pub fn get_all(&self) -> Vec<&T> {
        self.by_contig
            .values()
            .flat_map(|c| c.entries.iter().map(|e| &e.object))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The indexed query against the scan it replaced: same objects, same order, with a long
    /// interval in the mix, queries added between, and both buffer signs.
    #[test]
    fn the_index_answers_as_the_scan_did() {
        let mut seed = 7u64;
        let mut next = move |bound: u64| {
            seed = seed
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            ((seed >> 33) % bound) as i32
        };
        for (lhs, rhs) in [(0, 0), (3, -2), (-4, 5)] {
            let mut d = OverlapDetector::new(lhs, rhs);
            let mut stored: Vec<(&str, i32, i32, usize)> = Vec::new();
            for i in 0..2000 {
                let contig = ["a", "b"][next(2) as usize];
                let start = next(100_000);
                let length = if i == 500 { 60_000 } else { next(300) };
                d.add(contig, start, start + length, i);
                if start + lhs <= start + length - lhs {
                    stored.push((contig, start + lhs, start + length - lhs, i));
                }
                if i % 100 == 0 || i > 1900 {
                    let (c, s) = (["a", "b", "c"][next(3) as usize], next(100_000));
                    let e = s + next(2000);
                    let (qs, qe) = (s + rhs, e - rhs);
                    let expected: Vec<usize> = stored
                        .iter()
                        .filter(|x| qs <= qe && x.0 == c && x.1 <= qe && x.2 >= qs)
                        .map(|x| x.3)
                        .collect();
                    let got: Vec<usize> = d.get_overlaps(c, s, e).into_iter().copied().collect();
                    assert_eq!(got, expected, "query {c}:{s}-{e} after {i} adds");
                }
            }
        }
    }

    fn detector(intervals: &[(&str, i32, i32, &str)]) -> OverlapDetector<String> {
        let mut d = OverlapDetector::create();
        for (contig, start, end, name) in intervals {
            d.add(contig, *start, *end, name.to_string());
        }
        d
    }

    fn names(mut v: Vec<&String>) -> Vec<String> {
        v.sort();
        v.into_iter().cloned().collect()
    }

    #[test]
    fn overlap_is_the_closed_interval_test() {
        let d = detector(&[("chr1", 10, 20, "a")]);
        // Touching at a single base counts, both ends.
        assert_eq!(names(d.get_overlaps("chr1", 20, 30)), ["a"]);
        assert_eq!(names(d.get_overlaps("chr1", 1, 10)), ["a"]);
        // One base short of touching does not.
        assert!(d.get_overlaps("chr1", 21, 30).is_empty());
        assert!(d.get_overlaps("chr1", 1, 9).is_empty());
    }

    #[test]
    fn every_overlapping_object_is_returned() {
        let d = detector(&[
            ("chr1", 1000, 2000, "a"),
            ("chr1", 1500, 2500, "b"),
            ("chr1", 3000, 4000, "c"),
        ]);
        assert_eq!(names(d.get_overlaps("chr1", 1800, 1900)), ["a", "b"]);
        assert_eq!(names(d.get_overlaps("chr1", 3500, 3600)), ["c"]);
    }

    #[test]
    fn a_different_contig_never_overlaps() {
        let d = detector(&[("chr1", 10, 20, "a")]);
        assert!(d.get_overlaps("chr2", 10, 20).is_empty());
    }

    /// A positive buffer shrinks the interval. With `lhsBuffer = 5`, an interval stored as
    /// `[10, 20]` becomes `[15, 15]`, so a query at 12 no longer overlaps it.
    #[test]
    fn a_positive_buffer_shrinks_the_stored_interval() {
        let mut d = OverlapDetector::new(5, 0);
        d.add("chr1", 10, 20, "a".to_string());
        assert!(d.get_overlaps("chr1", 12, 12).is_empty(), "shrunk past 12");
        assert_eq!(names(d.get_overlaps("chr1", 15, 15)), ["a"]);
    }

    /// A negative buffer grows it, the counter-intuitive direction.
    #[test]
    fn a_negative_buffer_grows_the_query() {
        let mut d = OverlapDetector::new(0, -5);
        d.add("chr1", 100, 100, "a".to_string());
        // Query [108, 108] grows to [103, 113], which reaches the stored point at 100? No: 103 > 100.
        assert!(d.get_overlaps("chr1", 108, 108).is_empty());
        // Query [102,102] grows to [97,107], which contains 100.
        assert_eq!(names(d.get_overlaps("chr1", 102, 102)), ["a"]);
    }

    /// An interval with no overlappable bases after the buffer adjustment is dropped on add.
    #[test]
    fn an_interval_that_buffers_to_nothing_is_not_stored() {
        let mut d = OverlapDetector::new(10, 0);
        d.add("chr1", 10, 20, "a".to_string()); // becomes [20, 10], start > end, dropped
        assert!(d.get_all().is_empty());
    }

    /// A query that buffers to `start > end` returns nothing rather than scanning.
    #[test]
    fn a_query_that_buffers_to_nothing_returns_empty() {
        let mut d = OverlapDetector::new(0, 10);
        d.add("chr1", 1, 1000, "a".to_string());
        assert!(
            d.get_overlaps("chr1", 10, 20).is_empty(),
            "[20,10] after buffer"
        );
    }

    #[test]
    fn get_all_returns_every_object() {
        let d = detector(&[("chr1", 1, 2, "a"), ("chr2", 5, 6, "b")]);
        assert_eq!(names(d.get_all()), ["a", "b"]);
    }

    /// Two objects sharing the same interval are both returned, which is htsjdk's set-union at a
    /// coincident key.
    #[test]
    fn coincident_intervals_both_overlap() {
        let d = detector(&[("chr1", 100, 200, "a"), ("chr1", 100, 200, "b")]);
        assert_eq!(names(d.get_overlaps("chr1", 150, 150)), ["a", "b"]);
    }
}
