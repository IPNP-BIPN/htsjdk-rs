//! Intervals and interval lists.
//!
//! Ported from `htsjdk.samtools.util.Interval`, `IntervalCoordinateComparator`, `IntervalList`
//! and `IntervalListWriter` at htsjdk 4.2.0. Shared infrastructure: 16 of Picard's 44 metrics
//! tools reference it, and every GATK walker that takes `-L` does.
//!
//! ## Two orderings on one class, and a fifth rule for the tally
//!
//! `Interval` implements `Comparable` and sorts contigs by **`String.compareTo`**:
//!
//! ```java
//! int result = this.getContig().compareTo(that.getContig());
//! ```
//!
//! `IntervalCoordinateComparator`, which is what `IntervalList.sorted()` and therefore the file
//! writer actually use, sorts contigs by their **index in the sequence dictionary**:
//!
//! ```java
//! final int lhsIndex = this.header.getSequenceIndex(lhs.getContig());
//! ```
//!
//! The two disagree on any dictionary that is not in lexicographic order, which is every real
//! one: `chr10` precedes `chr2` under the natural ordering and follows it under the file's. So
//! the same list, sorted the two ways htsjdk offers, gives two different files. Decision 0018
//! counted four ordering rules in this library; this is a fifth, and the first where **one class
//! carries two of them**.
//!
//! The port exposes both under names that say which is which, because `impl Ord for Interval`
//! would silently pick the wrong one for file output.
//!
//! ## Two renderings, and only one of them is the format
//!
//! `Interval.toString()` produces `contig:start-end<TAB>strand<TAB>name`. `IntervalListWriter`
//! produces `contig<TAB>start<TAB>end<TAB>strand<TAB>name`. Only the second is the
//! `.interval_list` format; the first appears in exception messages. A port that implemented
//! `Display` and wrote it to the file would produce something that looks close enough to read.

use std::cmp::Ordering;

/// `Interval`. Coordinates are 1-based and closed, as htsjdk's are.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Interval {
    pub contig: String,
    pub start: i32,
    pub end: i32,
    pub negative_strand: bool,
    /// `null` in Java, which the writer renders as `.`.
    pub name: Option<String>,
}

impl Interval {
    pub fn new(contig: &str, start: i32, end: i32) -> Self {
        Interval {
            contig: contig.to_string(),
            start,
            end,
            negative_strand: false,
            name: None,
        }
    }

    pub fn with_strand_and_name(
        contig: &str,
        start: i32,
        end: i32,
        negative_strand: bool,
        name: Option<&str>,
    ) -> Self {
        Interval {
            contig: contig.to_string(),
            start,
            end,
            negative_strand,
            name: name.map(str::to_string),
        }
    }

    /// `Strand.encode()`.
    pub fn strand(&self) -> char {
        if self.negative_strand {
            '-'
        } else {
            '+'
        }
    }

    pub fn length(&self) -> i32 {
        self.end - self.start + 1
    }

    /// `intersects`: same contig and overlapping coordinates. Strand is not consulted.
    pub fn intersects(&self, other: &Interval) -> bool {
        self.contig == other.contig && self.start <= other.end && self.end >= other.start
    }

    /// `withinDistanceOf(other, distance)`, used by `uniqued` to combine abutting intervals.
    pub fn within_distance_of(&self, other: &Interval, distance: i32) -> bool {
        self.contig == other.contig
            && self.start - distance <= other.end
            && self.end + distance >= other.start
    }

    /// `Interval.compareTo`: contig by **string**, then start, then end, then strand, then name.
    ///
    /// This is the ordering a `TreeSet<Interval>` or `Collections.sort` gets, and it is **not**
    /// the ordering the file writer uses. See the module note.
    pub fn compare_natural(&self, other: &Interval) -> Ordering {
        self.contig
            .cmp(&other.contig)
            .then_with(|| self.start.cmp(&other.start))
            .then_with(|| self.end.cmp(&other.end))
            .then_with(|| compare_strand_then_name(self, other))
    }

    /// `Interval.toString()`, which is **not** the file format. Exception messages only.
    pub fn to_display_string(&self) -> String {
        format!(
            "{}:{}-{}\t{}\t{}",
            self.contig,
            self.start,
            self.end,
            self.strand(),
            self.name.as_deref().unwrap_or(".")
        )
    }

    /// One line as `IntervalListWriter.write` emits it, without the line terminator.
    pub fn to_file_line(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}",
            self.contig,
            self.start,
            self.end,
            self.strand(),
            self.name.as_deref().unwrap_or(".")
        )
    }
}

/// The tail shared by both comparators: positive strand first, then name with nulls first.
fn compare_strand_then_name(lhs: &Interval, rhs: &Interval) -> Ordering {
    if !lhs.negative_strand && rhs.negative_strand {
        return Ordering::Less;
    }
    if lhs.negative_strand && !rhs.negative_strand {
        return Ordering::Greater;
    }
    match (&lhs.name, &rhs.name) {
        (None, None) => Ordering::Equal,
        (None, Some(_)) => Ordering::Less,
        (Some(_), None) => Ordering::Greater,
        (Some(a), Some(b)) => a.cmp(b),
    }
}

/// `IntervalCoordinateComparator`: contig by **dictionary index**, then start, end, strand, name.
///
/// `getSequenceIndex` returns `-1` for a contig the dictionary does not have, so an unknown
/// contig sorts before every known one rather than raising. Reproduced: a list referring to a
/// contig outside its own header sorts to the front and writes without complaint.
pub fn compare_coordinate(dictionary: &[String], lhs: &Interval, rhs: &Interval) -> Ordering {
    let index = |c: &str| -> i32 {
        dictionary
            .iter()
            .position(|s| s == c)
            .map_or(-1, |i| i as i32)
    };
    index(&lhs.contig)
        .cmp(&index(&rhs.contig))
        .then_with(|| lhs.start.cmp(&rhs.start))
        .then_with(|| lhs.end.cmp(&rhs.end))
        .then_with(|| compare_strand_then_name(lhs, rhs))
}

/// `IntervalList`: a SAM header plus intervals.
#[derive(Debug, Clone, PartialEq)]
pub struct IntervalList {
    /// The sequence dictionary's contig names, in order. The index into this is what
    /// [`compare_coordinate`] sorts by.
    pub dictionary: Vec<String>,
    pub intervals: Vec<Interval>,
}

impl IntervalList {
    pub fn new(dictionary: Vec<String>) -> Self {
        IntervalList {
            dictionary,
            intervals: Vec::new(),
        }
    }

    /// `sorted()`: an independent list in coordinate order.
    pub fn sorted(&self) -> IntervalList {
        let mut out = self.clone();
        out.intervals
            .sort_by(|a, b| compare_coordinate(&self.dictionary, a, b));
        out
    }

    /// `uniqued()`: sorted, then overlapping and abutting intervals merged.
    ///
    /// Two defaults worth spelling out, because the no-argument Java form hides both.
    /// `IntervalList.uniqued()` is `uniqued(true)`, so the **default concatenates** the names of
    /// a merged run with `|` rather than keeping the first; `concatenate_names = false` is the
    /// path you have to ask for. And `combineAbuttingIntervals` is on, so `1-10` and `11-20`
    /// become `1-20` even though they do not overlap.
    pub fn uniqued(&self, concatenate_names: bool) -> IntervalList {
        let sorted = self.sorted();
        let mut out = IntervalList::new(self.dictionary.clone());
        let mut current: Option<Interval> = None;
        let mut names: Vec<String> = Vec::new();

        for next in &sorted.intervals {
            match &mut current {
                None => {
                    current = Some(next.clone());
                    names.clear();
                    if let Some(n) = &next.name {
                        names.push(n.clone());
                    }
                }
                Some(cur) if cur.intersects(next) || cur.within_distance_of(next, 1) => {
                    cur.end = cur.end.max(next.end);
                    if let Some(n) = &next.name {
                        names.push(n.clone());
                    }
                }
                Some(cur) => {
                    out.intervals.push(finish(cur, &names, concatenate_names));
                    current = Some(next.clone());
                    names.clear();
                    if let Some(n) = &next.name {
                        names.push(n.clone());
                    }
                }
            }
        }
        if let Some(cur) = &current {
            out.intervals.push(finish(cur, &names, concatenate_names));
        }
        out
    }

    /// `padded(before, after)`, clamped to the contig's bounds by the caller's dictionary.
    ///
    /// htsjdk clamps the start at 1 and the end at the sequence length. The length is not
    /// carried here, so only the start is clamped and the caller supplies the end bound.
    pub fn padded(
        &self,
        before: i32,
        after: i32,
        contig_length: impl Fn(&str) -> i32,
    ) -> IntervalList {
        let mut out = IntervalList::new(self.dictionary.clone());
        for i in &self.intervals {
            let start = (i.start - before).max(1);
            let end = (i.end + after).min(contig_length(&i.contig));
            if start <= end {
                out.intervals.push(Interval {
                    contig: i.contig.clone(),
                    start,
                    end,
                    negative_strand: i.negative_strand,
                    name: i.name.clone(),
                });
            }
        }
        out
    }

    /// The interval lines as `IntervalListWriter` emits them, in list order.
    ///
    /// The SAM header comes first in a real `.interval_list` and is written by the SAM header
    /// codec, which is already ported; this returns the body so the two can be joined without
    /// this module depending on the header writer.
    pub fn write_body(&self) -> String {
        let mut out = String::new();
        for i in &self.intervals {
            out.push_str(&i.to_file_line());
            // BufferedWriter.newLine() writes the *platform* line separator. On the pinned
            // linux/amd64 oracle that is "\n"; a Windows JVM would write "\r\n" and produce a
            // different file from the same list. The port writes "\n" and says so rather than
            // reproducing a host-dependent choice.
            out.push('\n');
        }
        out
    }

    /// Parses the interval lines of an `.interval_list`, skipping `@` header lines.
    pub fn parse_body(dictionary: Vec<String>, text: &str) -> Result<IntervalList, ParseError> {
        let mut out = IntervalList::new(dictionary);
        for line in text.lines() {
            if line.starts_with('@') || line.trim().is_empty() {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 5 {
                return Err(ParseError::TooFewFields(line.to_string()));
            }
            out.intervals.push(Interval {
                contig: f[0].to_string(),
                start: f[1]
                    .parse()
                    .map_err(|_| ParseError::BadNumber(f[1].into()))?,
                end: f[2]
                    .parse()
                    .map_err(|_| ParseError::BadNumber(f[2].into()))?,
                negative_strand: f[3] == "-",
                // The writer renders a null name as ".", so "." reads back as null. A list
                // whose interval is genuinely named "." cannot round-trip, and htsjdk has the
                // same hole.
                name: if f[4] == "." {
                    None
                } else {
                    Some(f[4].to_string())
                },
            });
        }
        Ok(out)
    }
}

fn finish(current: &Interval, names: &[String], concatenate_names: bool) -> Interval {
    Interval {
        contig: current.contig.clone(),
        start: current.start,
        end: current.end,
        negative_strand: current.negative_strand,
        name: if concatenate_names {
            if names.is_empty() {
                None
            } else {
                Some(names.join("|"))
            }
        } else {
            current.name.clone()
        },
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ParseError {
    TooFewFields(String),
    BadNumber(String),
}

#[cfg(test)]
mod tests {
    use super::*;

    fn iv(contig: &str, start: i32, end: i32) -> Interval {
        Interval::new(contig, start, end)
    }

    /// The finding: the two orderings htsjdk offers on `Interval` disagree on any dictionary
    /// that is not lexicographic, which is every real one.
    #[test]
    fn the_two_orderings_disagree_on_a_real_dictionary() {
        let dict: Vec<String> = ["chr1", "chr2", "chr10"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let a = iv("chr10", 1, 100);
        let b = iv("chr2", 1, 100);

        assert_eq!(
            a.compare_natural(&b),
            Ordering::Less,
            "chr10 < chr2 as strings"
        );
        assert_eq!(
            compare_coordinate(&dict, &a, &b),
            Ordering::Greater,
            "chr10 > chr2 by dictionary index"
        );
    }

    #[test]
    fn sorting_uses_the_dictionary_not_the_alphabet() {
        let dict: Vec<String> = ["chr1", "chr2", "chr10"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let mut list = IntervalList::new(dict);
        list.intervals = vec![iv("chr10", 1, 10), iv("chr1", 1, 10), iv("chr2", 1, 10)];
        let sorted = list.sorted();
        let names: Vec<&str> = sorted.intervals.iter().map(|i| i.contig.as_str()).collect();
        assert_eq!(names, ["chr1", "chr2", "chr10"]);
    }

    /// An unknown contig gets index -1 and sorts before everything, rather than raising.
    #[test]
    fn a_contig_outside_the_dictionary_sorts_first() {
        let dict = vec!["chr1".to_string()];
        let mut list = IntervalList::new(dict);
        list.intervals = vec![iv("chr1", 1, 10), iv("chrUnknown", 1, 10)];
        let sorted = list.sorted();
        assert_eq!(sorted.intervals[0].contig, "chrUnknown");
    }

    /// The two renderings differ, and only one is the file format.
    #[test]
    fn the_display_string_is_not_the_file_line() {
        let i = Interval::with_strand_and_name("chr1", 10, 20, true, Some("exon1"));
        assert_eq!(i.to_display_string(), "chr1:10-20\t-\texon1");
        assert_eq!(i.to_file_line(), "chr1\t10\t20\t-\texon1");
    }

    #[test]
    fn a_missing_name_is_written_as_a_dot_and_read_back_as_missing() {
        let i = iv("chr1", 1, 10);
        assert_eq!(i.to_file_line(), "chr1\t1\t10\t+\t.");
        let back = IntervalList::parse_body(vec!["chr1".to_string()], &i.to_file_line()).unwrap();
        assert_eq!(back.intervals[0].name, None);
    }

    /// `uniqued` combines abutting intervals, not only overlapping ones.
    #[test]
    fn abutting_intervals_are_merged() {
        let dict = vec!["chr1".to_string()];
        let mut list = IntervalList::new(dict);
        list.intervals = vec![iv("chr1", 1, 10), iv("chr1", 11, 20)];
        let u = list.uniqued(false);
        assert_eq!(u.intervals.len(), 1);
        assert_eq!((u.intervals[0].start, u.intervals[0].end), (1, 20));
    }

    #[test]
    fn separated_intervals_are_not_merged() {
        let dict = vec!["chr1".to_string()];
        let mut list = IntervalList::new(dict);
        list.intervals = vec![iv("chr1", 1, 10), iv("chr1", 12, 20)];
        assert_eq!(list.uniqued(false).intervals.len(), 2);
    }

    /// The merged interval keeps the first name of the run, not the last and not a join.
    #[test]
    fn a_merged_interval_keeps_the_first_name() {
        let dict = vec!["chr1".to_string()];
        let mut list = IntervalList::new(dict);
        list.intervals = vec![
            Interval::with_strand_and_name("chr1", 1, 10, false, Some("first")),
            Interval::with_strand_and_name("chr1", 5, 20, false, Some("second")),
        ];
        let kept = list.uniqued(false);
        assert_eq!(kept.intervals[0].name.as_deref(), Some("first"));
        let joined = list.uniqued(true);
        assert_eq!(joined.intervals[0].name.as_deref(), Some("first|second"));
    }

    #[test]
    fn padding_clamps_at_one() {
        let dict = vec!["chr1".to_string()];
        let mut list = IntervalList::new(dict);
        list.intervals = vec![iv("chr1", 5, 10)];
        let p = list.padded(100, 5, |_| 1000);
        assert_eq!((p.intervals[0].start, p.intervals[0].end), (1, 15));
    }

    #[test]
    fn a_body_round_trips() {
        let dict: Vec<String> = ["chr1", "chr2"].iter().map(|s| s.to_string()).collect();
        let mut list = IntervalList::new(dict.clone());
        list.intervals = vec![
            Interval::with_strand_and_name("chr1", 1, 10, false, Some("a")),
            Interval::with_strand_and_name("chr2", 5, 20, true, None),
        ];
        let text = list.write_body();
        let back = IntervalList::parse_body(dict, &text).unwrap();
        assert_eq!(back.intervals, list.intervals);
        assert_eq!(back.write_body(), text);
    }

    #[test]
    fn header_lines_are_skipped_by_the_body_parser() {
        let text = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\nchr1\t1\t10\t+\t.\n";
        let list = IntervalList::parse_body(vec!["chr1".to_string()], text).unwrap();
        assert_eq!(list.intervals.len(), 1);
    }

    #[test]
    fn a_short_line_is_rejected() {
        assert!(matches!(
            IntervalList::parse_body(vec![], "chr1\t1\t10\t+\n"),
            Err(ParseError::TooFewFields(_))
        ));
    }
}
