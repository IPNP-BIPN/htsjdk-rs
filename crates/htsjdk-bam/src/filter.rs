//! `htsjdk.samtools.filter`: the record predicates a tool composes rather than writes.
//!
//! Each one answers `filterOut`, which is **true when the record is dropped**. The naming is
//! htsjdk's and it is worth keeping, because inverting it is the easiest possible mistake: a filter
//! called `AlignedFilter(true)` keeps aligned reads by *filtering out* everything else.
//!
//! Every filter here has a **pair** form as well, and the pair form is not the single form applied
//! twice. Each is asymmetric in its own way, and each asymmetry is a decision the reference made:
//!
//! * [`AlignedFilter`] with `include_aligned` keeps a pair only when **both** ends are mapped, and
//!   without it keeps a pair when **either** end is unmapped.
//! * [`ReadNameFilter`] with `include_reads` keeps a pair only when **both** names are listed, and
//!   without it keeps a pair only when **neither** is.
//! * [`TagFilter`] with `include_reads` keeps a pair when **either** end carries a listed value,
//!   and without it drops a pair only when **both** do.
//!
//! `picard-rs` ported these inside `FilterSamReads`, which is the tool that composes them; they are
//! htsjdk's, and GATK's read filters reach the same package.

use crate::record::BamRecord;
use crate::tag::{Tag, TagValue};

const READ_UNMAPPED: u16 = 0x4;

fn is_unmapped(record: &BamRecord) -> bool {
    record.flags & READ_UNMAPPED != 0
}

/// `SamRecordFilter`: `filterOut` is true when the record is **dropped**.
pub trait SamRecordFilter {
    fn filter_out(&self, record: &BamRecord) -> bool;
    fn filter_out_pair(&self, first: &BamRecord, second: &BamRecord) -> bool;
}

/// `htsjdk.samtools.filter.AlignedFilter`.
pub struct AlignedFilter {
    pub include_aligned: bool,
}

impl SamRecordFilter for AlignedFilter {
    fn filter_out(&self, record: &BamRecord) -> bool {
        if self.include_aligned {
            is_unmapped(record)
        } else {
            !is_unmapped(record)
        }
    }

    fn filter_out_pair(&self, first: &BamRecord, second: &BamRecord) -> bool {
        if self.include_aligned {
            // Both ends must be mapped for the pair to survive.
            is_unmapped(first) || is_unmapped(second)
        } else {
            // Either end unmapped is enough to keep it, which is not the negation of the above.
            !(is_unmapped(first) || is_unmapped(second))
        }
    }
}

/// `htsjdk.samtools.filter.ReadNameFilter`.
pub struct ReadNameFilter {
    pub names: std::collections::HashSet<String>,
    pub include_reads: bool,
}

impl ReadNameFilter {
    /// The file constructor's parsing: a blank line is skipped, and a line's name is everything
    /// before its first run of whitespace (`line.split("\\s+")[0]`), so a name-plus-comment file
    /// reads as names.
    pub fn from_lines(text: &str, include_reads: bool) -> Self {
        let mut names = std::collections::HashSet::new();
        for line in text.lines() {
            if line.trim().is_empty() {
                continue;
            }
            names.insert(line.split_whitespace().next().unwrap_or("").to_string());
        }
        ReadNameFilter {
            names,
            include_reads,
        }
    }
}

impl SamRecordFilter for ReadNameFilter {
    fn filter_out(&self, record: &BamRecord) -> bool {
        self.names.contains(&record.read_name) != self.include_reads
    }

    fn filter_out_pair(&self, first: &BamRecord, second: &BamRecord) -> bool {
        let has_first = self.names.contains(&first.read_name);
        let has_second = self.names.contains(&second.read_name);
        if self.include_reads {
            !(has_first && has_second)
        } else {
            has_first || has_second
        }
    }
}

/// `htsjdk.samtools.filter.TagFilter`.
///
/// The comparison is `values.contains(record.getAttribute(tag))`, which is Java equality between
/// the boxed attribute and the listed value: a record with no such tag matches nothing, and a
/// numeric tag matches only a numerically equal listed value of the same box.
pub struct TagFilter {
    pub tag: Tag,
    pub values: Vec<TagValue>,
    pub include_reads: bool,
}

impl TagFilter {
    fn matches(&self, record: &BamRecord) -> bool {
        match record.tags.get(self.tag) {
            Some(actual) => self.values.iter().any(|value| value == actual),
            None => false,
        }
    }
}

impl SamRecordFilter for TagFilter {
    fn filter_out(&self, record: &BamRecord) -> bool {
        self.matches(record) != self.include_reads
    }

    fn filter_out_pair(&self, first: &BamRecord, second: &BamRecord) -> bool {
        if self.include_reads {
            // Any pair carrying the value survives.
            !(self.matches(first) || self.matches(second))
        } else {
            // Only a pair where BOTH ends carry it is dropped.
            self.matches(first) && self.matches(second)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(name: &str, unmapped: bool, rg: Option<&str>) -> BamRecord {
        let mut record = BamRecord {
            read_name: name.to_string(),
            flags: if unmapped { READ_UNMAPPED } else { 0 },
            ..Default::default()
        };
        if let Some(rg) = rg {
            record
                .tags
                .insert(Tag::new(b"RG"), TagValue::Str(rg.to_string()));
        }
        record
    }

    #[test]
    fn aligned_filter_keeps_what_it_says_it_keeps() {
        let include = AlignedFilter {
            include_aligned: true,
        };
        assert!(!include.filter_out(&record("a", false, None)));
        assert!(include.filter_out(&record("a", true, None)));
        let exclude = AlignedFilter {
            include_aligned: false,
        };
        assert!(exclude.filter_out(&record("a", false, None)));
        assert!(!exclude.filter_out(&record("a", true, None)));
    }

    #[test]
    fn the_pair_form_is_not_the_single_form_twice() {
        let mapped = record("a", false, None);
        let unmapped = record("a", true, None);
        // include_aligned: both ends must be mapped.
        let include = AlignedFilter {
            include_aligned: true,
        };
        assert!(include.filter_out_pair(&mapped, &unmapped));
        // exclude_aligned: either end unmapped keeps the pair, so the same pair survives here.
        let exclude = AlignedFilter {
            include_aligned: false,
        };
        assert!(!exclude.filter_out_pair(&mapped, &unmapped));
    }

    #[test]
    fn a_read_name_file_takes_the_first_field_of_a_line() {
        let filter = ReadNameFilter::from_lines("read1\nread2 a comment\n\n  \nread3\t2\n", true);
        assert_eq!(filter.names.len(), 3);
        assert!(filter.names.contains("read2"));
        assert!(!filter.filter_out(&record("read2", false, None)));
        assert!(filter.filter_out(&record("read9", false, None)));
    }

    #[test]
    fn a_missing_tag_matches_nothing() {
        let filter = TagFilter {
            tag: Tag::new(b"RG"),
            values: vec![TagValue::Str("rg1".to_string())],
            include_reads: true,
        };
        // include_reads: a record without the tag does not match, so it is filtered out.
        assert!(filter.filter_out(&record("a", false, None)));
        assert!(!filter.filter_out(&record("a", false, Some("rg1"))));
    }

    #[test]
    fn the_tag_pair_form_drops_only_when_both_carry_it() {
        let filter = TagFilter {
            tag: Tag::new(b"RG"),
            values: vec![TagValue::Str("rg1".to_string())],
            include_reads: false,
        };
        let tagged = record("a", false, Some("rg1"));
        let other = record("a", false, Some("rg2"));
        assert!(filter.filter_out_pair(&tagged, &tagged));
        assert!(!filter.filter_out_pair(&tagged, &other));
    }
}
