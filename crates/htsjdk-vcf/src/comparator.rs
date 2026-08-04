//! Ported from `htsjdk.variant.variantcontext.VariantContextComparator` (htsjdk 4.2.0).
//!
//! Orders variants by the position of their contig in a list, then by start. One of the two things
//! gatk-rs's `MultiVariantDataSource` needs, and therefore one of the two G1.6 handed over: the
//! multi-input walkers cannot be ported until this and `VCFUtils.smartMergeHeaders` exist.
//!
//! # The two constructors do not check the same things
//!
//! Both refuse duplicates, and they refuse different duplicates.
//!
//! From a **contig list** the index is the position in the list, so the only thing that can go
//! wrong is the same name twice, caught by comparing the map's size to the list's:
//!
//! ```java
//! if (protoContigIndexLookup.size() != contigs.size())
//!     throw new IllegalArgumentException("There are duplicate contigs/chromosomes in the input contig list.");
//! ```
//!
//! From **header lines** the index is carried by the line, so two things can go wrong and both are
//! checked, with different messages:
//!
//! ```java
//! if (protoContigIndexLookup.size() != headerLines.size())
//!     throw new IllegalArgumentException("There are duplicate contigs/chromosomes in the input header line collection.");
//! final Set<Integer> protoIndexValues = new HashSet<>(protoContigIndexLookup.values());
//! if (protoIndexValues.size() != headerLines.size())
//!     throw new IllegalArgumentException("One or more contigs share the same index number.");
//! ```
//!
//! So two contigs sharing an index is an error from header lines and cannot be expressed from a
//! list. A port that modelled one constructor and derived the other would lose that — and would
//! also miss that the two **word the empty case differently**: "One or more contigs must be in the
//! contig list." against "One or more header lines must be in the header line collection.". That
//! one was found by the oracle rather than by reading, which is what the dump is for.
//!
//! # An unknown contig is a `NullPointerException`, deliberately
//!
//! ```java
//! // Will throw NullPointerException -- happily -- if either of the chromosomes/contigs aren't
//! // present. This error checking should already have been done in the constructor but it's left
//! // in as defence anyway.
//! ```
//!
//! The comment is htsjdk's. It is not a lapse to be tidied into a `None`: a caller that sorts a
//! variant whose contig is not in the dictionary gets an exception rather than an arbitrary order,
//! and [`VariantContextComparator::compare`] returns an error for the same reason.
//!
//! # The subtraction is a subtraction, not a three-way compare
//!
//! `compare` returns `indexA - indexB`, then `startA - startB`, rather than `Integer.compare`.
//! For contig indexes that is harmless. For starts it is worth checking rather than assuming:
//! a VCF start is positive and at most `Integer.MAX_VALUE`, so the difference of two of them lies
//! in `[-(2^31-1), 2^31-1]` and cannot overflow an `i32`. The port keeps the subtraction so that
//! the *magnitude* of the returned value matches too, which a caller sorting with a comparator
//! that inspects it would see.

use std::collections::HashMap;

use crate::header::HeaderLine;
use crate::variant::VariantContext;

/// The `IllegalArgumentException`s the constructors throw, kept apart because their messages are
/// different and a caller reporting them should say which one it hit.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ComparatorError {
    /// The contig list was empty.
    NoContigs,
    /// The header-line collection was empty. A **different message** from [`Self::NoContigs`],
    /// which is the sort of thing only the oracle tells you: the two constructors both refuse an
    /// empty input and word it differently.
    NoHeaderLines,
    /// The same contig name twice, from a list.
    DuplicateInList,
    /// The same contig name twice, from header lines.
    DuplicateInHeaderLines,
    /// Two contigs carrying the same index. Only reachable from header lines.
    SharedIndex,
}

impl ComparatorError {
    /// The message htsjdk throws, which is what a caller prints.
    pub fn message(&self) -> &'static str {
        match self {
            ComparatorError::NoContigs => "One or more contigs must be in the contig list.",
            ComparatorError::NoHeaderLines => {
                "One or more header lines must be in the header line collection."
            }
            ComparatorError::DuplicateInList => {
                "There are duplicate contigs/chromosomes in the input contig list."
            }
            ComparatorError::DuplicateInHeaderLines => {
                "There are duplicate contigs/chromosomes in the input header line collection."
            }
            ComparatorError::SharedIndex => "One or more contigs share the same index number.",
        }
    }

    pub fn class(&self) -> &'static str {
        "java.lang.IllegalArgumentException"
    }
}

/// A contig this comparator was never told about.
///
/// htsjdk throws `NullPointerException` here on purpose; this is that, as a value.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UnknownContig {
    pub contig: String,
}

/// `VariantContextComparator`.
#[derive(Debug, Clone)]
pub struct VariantContextComparator {
    contig_index_lookup: HashMap<String, i32>,
}

impl VariantContextComparator {
    /// `VariantContextComparator(List<String> contigs)`: the index is the position in the list.
    pub fn from_contigs(contigs: &[String]) -> Result<Self, ComparatorError> {
        if contigs.is_empty() {
            return Err(ComparatorError::NoContigs);
        }
        let mut lookup: HashMap<String, i32> = HashMap::new();
        for (index, contig) in contigs.iter().enumerate() {
            // A repeat overwrites rather than being rejected here; the size check below is what
            // catches it, exactly as the Java does.
            lookup.insert(contig.clone(), index as i32);
        }
        if lookup.len() != contigs.len() {
            return Err(ComparatorError::DuplicateInList);
        }
        Ok(Self {
            contig_index_lookup: lookup,
        })
    }

    /// `VariantContextComparator(Collection<VCFContigHeaderLine> headerLines)`: the index is the
    /// line's own, so it can repeat, and that is a separate error.
    ///
    /// Lines that are not contig lines are not filtered out here because the Java signature cannot
    /// receive them; a caller holding a mixed header selects the contig lines first.
    pub fn from_header_lines(lines: &[HeaderLine]) -> Result<Self, ComparatorError> {
        if lines.is_empty() {
            return Err(ComparatorError::NoHeaderLines);
        }
        let mut lookup: HashMap<String, i32> = HashMap::new();
        for line in lines {
            let HeaderLine::Contig { index, fields } = line else {
                continue;
            };
            let id = fields
                .iter()
                .find(|(key, _)| key == "ID")
                .map(|(_, value)| value.clone())
                .unwrap_or_default();
            lookup.insert(id, *index);
        }
        if lookup.len() != lines.len() {
            return Err(ComparatorError::DuplicateInHeaderLines);
        }
        // The index check is the one a contig list cannot express.
        let mut indexes: Vec<i32> = lookup.values().copied().collect();
        indexes.sort_unstable();
        indexes.dedup();
        if indexes.len() != lines.len() {
            return Err(ComparatorError::SharedIndex);
        }
        Ok(Self {
            contig_index_lookup: lookup,
        })
    }

    /// `compare(VariantContext, VariantContext)`.
    ///
    /// The subtraction is htsjdk's, kept so the magnitude matches and not only the sign.
    pub fn compare(
        &self,
        first: &VariantContext,
        second: &VariantContext,
    ) -> Result<i32, UnknownContig> {
        let first_index = self.index_of(&first.contig)?;
        let second_index = self.index_of(&second.contig)?;
        let contig_compare = first_index - second_index;
        if contig_compare != 0 {
            return Ok(contig_compare);
        }
        // Starts are positive and fit in an i32, so this difference cannot overflow one.
        Ok((first.start - second.start) as i32)
    }

    fn index_of(&self, contig: &str) -> Result<i32, UnknownContig> {
        self.contig_index_lookup
            .get(contig)
            .copied()
            .ok_or_else(|| UnknownContig {
                contig: contig.to_string(),
            })
    }

    /// `isCompatible(Collection<VCFContigHeaderLine>)`: every line's contig must be known **and**
    /// carry the same index.
    pub fn is_compatible(&self, lines: &[HeaderLine]) -> bool {
        lines.iter().all(|line| {
            let HeaderLine::Contig { index, fields } = line else {
                return true;
            };
            let id = fields
                .iter()
                .find(|(key, _)| key == "ID")
                .map(|(_, value)| value.as_str())
                .unwrap_or("");
            self.contig_index_lookup.get(id) == Some(index)
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn contigs(names: &[&str]) -> Vec<String> {
        names.iter().map(|n| (*n).to_string()).collect()
    }

    fn contig_line(id: &str, index: i32) -> HeaderLine {
        HeaderLine::Contig {
            index,
            fields: vec![("ID".to_string(), id.to_string())],
        }
    }

    fn variant(contig: &str, start: i64) -> VariantContext {
        VariantContext::new(contig, start, Vec::new())
    }

    #[test]
    fn contig_order_then_start() {
        let comparator = VariantContextComparator::from_contigs(&contigs(&["chr1", "chr2"]))
            .expect("two contigs");
        assert!(
            comparator
                .compare(&variant("chr1", 10), &variant("chr2", 1))
                .unwrap()
                < 0
        );
        assert!(
            comparator
                .compare(&variant("chr2", 1), &variant("chr1", 10))
                .unwrap()
                > 0
        );
        assert_eq!(
            comparator.compare(&variant("chr1", 10), &variant("chr1", 10)),
            Ok(0)
        );
    }

    /// The value is a subtraction, not a normalised -1/0/1, and a caller can see the difference.
    #[test]
    fn the_returned_value_is_the_difference_itself() {
        let comparator =
            VariantContextComparator::from_contigs(&contigs(&["chr1", "chr2", "chr3", "chr4"]))
                .expect("four contigs");
        assert_eq!(
            comparator.compare(&variant("chr1", 1), &variant("chr4", 1)),
            Ok(-3)
        );
        assert_eq!(
            comparator.compare(&variant("chr1", 500), &variant("chr1", 100)),
            Ok(400)
        );
    }

    /// The two constructors refuse different things, which is why both are modelled.
    #[test]
    fn the_two_constructors_refuse_different_things() {
        assert_eq!(
            VariantContextComparator::from_contigs(&[]).unwrap_err(),
            ComparatorError::NoContigs
        );
        // The same situation, a different sentence. Measured, not guessed.
        assert_eq!(
            VariantContextComparator::from_header_lines(&[]).unwrap_err(),
            ComparatorError::NoHeaderLines
        );
        assert_ne!(
            ComparatorError::NoContigs.message(),
            ComparatorError::NoHeaderLines.message()
        );
        assert_eq!(
            VariantContextComparator::from_contigs(&contigs(&["chr1", "chr1"])).unwrap_err(),
            ComparatorError::DuplicateInList
        );
        assert_eq!(
            VariantContextComparator::from_header_lines(&[
                contig_line("chr1", 0),
                contig_line("chr1", 1)
            ])
            .unwrap_err(),
            ComparatorError::DuplicateInHeaderLines
        );
        // Two names sharing one index: an error a contig list cannot even express.
        assert_eq!(
            VariantContextComparator::from_header_lines(&[
                contig_line("chr1", 0),
                contig_line("chr2", 0)
            ])
            .unwrap_err(),
            ComparatorError::SharedIndex
        );
    }

    /// An unknown contig is an exception in htsjdk, on purpose, and stays one here.
    #[test]
    fn an_unknown_contig_is_refused_rather_than_ordered() {
        let comparator =
            VariantContextComparator::from_contigs(&contigs(&["chr1"])).expect("one contig");
        assert_eq!(
            comparator.compare(&variant("chr1", 1), &variant("chrX", 1)),
            Err(UnknownContig {
                contig: "chrX".to_string()
            })
        );
    }

    #[test]
    fn compatibility_needs_the_same_index_not_just_the_same_name() {
        let comparator = VariantContextComparator::from_header_lines(&[
            contig_line("chr1", 0),
            contig_line("chr2", 1),
        ])
        .expect("two lines");
        assert!(comparator.is_compatible(&[contig_line("chr1", 0)]));
        assert!(!comparator.is_compatible(&[contig_line("chr1", 1)]));
        assert!(!comparator.is_compatible(&[contig_line("chrX", 0)]));
    }
}
