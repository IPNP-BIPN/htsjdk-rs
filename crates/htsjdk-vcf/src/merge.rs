//! Ported from `htsjdk.variant.vcf.VCFUtils.smartMergeHeaders` (htsjdk 4.2.0).
//!
//! The second of the two things gatk-rs's `MultiVariantDataSource` needs, alongside
//! [`crate::comparator`]. With both, the multi-input walkers G1.6 handed over become portable.
//!
//! Eighty lines of Java, and most of what matters in them is not the merge.
//!
//! # The key is the key plus the ID, and that makes one branch dead
//!
//! ```java
//! String key = line.getKey();
//! if (line instanceof VCFIDHeaderLine)
//!     key = key + "-" + ((VCFIDHeaderLine) line).getID();
//! ```
//!
//! So `INFO`, `FORMAT`, `FILTER` and `contig` key on `KEY-ID`, and an unstructured line keys on its
//! bare key. One consequence is worth naming rather than silently dropping: the `VCFFilterHeaderLine`
//! arm throws when the two lines' IDs differ, but two lines can only collide in the map when their
//! IDs already match. **That arm is unreachable**, and this port says so at the place it would have
//! been rather than leaving a reader to wonder where it went.
//!
//! # Every merge output carries a `fileformat` line that no source wrote
//!
//! ```java
//! private Set<VCFHeaderLine> makeGetMetaDataSet(final Set<VCFHeaderLine> headerLinesInSomeOrder) {
//!     final Set<VCFHeaderLine> lines = new LinkedHashSet<>();
//!     if (vcfHeaderVersion != null && vcfHeaderVersion.isAtLeastAsRecentAs(VCFHeaderVersion.VCF4_3)) {
//!         lines.add(new VCFHeaderLine(VCF4_3.getFormatString(), VCF4_3.getVersionString()));
//!     } else {
//!         lines.add(new VCFHeaderLine(VCF4_2.getFormatString(), VCF4_2.getVersionString()));
//!     }
//! ```
//!
//! `getMetaDataInSortedOrder` **prepends** a version line before the sorted metadata, so the merge
//! sees it as the first line of every source and it wins the `fileformat` key from the first one.
//! A merge of two headers carrying one `INFO` line each therefore returns **two** lines, not one,
//! and the extra one was written by neither source. Measured; it is not the sort of thing reading
//! the merge alone would reveal.
//!
//! It is `VCFv4.2` unless the header's **version field** is 4.3 or later — a field, not a
//! `##fileformat` line among the metadata. That distinction matters more than it looks: a header
//! assembled in memory from a set of lines has a null version whatever `fileformat` line the set
//! contains, so it prepends `VCFv4.2` and **the version policy below never fires for it**.
//!
//! # First one wins, and the two promotion arms are both no-ops
//!
//! ```java
//! } else if (compLine.getType() == VCFHeaderLineType.Integer && compOther.getType() == VCFHeaderLineType.Float) {
//!     conflictWarner.warn(line, "Promoting Integer to Float in header: " + compOther);
//!     map.put(key, compOther);
//! } else if (compLine.getType() == VCFHeaderLineType.Float && compOther.getType() == VCFHeaderLineType.Integer) {
//!     conflictWarner.warn(line, "Promoting Integer to Float in header: " + compOther);
//! }
//! ```
//!
//! The comment says "promote key to Float" in both arms and **neither arm does it**. `compOther` *is*
//! `map.get(key)`, so the `put` writes back what is already there. Measured rather than assumed,
//! because the reading is easy to get wrong: an Integer seen first stays Integer, and a Float seen
//! first stays Float. The first line always wins, and the two arms differ only in emitting the same
//! message naming the stored line.
//!
//! The one branch that does change the map is `setNumberToUnbounded()`, which **mutates the line
//! already in the map in place**
//! rather than replacing it. Rust will not alias like that, so the entry is rebuilt; the resulting
//! map is the same.
//!
//! # The input order is the sorted order, not the file order
//!
//! The loop reads `source.getMetaDataInSortedOrder()`, which is a `TreeSet`. So the output's
//! `LinkedHashMap` order is first-seen order **across** sources and **sorted** order within one —
//! not the order the lines appeared in the file. Getting this wrong would scramble the contig lines,
//! which is precisely what the method's own comment says it exists to avoid.
//!
//! # The version policy throws a different exception from everything else
//!
//! `enforceHeaderVersionMergePolicy` throws `IllegalArgumentException` where the rest of the method
//! throws `IllegalStateException`, and only when a **4.3** header meets any other version. A caller
//! that catches one and not the other behaves differently, so the two are separate variants here.

use std::collections::BTreeSet;

use crate::header::{Cardinality, HeaderLine, LineType, VcfHeader};

/// The two exceptions the merge throws, kept apart because they are different Java classes and a
/// caller catching one and not the other behaves differently.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MergeError {
    /// `IllegalStateException`. Two lines under one key that cannot be reconciled.
    Incompatible { message: String },
    /// `IllegalArgumentException`, from the version policy alone.
    IncompatibleVersion { message: String },
}

impl MergeError {
    pub fn class(&self) -> &'static str {
        match self {
            MergeError::Incompatible { .. } => "java.lang.IllegalStateException",
            MergeError::IncompatibleVersion { .. } => "java.lang.IllegalArgumentException",
        }
    }

    pub fn message(&self) -> &str {
        match self {
            MergeError::Incompatible { message } | MergeError::IncompatibleVersion { message } => {
                message
            }
        }
    }
}

/// What `HeaderConflictWarner.warn` would have printed.
///
/// Returned rather than logged: the caller decides, and a test can assert on them. htsjdk's warner
/// also suppresses repeats of the same line, which is a property of the logger rather than of the
/// merge, so it is left to the caller.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Warning {
    /// The line that triggered it, rendered.
    pub line: String,
    pub message: String,
}

/// The merge key: `getKey()`, plus `-ID` for a line that has one.
fn merge_key(line: &HeaderLine) -> String {
    match line {
        HeaderLine::Compound { key, id, .. } => format!("{key}-{id}"),
        HeaderLine::Filter { id, .. } => format!("FILTER-{id}"),
        HeaderLine::Contig { fields, .. } => {
            let id = fields
                .iter()
                .find(|(k, _)| k == "ID")
                .map(|(_, v)| v.as_str())
                .unwrap_or("");
            format!("contig-{id}")
        }
        // A structured line that is not one of the above is a `VCFSimpleHeaderLine`, which *is* a
        // `VCFIDHeaderLine`, so it keys on its ID too.
        HeaderLine::Structured { key, fields } => {
            let id = fields
                .iter()
                .find(|(k, _)| k == "ID")
                .map(|(_, v)| v.as_str())
                .unwrap_or("");
            format!("{key}-{id}")
        }
        HeaderLine::Unstructured { key, .. } => key.clone(),
    }
}

/// `VCFCompoundHeaderLine.equalsExcludingDescription`: name, count, count type and type.
///
/// Not the trailing extra fields, which htsjdk does not compare here either.
fn equals_excluding_description(
    first: (&str, &str, Cardinality, LineType),
    second: (&str, &str, Cardinality, LineType),
) -> bool {
    first.0 == second.0 && first.1 == second.1 && first.2 == second.2 && first.3 == second.3
}

/// The lines of one header in `getMetaDataInSortedOrder`, which is a `TreeSet` and therefore
/// sorted rather than in file order.
fn sorted_lines(source: &Source) -> Vec<HeaderLine> {
    let mut lines: Vec<HeaderLine> = source.header.lines.clone();
    lines.sort_by_key(HeaderLine::sort_key);
    // Prepended, not sorted in: `makeGetMetaDataSet` adds it to a fresh `LinkedHashSet` before the
    // sorted ones, so it is first however it would have sorted.
    let mut all = vec![prepended_version_line(source.version)];
    all.extend(lines);
    all
}

/// One source of the merge: a header and the version htsjdk keeps beside it.
///
/// The version is a **field** on `VCFHeader`, set when a file is parsed, and not derived from any
/// `##fileformat` line among the metadata. This port's `VcfHeader` has no such field, so it is
/// passed in — which also makes visible the thing the oracle showed: a header assembled in memory
/// has no version, so it prepends `VCFv4.2` and never reaches the version policy however its lines
/// are labelled.
#[derive(Debug, Clone, Copy)]
pub struct Source<'a> {
    pub header: &'a VcfHeader,
    /// `getVCFHeaderVersion()`, `None` for a header that was not parsed from a file.
    pub version: Option<&'a str>,
}

/// The version string htsjdk special-cases, and the only one it special-cases.
pub const VCF4_3: &str = "VCFv4.3";

/// The version every other header's prepended line carries.
pub const VCF4_2: &str = "VCFv4.2";

/// The line `makeGetMetaDataSet` puts in front of the sorted metadata.
fn prepended_version_line(version: Option<&str>) -> HeaderLine {
    // `isAtLeastAsRecentAs(VCF4_3)`: 4.3 and anything later keep 4.3, everything else, null
    // included, gets 4.2.
    let value = match version {
        Some(v) if v >= VCF4_3 => VCF4_3,
        _ => VCF4_2,
    };
    HeaderLine::Unstructured {
        key: "fileformat".to_string(),
        value: value.to_string(),
    }
}

/// `enforceHeaderVersionMergePolicy`: a 4.3 header may not be merged with any other version.
fn enforce_version_policy(
    versions: &mut BTreeSet<String>,
    candidate: Option<&str>,
) -> Result<(), MergeError> {
    let Some(candidate) = candidate else {
        return Ok(());
    };
    versions.insert(candidate.to_string());
    if versions.len() > 1 && versions.contains(VCF4_3) {
        let others: Vec<&str> = versions
            .iter()
            .filter(|v| v.as_str() != VCF4_3)
            .map(String::as_str)
            .collect();
        return Err(MergeError::IncompatibleVersion {
            message: format!(
                "Attempt to merge version {VCF4_3} header with incompatible header version {}",
                others.join(" ")
            ),
        });
    }
    Ok(())
}

/// `smartMergeHeaders(headers, emitWarnings)`.
///
/// Returns the merged lines in the `LinkedHashMap`'s order and the warnings that would have been
/// printed. `emit_warnings` gates only whether they are collected, exactly as it gates only whether
/// htsjdk logs them: the merge itself is the same either way.
pub fn smart_merge_headers(
    sources: &[Source],
    emit_warnings: bool,
) -> Result<(Vec<HeaderLine>, Vec<Warning>), MergeError> {
    // Insertion-ordered, which is the whole reason htsjdk uses a LinkedHashMap here.
    let mut order: Vec<String> = Vec::new();
    let mut map: std::collections::HashMap<String, HeaderLine> = std::collections::HashMap::new();
    let mut warnings: Vec<Warning> = Vec::new();
    let mut versions: BTreeSet<String> = BTreeSet::new();

    let mut warn = |line: &HeaderLine, message: String| {
        if emit_warnings {
            warnings.push(Warning {
                line: line.render(),
                message,
            });
        }
    };

    for source in sources {
        for line in sorted_lines(source) {
            let line = &line;
            // Inside the per-line loop, as the Java has it. Every header has at least the
            // prepended version line, so unlike the Java's `todo` suggests, the policy is reached
            // for any header at all — but only when a version was actually set.
            enforce_version_policy(&mut versions, source.version)?;

            let key = merge_key(line);
            let Some(other) = map.get(&key).cloned() else {
                order.push(key.clone());
                map.insert(key, line.clone());
                continue;
            };

            if *line == other {
                continue;
            }
            if std::mem::discriminant(line) != std::mem::discriminant(&other) {
                return Err(MergeError::Incompatible {
                    message: format!(
                        "Incompatible header types: {} {}",
                        line.render(),
                        other.render()
                    ),
                });
            }

            match (line, &other) {
                // The `VCFFilterHeaderLine` arm of the Java lives here and is UNREACHABLE: it
                // throws when the two IDs differ, and two lines only collide when their IDs are
                // already equal because the ID is part of the key. Left as a comment rather than
                // as dead code, so the shape of the original is still legible.
                (HeaderLine::Filter { .. }, HeaderLine::Filter { .. }) => {}
                (
                    HeaderLine::Compound {
                        key: line_key,
                        id: line_id,
                        number: line_number,
                        line_type: line_type_,
                        description: line_description,
                        ..
                    },
                    HeaderLine::Compound {
                        key: other_key,
                        id: other_id,
                        number: other_number,
                        line_type: other_type,
                        description: other_description,
                        ..
                    },
                ) => {
                    if !equals_excluding_description(
                        (line_key, line_id, *line_number, *line_type_),
                        (other_key, other_id, *other_number, *other_type),
                    ) {
                        if line_type_ == other_type {
                            warn(
                                line,
                                format!(
                                    "Promoting header field Number to . due to number differences \
                                     in header lines: {} {}",
                                    line.render(),
                                    other.render()
                                ),
                            );
                            // `setNumberToUnbounded()` mutates the stored line in place; the entry
                            // is rebuilt here to the same effect.
                            if let Some(HeaderLine::Compound { number, .. }) = map.get_mut(&key) {
                                *number = Cardinality::Unbounded;
                            }
                        } else if *line_type_ == LineType::Integer && *other_type == LineType::Float
                        {
                            // The Java's `map.put(key, compOther)` here is a no-op: compOther is
                            // already the stored value. Both arms keep the Float.
                            warn(
                                line,
                                format!("Promoting Integer to Float in header: {}", other.render()),
                            );
                        } else if *line_type_ == LineType::Float && *other_type == LineType::Integer
                        {
                            // Same message, naming the Integer line, which is what the Java does.
                            warn(
                                line,
                                format!("Promoting Integer to Float in header: {}", other.render()),
                            );
                        } else {
                            return Err(MergeError::Incompatible {
                                message: format!(
                                    "Incompatible header types, collision between these two types: \
                                     {} {}",
                                    line.render(),
                                    other.render()
                                ),
                            });
                        }
                    }
                    if line_description != other_description {
                        warn(
                            line,
                            format!(
                                "Allowing unequal description fields through: keeping {} excluding \
                                 {}",
                                other.render(),
                                line.render()
                            ),
                        );
                    }
                }
                _ => {
                    warn(
                        line,
                        format!(
                            "Ignoring header line already in map: this header line = {} already \
                             present header = {}",
                            line.render(),
                            other.render()
                        ),
                    );
                }
            }
        }
    }

    let merged = order
        .into_iter()
        .filter_map(|key| map.remove(&key))
        .collect();
    Ok((merged, warnings))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn compound(
        key: &str,
        id: &str,
        number: Cardinality,
        line_type: LineType,
        desc: &str,
    ) -> HeaderLine {
        HeaderLine::Compound {
            key: key.to_string(),
            id: id.to_string(),
            number,
            line_type,
            description: desc.to_string(),
            extra: Vec::new(),
        }
    }

    fn header(lines: Vec<HeaderLine>) -> VcfHeader {
        VcfHeader {
            lines,
            samples: Vec::new(),
        }
    }

    /// Two in-memory headers, which is the shape the oracle showed has **no version**: assembling
    /// a header from lines leaves the field null however the lines are labelled.
    fn merge_two(
        first: &VcfHeader,
        second: &VcfHeader,
    ) -> Result<(Vec<HeaderLine>, Vec<Warning>), MergeError> {
        smart_merge_headers(
            &[
                Source {
                    header: first,
                    version: None,
                },
                Source {
                    header: second,
                    version: None,
                },
            ],
            true,
        )
    }

    /// The line every merge carries whether a source wrote one or not.
    fn version_line(value: &str) -> HeaderLine {
        HeaderLine::Unstructured {
            key: "fileformat".to_string(),
            value: value.to_string(),
        }
    }

    #[test]
    fn every_merge_prepends_a_version_line_no_source_wrote() {
        let line = compound(
            "INFO",
            "DP",
            Cardinality::Fixed(1),
            LineType::Integer,
            "depth",
        );
        let (merged, warnings) =
            merge_two(&header(vec![line.clone()]), &header(vec![line])).expect("merges");
        // Two lines out of two sources carrying one apiece: the extra one is the prepended
        // version, and it comes first.
        assert_eq!(merged.len(), 2);
        assert_eq!(merged[0], version_line(VCF4_2));
        assert!(warnings.is_empty());
    }

    /// A header whose version field is 4.3 or later prepends 4.3 instead.
    #[test]
    fn a_four_three_header_prepends_its_own_version() {
        let empty = header(vec![]);
        let (merged, _) = smart_merge_headers(
            &[Source {
                header: &empty,
                version: Some(VCF4_3),
            }],
            true,
        )
        .expect("merges");
        assert_eq!(merged, vec![version_line(VCF4_3)]);
    }

    #[test]
    fn a_number_difference_promotes_the_stored_line_to_unbounded() {
        let first = compound("INFO", "AF", Cardinality::Fixed(1), LineType::Float, "af");
        let second = compound("INFO", "AF", Cardinality::Fixed(2), LineType::Float, "af");
        let (merged, warnings) =
            merge_two(&header(vec![first]), &header(vec![second])).expect("merges");
        let HeaderLine::Compound { number, .. } = &merged[1] else {
            panic!("a compound line after the version line");
        };
        assert_eq!(*number, Cardinality::Unbounded);
        assert!(warnings[0]
            .message
            .starts_with("Promoting header field Number to ."));
    }

    /// The finding the oracle corrected: **the first line wins in both directions**. The Java's
    /// two arms both claim to promote to Float and neither does, because the `put` writes back
    /// what the map already holds.
    #[test]
    fn integer_against_float_keeps_whichever_came_first() {
        let integer = compound("INFO", "X", Cardinality::Fixed(1), LineType::Integer, "x");
        let float = compound("INFO", "X", Cardinality::Fixed(1), LineType::Float, "x");
        for [first, second] in [
            [integer.clone(), float.clone()],
            [float.clone(), integer.clone()],
        ] {
            let (merged, warnings) =
                merge_two(&header(vec![first.clone()]), &header(vec![second])).expect("merges");
            assert_eq!(merged[1], first, "the first line keeps the slot");
            assert_eq!(
                warnings[0].message,
                format!("Promoting Integer to Float in header: {}", first.render()),
                "the message names the stored line, whichever type it is"
            );
        }
    }

    #[test]
    fn an_unpromotable_type_collision_is_refused() {
        let integer = compound("INFO", "X", Cardinality::Fixed(1), LineType::Integer, "x");
        let string = compound("INFO", "X", Cardinality::Fixed(1), LineType::String, "x");
        let error = merge_two(&header(vec![integer]), &header(vec![string])).expect_err("refuses");
        assert_eq!(error.class(), "java.lang.IllegalStateException");
        assert!(error
            .message()
            .starts_with("Incompatible header types, collision between these two types:"));
    }

    #[test]
    fn a_description_difference_is_a_warning_and_the_stored_one_wins() {
        let first = compound(
            "INFO",
            "DP",
            Cardinality::Fixed(1),
            LineType::Integer,
            "one",
        );
        let second = compound(
            "INFO",
            "DP",
            Cardinality::Fixed(1),
            LineType::Integer,
            "two",
        );
        let (merged, warnings) =
            merge_two(&header(vec![first.clone()]), &header(vec![second])).expect("merges");
        assert_eq!(merged[1], first);
        assert!(warnings[0]
            .message
            .starts_with("Allowing unequal description fields through: keeping"));
    }

    /// The version policy is the only `IllegalArgumentException` in the method — and it fires only
    /// when the version was actually set, which an in-memory header never has.
    #[test]
    fn a_four_three_version_may_not_meet_another_one() {
        let line = compound("INFO", "DP", Cardinality::Fixed(1), LineType::Integer, "d");
        let first = header(vec![line.clone()]);
        let second = header(vec![line]);
        let error = smart_merge_headers(
            &[
                Source {
                    header: &first,
                    version: Some(VCF4_3),
                },
                Source {
                    header: &second,
                    version: Some(VCF4_2),
                },
            ],
            true,
        )
        .expect_err("refuses");
        assert_eq!(error.class(), "java.lang.IllegalArgumentException");
        assert_eq!(
            error.message(),
            "Attempt to merge version VCFv4.3 header with incompatible header version VCFv4.2"
        );
        // The same two headers with no version set merge happily, which is what an in-memory
        // header does and what the oracle showed.
        assert!(merge_two(&first, &second).is_ok());
    }

    #[test]
    fn the_output_is_first_seen_order_across_sources_and_sorted_within_one() {
        let first = header(vec![
            compound("INFO", "B", Cardinality::Fixed(1), LineType::Integer, "b"),
            compound("INFO", "A", Cardinality::Fixed(1), LineType::Integer, "a"),
        ]);
        let second = header(vec![compound(
            "INFO",
            "C",
            Cardinality::Fixed(1),
            LineType::Integer,
            "c",
        )]);
        let (merged, _) = merge_two(&first, &second).expect("merges");
        let ids: Vec<String> = merged
            .iter()
            .skip(1)
            .map(|line| match line {
                HeaderLine::Compound { id, .. } => id.clone(),
                other => other.render(),
            })
            .collect();
        // A before B because the first source is read in *sorted* order, not the order given.
        assert_eq!(ids, vec!["A", "B", "C"]);
    }
}
