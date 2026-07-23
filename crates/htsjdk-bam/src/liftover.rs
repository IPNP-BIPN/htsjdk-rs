//! UCSC-style liftover.
//!
//! Ports `htsjdk.samtools.liftover.LiftOver` and `htsjdk.samtools.liftover.Chain` at htsjdk 4.2.0:
//! a UCSC chain file describes a base-by-base correspondence between two reference builds, and
//! [`LiftOver::lift_over`] maps an [`Interval`] from the "from" build to the "to" build.
//!
//! A chain is a header line plus a list of *continuous blocks*, each a range of the "from" build
//! that lines up with an equal-length range of the "to" build; the gaps between blocks are
//! un-liftable. Chain coordinates are 0-based half-open, but the public API is Picard's 1-based
//! inclusive [`Interval`], as in htsjdk.
//!
//! Only the surface that reaches a tool's output is ported. `Chain.write`, `diagnosticLiftover`
//! and `PartialLiftover` exist in htsjdk purely to serialise a chain back to disk and to log why
//! an interval failed; neither reaches the interval-list a consumer such as
//! `LiftOverIntervalList` writes, so they are omitted. The failed-below-threshold tally
//! (`getFailedIntervalsBelowThreshold`) is kept because it is part of the observable API, but it
//! never affects a lifted interval.
//!
//! ## What is faithfully reproduced
//!
//! * **Chain parsing** including UCSC's whitespace splitting (`Pattern.compile("\\s")`, one
//!   whitespace character per split, trailing empties dropped), the 13-field header, the
//!   1-or-3-field block lines, the terminal single-field block, and every `validate()` check.
//! * **`targetIntersection`**: the walk that accumulates how many bases of the interval fall
//!   inside the chain's blocks, and the offsets into the first and last hit block.
//! * **`createToInterval`**: the 0-based to 1-based conversion, and the strand flip when the
//!   chain maps to the opposite strand (`toStart = toSequenceSize - toEnd`).
//! * **The single-hit rule**: a basic liftover that intersects two chains above the threshold
//!   maps to neither and returns `None`. Because a second qualifying hit forces `None`, the order
//!   [`OverlapDetector`] returns the candidate chains in does not change the result, so the
//!   unordered overlap set this port shares with htsjdk is sufficient.

use std::cell::Cell;
use std::collections::{HashMap, HashSet};
use std::fmt;

use crate::interval::Interval;
use crate::overlap::OverlapDetector;

/// `LiftOver.DEFAULT_LIFTOVER_MINMATCH`: the default minimum fraction of bases that must map.
pub const DEFAULT_LIFTOVER_MINMATCH: f64 = 0.95;

/// A UCSC chain file could not be parsed, or a chain failed `validate()`.
///
/// The message mirrors htsjdk's `SAMException` text; it is not part of any tool's byte output, so
/// only its presence is contractual, not its wording.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChainParseError(pub String);

impl fmt::Display for ChainParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for ChainParseError {}

/// A "to" sequence named by a chain is absent from the target sequence dictionary.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MissingToSequence(pub String);

impl fmt::Display for MissingToSequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Sequence {} from chain file is not found in sequence dictionary.",
            self.0
        )
    }
}

impl std::error::Error for MissingToSequence {}

/// `Chain.ContinuousBlock`: a range that lines up between the two builds. 0-based, half-open.
#[derive(Debug, Clone, PartialEq, Eq)]
struct ContinuousBlock {
    from_start: i32,
    to_start: i32,
    block_length: i32,
}

// `from_end`/`to_end` name the "from"/"to" build ends, mirroring htsjdk's `getFromEnd`/`getToEnd`;
// the `from_` prefix is a domain term, not the constructor convention clippy assumes.
#[allow(clippy::wrong_self_convention)]
impl ContinuousBlock {
    /// 0-based, half-open end of the region in "from".
    fn from_end(&self) -> i32 {
        self.from_start + self.block_length
    }

    /// 0-based, half-open end of the region in "to".
    fn to_end(&self) -> i32 {
        self.to_start + self.block_length
    }
}

/// A single chain from a UCSC chain file: the header line plus its continuous blocks.
#[derive(Debug, Clone)]
pub struct Chain {
    #[allow(dead_code)] // parsed and validated, but only used by Chain.write, which is not ported.
    score: f64,
    from_sequence_name: String,
    from_sequence_size: i32,
    from_chain_start: i32,
    from_chain_end: i32,
    to_sequence_name: String,
    to_sequence_size: i32,
    to_opposite_strand: bool,
    to_chain_start: i32,
    to_chain_end: i32,
    id: i32,
    blocks: Vec<ContinuousBlock>,
}

impl Chain {
    /// The chain's "from" range as a 1-based inclusive interval, the key it is stored under in the
    /// [`OverlapDetector`].
    fn interval(&self) -> Interval {
        Interval::new(
            &self.from_sequence_name,
            self.from_chain_start + 1,
            self.from_chain_end,
        )
    }
}

/// The bases of an interval that fall inside a chain's blocks, plus where they start and end.
///
/// Ports `LiftOver.TargetIntersection`.
struct TargetIntersection<'a> {
    chain: &'a Chain,
    intersection_length: i32,
    start_offset: i32,
    offset_from_end: i32,
    first_block_index: usize,
    last_block_index: usize,
}

/// `LiftOver`: a loaded chain file that maps intervals between two reference builds.
pub struct LiftOver {
    chains: OverlapDetector<Chain>,
    contig_map: HashMap<String, HashSet<String>>,
    lift_over_min_match: f64,
    log_failed_intervals: bool,
    total_failed_below_threshold: Cell<u64>,
}

impl LiftOver {
    /// `new LiftOver(...)`: build from the text of a UCSC chain file.
    pub fn load(chain_file_text: &str) -> Result<LiftOver, ChainParseError> {
        let chains = load_chains(chain_file_text)?;
        let mut contig_map: HashMap<String, HashSet<String>> = HashMap::new();
        for chain in chains.get_all() {
            contig_map
                .entry(chain.from_sequence_name.clone())
                .or_default()
                .insert(chain.to_sequence_name.clone());
        }
        Ok(LiftOver {
            chains,
            contig_map,
            lift_over_min_match: DEFAULT_LIFTOVER_MINMATCH,
            log_failed_intervals: true,
            total_failed_below_threshold: Cell::new(0),
        })
    }

    /// `getLiftOverMinMatch`.
    pub fn lift_over_min_match(&self) -> f64 {
        self.lift_over_min_match
    }

    /// `setLiftOverMinMatch`.
    pub fn set_lift_over_min_match(&mut self, value: f64) {
        self.lift_over_min_match = value;
    }

    /// `setShouldLogFailedIntervalsBelowThreshold`. Kept for API parity; this port emits no log,
    /// but the flag still gates whether a below-threshold overlap is counted the way htsjdk does.
    pub fn set_should_log_failed_intervals_below_threshold(&mut self, value: bool) {
        self.log_failed_intervals = value;
    }

    /// `getFailedIntervalsBelowThreshold`.
    pub fn failed_intervals_below_threshold(&self) -> u64 {
        self.total_failed_below_threshold.get()
    }

    /// `resetFailedIntervalsBelowThresholdCounter`.
    pub fn reset_failed_intervals_below_threshold_counter(&self) {
        self.total_failed_below_threshold.set(0);
    }

    /// `getContigMap`: the set of destination contigs for each source contig in the chain file.
    pub fn contig_map(&self) -> &HashMap<String, HashSet<String>> {
        &self.contig_map
    }

    /// `validateToSequences`: every chain's "to" sequence must be present in `dictionary` (the
    /// ordered list of sequence names of the target reference), else the first missing one is
    /// returned.
    pub fn validate_to_sequences(&self, dictionary: &[String]) -> Result<(), MissingToSequence> {
        for chain in self.chains.get_all() {
            if !dictionary.contains(&chain.to_sequence_name) {
                return Err(MissingToSequence(chain.to_sequence_name.clone()));
            }
        }
        Ok(())
    }

    /// `liftOver(interval)`: lift using this object's configured minimum match.
    pub fn lift_over(&self, interval: &Interval) -> Option<Interval> {
        self.lift_over_with_min_match(interval, self.lift_over_min_match)
    }

    /// `liftOver(interval, liftOverMinMatch)`: map `interval` to the "to" build, or `None` if it
    /// cannot be lifted (no qualifying chain, or more than one).
    ///
    /// # Panics
    ///
    /// Mirrors htsjdk's unchecked exceptions: a zero-length input, or a chain whose blocks produce
    /// an impossible target range, are data/invariant errors and panic exactly where htsjdk throws
    /// `IllegalArgumentException` / `SAMException`.
    pub fn lift_over_with_min_match(
        &self,
        interval: &Interval,
        lift_over_min_match: f64,
    ) -> Option<Interval> {
        if interval.length() == 0 {
            panic!(
                "Zero-length interval cannot be lifted over.  Interval: {}",
                interval.name.as_deref().unwrap_or("null")
            );
        }

        // Number of bases in interval that can be lifted over must be >= this.
        let min_match_size = lift_over_min_match * interval.length() as f64;

        // In basic liftOver at most one chain may qualify, so a single owned slot is enough: a
        // second qualifying hit returns `None` outright, which is why the overlap set's order does
        // not matter.
        let mut hit: Option<TargetIntersection> = None;
        let mut has_overlap_below_threshold = false;

        for chain in self
            .chains
            .get_overlaps(&interval.contig, interval.start, interval.end)
        {
            if let Some(candidate) = target_intersection(chain, interval) {
                if candidate.intersection_length as f64 >= min_match_size {
                    if hit.is_some() {
                        // In basic liftOver, multiple hits are not allowed.
                        return None;
                    }
                    hit = Some(candidate);
                } else {
                    has_overlap_below_threshold = true;
                }
            }
        }

        match hit {
            None => {
                if has_overlap_below_threshold && self.log_failed_intervals {
                    self.total_failed_below_threshold
                        .set(self.total_failed_below_threshold.get() + 1);
                }
                None
            }
            Some(ti) => Some(create_to_interval(
                interval.name.as_deref(),
                interval.negative_strand,
                &ti,
            )),
        }
    }
}

/// `LiftOver.createToInterval`.
fn create_to_interval(
    interval_name: Option<&str>,
    source_negative_strand: bool,
    ti: &TargetIntersection,
) -> Interval {
    let first = &ti.chain.blocks[ti.first_block_index];
    let last = &ti.chain.blocks[ti.last_block_index];
    let mut to_start = first.to_start + ti.start_offset;
    let mut to_end = last.to_end() - ti.offset_from_end;
    if to_end <= to_start || to_start < 0 {
        panic!(
            "Something strange lifting over interval {}",
            interval_name.unwrap_or("null")
        );
    }

    if ti.chain.to_opposite_strand {
        // Flip if query is negative.
        let negative_start = ti.chain.to_sequence_size - to_end;
        let negative_end = ti.chain.to_sequence_size - to_start;
        to_start = negative_start;
        to_end = negative_end;
    }
    // Convert to 1-based, inclusive.
    let negative_strand = if ti.chain.to_opposite_strand {
        !source_negative_strand
    } else {
        source_negative_strand
    };
    Interval::with_strand_and_name(
        &ti.chain.to_sequence_name,
        to_start + 1,
        to_end,
        negative_strand,
        interval_name,
    )
}

/// `LiftOver.targetIntersection`: add up the overlap between a chain's blocks and the interval.
fn target_intersection<'a>(
    chain: &'a Chain,
    interval: &Interval,
) -> Option<TargetIntersection<'a>> {
    let mut intersection_length = 0;
    // Convert interval to 0-based, half-open.
    let start = interval.start - 1;
    let end = interval.end;
    let mut first_block_index: i32 = -1;
    let mut last_block_index: i32 = -1;
    let mut start_offset: i32 = -1;
    let mut offset_from_end: i32 = -1;

    for (i, block) in chain.blocks.iter().enumerate() {
        if block.from_start >= end {
            break;
        } else if block.from_end() <= start {
            continue;
        }
        if first_block_index == -1 {
            first_block_index = i as i32;
            start_offset = if start > block.from_start {
                start - block.from_start
            } else {
                0
            };
        }
        last_block_index = i as i32;
        offset_from_end = if block.from_end() > end {
            block.from_end() - end
        } else {
            0
        };
        let this_intersection = end.min(block.from_end()) - start.max(block.from_start);
        if this_intersection <= 0 {
            panic!("Should have been some intersection.");
        }
        intersection_length += this_intersection;
    }
    if intersection_length == 0 {
        return None;
    }
    Some(TargetIntersection {
        chain,
        intersection_length,
        start_offset,
        offset_from_end,
        first_block_index: first_block_index as usize,
        last_block_index: last_block_index as usize,
    })
}

/// Split a line the way `Pattern.compile("\\s").split(line)` does with the default limit: split on
/// each single Java-whitespace character (interior empties kept), then drop trailing empties.
fn split_ws(line: &str) -> Vec<&str> {
    fn is_java_ws(c: char) -> bool {
        matches!(c, ' ' | '\t' | '\n' | '\u{000B}' | '\u{000C}' | '\r')
    }
    let mut parts: Vec<&str> = Vec::new();
    let mut start = 0;
    for (i, c) in line.char_indices() {
        if is_java_ws(c) {
            parts.push(&line[start..i]);
            start = i + c.len_utf8();
        }
    }
    parts.push(&line[start..]);
    while matches!(parts.last(), Some(last) if last.is_empty()) {
        parts.pop();
    }
    parts
}

/// A cursor over the chain file's lines, mirroring `BufferedLineReader`. `read_line` returns `None`
/// at end of input, matching htsjdk's `readLine() == null`; a blank line is a real `Some("")`.
struct LineCursor<'a> {
    lines: Vec<&'a str>,
    idx: usize,
}

impl<'a> LineCursor<'a> {
    fn read_line(&mut self) -> Option<&'a str> {
        let line = self.lines.get(self.idx).copied();
        if line.is_some() {
            self.idx += 1;
        }
        line
    }

    /// 1-based number of the line most recently returned, for error messages.
    fn line_number(&self) -> usize {
        self.idx
    }
}

/// `Chain.loadChains`: read every chain into an [`OverlapDetector`].
fn load_chains(text: &str) -> Result<OverlapDetector<Chain>, ChainParseError> {
    let mut cursor = LineCursor {
        lines: text.lines().collect(),
        idx: 0,
    };
    let mut detector = OverlapDetector::new(0, 0);
    while let Some(chain) = load_chain(&mut cursor)? {
        let interval = chain.interval();
        detector.add(&interval.contig, interval.start, interval.end, chain);
    }
    Ok(detector)
}

fn parse_err(message: &str, line_number: usize) -> ChainParseError {
    ChainParseError(format!("{message} in chain file at line {line_number}"))
}

/// `Chain.loadChain`: read one chain, or `None` at end of file.
fn load_chain(cursor: &mut LineCursor) -> Result<Option<Chain>, ChainParseError> {
    // Skip comment lines to reach the header line.
    let header = loop {
        match cursor.read_line() {
            None => return Ok(None),
            Some(line) => {
                if !line.starts_with('#') {
                    break line;
                }
            }
        }
    };

    let fields = split_ws(header);
    if fields.len() != 13 {
        return Err(parse_err(
            "chain line has wrong number of fields",
            cursor.line_number(),
        ));
    }
    if fields[0] != "chain" {
        return Err(parse_err(
            "chain line does not start with 'chain'",
            cursor.line_number(),
        ));
    }

    let invalid = || parse_err("Invalid field", cursor.line_number());
    let score: f64 = fields[1].parse().map_err(|_| invalid())?;
    let from_sequence_name = fields[2].to_string();
    let from_sequence_size: i32 = fields[3].parse().map_err(|_| invalid())?;
    // Field 4 (strand) is ignored because it is always +.
    let from_chain_start: i32 = fields[5].parse().map_err(|_| invalid())?;
    let from_chain_end: i32 = fields[6].parse().map_err(|_| invalid())?;
    let to_sequence_name = fields[7].to_string();
    let to_sequence_size: i32 = fields[8].parse().map_err(|_| invalid())?;
    let to_opposite_strand = fields[9] == "-";
    let to_chain_start: i32 = fields[10].parse().map_err(|_| invalid())?;
    let to_chain_end: i32 = fields[11].parse().map_err(|_| invalid())?;
    let id: i32 = fields[12].parse().map_err(|_| invalid())?;

    let mut chain = Chain {
        score,
        from_sequence_name,
        from_sequence_size,
        from_chain_start,
        from_chain_end,
        to_sequence_name,
        to_sequence_size,
        to_opposite_strand,
        to_chain_start,
        to_chain_end,
        id,
        blocks: Vec::new(),
    };

    let mut to_block_start = chain.to_chain_start;
    let mut from_block_start = chain.from_chain_start;
    let mut saw_last_line = false;
    loop {
        let line = cursor.read_line();
        match line {
            None | Some("") => {
                if !saw_last_line {
                    return Err(parse_err(
                        "Reached end of chain without seeing terminal block",
                        cursor.line_number(),
                    ));
                }
                break;
            }
            Some(line) => {
                if saw_last_line {
                    return Err(parse_err(
                        "Terminal block seen before end of chain",
                        cursor.line_number(),
                    ));
                }
                let block_fields = split_ws(line);
                if block_fields.len() == 1 {
                    saw_last_line = true;
                } else if block_fields.len() != 3 {
                    return Err(parse_err(
                        "Block line has unexpected number of fields",
                        cursor.line_number(),
                    ));
                }
                let size: i32 = block_fields[0]
                    .parse()
                    .map_err(|_| parse_err("Invalid field", cursor.line_number()))?;
                chain.blocks.push(ContinuousBlock {
                    from_start: from_block_start,
                    to_start: to_block_start,
                    block_length: size,
                });
                if !saw_last_line {
                    let from_gap: i32 = block_fields[1]
                        .parse()
                        .map_err(|_| parse_err("Invalid field", cursor.line_number()))?;
                    let to_gap: i32 = block_fields[2]
                        .parse()
                        .map_err(|_| parse_err("Invalid field", cursor.line_number()))?;
                    from_block_start += from_gap + size;
                    to_block_start += to_gap + size;
                }
            }
        }
    }

    validate(&chain)?;
    Ok(Some(chain))
}

/// `Chain.validate`: throw if the chain looks malformed.
fn validate(chain: &Chain) -> Result<(), ChainParseError> {
    let positive = |name: &str, value: i32| -> Result<(), ChainParseError> {
        if value <= 0 {
            Err(ChainParseError(format!(
                "{name} is not positive: {value} for chain {}",
                chain.id
            )))
        } else {
            Ok(())
        }
    };
    let non_negative = |name: &str, value: i32| -> Result<(), ChainParseError> {
        if value < 0 {
            Err(ChainParseError(format!(
                "{name} is negative: {value} for chain {}",
                chain.id
            )))
        } else {
            Ok(())
        }
    };

    positive("fromSequenceSize", chain.from_sequence_size)?;
    non_negative("fromChainStart", chain.from_chain_start)?;
    non_negative("fromChainEnd", chain.from_chain_end)?;
    positive("toSequenceSize", chain.to_sequence_size)?;
    non_negative("toChainStart", chain.to_chain_start)?;
    non_negative("toChainEnd", chain.to_chain_end)?;
    let from_length = chain.from_chain_end - chain.from_chain_start;
    positive("from length", from_length)?;
    let to_length = chain.to_chain_end - chain.to_chain_start;
    positive("to length", to_length)?;
    if from_length > chain.from_sequence_size {
        return Err(ChainParseError(format!(
            "From chain length ({from_length}) < from sequence length ({}) for chain {}",
            chain.from_sequence_size, chain.id
        )));
    }
    if to_length > chain.to_sequence_size {
        return Err(ChainParseError(format!(
            "To chain length ({to_length}) < to sequence length ({}) for chain {}",
            chain.to_sequence_size, chain.id
        )));
    }
    if chain.from_sequence_name.is_empty() {
        return Err(ChainParseError(format!(
            "Chain {}has empty from sequence name.",
            chain.id
        )));
    }
    if chain.to_sequence_name.is_empty() {
        return Err(ChainParseError(format!(
            "Chain {}has empty to sequence name.",
            chain.id
        )));
    }
    if chain.blocks.is_empty() {
        return Err(ChainParseError(format!(
            "Chain {} has empty block list.",
            chain.id
        )));
    }
    let first = &chain.blocks[0];
    if first.from_start != chain.from_chain_start {
        return Err(ChainParseError(format!(
            "First block from start != chain from start for chain {}",
            chain.id
        )));
    }
    if first.to_start != chain.to_chain_start {
        return Err(ChainParseError(format!(
            "First block to start != chain to start for chain {}",
            chain.id
        )));
    }
    let last = &chain.blocks[chain.blocks.len() - 1];
    if last.from_end() != chain.from_chain_end {
        return Err(ChainParseError(format!(
            "Last block from end != chain from end for chain {}",
            chain.id
        )));
    }
    if last.to_end() != chain.to_chain_end {
        return Err(ChainParseError(format!(
            "Last block to end < chain to end for chain {}",
            chain.id
        )));
    }
    for i in 1..chain.blocks.len() {
        let this_block = &chain.blocks[i];
        let prev_block = &chain.blocks[i - 1];
        if this_block.from_start < prev_block.from_end() {
            return Err(ChainParseError(format!(
                "Continuous block {i} from starts before previous block ends for chain {}",
                chain.id
            )));
        }
        if this_block.to_start < prev_block.to_end() {
            return Err(ChainParseError(format!(
                "Continuous block {i} to starts before previous block ends for chain {}",
                chain.id
            )));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A chain covering `chr1:0-100` (from) mapping to `chr2:1000-1100` (to), same strand, with a
    /// single continuous block of length 100 (no gaps). 0-based half-open in the file.
    const SIMPLE_CHAIN: &str =
        "chain 4900 chr1 249250621 + 0 100 chr2 243199373 + 1000 1100 1\n100\n\n";

    #[test]
    fn a_whole_interval_inside_one_block_lifts_by_the_block_offset() {
        let lift = LiftOver::load(SIMPLE_CHAIN).unwrap();
        // from 1-based chr1:11-20 -> 0-based 10..20; block from_start 0, to_start 1000, so
        // to 0-based 1010..1020 -> 1-based chr2:1011-1020.
        let out = lift
            .lift_over(&Interval::with_strand_and_name(
                "chr1",
                11,
                20,
                false,
                Some("x"),
            ))
            .unwrap();
        assert_eq!(out.contig, "chr2");
        assert_eq!(out.start, 1011);
        assert_eq!(out.end, 1020);
        assert!(!out.negative_strand);
        assert_eq!(out.name.as_deref(), Some("x"));
    }

    #[test]
    fn an_interval_outside_every_block_does_not_lift() {
        let lift = LiftOver::load(SIMPLE_CHAIN).unwrap();
        assert!(lift.lift_over(&Interval::new("chr1", 200, 210)).is_none());
    }

    #[test]
    fn a_two_block_chain_gaps_the_dropped_bases() {
        // Two blocks: from [0,10) -> to [1000,1010), then after a 5-base "from" gap and no "to"
        // gap, from [15,25) -> to [1010,1020). Interval chr1:1-10 sits entirely in block 0.
        let chain = "chain 100 chr1 1000 + 0 25 chr2 2000 + 1000 1020 7\n10\t5\t0\n10\n\n";
        let lift = LiftOver::load(chain).unwrap();
        let out = lift.lift_over(&Interval::new("chr1", 1, 10)).unwrap();
        assert_eq!(
            (out.contig.as_str(), out.start, out.end),
            ("chr2", 1001, 1010)
        );
    }

    #[test]
    fn a_negative_to_strand_flips_the_coordinates() {
        // to strand '-', toSequenceSize 2000, single block from [0,100) to [500,600).
        let chain = "chain 100 chr1 1000 + 0 100 chr2 2000 - 500 600 3\n100\n\n";
        let lift = LiftOver::load(chain).unwrap();
        // from 1-based chr1:11-20 -> 0-based [10,20); to 0-based [510,520) before flip.
        // flipped: toStart = 2000-520 = 1480, toEnd = 2000-510 = 1490 -> 1-based chr2:1481-1490.
        let out = lift
            .lift_over(&Interval::with_strand_and_name(
                "chr1",
                11,
                20,
                false,
                Some("y"),
            ))
            .unwrap();
        assert_eq!(
            (out.contig.as_str(), out.start, out.end),
            ("chr2", 1481, 1490)
        );
        assert!(out.negative_strand, "source + on a - chain becomes -");
    }

    #[test]
    fn below_the_min_match_threshold_does_not_lift() {
        // Block covers from [0,10); interval chr1:1-20 (20 bases) only overlaps 10 of them = 50%.
        let chain = "chain 100 chr1 1000 + 0 10 chr2 2000 + 0 10 1\n10\n\n";
        let lift = LiftOver::load(chain).unwrap();
        // Default min match 0.95 -> fails; the below-threshold counter ticks.
        assert!(lift.lift_over(&Interval::new("chr1", 1, 20)).is_none());
        assert_eq!(lift.failed_intervals_below_threshold(), 1);
        // A lower threshold (0.4) accepts the 50% overlap.
        let out = lift
            .lift_over_with_min_match(&Interval::new("chr1", 1, 20), 0.4)
            .unwrap();
        assert_eq!((out.start, out.end), (1, 10));
    }

    #[test]
    fn validate_to_sequences_reports_the_missing_contig() {
        let lift = LiftOver::load(SIMPLE_CHAIN).unwrap();
        assert!(lift.validate_to_sequences(&["chr2".to_string()]).is_ok());
        assert_eq!(
            lift.validate_to_sequences(&["chrX".to_string()]),
            Err(MissingToSequence("chr2".to_string()))
        );
    }

    #[test]
    fn contig_map_lists_destination_contigs() {
        let lift = LiftOver::load(SIMPLE_CHAIN).unwrap();
        let dests = lift.contig_map().get("chr1").unwrap();
        assert!(dests.contains("chr2"));
    }

    #[test]
    fn comment_lines_are_skipped() {
        let chain = format!("# a comment\n{SIMPLE_CHAIN}");
        let lift = LiftOver::load(&chain).unwrap();
        assert!(lift.lift_over(&Interval::new("chr1", 11, 20)).is_some());
    }

    #[test]
    fn a_header_with_the_wrong_field_count_is_an_error() {
        let err = LiftOver::load("chain 100 chr1\n10\n\n").err().unwrap();
        assert!(err.0.contains("wrong number of fields"), "{}", err.0);
    }

    #[test]
    fn a_chain_without_a_terminal_block_is_an_error() {
        // Block line has 3 fields and then EOF, never a single-field terminal line.
        let err = LiftOver::load("chain 100 chr1 1000 + 0 20 chr2 2000 + 0 20 1\n10\t0\t0\n")
            .err()
            .unwrap();
        assert!(err.0.contains("without seeing terminal block"), "{}", err.0);
    }

    #[test]
    fn split_ws_drops_trailing_empties_but_keeps_interior() {
        assert_eq!(split_ws("a b"), vec!["a", "b"]);
        assert_eq!(split_ws("a\tb"), vec!["a", "b"]);
        assert_eq!(split_ws("a  b"), vec!["a", "", "b"]);
        assert_eq!(split_ws("a "), vec!["a"]);
        assert_eq!(split_ws(" a"), vec!["", "a"]);
    }
}
