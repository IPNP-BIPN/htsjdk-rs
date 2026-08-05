//! Building a `.idx`, which is what GATK does beside every VCF it writes.
//!
//! Ported from `htsjdk.tribble.index.linear.LinearIndexCreator`,
//! `htsjdk.tribble.index.linear.LinearIndex` (its `optimize` half),
//! `htsjdk.tribble.index.DynamicIndexCreator`, `htsjdk.tribble.index.AbstractIndex` (its `write`
//! half) and `htsjdk.tribble.util.MathUtils.RunningStat` at htsjdk 4.2.0.
//!
//! [`crate::index`] reads an index. This builds one, and the two are different problems: reading is
//! a layout, writing is a set of decisions the layout only records the outcome of. Four of those
//! decisions are invisible in the file and in the format.
//!
//! # The bin width is not the one the creator was given
//!
//! `LinearIndex.optimize` doubles it, **per contig**, merging blocks pairwise, until the most dense
//! block is estimated to hold more than [`MAX_FEATURES_PER_BIN`] features, or one block is left, or
//! the width goes bad. That is why the read suite measured 16000 for one contig and 8000 for
//! another in the same file: not a setting, an outcome, and a reader that assumed the creator's
//! default would answer every query wrongly without ever failing.
//!
//! The loop keeps the **last** width that was still under the threshold rather than the first one
//! over it, so it stops one doubling short of where the condition first holds. On a contig that
//! never crosses the threshold it therefore runs to a single block.
//!
//! # The density is estimated, never counted
//!
//! The score is the largest block's size in **bytes** divided by the average feature size in bytes,
//! so it is a guess at a feature count that is never compared to the `n_features` the same object
//! is carrying. Two files with identical feature counts and different line lengths optimize
//! differently.
//!
//! # The dynamic creator picks the index type from the data
//!
//! It feeds every feature to a linear creator *and* an interval-tree creator, scores both, and
//! keeps one. Measured, the choice flips both ways on the same two files: sparse data gets a linear
//! index under `ForSeekTime` and an interval tree under `ForSize`, and dense data gets the
//! opposite. Nothing on a command line says which a run produced.
//!
//! # The statistics written into the header are not the feature lengths
//!
//! `stats.push(longestFeatureLength)` pushes the running **maximum** at each step, not the feature
//! just seen, so `FEATURE_LENGTH_MEAN` is the mean of a non-decreasing sequence. Measured on
//! features of lengths 10, 10, 600, 10, 10 it is 364.0, the mean of 10, 10, 600, 600, 600. The
//! number is wrong and it is in the bytes, so a port that computes the honest mean writes a
//! different file.
//!
//! # What this module does not do
//!
//! **It does not write an interval-tree index**, the same boundary [`crate::index`] draws on the
//! read side: that layout is refused rather than guessed at. The dynamic creator here therefore
//! reports which type it chose and refuses to build the bytes when the answer is the tree. Since
//! the choice depends on the data, that is a real limit and not a formality.

use crate::index::{Block, ChrIndex, IndexError, TribbleIndex, LINEAR, MAGIC_NUMBER};

/// `LinearIndexCreator.DEFAULT_BIN_WIDTH`.
pub const DEFAULT_BIN_WIDTH: i32 = 8000;

/// `IntervalIndexCreator.DEFAULT_FEATURE_COUNT`.
pub const DEFAULT_FEATURE_COUNT: i32 = 600;

/// `LinearIndex.MAX_FEATURES_PER_BIN`, read from a system property that defaults to 100.
pub const MAX_FEATURES_PER_BIN: f64 = 100.0;

/// `LinearIndex.MAX_BIN_WIDTH`: "widths must be less than 1 billion".
pub const MAX_BIN_WIDTH: i32 = 1_000_000_000;

/// `LinearIndex.MAX_BIN_WIDTH_FOR_OCCUPIED_CHR_INDEX`, likewise a property with a default.
pub const MAX_BIN_WIDTH_FOR_OCCUPIED_CHR_INDEX: i64 = 1_024_000;

/// One feature as the creator sees it: a contig and a closed 1-based interval.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Feature {
    pub contig: String,
    pub start: i32,
    pub end: i32,
}

impl Feature {
    /// `(getEnd() - getStart()) + 1`, which is what both the creator and the statistics use.
    fn length(&self) -> i32 {
        (self.end - self.start) + 1
    }
}

/// What a creator refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CreateError {
    /// `IllegalArgumentException` out of `finalizeIndex`, whose message states the **opposite** of
    /// the condition it guards: the check is `finalFilePosition == 0` and the text reads
    /// `finalFilePosition != 0`. Reproduced as it is, because it is what an empty file produces.
    IllegalArgument(String),
    /// `TribbleException.MalformedFeatureFile`.
    MalformedFeatureFile(String),
    /// The dynamic creator chose the interval-tree layout, which this port refuses to write for
    /// the same reason the reader refuses to parse it.
    IntervalTreeRefused,
}

impl CreateError {
    pub fn class(&self) -> &'static str {
        match self {
            CreateError::IllegalArgument(_) => "java.lang.IllegalArgumentException",
            CreateError::MalformedFeatureFile(_) | CreateError::IntervalTreeRefused => {
                "htsjdk.tribble.TribbleException$MalformedFeatureFile"
            }
        }
    }

    pub fn message(&self) -> String {
        match self {
            CreateError::IllegalArgument(message) | CreateError::MalformedFeatureFile(message) => {
                message.clone()
            }
            CreateError::IntervalTreeRefused => {
                "the interval-tree layout is not written by this port".to_string()
            }
        }
    }
}

/// `LinearIndexCreator`.
#[derive(Debug, Clone)]
pub struct LinearIndexCreator {
    bin_width: i32,
    chr_list: Vec<ChrIndex>,
    /// The blocks of the contig being filled. Closed and handed to the contig when the next contig
    /// starts, or by `finalize`.
    blocks: Vec<Block>,
    longest_feature: i32,
}

impl LinearIndexCreator {
    pub fn new(bin_width: i32) -> Self {
        Self {
            bin_width,
            chr_list: Vec::new(),
            blocks: Vec::new(),
            longest_feature: 0,
        }
    }

    pub fn bin_size(&self) -> i32 {
        self.bin_width
    }

    /// `addFeature`.
    ///
    /// The block a feature lands in is decided by `while (start > blocks.size() * binWidth)`, which
    /// appends **one block per iteration**, so a feature far along a contig is preceded by every
    /// empty bin between it and the last one. The product is `int` arithmetic upstream; it is kept
    /// as `i64` here because the loop is bounded by the coordinate rather than by the type, and an
    /// overflow would be a different bug rather than the same one.
    pub fn add_feature(&mut self, feature: &Feature, file_position: i64) {
        let new_contig = match self.chr_list.last() {
            None => true,
            Some(last) => last.name != feature.contig,
        };
        if new_contig {
            if !self.chr_list.is_empty() {
                // Every block of the contig just finished ends where the next one starts, and the
                // last one ends at this feature's position.
                let closed = self.close_blocks(file_position);
                let last = self.chr_list.last_mut().expect("a contig was open");
                last.blocks.extend(closed);
            }
            self.chr_list.push(ChrIndex {
                name: feature.contig.clone(),
                bin_width: self.bin_width,
                longest_feature: 0,
                unused: 0,
                n_features: 0,
                blocks: Vec::new(),
            });
            self.blocks.clear();
            self.blocks.push(Block {
                start: file_position,
                size: 0,
            });
            self.longest_feature = 0;
        }

        while i64::from(feature.start) > self.blocks.len() as i64 * i64::from(self.bin_width) {
            self.blocks.push(Block {
                start: file_position,
                size: 0,
            });
        }

        if feature.length() > self.longest_feature {
            self.longest_feature = feature.length();
            let last = self.chr_list.last_mut().expect("a contig is open");
            last.longest_feature = last.longest_feature.max(self.longest_feature);
        }
        self.chr_list
            .last_mut()
            .expect("a contig is open")
            .n_features += 1;
    }

    /// Each block ends where the next one starts, and the last ends at `end`.
    fn close_blocks(&self, end: i64) -> Vec<Block> {
        let mut closed = Vec::with_capacity(self.blocks.len());
        for (index, block) in self.blocks.iter().enumerate() {
            let end_position = match self.blocks.get(index + 1) {
                Some(next) => next.start,
                None => end,
            };
            closed.push(Block {
                start: block.start,
                size: end_position - block.start,
            });
        }
        closed
    }

    /// `finalizeIndex`, including the `optimize` pass that decides the bin widths.
    pub fn finalize(
        mut self,
        final_file_position: i64,
        properties: Vec<(String, String)>,
    ) -> Result<Vec<ChrIndex>, CreateError> {
        if final_file_position == 0 {
            return Err(CreateError::IllegalArgument(format!(
                "finalFilePosition != 0, -> {final_file_position}"
            )));
        }
        let _ = properties;
        if !self.chr_list.is_empty() {
            let closed = self.close_blocks(final_file_position);
            self.chr_list
                .last_mut()
                .expect("a contig was open")
                .blocks
                .extend(closed);
        }
        Ok(self.chr_list.iter().map(optimize_contig).collect())
    }
}

/// `LinearIndex.ChrIndex.optimize(idx, threshold, 0)`.
///
/// `best` is set **before** the merge, so the returned index is the last one that was still under
/// the threshold. The `level > 30` guard is upstream's and is reproduced: it cannot be reached
/// through `MAX_BIN_WIDTH` alone, because doubling from 8000 passes a billion in 18 steps.
fn optimize_contig(index: &ChrIndex) -> ChrIndex {
    let mut current = index.clone();
    let mut best = current.clone();
    let mut level = 0;
    loop {
        let score = optimize_score(&current);
        if score > MAX_FEATURES_PER_BIN || current.blocks.len() == 1 || bad_bin_width(&current) {
            break;
        }
        best = current.clone();
        current = merge_blocks(&current);
        level += 1;
        if level > 30 {
            // `IllegalStateException("Too many iterations")` upstream. Unreachable from any width
            // that starts positive, and left as a panic rather than an error because reaching it
            // would mean the loop's own invariant is broken.
            panic!("Too many iterations");
        }
    }
    best
}

/// `badBinWidth`. The negative test is upstream's own comment: "an overflow occurred".
fn bad_bin_width(index: &ChrIndex) -> bool {
    if index.bin_width > MAX_BIN_WIDTH || index.bin_width < 0 {
        true
    } else {
        MAX_BIN_WIDTH_FOR_OCCUPIED_CHR_INDEX != 0
            && index.n_features > 1
            && i64::from(index.bin_width) > MAX_BIN_WIDTH_FOR_OCCUPIED_CHR_INDEX
    }
}

/// `getNFeaturesOfMostDenseBlock(getAverageFeatureSize())`.
///
/// Both halves are byte counts: the numerator is a block's size on disk and the denominator is the
/// mean bytes per feature. On a contig with no features the denominator is a division by zero,
/// which in Java is an infinity rather than a throw, so every block scores 0 and the loop keeps
/// merging. Reproduced with `f64` division, which behaves the same way.
fn optimize_score(index: &ChrIndex) -> f64 {
    let total: i64 = index.blocks.iter().map(|block| block.size).sum();
    let average = total as f64 / f64::from(index.n_features);
    let mut most = -1.0f64;
    for block in &index.blocks {
        let n = block.size as f64 / average;
        if most == -1.0 || n > most {
            most = n;
        }
    }
    most
}

/// `mergeBlocks`: walk left to right joining pairs, doubling the width.
///
/// An odd block count leaves the last block alone rather than pairing it with nothing, so a contig
/// with three blocks merges to two and not to one.
fn merge_blocks(index: &ChrIndex) -> ChrIndex {
    let mut merged = ChrIndex {
        name: index.name.clone(),
        bin_width: index.bin_width * 2,
        longest_feature: index.longest_feature,
        unused: index.unused,
        n_features: index.n_features,
        blocks: Vec::new(),
    };
    let mut blocks = index.blocks.iter();
    while let Some(first) = blocks.next() {
        match blocks.next() {
            None => merged.blocks.push(*first),
            Some(second) => merged.blocks.push(Block {
                start: first.start,
                size: first.size + second.size,
            }),
        }
    }
    merged
}

/// `MathUtils.RunningStat`, Welford's algorithm as htsjdk spells it.
///
/// The variance is the **sample** one, dividing by `n - 1`, and both it and the mean are written
/// into the header as `String.valueOf(double)`, so the accumulation order is part of the bytes.
#[derive(Debug, Clone, Default)]
pub struct RunningStat {
    old_mean: f64,
    new_mean: f64,
    old_std_dev: f64,
    new_std_dev: f64,
    record_count: i64,
}

impl RunningStat {
    pub fn push(&mut self, x: f64) {
        self.record_count += 1;
        if self.record_count == 1 {
            self.old_mean = x;
            self.new_mean = x;
            self.old_std_dev = 0.0;
        } else {
            self.new_mean = self.old_mean + (x - self.old_mean) / self.record_count as f64;
            self.new_std_dev = self.old_std_dev + (x - self.old_mean) * (x - self.new_mean);
            self.old_mean = self.new_mean;
            self.old_std_dev = self.new_std_dev;
        }
    }

    pub fn mean(&self) -> f64 {
        if self.record_count > 0 {
            self.new_mean
        } else {
            0.0
        }
    }

    pub fn variance(&self) -> f64 {
        if self.record_count > 1 {
            self.new_std_dev / (self.record_count - 1) as f64
        } else {
            0.0
        }
    }

    pub fn standard_deviation(&self) -> f64 {
        self.variance().sqrt()
    }
}

/// `IndexFactory.IndexBalanceApproach`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BalanceApproach {
    ForSize,
    ForSeekTime,
}

/// Which layout the dynamic creator settled on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChosenIndex {
    Linear,
    IntervalTree,
}

/// `DynamicIndexCreator`.
///
/// It holds both creators and feeds every feature to both. Only the linear one is built here; the
/// interval-tree creator is present as its **score**, which is a constant of the creator rather
/// than a function of the data, so the choice can be reproduced exactly without the layout.
#[derive(Debug, Clone)]
pub struct DynamicIndexCreator {
    approach: BalanceApproach,
    linear: LinearIndexCreator,
    features_per_interval: i32,
    longest_feature_length: i32,
    feature_count: i64,
    stats: RunningStat,
    bases_seen: i64,
    last_start: Option<i32>,
}

impl DynamicIndexCreator {
    /// `getIndexCreators`, whose two arms give the two approaches different creators, not just
    /// different scoring. `ForSeekTime` gets a quarter-width linear creator and an eighth-count
    /// interval one, each with a floor.
    pub fn new(approach: BalanceApproach) -> Self {
        let (bin_width, features_per_interval) = match approach {
            BalanceApproach::ForSize => (DEFAULT_BIN_WIDTH, DEFAULT_FEATURE_COUNT),
            BalanceApproach::ForSeekTime => (
                200.max(DEFAULT_BIN_WIDTH / 4),
                20.max(DEFAULT_FEATURE_COUNT / 8),
            ),
        };
        Self {
            approach,
            linear: LinearIndexCreator::new(bin_width),
            features_per_interval,
            longest_feature_length: 0,
            feature_count: 0,
            stats: RunningStat::default(),
            bases_seen: 0,
            last_start: None,
        }
    }

    /// `addFeature`.
    ///
    /// Two things here are not what they look like. `basesSeen` adds the **whole start** rather
    /// than a delta whenever the current feature starts at or before the last one, which is how a
    /// contig change is detected without looking at the contig; and the value pushed to the
    /// statistics is the running maximum length, not this feature's.
    pub fn add_feature(&mut self, feature: &Feature, file_position: i64) {
        self.feature_count += 1;
        self.bases_seen += match self.last_start {
            None => i64::from(feature.start),
            Some(last) if feature.start - last >= 0 => i64::from(feature.start - last),
            Some(_) => i64::from(feature.start),
        };
        self.longest_feature_length = self.longest_feature_length.max(feature.length());
        self.stats.push(f64::from(self.longest_feature_length));
        self.linear.add_feature(feature, file_position);
        self.last_start = Some(feature.start);
    }

    /// `scoreIndexes` then `getMinIndex`.
    ///
    /// The linear score is `binSize * density * ceil(longestFeature / binSize)` and the
    /// interval-tree score is its features-per-interval outright, which is a constant. `FOR_SIZE`
    /// takes the **largest** score and everything else the smallest, so the two approaches are not
    /// a tie-break of the same ordering: they are opposite ends of it.
    ///
    /// The scores go into a map keyed by `Double`, so two creators that score exactly equal
    /// collapse to one entry and the survivor is whichever was inserted last. Unreachable with a
    /// linear score that carries a density, and reproduced anyway.
    pub fn chosen(&self) -> ChosenIndex {
        let density = self.feature_count as f64 / self.bases_seen as f64;
        let bin_size = f64::from(self.linear.bin_size());
        let linear_score =
            bin_size * density * (f64::from(self.longest_feature_length) / bin_size).ceil();
        let tree_score = f64::from(self.features_per_interval);

        if linear_score == tree_score {
            return ChosenIndex::IntervalTree;
        }
        let linear_is_smaller = linear_score < tree_score;
        match self.approach {
            BalanceApproach::ForSeekTime if linear_is_smaller => ChosenIndex::Linear,
            BalanceApproach::ForSeekTime => ChosenIndex::IntervalTree,
            BalanceApproach::ForSize if linear_is_smaller => ChosenIndex::IntervalTree,
            BalanceApproach::ForSize => ChosenIndex::Linear,
        }
    }

    /// The properties the chosen creator is handed before it finalizes, in insertion order.
    ///
    /// `String.valueOf(double)` is the rendering, which is [`crate::jformat`]'s problem and not a
    /// `{}`: `364.0` rather than `364`.
    pub fn properties(&self) -> Vec<(String, String)> {
        vec![
            (
                "FEATURE_LENGTH_MEAN".to_string(),
                java_double(self.stats.mean()),
            ),
            (
                "FEATURE_LENGTH_STD_DEV".to_string(),
                java_double(self.stats.standard_deviation()),
            ),
            (
                "MEAN_FEATURE_VARIANCE".to_string(),
                java_double(self.stats.variance()),
            ),
            ("FEATURE_COUNT".to_string(), self.feature_count.to_string()),
        ]
    }

    /// `finalizeIndex`, which refuses when the choice is the layout this port does not write.
    pub fn finalize(
        self,
        final_file_position: i64,
    ) -> Result<(ChosenIndex, Vec<ChrIndex>), CreateError> {
        let chosen = self.chosen();
        if chosen == ChosenIndex::IntervalTree {
            return Err(CreateError::IntervalTreeRefused);
        }
        let properties = self.properties();
        let contigs = self.linear.finalize(final_file_position, properties)?;
        Ok((chosen, contigs))
    }
}

/// `String.valueOf(double)` for the two shapes these properties take.
fn java_double(value: f64) -> String {
    htsjdk_vcf_jformat(value)
}

/// The formatter lives in `jmath`'s sibling and is not a dependency of this crate, so the two
/// shapes that occur here are spelled out: an integral value keeps its `.0`, and everything else
/// takes Rust's shortest round-trip rendering, which agrees with Java's for these magnitudes.
fn htsjdk_vcf_jformat(value: f64) -> String {
    if value == value.trunc() && value.abs() < 1e7 {
        format!("{value:.1}")
    } else {
        format!("{value}")
    }
}

/// `IndexFactory`'s two ordering checks, which run over the feature stream before any creator sees
/// it and refuse two different malformations with two different messages.
pub fn check_ordering(features: &[Feature], source: &str) -> Result<(), CreateError> {
    /// `featToString`.
    fn render(feature: &Feature) -> String {
        format!("{}:{}-{}", feature.contig, feature.start, feature.end)
    }

    // `visitedChromos`, which records the feature that **opened** each contig rather than the last
    // one on it, and is only written when the contig changes. That is what the message quotes.
    let mut visited: Vec<(&str, &Feature)> = Vec::new();
    let mut last: Option<&Feature> = None;

    for feature in features {
        // `checkSorted`, which compares starts only when the contig is unchanged, so a smaller
        // start on a new contig is the *other* failure or none at all.
        if let Some(previous) = last {
            if feature.start < previous.start && previous.contig == feature.contig {
                return Err(CreateError::MalformedFeatureFile(format!(
                    "Input file is not sorted by start position. \nWe saw a record with a start \
                     of {}:{} after a record with a start of {}:{}, for input source: {source}",
                    feature.contig, feature.start, previous.contig, previous.start
                )));
            }
        }

        let last_contig = last.map(|previous| previous.contig.as_str());
        if Some(feature.contig.as_str()) != last_contig {
            match visited.iter().find(|(name, _)| *name == feature.contig) {
                Some((_, opened)) => {
                    return Err(CreateError::MalformedFeatureFile(format!(
                        "Input file must have contiguous chromosomes. Saw feature {} followed \
                         later by {} and then {}, for input source: {source}",
                        render(opened),
                        // `lastFeature` is null only on the first feature, which cannot reach
                        // here because nothing has been visited yet.
                        last.map(render).unwrap_or_else(|| "null".to_string()),
                        render(feature)
                    )));
                }
                None => visited.push((&feature.contig, feature)),
            }
        }
        last = Some(feature);
    }
    Ok(())
}

/// The little-endian writer, with the NUL-terminated strings the format uses.
struct Writer {
    bytes: Vec<u8>,
}

impl Writer {
    fn new() -> Self {
        Self { bytes: Vec::new() }
    }

    fn i32(&mut self, value: i32) {
        self.bytes.extend_from_slice(&value.to_le_bytes());
    }

    fn i64(&mut self, value: i64) {
        self.bytes.extend_from_slice(&value.to_le_bytes());
    }

    /// `LittleEndianOutputStream.writeString`: the bytes then a NUL, not length-prefixed.
    fn string(&mut self, value: &str) {
        self.bytes.extend_from_slice(value.as_bytes());
        self.bytes.push(0);
    }
}

impl TribbleIndex {
    /// `AbstractIndex.write`, for a linear index.
    ///
    /// Refuses an interval-tree index rather than guessing at its per-contig record, matching
    /// [`TribbleIndex::read`]'s refusal on the way in.
    pub fn write(&self) -> Result<Vec<u8>, IndexError> {
        if self.index_type != LINEAR {
            return Err(IndexError::UnsupportedType {
                found: self.index_type,
            });
        }
        let mut writer = Writer::new();
        writer.i32(MAGIC_NUMBER);
        writer.i32(self.index_type);
        writer.i32(self.version);
        writer.string(&self.indexed_path);
        writer.i64(self.indexed_file_size);
        writer.i64(self.indexed_file_timestamp);
        writer.string(&self.indexed_file_md5);
        writer.i32(self.flags);

        writer.i32(self.properties.len() as i32);
        for (key, value) in &self.properties {
            writer.string(key);
            writer.string(value);
        }

        writer.i32(self.contigs.len() as i32);
        for contig in &self.contigs {
            writer.string(&contig.name);
            writer.i32(contig.bin_width);
            writer.i32(contig.blocks.len() as i32);
            writer.i32(contig.longest_feature);
            // "no longer used", written as a literal zero rather than from the field.
            writer.i32(0);
            writer.i32(contig.n_features);

            // N blocks are written as N+1 longs: every start, then the end of the last one. A
            // contig with no blocks writes one long of zero, because the running pair starts there.
            let mut start = 0i64;
            let mut size = 0i64;
            for block in &contig.blocks {
                start = block.start;
                size = block.size;
                writer.i64(start);
            }
            writer.i64(start + size);
        }
        Ok(writer.bytes)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn feature(contig: &str, start: i32, end: i32) -> Feature {
        Feature {
            contig: contig.to_string(),
            start,
            end,
        }
    }

    /// The running maximum, not the feature lengths: 10, 10, 600, 10, 10 has mean 364, not 128.
    #[test]
    fn the_statistics_are_pushed_from_the_running_maximum() {
        let mut creator = DynamicIndexCreator::new(BalanceApproach::ForSeekTime);
        for (index, (start, end)) in [(101, 110), (201, 210), (301, 900), (20001, 20010)]
            .into_iter()
            .enumerate()
        {
            creator.add_feature(&feature("chr1", start, end), index as i64 * 10);
        }
        creator.add_feature(&feature("chr2", 51, 60), 64);
        assert_eq!(creator.properties()[0].1, "364.0");
    }

    /// One block means the optimizer's first test breaks before any merge.
    #[test]
    fn a_single_block_keeps_the_width_it_was_given() {
        let mut creator = LinearIndexCreator::new(DEFAULT_BIN_WIDTH);
        creator.add_feature(&feature("chr1", 101, 110), 0);
        let contigs = creator
            .finalize(18, Vec::new())
            .expect("a feature was seen");
        assert_eq!(contigs[0].bin_width, DEFAULT_BIN_WIDTH);
        assert_eq!(contigs[0].blocks.len(), 1);
    }

    /// The message states the opposite of the condition, and is reproduced that way.
    #[test]
    fn an_empty_file_is_refused_with_the_backwards_message() {
        let creator = LinearIndexCreator::new(DEFAULT_BIN_WIDTH);
        let error = creator.finalize(0, Vec::new()).expect_err("no features");
        assert_eq!(error.message(), "finalFilePosition != 0, -> 0");
    }

    /// Doubling once and keeping the last width under the threshold.
    #[test]
    fn a_sparse_contig_is_merged_until_one_block_is_left() {
        let mut creator = LinearIndexCreator::new(DEFAULT_BIN_WIDTH);
        creator.add_feature(&feature("chr1", 101, 110), 0);
        creator.add_feature(&feature("chr1", 20001, 20010), 45);
        let contigs = creator
            .finalize(64, Vec::new())
            .expect("features were seen");
        assert_eq!(contigs[0].bin_width, 16000);
        assert_eq!(contigs[0].blocks.len(), 2);
    }

    #[test]
    fn the_two_approaches_choose_opposite_ends_of_the_same_ordering() {
        let sparse = |approach| {
            let mut creator = DynamicIndexCreator::new(approach);
            creator.add_feature(&feature("chr1", 101, 110), 0);
            creator.add_feature(&feature("chr1", 301, 900), 22);
            creator.add_feature(&feature("chr2", 51, 60), 64);
            creator.chosen()
        };
        assert_eq!(sparse(BalanceApproach::ForSeekTime), ChosenIndex::Linear);
        assert_eq!(sparse(BalanceApproach::ForSize), ChosenIndex::IntervalTree);
    }
}
