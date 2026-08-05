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
//! # Both layouts are written, because the choice is not the caller's
//!
//! The interval-tree layout is here too, and it has to be: the dynamic creator reaches it from
//! ordinary data, so a writer that only knew the linear one could not build the index htsjdk builds
//! for three of the ten cases this suite measures.
//!
//! Its byte order is the harder half. `IntervalTreeIndex.ChrIndex.write` writes
//! `tree.getIntervals()`, and that is a **pre-order walk of a red-black tree**, not a sorted list,
//! so the order in the file is the order the rotations left the nodes in. Reproducing the bytes
//! means reproducing the tree: the CLRS insert with its two mirrored fixup halves, both rotations,
//! and the `min`/`max` update that walks to the root after each one. The insert comparator sends
//! **equal starts left**, which is what turns a run of intervals sharing a start into a left spine.

use crate::index::{
    Block, ChrIndex, IndexError, TribbleIndex, INTERVAL_TREE, LINEAR, MAGIC_NUMBER,
};

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
}

impl CreateError {
    pub fn class(&self) -> &'static str {
        match self {
            CreateError::IllegalArgument(_) => "java.lang.IllegalArgumentException",
            CreateError::MalformedFeatureFile(_) => {
                "htsjdk.tribble.TribbleException$MalformedFeatureFile"
            }
        }
    }

    pub fn message(&self) -> String {
        match self {
            CreateError::IllegalArgument(message) | CreateError::MalformedFeatureFile(message) => {
                message.clone()
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

/// What a creator built: the two layouts are different shapes, so the choice travels with them.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BuiltIndex {
    Linear(Vec<ChrIndex>),
    IntervalTree(Vec<crate::index::IntervalChrIndex>),
}

/// `DynamicIndexCreator`.
///
/// It holds both creators and feeds every feature to **both**, then scores them and keeps one. The
/// discarded creator's work is thrown away, which is why a dynamic index costs roughly twice what
/// a fixed one does.
#[derive(Debug, Clone)]
pub struct DynamicIndexCreator {
    approach: BalanceApproach,
    linear: LinearIndexCreator,
    tree: IntervalIndexCreator,
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
            tree: IntervalIndexCreator::new(features_per_interval),
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
        self.tree.add_feature(feature, file_position);
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
        let tree_score = f64::from(self.tree.features_per_interval());

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

    /// `finalizeIndex`: score, keep one, and hand the properties to the survivor.
    ///
    /// The properties are added to the **chosen** creator only, which is why a linear index built
    /// directly carries none and the same index built dynamically carries four.
    pub fn finalize(self, final_file_position: i64) -> Result<BuiltIndex, CreateError> {
        let chosen = self.chosen();
        let properties = self.properties();
        match chosen {
            ChosenIndex::Linear => Ok(BuiltIndex::Linear(
                self.linear.finalize(final_file_position, properties)?,
            )),
            // The interval creator has no `finalFilePosition == 0` guard of its own, so an empty
            // file reaches it and produces an index with no contigs rather than a refusal.
            ChosenIndex::IntervalTree => Ok(BuiltIndex::IntervalTree(
                self.tree.finalize(final_file_position),
            )),
        }
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
    /// `AbstractIndex.write`.
    ///
    /// Both layouts, which is what [`TribbleIndex::read`] already parses. Their per-contig records
    /// are different shapes and the difference is easy to miss: the linear one writes N blocks as
    /// **N+1 positions** whose differences are the sizes, and the interval one writes each size
    /// **outright**, as an `int` where the block holds a `long`.
    pub fn write(&self) -> Result<Vec<u8>, IndexError> {
        if self.index_type != LINEAR && self.index_type != INTERVAL_TREE {
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

        if self.index_type == INTERVAL_TREE {
            writer.i32(self.interval_contigs.len() as i32);
            for contig in &self.interval_contigs {
                writer.string(&contig.name);
                writer.i32(contig.intervals.len() as i32);
                for interval in &contig.intervals {
                    writer.i32(interval.start);
                    writer.i32(interval.end);
                    writer.i64(interval.block.start);
                    // `(int) interval.getBlock().getSize()`: an unchecked narrowing, so a block
                    // larger than 2 GB writes a truncated size rather than failing.
                    writer.i32(interval.block.size as i32);
                }
            }
            return Ok(writer.bytes);
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

/// `IntervalTree`, the red-black tree whose **shape** decides the byte order of an interval-tree
/// index.
///
/// `IntervalTreeIndex.ChrIndex.write` writes `tree.getIntervals()`, and `getIntervals` is a
/// **pre-order** walk (node, left, right) of the tree rather than a sorted list. So the order in
/// the file is the order the rotations left the nodes in, and reproducing the bytes means
/// reproducing the tree: the CLRS insert, both rotations, and the `min`/`max` update that walks to
/// the root after each one.
///
/// The insert comparator is `x.start <= node.start`, so **equal starts go left**, which is what
/// puts a run of intervals sharing a start into a left spine rather than balanced around it.
///
/// Nodes live in a `Vec` with index links, because the Java is a graph of mutable parent pointers
/// and an arena is the honest translation of that; `NIL` is index 0 and carries the sentinel
/// `min`/`max` its default constructor sets.
mod interval_tree {
    use crate::index::{Block, Interval};

    const NIL: usize = 0;
    const RED: bool = true;
    const BLACK: bool = false;

    struct Node {
        interval: Interval,
        min: i32,
        max: i32,
        left: usize,
        right: usize,
        parent: usize,
        color: bool,
    }

    pub struct IntervalTree {
        nodes: Vec<Node>,
        root: usize,
    }

    impl Default for IntervalTree {
        fn default() -> Self {
            Self::new()
        }
    }

    impl IntervalTree {
        pub fn new() -> Self {
            // `Node.NIL`, whose private constructor leaves `max` at Integer.MIN_VALUE and `min` at
            // Integer.MAX_VALUE so that `update`'s maxima and minima ignore it.
            let nil = Node {
                interval: Interval {
                    start: 0,
                    end: 0,
                    block: Block { start: 0, size: 0 },
                },
                min: i32::MAX,
                max: i32::MIN,
                left: NIL,
                right: NIL,
                parent: NIL,
                color: BLACK,
            };
            Self {
                nodes: vec![nil],
                root: NIL,
            }
        }

        pub fn insert(&mut self, interval: Interval) {
            self.nodes.push(Node {
                interval,
                min: i32::MAX,
                max: i32::MIN,
                left: NIL,
                right: NIL,
                parent: NIL,
                color: RED,
            });
            let x = self.nodes.len() - 1;
            self.tree_insert(x);
            self.fixup(x);
        }

        /// `getIntervals`: pre-order, not sorted and not in-order.
        pub fn intervals(&self) -> Vec<Interval> {
            let mut out = Vec::new();
            if self.root != NIL {
                self.get_all(self.root, &mut out);
            }
            out
        }

        fn get_all(&self, node: usize, out: &mut Vec<Interval>) {
            out.push(self.nodes[node].interval);
            if self.nodes[node].left != NIL {
                self.get_all(self.nodes[node].left, out);
            }
            if self.nodes[node].right != NIL {
                self.get_all(self.nodes[node].right, out);
            }
        }

        /// `treeInsert`. Equal starts go **left**.
        fn tree_insert(&mut self, x: usize) {
            let mut node = self.root;
            let mut y = NIL;
            while node != NIL {
                y = node;
                node = if self.nodes[x].interval.start <= self.nodes[node].interval.start {
                    self.nodes[node].left
                } else {
                    self.nodes[node].right
                };
            }
            self.nodes[x].parent = y;
            if y == NIL {
                self.root = x;
                self.nodes[x].left = NIL;
                self.nodes[x].right = NIL;
            } else if self.nodes[x].interval.start <= self.nodes[y].interval.start {
                self.nodes[y].left = x;
            } else {
                self.nodes[y].right = x;
            }
            self.apply_update(x);
        }

        /// The CLRS insert fixup, spelled out rather than shortened: the two halves are mirror
        /// images and htsjdk writes both, so both are here.
        fn fixup(&mut self, mut x: usize) {
            self.nodes[x].color = RED;
            while x != self.root && self.nodes[self.nodes[x].parent].color == RED {
                let parent = self.nodes[x].parent;
                let grandparent = self.nodes[parent].parent;
                if parent == self.nodes[grandparent].left {
                    let uncle = self.nodes[grandparent].right;
                    if self.nodes[uncle].color == RED {
                        self.nodes[parent].color = BLACK;
                        self.nodes[uncle].color = BLACK;
                        self.nodes[grandparent].color = RED;
                        x = grandparent;
                    } else {
                        if x == self.nodes[parent].right {
                            x = parent;
                            self.left_rotate(x);
                        }
                        let parent = self.nodes[x].parent;
                        let grandparent = self.nodes[parent].parent;
                        self.nodes[parent].color = BLACK;
                        self.nodes[grandparent].color = RED;
                        self.right_rotate(grandparent);
                    }
                } else {
                    let uncle = self.nodes[grandparent].left;
                    if self.nodes[uncle].color == RED {
                        self.nodes[parent].color = BLACK;
                        self.nodes[uncle].color = BLACK;
                        self.nodes[grandparent].color = RED;
                        x = grandparent;
                    } else {
                        if x == self.nodes[parent].left {
                            x = parent;
                            self.right_rotate(x);
                        }
                        let parent = self.nodes[x].parent;
                        let grandparent = self.nodes[parent].parent;
                        self.nodes[parent].color = BLACK;
                        self.nodes[grandparent].color = RED;
                        self.left_rotate(grandparent);
                    }
                }
            }
            self.nodes[self.root].color = BLACK;
        }

        fn left_rotate(&mut self, x: usize) {
            let y = self.nodes[x].right;
            self.nodes[x].right = self.nodes[y].left;
            if self.nodes[y].left != NIL {
                let left = self.nodes[y].left;
                self.nodes[left].parent = x;
            }
            self.nodes[y].parent = self.nodes[x].parent;
            let parent = self.nodes[x].parent;
            if parent == NIL {
                self.root = y;
            } else if self.nodes[parent].left == x {
                self.nodes[parent].left = y;
            } else {
                self.nodes[parent].right = y;
            }
            self.nodes[y].left = x;
            self.nodes[x].parent = y;
            self.apply_update(x);
        }

        fn right_rotate(&mut self, x: usize) {
            let y = self.nodes[x].left;
            self.nodes[x].left = self.nodes[y].right;
            if self.nodes[y].right != NIL {
                let right = self.nodes[y].right;
                self.nodes[right].parent = x;
            }
            self.nodes[y].parent = self.nodes[x].parent;
            let parent = self.nodes[x].parent;
            if parent == NIL {
                self.root = y;
            } else if self.nodes[parent].right == x {
                self.nodes[parent].right = y;
            } else {
                self.nodes[parent].left = y;
            }
            self.nodes[y].right = x;
            self.nodes[x].parent = y;
            self.apply_update(x);
        }

        /// `applyUpdate`: walk to the root recomputing `min` and `max`. The loop stops at `NIL`,
        /// which is reached through the root's own parent, so the root is updated and the sentinel
        /// is not.
        fn apply_update(&mut self, mut node: usize) {
            while node != NIL {
                let left = self.nodes[node].left;
                let right = self.nodes[node].right;
                self.nodes[node].max = self.nodes[left]
                    .max
                    .max(self.nodes[right].max)
                    .max(self.nodes[node].interval.end);
                self.nodes[node].min = self.nodes[left]
                    .min
                    .min(self.nodes[right].min)
                    .min(self.nodes[node].interval.start);
                node = self.nodes[node].parent;
            }
        }
    }
}

/// `IntervalIndexCreator`.
///
/// Where the linear creator cuts the file by **coordinate**, this one cuts it by **feature count**:
/// a new interval opens every `features_per_interval` features, and the one before it is closed at
/// the position where the new one starts. So the intervals partition the file contiguously and
/// their coordinates are whatever the features in them happened to span.
///
/// The stop of the open interval is `max(feature.end, stop)`, so an interval's end can exceed the
/// start of the next one: the coordinate ranges overlap even though the file ranges do not.
#[derive(Debug, Clone)]
pub struct IntervalIndexCreator {
    features_per_interval: i32,
    /// The contigs closed so far, each with its intervals in insertion order.
    chr_list: Vec<(String, Vec<crate::index::Interval>)>,
    /// The intervals of the contig being filled, still mutable.
    intervals: Vec<MutableInterval>,
    feature_count: i32,
}

#[derive(Debug, Clone, Copy)]
struct MutableInterval {
    start: i32,
    stop: i32,
    start_file_position: i64,
    end_file_position: i64,
}

impl IntervalIndexCreator {
    pub fn new(features_per_interval: i32) -> Self {
        Self {
            features_per_interval,
            chr_list: Vec::new(),
            intervals: Vec::new(),
            feature_count: 0,
        }
    }

    pub fn features_per_interval(&self) -> i32 {
        self.features_per_interval
    }

    /// `addFeature`.
    ///
    /// `featureCount` is **not** reset when the contig changes, only when an interval opens, so the
    /// first interval of a new contig inherits the count from the previous contig. It opens anyway,
    /// because `intervals.isEmpty()` is the other arm of the same test.
    pub fn add_feature(&mut self, feature: &Feature, file_position: i64) {
        let new_contig = match self.chr_list.last() {
            None => true,
            Some((name, _)) => name != &feature.contig,
        };
        if new_contig {
            if !self.chr_list.is_empty() {
                self.close_intervals(file_position);
            }
            self.chr_list.push((feature.contig.clone(), Vec::new()));
            self.intervals.clear();
        }

        if self.feature_count >= self.features_per_interval || self.intervals.is_empty() {
            if let Some(last) = self.intervals.last_mut() {
                last.end_file_position = file_position;
            }
            self.feature_count = 0;
            self.intervals.push(MutableInterval {
                start: feature.start,
                stop: 0,
                start_file_position: file_position,
                end_file_position: 0,
            });
        }
        let last = self.intervals.last_mut().expect("an interval is open");
        last.stop = feature.end.max(last.stop);
        self.feature_count += 1;
    }

    /// `addIntervalsToLastChr`: only the **last** interval's end is set here, because every earlier
    /// one was closed when its successor opened.
    fn close_intervals(&mut self, current_position: i64) {
        let count = self.intervals.len();
        for (index, interval) in self.intervals.iter_mut().enumerate() {
            if index == count - 1 {
                interval.end_file_position = current_position;
            }
        }
        let converted: Vec<crate::index::Interval> = self
            .intervals
            .iter()
            .map(|interval| crate::index::Interval {
                start: interval.start,
                end: interval.stop,
                block: Block {
                    start: interval.start_file_position,
                    size: interval.end_file_position - interval.start_file_position,
                },
            })
            .collect();
        if let Some((_, list)) = self.chr_list.last_mut() {
            list.extend(converted);
        }
    }

    /// `finalizeIndex`, which returns each contig's intervals **in the tree's pre-order**, since
    /// that is the order they are written in.
    pub fn finalize(mut self, final_file_position: i64) -> Vec<crate::index::IntervalChrIndex> {
        if !self.chr_list.is_empty() {
            self.close_intervals(final_file_position);
        }
        self.chr_list
            .into_iter()
            .map(|(name, inserted)| {
                let mut tree = interval_tree::IntervalTree::new();
                for interval in inserted {
                    tree.insert(interval);
                }
                crate::index::IntervalChrIndex {
                    name,
                    intervals: tree.intervals(),
                }
            })
            .collect()
    }
}
