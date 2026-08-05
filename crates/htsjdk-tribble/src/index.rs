//! Ported from `htsjdk.tribble.index.AbstractIndex` and `htsjdk.tribble.index.linear.LinearIndex`
//! (htsjdk 4.2.0).
//!
//! The `.idx` file, and the blocks a query resolves to. This is what turns a Feature file into a
//! random-access source rather than a linear read, and it is the item gatk-rs G1.3 named when it
//! closed.
//!
//! # Written against measured bytes, not against a reading of the format
//!
//! The dump came first here, and it found three things this file would otherwise have got wrong —
//! two of them without ever failing loudly:
//!
//!  * **the type identifiers cannot be read out of the Java at all.** `LinearIndex.INDEX_TYPE`
//!    reads a field of the `IndexType` enum whose own constructor is handed
//!    `LinearIndex.INDEX_TYPE`. Measured: `1` for linear, `2` for interval tree;
//!  * **the bin width is per contig**, not a constant. One file, two contigs, 16000 and 8000. A
//!    reader that assumed the creator's default would answer every query wrongly and never fail;
//!  * **the header carries a timestamp**, so two `.idx` files built from the same input at
//!    different moments differ in eight bytes. Nothing here depends on it, and the conformance
//!    suite masks it rather than pretending otherwise.
//!
//! # The layout
//!
//! Little-endian throughout, and strings are **NUL-terminated rather than length-prefixed**:
//!
//! ```text
//! magic i32 = 1480870228 ("TIDX")   type i32   version i32
//! path  NUL-terminated               size i64   timestamp i64
//! md5   NUL-terminated               flags i32  nProperties i32, then that many key/value pairs
//! nContigs i32, then per contig:
//!   name NUL   binWidth i32   nBins i32   longestFeature i32   unused i32   nFeatures i32
//!   then nBins + 1 longs
//! ```
//!
//! That last line is the one to be careful with: **N blocks are stored as N+1 positions**, and a
//! block's size is the difference to the next. The final position is the end of the last block.
//!
//! # The query is three decisions, not one lookup
//!
//! ```java
//! final int adjustedPosition = Math.max(start - longestFeature, 0);
//! final int startBinNumber = adjustedPosition / binWidth;
//! if (startBinNumber >= blocks.size()) return Collections.emptyList();
//! final int endBinNumber = Math.min((end - 1) / binWidth, blocks.size() - 1);
//! ```
//!
//!  * **the longest feature is subtracted from the start**, so a query reaches features that begin
//!    before the interval and extend into it. That is the only reason the index records it;
//!  * **the answer is one merged block or none**, never a list: linear blocks are adjacent by
//!    definition, so everything from the first bin to the last is contiguous on disk;
//!  * **an unknown contig throws.** `IllegalArgumentException`, not an empty list, so a caller that
//!    mistypes a contig learns rather than silently reading nothing.

/// `AbstractIndex.MAGIC_NUMBER`, which is the bytes `TIDX`.
pub const MAGIC_NUMBER: i32 = 1_480_870_228;

/// `AbstractIndex.VERSION`.
pub const VERSION: i32 = 3;

/// The type identifiers, measured rather than read: the Java defines them through a circular
/// reference between `LinearIndex.INDEX_TYPE` and the `IndexType` enum.
pub const LINEAR: i32 = 1;
pub const INTERVAL_TREE: i32 = 2;

/// What can go wrong reading or querying an index.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IndexError {
    /// `TribbleException`: the file does not start with `TIDX`.
    BadMagic { found: i32 },
    /// A type this reader does not parse. The interval-tree chromosome record is a different shape
    /// and is refused rather than mis-read.
    UnsupportedType { found: i32 },
    /// The stream ended mid-record.
    Truncated,
    /// `IllegalArgumentException` from a query naming a contig the index does not hold.
    UnknownContig { contig: String },
}

impl IndexError {
    pub fn class(&self) -> &'static str {
        match self {
            IndexError::UnknownContig { .. } => "java.lang.IllegalArgumentException",
            _ => "htsjdk.tribble.TribbleException",
        }
    }
}

/// One contiguous run of the indexed file.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Block {
    pub start: i64,
    pub size: i64,
}

/// One interval of an `IntervalTreeIndex.ChrIndex`, with the block it points at.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Interval {
    pub start: i32,
    pub end: i32,
    pub block: Block,
}

/// `LinearIndex.ChrIndex`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChrIndex {
    pub name: String,
    /// Chosen per contig by the creator, so it differs between contigs of one file.
    pub bin_width: i32,
    /// The longest feature on this contig, which the query subtracts from its start.
    pub longest_feature: i32,
    /// `OLD_V3_INDEX` is derived from this being positive; it is written as zero and otherwise
    /// unused, and is kept because the layout has a slot for it.
    pub unused: i32,
    pub n_features: i32,
    pub blocks: Vec<Block>,
}

/// `IntervalTreeIndex.ChrIndex`, whose per-contig record is a different shape from the linear one:
/// a flat list of intervals rather than a bin width and a run of positions.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IntervalChrIndex {
    pub name: String,
    /// In the order the file holds them, which is the tree's own traversal order rather than
    /// anything sorted.
    pub intervals: Vec<Interval>,
}

/// The gap below which two blocks are merged into one read.
///
/// `block.getStartPosition() < lastBlock.getEndPosition() + 1000`. Not a tunable: the constant is
/// written into the method.
pub const CONSOLIDATION_GAP: i64 = 1000;

/// A parsed `.idx`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TribbleIndex {
    pub index_type: i32,
    pub version: i32,
    /// The URI of the file this index was built from, as the header records it.
    pub indexed_path: String,
    pub indexed_file_size: i64,
    /// The source file's modification time. Nothing reads it here, and it is why two indexes over
    /// the same input differ byte for byte.
    pub indexed_file_timestamp: i64,
    pub indexed_file_md5: String,
    pub flags: i32,
    pub properties: Vec<(String, String)>,
    /// In file order, which is the order the creator wrote the contigs in. Empty for an
    /// interval-tree index.
    pub contigs: Vec<ChrIndex>,
    /// The interval-tree contigs, likewise in file order. Empty for a linear index.
    pub interval_contigs: Vec<IntervalChrIndex>,
}

/// A cursor over the little-endian stream, with the NUL-terminated strings the format uses.
struct Reader<'a> {
    bytes: &'a [u8],
    at: usize,
}

impl<'a> Reader<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, at: 0 }
    }

    fn i32(&mut self) -> Result<i32, IndexError> {
        let end = self.at + 4;
        let slice = self.bytes.get(self.at..end).ok_or(IndexError::Truncated)?;
        self.at = end;
        Ok(i32::from_le_bytes(slice.try_into().expect("four bytes")))
    }

    fn i64(&mut self) -> Result<i64, IndexError> {
        let end = self.at + 8;
        let slice = self.bytes.get(self.at..end).ok_or(IndexError::Truncated)?;
        self.at = end;
        Ok(i64::from_le_bytes(slice.try_into().expect("eight bytes")))
    }

    /// `LittleEndianInputStream.readString`: bytes up to a NUL, which is written by
    /// `writeString` as `writeBytes(s)` followed by a zero. Not length-prefixed.
    fn string(&mut self) -> Result<String, IndexError> {
        let start = self.at;
        while *self.bytes.get(self.at).ok_or(IndexError::Truncated)? != 0 {
            self.at += 1;
        }
        let text = String::from_utf8_lossy(&self.bytes[start..self.at]).into_owned();
        self.at += 1;
        Ok(text)
    }
}

impl TribbleIndex {
    /// Parse an `.idx`.
    pub fn read(bytes: &[u8]) -> Result<Self, IndexError> {
        let mut reader = Reader::new(bytes);
        let magic = reader.i32()?;
        if magic != MAGIC_NUMBER {
            return Err(IndexError::BadMagic { found: magic });
        }
        let index_type = reader.i32()?;
        if index_type != LINEAR && index_type != INTERVAL_TREE {
            return Err(IndexError::UnsupportedType { found: index_type });
        }
        let version = reader.i32()?;
        let indexed_path = reader.string()?;
        let indexed_file_size = reader.i64()?;
        let indexed_file_timestamp = reader.i64()?;
        let indexed_file_md5 = reader.string()?;
        let flags = reader.i32()?;

        let mut properties = Vec::new();
        let mut remaining = reader.i32()?;
        while remaining > 0 {
            let key = reader.string()?;
            let value = reader.string()?;
            properties.push((key, value));
            remaining -= 1;
        }

        let mut contigs = Vec::new();
        let mut interval_contigs = Vec::new();
        let mut remaining = reader.i32()?;
        while remaining > 0 && index_type == INTERVAL_TREE {
            // `IntervalTreeIndex.ChrIndex.read`: a name, a count, then that many
            // (start, end, position, size) records. Sizes are stored, not derived, which is the
            // opposite of the linear layout and the thing to be careful about when reading both.
            let name = reader.string()?;
            let mut count = reader.i32()?;
            let mut intervals = Vec::new();
            while count > 0 {
                let start = reader.i32()?;
                let end = reader.i32()?;
                let position = reader.i64()?;
                let size = i64::from(reader.i32()?);
                intervals.push(Interval {
                    start,
                    end,
                    block: Block {
                        start: position,
                        size,
                    },
                });
                count -= 1;
            }
            interval_contigs.push(IntervalChrIndex { name, intervals });
            remaining -= 1;
        }
        while remaining > 0 {
            let name = reader.string()?;
            let bin_width = reader.i32()?;
            let n_bins = reader.i32()?;
            let longest_feature = reader.i32()?;
            let unused = reader.i32()?;
            let n_features = reader.i32()?;

            // N blocks, N+1 positions: a block's size is the distance to the next one.
            let mut position = reader.i64()?;
            let mut blocks = Vec::with_capacity(n_bins.max(0) as usize);
            for _ in 0..n_bins {
                let next = reader.i64()?;
                blocks.push(Block {
                    start: position,
                    size: next - position,
                });
                position = next;
            }

            contigs.push(ChrIndex {
                name,
                bin_width,
                longest_feature,
                unused,
                n_features,
                blocks,
            });
            remaining -= 1;
        }

        Ok(Self {
            index_type,
            version,
            indexed_path,
            indexed_file_size,
            indexed_file_timestamp,
            indexed_file_md5,
            flags,
            properties,
            contigs,
            interval_contigs,
        })
    }

    /// `AbstractIndex.getBlocks(chr, start, end)`.
    ///
    /// An unknown contig is refused rather than answered with nothing, which is htsjdk's choice
    /// and a useful one: a mistyped contig otherwise reads as a region with no features.
    pub fn blocks(&self, contig: &str, start: i32, end: i32) -> Result<Vec<Block>, IndexError> {
        if let Some(chr) = self.contigs.iter().find(|c| c.name == contig) {
            return Ok(chr.blocks_for(start, end));
        }
        if let Some(chr) = self.interval_contigs.iter().find(|c| c.name == contig) {
            return Ok(chr.blocks_for(start, end));
        }
        Err(IndexError::UnknownContig {
            contig: contig.to_string(),
        })
    }

    /// `getSequenceNames()`, in file order.
    pub fn sequence_names(&self) -> Vec<&str> {
        self.contigs
            .iter()
            .map(|c| c.name.as_str())
            .chain(self.interval_contigs.iter().map(|c| c.name.as_str()))
            .collect()
    }
}

impl IntervalChrIndex {
    /// `IntervalTreeIndex.ChrIndex.getBlocks(start, end)`.
    ///
    /// Three steps, and the middle one is where a port has to say what it cannot reproduce.
    ///
    /// **Overlap** is the tree's `findOverlapping`, which is inclusive at both ends.
    ///
    /// **The sort** is by start position, through a comparator htsjdk itself calls "a little
    /// cryptic":
    ///
    /// ```java
    /// return b1.getStartPosition() - b2.getStartPosition() < 1 ? -1
    ///      : (b1.getStartPosition() - b2.getStartPosition() > 1 ? 1 : 0);
    /// ```
    ///
    /// It is **not a consistent comparator**. `compare(a, a)` is `-1` rather than `0`, and two
    /// blocks one byte apart compare *equal*. For any pair whose starts differ by two or more it
    /// is an ordinary ascending sort, which is what this does.
    ///
    /// The other two cases cannot arise, and that is a measurement rather than an absence of one.
    /// An index was built over a 3.7 MB file and read back **before any query touched it**: seven
    /// intervals, and the closest two block starts once sorted were **557,062 bytes** apart. The
    /// creator emits one interval per run of features and the runs partition the file, so two
    /// blocks starting within one byte of each other would need two intervals pointing at the same
    /// offset — something the creator never produces. The sort is therefore an ordinary ascending
    /// one on every index the creator can make, and the port reproduces it exactly.
    ///
    /// Worth knowing why the sort is there at all: the intervals come out of the tree in **tree
    /// order, not file order**. The same measurement showed consecutive stored intervals with
    /// block starts *decreasing* by 559,200 bytes.
    ///
    /// **The consolidation** merges anything starting within [`CONSOLIDATION_GAP`] of the previous
    /// block's end, so the answer is a list of reads rather than one block per interval. In the
    /// reference this is done by **mutating the stored block in place**, so a query widens the
    /// index it queried; this returns new blocks instead.
    ///
    /// That difference cannot reach an answer, and again the reason is measured. Sorted, the blocks
    /// **partition the file contiguously**: seven blocks over 3.7 MB with starts 557,062 apart and
    /// sizes to match, so the gap between one block's end and the next block's start is nothing
    /// like the 1000-byte threshold. Every query in a contiguous range therefore consolidates to a
    /// single block, and a widened block already covers what it merged, so a later and wider query
    /// re-merges the same intervals to the same union. Probed directly as well — a narrow query,
    /// then a wide one, against a fresh index asked the wide one straight away, on both a 9 KB and
    /// a 3.7 MB corpus — and all of them agreed.
    pub fn blocks_for(&self, start: i32, end: i32) -> Vec<Block> {
        let mut blocks: Vec<Block> = self
            .intervals
            .iter()
            .filter(|interval| interval.start <= end && start <= interval.end)
            .map(|interval| interval.block)
            .collect();
        if blocks.is_empty() {
            return Vec::new();
        }
        blocks.sort_by_key(|block| block.start);

        let mut consolidated: Vec<Block> = Vec::with_capacity(blocks.len());
        consolidated.push(blocks[0]);
        for block in &blocks[1..] {
            let last = consolidated.last_mut().expect("pushed one");
            if block.start < last.start + last.size + CONSOLIDATION_GAP {
                // `lastBlock.setEndPosition(block.getEndPosition())`: the end moves, the start
                // does not, so the size is the distance from the original start.
                let end_position = block.start + block.size;
                last.size = end_position - last.start;
            } else {
                consolidated.push(*block);
            }
        }
        consolidated
    }
}

impl ChrIndex {
    /// `LinearIndex.ChrIndex.getBlocks(start, end)`.
    pub fn blocks_for(&self, start: i32, end: i32) -> Vec<Block> {
        if self.blocks.is_empty() {
            return Vec::new();
        }
        // The back-off: a feature beginning before the interval can still reach into it, and the
        // longest one on this contig bounds how far back that can be.
        let adjusted = (start - self.longest_feature).max(0);
        let start_bin = (adjusted / self.bin_width) as usize;
        if start_bin >= self.blocks.len() {
            return Vec::new();
        }
        let end_bin = (((end - 1) / self.bin_width) as usize).min(self.blocks.len() - 1);

        // Blocks are adjacent by definition, so the answer is the single run from the first bin's
        // start to the last bin's end rather than a list of bins.
        let start_position = self.blocks[start_bin].start;
        let end_position = self.blocks[end_bin].start + self.blocks[end_bin].size;
        let size = end_position - start_position;
        if size == 0 {
            return Vec::new();
        }
        vec![Block {
            start: start_position,
            size,
        }]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn chr(name: &str, bin_width: i32, longest: i32, positions: &[i64]) -> ChrIndex {
        let blocks = positions
            .windows(2)
            .map(|pair| Block {
                start: pair[0],
                size: pair[1] - pair[0],
            })
            .collect();
        ChrIndex {
            name: name.to_string(),
            bin_width,
            longest_feature: longest,
            unused: 0,
            n_features: 0,
            blocks,
        }
    }

    /// The numbers are the oracle's, from `chr linear-bed chr1 16000 2 600 0 4 0,59,80`.
    #[test]
    fn the_answer_is_one_merged_block() {
        let chr1 = chr("chr1", 16000, 600, &[0, 59, 80]);
        assert_eq!(
            chr1.blocks_for(100, 120),
            vec![Block { start: 0, size: 59 }]
        );
        // Both bins, merged into one run rather than returned as two.
        assert_eq!(
            chr1.blocks_for(100, 32001),
            vec![Block { start: 0, size: 80 }]
        );
    }

    /// The back-off is why the index records the longest feature at all.
    #[test]
    fn the_longest_feature_is_subtracted_from_the_start() {
        // One bin per 100 bases, and a 600-base longest feature: a query at 700 reaches back to
        // bin 1 rather than starting at bin 7.
        let chr1 = chr("chr1", 100, 600, &[0, 10, 20, 30, 40, 50, 60, 70, 80]);
        assert_eq!(
            chr1.blocks_for(700, 710),
            vec![Block {
                start: 10,
                size: 70
            }]
        );
        // Without the back-off it would have started at bin 7.
        let no_backoff = chr("chr1", 100, 0, &[0, 10, 20, 30, 40, 50, 60, 70, 80]);
        assert_eq!(
            no_backoff.blocks_for(700, 710),
            vec![Block {
                start: 70,
                size: 10
            }]
        );
    }

    #[test]
    fn a_query_past_the_last_bin_is_empty_and_an_unknown_contig_is_refused() {
        let index = TribbleIndex {
            index_type: LINEAR,
            version: VERSION,
            indexed_path: String::new(),
            indexed_file_size: 0,
            indexed_file_timestamp: 0,
            indexed_file_md5: String::new(),
            flags: 0,
            properties: Vec::new(),
            contigs: vec![chr("chr1", 16000, 600, &[0, 59, 80])],
            interval_contigs: Vec::new(),
        };
        assert_eq!(index.blocks("chr1", 1_000_000, 1_000_010), Ok(Vec::new()));
        let error = index.blocks("chrX", 1, 10).expect_err("refused");
        assert_eq!(
            error,
            IndexError::UnknownContig {
                contig: "chrX".to_string()
            }
        );
        assert_eq!(error.class(), "java.lang.IllegalArgumentException");
    }

    /// A zero-width merge is the second way to get nothing, and it is not the same as being off
    /// the end.
    #[test]
    fn a_zero_sized_merge_is_empty_too() {
        let chr1 = chr("chr1", 100, 0, &[40, 40, 40]);
        assert!(chr1.blocks_for(1, 10).is_empty());
    }

    #[test]
    fn a_file_that_is_not_an_index_is_refused_by_its_magic() {
        let error = TribbleIndex::read(b"not an index at all").expect_err("refused");
        assert!(matches!(error, IndexError::BadMagic { .. }));
        // A type that is neither linear nor interval tree is refused rather than guessed at.
        let mut bytes = MAGIC_NUMBER.to_le_bytes().to_vec();
        bytes.extend(99i32.to_le_bytes());
        assert_eq!(
            TribbleIndex::read(&bytes),
            Err(IndexError::UnsupportedType { found: 99 })
        );
    }

    fn interval(start: i32, end: i32, block_start: i64, size: i64) -> Interval {
        Interval {
            start,
            end,
            block: Block {
                start: block_start,
                size,
            },
        }
    }

    /// The consolidation is what makes an interval-tree answer a list of *reads* rather than one
    /// block per interval, and the gap is a constant written into the method.
    #[test]
    fn intervals_within_a_thousand_bytes_become_one_read() {
        let chr = IntervalChrIndex {
            name: "chr1".to_string(),
            intervals: vec![
                interval(100, 200, 0, 50),
                // Starts 40 bytes after the first one ends: inside the gap, so merged.
                interval(300, 400, 90, 10),
                // Starts 1000 past the previous end exactly, which is NOT less than the gap.
                interval(500, 600, 1100, 20),
            ],
        };
        assert_eq!(
            chr.blocks_for(1, 1000),
            vec![
                Block {
                    start: 0,
                    size: 100
                },
                Block {
                    start: 1100,
                    size: 20
                }
            ]
        );
    }

    /// Overlap is inclusive at both ends, so a query touching an interval's first or last base
    /// still finds it.
    #[test]
    fn overlap_is_inclusive_at_both_ends() {
        let chr = IntervalChrIndex {
            name: "chr1".to_string(),
            intervals: vec![interval(100, 200, 0, 50)],
        };
        assert_eq!(chr.blocks_for(200, 300).len(), 1, "touching the end");
        assert_eq!(chr.blocks_for(1, 100).len(), 1, "touching the start");
        assert!(chr.blocks_for(201, 300).is_empty(), "just past the end");
        assert!(chr.blocks_for(1, 99).is_empty(), "just before the start");
    }
}
