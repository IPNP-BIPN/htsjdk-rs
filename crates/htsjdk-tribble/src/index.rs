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
    /// In file order, which is the order the creator wrote the contigs in.
    pub contigs: Vec<ChrIndex>,
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
        if index_type != LINEAR {
            // The interval-tree chromosome record is a different shape. Refused rather than
            // mis-parsed, because a reader that guessed would produce plausible nonsense.
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
        let mut remaining = reader.i32()?;
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
        })
    }

    /// `AbstractIndex.getBlocks(chr, start, end)`.
    ///
    /// An unknown contig is refused rather than answered with nothing, which is htsjdk's choice
    /// and a useful one: a mistyped contig otherwise reads as a region with no features.
    pub fn blocks(&self, contig: &str, start: i32, end: i32) -> Result<Vec<Block>, IndexError> {
        let chr = self
            .contigs
            .iter()
            .find(|c| c.name == contig)
            .ok_or_else(|| IndexError::UnknownContig {
                contig: contig.to_string(),
            })?;
        Ok(chr.blocks_for(start, end))
    }

    /// `getSequenceNames()`, in file order.
    pub fn sequence_names(&self) -> Vec<&str> {
        self.contigs.iter().map(|c| c.name.as_str()).collect()
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
        // An interval-tree index is refused rather than mis-parsed.
        let mut bytes = MAGIC_NUMBER.to_le_bytes().to_vec();
        bytes.extend(INTERVAL_TREE.to_le_bytes());
        assert_eq!(
            TribbleIndex::read(&bytes),
            Err(IndexError::UnsupportedType {
                found: INTERVAL_TREE
            })
        );
    }
}
