//! The BAI index.
//!
//! Ported from `htsjdk.samtools.BinningIndexBuilder`, `Bin.addChunk`, `BAMIndexMetaData`,
//! `BAMIndexer` and `BinaryBAMIndexWriter`.
//!
//! The `.bai` is a pure side file: nothing in it changes what a BAM contains, and a wrong index
//! degrades to a slow query rather than to a visible error. That makes it the part of the
//! format where a divergence is least likely to be noticed and most likely to persist.
//!
//! Three of its choices are htsjdk's rather than the format's:
//!
//! - **The linear index is back-filled.** A 16 kb window with no read is written as the last
//!   non-empty window's offset, not as zero and not as "absent". htsjdk's own comment says this
//!   is unnecessary and that it does it because samtools does.
//! - **Chunks coalesce across adjacent blocks**, where "adjacent" means the block address plus
//!   one. Block addresses are byte offsets, so that second case essentially never fires; the
//!   rule is really "same block".
//! - **A pseudo-bin numbered 37450 carries statistics**, and it is counted in `n_bin` while
//!   being excluded from the loop that writes the real bins.

use crate::bin::MAX_BINS;

/// `LinearIndex.BAM_LIDX_SHIFT`: the linear index window is 16 kb.
pub const BAM_LIDX_SHIFT: u32 = 14;

/// `GenomicIndexUtil.MAX_LINEAR_INDEX_SIZE`.
pub const MAX_LINEAR_INDEX_SIZE: usize = 32770;

/// `BinningIndexBuilder.UNINITIALIZED_WINDOW`.
///
/// Not zero, and the distinction is load-bearing: 0 is a legitimate virtual file pointer, so a
/// zero-initialised array cannot tell "no read here" from "a read at the very start".
pub const UNINITIALIZED_WINDOW: i64 = -1;

/// `BAMFileConstants.BAM_INDEX_MAGIC`.
pub const BAM_INDEX_MAGIC: [u8; 4] = *b"BAI\x01";

/// `Chunk`: a half-open range of virtual file pointers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Chunk {
    pub start: u64,
    pub end: u64,
}

/// `Bin`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Bin {
    pub bin_number: i32,
    pub chunks: Vec<Chunk>,
}

impl Bin {
    /// `Bin.addChunk`.
    ///
    /// Coalesces onto the last chunk when the new one starts in the same or the next BGZF
    /// block. Only the *last* chunk is considered, so the result depends on arrival order; that
    /// is fine because records arrive in coordinate order, and it is why an index built from
    /// out-of-order input would differ.
    pub fn add_chunk(&mut self, new_chunk: Chunk) {
        match self.chunks.last_mut() {
            None => self.chunks.push(new_chunk),
            Some(last) => {
                if htsjdk_bgzf::vfp::are_in_same_or_adjacent_blocks(last.end, new_chunk.start) {
                    last.end = new_chunk.end;
                } else {
                    self.chunks.push(new_chunk);
                }
            }
        }
    }
}

/// `LinearIndex.convertToLinearIndexOffset`.
pub fn to_linear_index_offset(contig_pos: i32) -> usize {
    let index_pos = if contig_pos <= 0 { 0 } else { contig_pos - 1 };
    (index_pos >> BAM_LIDX_SHIFT) as usize
}

/// `BAMIndexMetaData`: the statistics written into pseudo-bin 37450.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IndexMetaData {
    pub first_offset: i64,
    pub last_offset: i64,
    pub aligned_records: i32,
    pub unaligned_records: i32,
}

impl Default for IndexMetaData {
    /// `firstOffset` starts at -1, `lastOffset` at 0. Not symmetric, and both values are
    /// written to the file as they stand for a reference with no reads.
    fn default() -> Self {
        IndexMetaData {
            first_offset: -1,
            last_offset: 0,
            aligned_records: 0,
            unaligned_records: 0,
        }
    }
}

/// `BlockCompressedFilePointerUtil.compare`, which treats the pointers as unsigned.
fn compare_vfp(a: i64, b: i64) -> i32 {
    if a == b {
        return 0;
    }
    // When treating as unsigned, a negative number is greater than a positive one.
    if a < 0 && b >= 0 {
        return 1;
    }
    if a >= 0 && b < 0 {
        return -1;
    }
    if a < b {
        -1
    } else {
        1
    }
}

impl IndexMetaData {
    /// `BAMIndexMetaData.recordMetaData`, for a record that has coordinates.
    pub fn record(&mut self, chunk: Chunk, read_unmapped: bool) {
        let (start, end) = (chunk.start as i64, chunk.end as i64);
        if read_unmapped {
            self.unaligned_records += 1;
        } else {
            self.aligned_records += 1;
        }
        // `compare(start, firstOffset) < 1` is `start <= firstOffset`. The `|| firstOffset == -1`
        // arm is what seeds the very first record, and it comes second, so a genuine pointer of
        // -1 would take the first arm anyway.
        if compare_vfp(start, self.first_offset) < 1 || self.first_offset == -1 {
            self.first_offset = start;
        }
        if compare_vfp(self.last_offset, end) < 1 {
            self.last_offset = end;
        }
    }
}

/// `BinningIndexBuilder`.
pub struct BinningIndexBuilder {
    reference_sequence: i32,
    bins: Vec<Option<Bin>>,
    bins_seen: usize,
    linear: Vec<i64>,
    largest_index_seen: i64,
    fill_in_uninitialized_values: bool,
}

/// `AbstractBAMFileIndex.getMaxBinNumberForSequenceLength`.
fn max_bin_number_for_sequence_length(sequence_length: i32) -> usize {
    // getFirstBinInLevel(numIndexLevels - 1) is 4681.
    4681 + (sequence_length >> 14) as usize
}

impl BinningIndexBuilder {
    /// `sequence_length <= 0` means unknown, which only affects how much is allocated.
    pub fn new(
        reference_sequence: i32,
        sequence_length: i32,
        fill_in_uninitialized_values: bool,
    ) -> Self {
        let num_bins = if sequence_length <= 0 {
            MAX_BINS as usize + 1
        } else {
            max_bin_number_for_sequence_length(sequence_length) + 1
        };
        BinningIndexBuilder {
            reference_sequence,
            bins: vec![None; num_bins],
            bins_seen: 0,
            linear: vec![UNINITIALIZED_WINDOW; MAX_LINEAR_INDEX_SIZE],
            largest_index_seen: -1,
            fill_in_uninitialized_values,
        }
    }

    /// `BinningIndexBuilder.processFeature`.
    ///
    /// `start` and `end` are 1-based inclusive, as htsjdk's `FeatureToBeIndexed` specifies.
    /// `bin_number` is the record's own computed bin, which `BAMIndexer` always supplies.
    pub fn process_feature(&mut self, start: i32, end: i32, bin_number: i32, chunk: Chunk) {
        let bin_num = bin_number as usize;
        if self.bins[bin_num].is_none() {
            self.bins[bin_num] = Some(Bin {
                bin_number,
                chunks: Vec::new(),
            });
            self.bins_seen += 1;
        }
        let chunk_start = chunk.start as i64;
        self.bins[bin_num].as_mut().unwrap().add_chunk(chunk);

        // The window range this feature touches.
        let (start_window, end_window) = if end == crate::bin::NO_ALIGNMENT_START {
            // "Next line for C (samtools index) compatibility. Differs only when on a window
            // boundary": the start is decremented before conversion, which shifts a feature
            // sitting exactly on a boundary into the previous window.
            let w = to_linear_index_offset(start - 1);
            (w, w)
        } else {
            (to_linear_index_offset(start), to_linear_index_offset(end))
        };

        if end_window as i64 > self.largest_index_seen {
            self.largest_index_seen = end_window as i64;
        }

        for win in start_window..=end_window {
            if self.linear[win] == UNINITIALIZED_WINDOW || chunk_start < self.linear[win] {
                self.linear[win] = chunk_start;
            }
        }
    }

    /// `BinningIndexBuilder.generateIndexContent`. `None` when the reference has no reads.
    pub fn generate(mut self) -> Option<IndexContent> {
        if self.bins_seen == 0 {
            return None;
        }

        let largest = self.largest_index_seen as usize;
        let mut entries = vec![0i64; largest + 1];
        // "C (samtools index) also fills in intermediate 0's with values. This seems
        // unnecessary, but safe." An empty window inherits the previous non-empty offset, so a
        // port that left it at -1, or at 0, produces a valid index with different bytes.
        let mut last_non_zero_offset = 0i64;
        let fill = self.fill_in_uninitialized_values;
        for (entry, slot) in entries.iter_mut().zip(self.linear.iter_mut()) {
            if *slot == UNINITIALIZED_WINDOW {
                if fill {
                    *slot = last_non_zero_offset;
                }
            } else {
                last_non_zero_offset = *slot;
            }
            *entry = *slot;
        }

        Some(IndexContent {
            reference_sequence: self.reference_sequence,
            bins: self.bins.into_iter().flatten().collect(),
            linear_index: entries,
        })
    }
}

/// `BinningIndexContent` plus the metadata `BAMIndexContent` carries alongside it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IndexContent {
    pub reference_sequence: i32,
    /// Non-empty bins, ascending by bin number because they come out of a dense array.
    pub bins: Vec<Bin>,
    pub linear_index: Vec<i64>,
}

/// `BAMIndexer`: drives one builder per reference and writes the result.
pub struct BamIndexer {
    sequence_lengths: Vec<i32>,
    current_reference: i32,
    builder: BinningIndexBuilder,
    meta: IndexMetaData,
    /// Per-reference results, in reference order.
    finished: Vec<Option<(IndexContent, IndexMetaData)>>,
    no_coordinate_records: i64,
}

impl BamIndexer {
    pub fn new(sequence_lengths: Vec<i32>) -> Self {
        let first_len = sequence_lengths.first().copied().unwrap_or(0);
        BamIndexer {
            sequence_lengths,
            current_reference: 0,
            builder: BinningIndexBuilder::new(0, first_len, true),
            meta: IndexMetaData::default(),
            finished: Vec::new(),
            no_coordinate_records: 0,
        }
    }

    /// `BAMIndexer.processAlignment`.
    ///
    /// `chunk` is the pair of virtual file pointers taken around the record's own bytes.
    pub fn process(
        &mut self,
        reference_index: i32,
        alignment_start: i32,
        alignment_end: i32,
        index_bin: i32,
        read_unmapped: bool,
        chunk: Chunk,
    ) {
        if reference_index != crate::bin::NO_ALIGNMENT_REFERENCE_INDEX
            && reference_index != self.current_reference
        {
            self.advance_to_reference(reference_index);
        }

        // Metadata is recorded first, and counts records the index itself then skips.
        if alignment_start == crate::bin::NO_ALIGNMENT_START {
            self.no_coordinate_records += 1;
            return;
        }
        self.meta.record(chunk, read_unmapped);
        self.builder
            .process_feature(alignment_start, alignment_end, index_bin, chunk);
    }

    fn advance_to_reference(&mut self, next_reference: i32) {
        while self.current_reference < next_reference {
            let len = self
                .sequence_lengths
                .get(self.current_reference as usize + 1)
                .copied()
                .unwrap_or(0);
            let done = std::mem::replace(
                &mut self.builder,
                BinningIndexBuilder::new(self.current_reference + 1, len, true),
            );
            let meta = std::mem::take(&mut self.meta);
            self.finished.push(done.generate().map(|c| (c, meta)));
            self.current_reference += 1;
        }
    }

    /// `BAMIndexer.finish`, then `BinaryBAMIndexWriter` over the whole result.
    pub fn finish(mut self) -> Vec<u8> {
        let n_ref = self.sequence_lengths.len() as i32;
        self.advance_to_reference(n_ref);

        let mut out = Vec::new();
        out.extend_from_slice(&BAM_INDEX_MAGIC);
        out.extend_from_slice(&n_ref.to_le_bytes());
        for entry in &self.finished {
            write_reference(&mut out, entry.as_ref());
        }
        // `writeNoCoordinateRecordCount`, always written, even as zero.
        out.extend_from_slice(&self.no_coordinate_records.to_le_bytes());
        out
    }
}

/// `BinaryBAMIndexWriter.writeReference`.
fn write_reference(out: &mut Vec<u8>, content: Option<&(IndexContent, IndexMetaData)>) {
    let Some((content, meta)) = content else {
        // `writeNullContent`: a single zero long, which is 0 bins followed by 0 intervals.
        out.extend_from_slice(&0i64.to_le_bytes());
        return;
    };
    if content.bins.is_empty() {
        out.extend_from_slice(&0i64.to_le_bytes());
        return;
    }

    // The pseudo-bin is counted here but written by writeChunkMetaData below, outside the loop.
    out.extend_from_slice(&((content.bins.len() + 1) as i32).to_le_bytes());
    for bin in &content.bins {
        if bin.bin_number == MAX_BINS {
            continue;
        }
        out.extend_from_slice(&bin.bin_number.to_le_bytes());
        out.extend_from_slice(&(bin.chunks.len() as i32).to_le_bytes());
        for c in &bin.chunks {
            out.extend_from_slice(&c.start.to_le_bytes());
            out.extend_from_slice(&c.end.to_le_bytes());
        }
    }

    // `writeChunkMetaData`: pseudo-bin 37450, declaring two chunks that are not chunks at all
    // but four statistics packed into the same eight-byte slots.
    out.extend_from_slice(&MAX_BINS.to_le_bytes());
    out.extend_from_slice(&2i32.to_le_bytes());
    out.extend_from_slice(&meta.first_offset.to_le_bytes());
    out.extend_from_slice(&meta.last_offset.to_le_bytes());
    out.extend_from_slice(&(meta.aligned_records as i64).to_le_bytes());
    out.extend_from_slice(&(meta.unaligned_records as i64).to_le_bytes());

    out.extend_from_slice(&(content.linear_index.len() as i32).to_le_bytes());
    for &e in &content.linear_index {
        out.extend_from_slice(&e.to_le_bytes());
    }
}

/// Per-reference index metadata: the aligned and unaligned record counts from the pseudo-bin, or
/// `None` when the reference has no content (`writeNullContent`, so no pseudo-bin).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReferenceMetadata {
    pub aligned: i64,
    pub unaligned: i64,
}

/// The statistics `BAMIndexMetaData.printIndexStats` reads back out of a `.bai`: per-reference
/// aligned/unaligned counts (in reference order) and the total no-coordinate record count.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IndexStats {
    pub references: Vec<Option<ReferenceMetadata>>,
    pub no_coordinate_records: i64,
}

/// Why a `.bai` could not be parsed.
#[derive(Debug)]
pub enum BaiParseError {
    NotABai,
    Truncated,
}

/// Reads the metadata out of a `.bai`, the way `BAMIndexMetaData(List<Chunk>)` does: the pseudo-bin
/// 37450 declares two "chunks" whose eight-byte slots are really four statistics; the second pair is
/// `(alignedRecords, unalignedRecords)`. A reference written as null content has no pseudo-bin and so
/// no metadata. The no-coordinate count is the trailing long.
pub fn parse_bai_metadata(bai: &[u8]) -> Result<IndexStats, BaiParseError> {
    let mut p = 0usize;
    let take = |p: &mut usize, n: usize| -> Result<&[u8], BaiParseError> {
        let s = bai.get(*p..*p + n).ok_or(BaiParseError::Truncated)?;
        *p += n;
        Ok(s)
    };
    let i32_at = |p: &mut usize| -> Result<i32, BaiParseError> {
        Ok(i32::from_le_bytes(take(p, 4)?.try_into().unwrap()))
    };
    let i64_at = |p: &mut usize| -> Result<i64, BaiParseError> {
        Ok(i64::from_le_bytes(take(p, 8)?.try_into().unwrap()))
    };

    if take(&mut p, 4)? != BAM_INDEX_MAGIC {
        return Err(BaiParseError::NotABai);
    }
    let n_ref = i32_at(&mut p)?;
    let mut references = Vec::with_capacity(n_ref.max(0) as usize);
    for _ in 0..n_ref {
        let n_bin = i32_at(&mut p)?;
        let mut meta = None;
        for _ in 0..n_bin {
            let bin_number = i32_at(&mut p)?;
            let n_chunk = i32_at(&mut p)?;
            if bin_number == MAX_BINS && n_chunk == 2 {
                let _first = i64_at(&mut p)?;
                let _last = i64_at(&mut p)?;
                let aligned = i64_at(&mut p)?;
                let unaligned = i64_at(&mut p)?;
                meta = Some(ReferenceMetadata { aligned, unaligned });
            } else {
                // Skip this bin's chunks (two eight-byte offsets each).
                take(&mut p, 16 * n_chunk.max(0) as usize)?;
            }
        }
        let n_intv = i32_at(&mut p)?;
        take(&mut p, 8 * n_intv.max(0) as usize)?;
        references.push(meta);
    }
    let no_coordinate_records = i64_at(&mut p)?;
    Ok(IndexStats {
        references,
        no_coordinate_records,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn chunk(start: u64, end: u64) -> Chunk {
        Chunk { start, end }
    }

    fn vfp(block: u64, offset: u32) -> u64 {
        htsjdk_bgzf::vfp::make_file_pointer(block, offset).unwrap()
    }

    #[test]
    fn the_linear_window_is_sixteen_kilobases() {
        assert_eq!(to_linear_index_offset(1), 0);
        assert_eq!(to_linear_index_offset(16_384), 0);
        assert_eq!(to_linear_index_offset(16_385), 1);
        // Non-positive positions clamp to window 0 rather than going negative.
        assert_eq!(to_linear_index_offset(0), 0);
        assert_eq!(to_linear_index_offset(-5), 0);
    }

    #[test]
    fn chunks_in_the_same_block_coalesce() {
        let mut bin = Bin {
            bin_number: 4681,
            chunks: Vec::new(),
        };
        bin.add_chunk(chunk(vfp(0, 100), vfp(0, 200)));
        bin.add_chunk(chunk(vfp(0, 200), vfp(0, 300)));
        assert_eq!(bin.chunks.len(), 1, "same block must coalesce");
        assert_eq!(bin.chunks[0], chunk(vfp(0, 100), vfp(0, 300)));
    }

    /// "Adjacent" is block address + 1. Block addresses are byte offsets, so two real blocks
    /// are adjacent only if one is a single byte long, which cannot happen. The rule is
    /// effectively "same block", and reproducing it literally costs nothing.
    #[test]
    fn adjacency_is_block_address_plus_one_not_the_next_real_block() {
        let mut bin = Bin {
            bin_number: 4681,
            chunks: Vec::new(),
        };
        bin.add_chunk(chunk(vfp(0, 100), vfp(0, 200)));
        // A realistic next block starts thousands of bytes later: no coalescing.
        bin.add_chunk(chunk(vfp(5000, 0), vfp(5000, 50)));
        assert_eq!(bin.chunks.len(), 2);

        let mut bin2 = Bin {
            bin_number: 4681,
            chunks: Vec::new(),
        };
        bin2.add_chunk(chunk(vfp(0, 100), vfp(0, 200)));
        // Block address literally one greater does coalesce, however unreachable that is.
        bin2.add_chunk(chunk(vfp(1, 0), vfp(1, 50)));
        assert_eq!(bin2.chunks.len(), 1);
    }

    /// The back-fill: empty windows carry the previous non-empty offset, not zero.
    #[test]
    fn empty_linear_windows_inherit_the_previous_offset() {
        let mut b = BinningIndexBuilder::new(0, 250_000_000, true);
        // A read in window 0, then one far away in window 5. Windows 1..4 are untouched.
        b.process_feature(1, 100, 4681, chunk(vfp(0, 0), vfp(0, 100)));
        let far = 5 * (1 << BAM_LIDX_SHIFT) + 1;
        b.process_feature(far, far + 99, 4686, chunk(vfp(0, 500), vfp(0, 600)));
        let content = b.generate().unwrap();

        assert_eq!(content.linear_index.len(), 6);
        assert_eq!(content.linear_index[0], vfp(0, 0) as i64);
        for w in 1..=4 {
            assert_eq!(
                content.linear_index[w],
                vfp(0, 0) as i64,
                "window {w} must inherit, not be zero or -1"
            );
        }
        assert_eq!(content.linear_index[5], vfp(0, 500) as i64);
    }

    /// Without the fill, empty windows stay at -1. That mode exists for index merging, and
    /// keeping both proves the fill is a deliberate choice rather than an artefact.
    #[test]
    fn the_fill_can_be_turned_off() {
        let mut b = BinningIndexBuilder::new(0, 250_000_000, false);
        b.process_feature(1, 100, 4681, chunk(vfp(0, 0), vfp(0, 100)));
        let far = 3 * (1 << BAM_LIDX_SHIFT) + 1;
        b.process_feature(far, far + 99, 4684, chunk(vfp(0, 500), vfp(0, 600)));
        let content = b.generate().unwrap();
        assert_eq!(content.linear_index[1], UNINITIALIZED_WINDOW);
    }

    /// A window whose only read starts at virtual offset 0 must not be confused with an empty
    /// one. This is why the sentinel is -1 and not 0.
    #[test]
    fn a_read_at_offset_zero_is_not_an_empty_window() {
        let mut b = BinningIndexBuilder::new(0, 250_000_000, true);
        b.process_feature(1, 100, 4681, chunk(0, 100));
        let content = b.generate().unwrap();
        assert_eq!(content.linear_index[0], 0);
    }

    #[test]
    fn a_reference_with_no_reads_generates_nothing() {
        let b = BinningIndexBuilder::new(0, 250_000_000, true);
        assert!(b.generate().is_none());
    }

    #[test]
    fn the_file_opens_with_the_bai_magic_and_reference_count() {
        let idx = BamIndexer::new(vec![250_000_000, 200_000_000]).finish();
        assert_eq!(&idx[0..4], b"BAI\x01");
        assert_eq!(i32::from_le_bytes(idx[4..8].try_into().unwrap()), 2);
        // Two empty references, each a single zero long, then the no-coordinate count.
        assert_eq!(idx.len(), 8 + 8 + 8 + 8);
        assert_eq!(
            i64::from_le_bytes(idx[idx.len() - 8..].try_into().unwrap()),
            0
        );
    }

    /// The pseudo-bin is counted in `n_bin` but skipped by the loop that writes real bins.
    #[test]
    fn the_pseudo_bin_is_counted_and_written_separately() {
        let mut ix = BamIndexer::new(vec![250_000_000]);
        ix.process(0, 100, 103, 4681, false, chunk(vfp(0, 0), vfp(0, 50)));
        let idx = ix.finish();

        let n_bin = i32::from_le_bytes(idx[8..12].try_into().unwrap());
        assert_eq!(n_bin, 2, "one real bin plus the pseudo-bin");

        // Real bin first.
        assert_eq!(i32::from_le_bytes(idx[12..16].try_into().unwrap()), 4681);
        assert_eq!(i32::from_le_bytes(idx[16..20].try_into().unwrap()), 1);
        // Then the pseudo-bin, declaring two chunks that are really four statistics.
        let p = 20 + 16;
        assert_eq!(
            i32::from_le_bytes(idx[p..p + 4].try_into().unwrap()),
            MAX_BINS
        );
        assert_eq!(i32::from_le_bytes(idx[p + 4..p + 8].try_into().unwrap()), 2);
        assert_eq!(
            i64::from_le_bytes(idx[p + 24..p + 32].try_into().unwrap()),
            1,
            "one aligned record"
        );
        assert_eq!(
            i64::from_le_bytes(idx[p + 32..p + 40].try_into().unwrap()),
            0,
            "no unaligned records"
        );
    }

    #[test]
    fn records_without_coordinates_are_counted_not_indexed() {
        let mut ix = BamIndexer::new(vec![250_000_000]);
        ix.process(0, 100, 103, 4681, false, chunk(vfp(0, 0), vfp(0, 50)));
        for _ in 0..7 {
            ix.process(-1, 0, 0, 0, true, chunk(vfp(0, 50), vfp(0, 60)));
        }
        let idx = ix.finish();
        assert_eq!(
            i64::from_le_bytes(idx[idx.len() - 8..].try_into().unwrap()),
            7
        );
    }

    /// An unmapped-but-placed read counts as unaligned in the metadata while still being
    /// indexed. The two counters answer different questions and are easy to conflate.
    #[test]
    fn a_placed_unmapped_read_is_indexed_but_counted_as_unaligned() {
        let mut meta = IndexMetaData::default();
        meta.record(chunk(vfp(0, 0), vfp(0, 50)), true);
        meta.record(chunk(vfp(0, 50), vfp(0, 90)), false);
        assert_eq!(meta.aligned_records, 1);
        assert_eq!(meta.unaligned_records, 1);
    }

    #[test]
    fn the_first_offset_starts_at_minus_one_and_the_last_at_zero() {
        let meta = IndexMetaData::default();
        assert_eq!((meta.first_offset, meta.last_offset), (-1, 0));
    }

    #[test]
    fn metadata_round_trips_through_a_built_index() {
        use crate::build_index::build_bam_index;
        use crate::cigar::{Cigar, CigarElement, Op};
        use crate::header::{SamHeader, SequenceRecord};
        use crate::record::{BamRecord, READ_UNMAPPED_FLAG};
        use crate::writer::BamWriter;

        let mut h = SamHeader::new();
        h.set_sort_order("coordinate");
        h.sequences.push(SequenceRecord::new("chr1", 100_000));
        h.sequences.push(SequenceRecord::new("chr2", 100_000));
        let m10 = || {
            Cigar::new(vec![CigarElement {
                length: 10,
                op: Op::M,
            }])
        };
        let mut w = BamWriter::new(Vec::new(), &h).unwrap();
        // chr1: two mapped reads. chr2: none. Then two unplaced (no-coordinate) reads.
        for start in [10, 500] {
            w.write(&BamRecord {
                read_name: "m".into(),
                reference_index: 0,
                alignment_start: start,
                mapping_quality: 60,
                cigar: m10(),
                ..BamRecord::default()
            })
            .unwrap();
        }
        for _ in 0..2 {
            w.write(&BamRecord {
                read_name: "u".into(),
                reference_index: -1,
                flags: READ_UNMAPPED_FLAG,
                ..BamRecord::default()
            })
            .unwrap();
        }
        let bam = w.finish().unwrap();

        let bai = build_bam_index(&bam).unwrap();
        let stats = parse_bai_metadata(&bai).unwrap();

        assert_eq!(stats.references.len(), 2);
        assert_eq!(
            stats.references[0],
            Some(ReferenceMetadata {
                aligned: 2,
                unaligned: 0
            })
        );
        // chr2 has no reads, so no pseudo-bin and no metadata.
        assert_eq!(stats.references[1], None);
        assert_eq!(stats.no_coordinate_records, 2);
    }
}
