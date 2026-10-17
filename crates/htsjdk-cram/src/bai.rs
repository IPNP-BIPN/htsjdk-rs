//! A BAI index of a CRAM, as `CRAMBAIIndexer` builds one while `CRAMFileWriter` writes: what
//! `SAMFileWriterFactory` gives a `.cram` output when `CREATE_INDEX` is set, as `<output>.bai`.
//!
//! Ported from `CRAMBAIIndexer`, `BAIEntry`, `Slice.getBAIEntries` and `BAMIndexMetaData` at
//! htsjdk 5.0.0.
//!
//! # A chunk is a slice, not a record
//!
//! Its virtual pointer is the container's byte offset shifted sixteen, with the slice's index in
//! the container where a BGZF pointer would hold the offset in the block, and it ends one past
//! that. The metadata's first offset is something else again: the slice's byte offset from the
//! end of the container header, which is what `recordMetaData` is handed for a CRAM.

/// `BAIEntry`: one slice's contribution, or one reference's share of a multi-reference slice.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BaiEntry {
    /// A sequence index, or -1 for the unmapped and unplaced.
    pub reference: i32,
    pub alignment_start: i32,
    pub alignment_span: i32,
    pub mapped: i32,
    pub unmapped: i32,
    pub unplaced: i32,
    pub container_offset: i64,
    pub slice_offset: i64,
    pub landmark_index: i32,
}

impl BaiEntry {
    /// `BAIEntry.compareTo`: references in order with the unmapped last, then start for placed
    /// entries, then the container and the slice.
    pub fn compare(&self, other: &Self) -> std::cmp::Ordering {
        use std::cmp::Ordering;
        if self.reference != other.reference {
            if self.reference == -1 {
                return Ordering::Greater;
            }
            if other.reference == -1 {
                return Ordering::Less;
            }
            return self.reference.cmp(&other.reference);
        }
        if self.reference != -1 && self.alignment_start != other.alignment_start {
            return self.alignment_start.cmp(&other.alignment_start);
        }
        self.container_offset
            .cmp(&other.container_offset)
            .then(self.slice_offset.cmp(&other.slice_offset))
    }
}

/// `AlignmentSpan`, for merging a multi-reference slice's records per reference.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Span {
    pub start: i32,
    pub span: i32,
    pub mapped: i32,
    pub unmapped: i32,
    pub unplaced: i32,
}

impl Span {
    /// `AlignmentSpan.combine`.
    pub fn combine(a: Span, b: Span) -> Span {
        let start = a.start.min(b.start);
        let span = if a.start == b.start {
            a.span.max(b.span)
        } else {
            (a.start + a.span).max(b.start + b.span) - start
        };
        Span {
            start,
            span,
            mapped: a.mapped + b.mapped,
            unmapped: a.unmapped + b.unmapped,
            unplaced: a.unplaced + b.unplaced,
        }
    }
}

const MAX_BINS: usize = 37450;
const MAX_LINEAR_INDEX_SIZE: usize = 32770;

#[derive(Clone, Copy)]
struct Chunk {
    start: u64,
    end: u64,
}

/// `BAMIndexMetaData`, as a CRAM's entries feed it.
#[derive(Clone, Copy)]
struct Meta {
    first_offset: i64,
    last_offset: i64,
    aligned: i64,
    unaligned: i64,
}

impl Meta {
    fn new() -> Self {
        Self {
            first_offset: -1,
            last_offset: 0,
            aligned: 0,
            unaligned: 0,
        }
    }
}

/// `BlockCompressedFilePointerUtil.compare`: the pointers as unsigned.
fn compare_unsigned(a: i64, b: i64) -> i32 {
    match (a as u64).cmp(&(b as u64)) {
        std::cmp::Ordering::Less => -1,
        std::cmp::Ordering::Equal => 0,
        std::cmp::Ordering::Greater => 1,
    }
}

/// `CRAMBAIIndexer` with its `CRAMBAIIndexBuilder`.
pub struct CramBaiIndexer {
    sequence_lengths: Vec<i32>,
    current: i32,
    bins: Option<Vec<Option<Vec<Chunk>>>>,
    bins_seen: usize,
    linear: Vec<u64>,
    largest_window: i64,
    meta: Meta,
    no_coordinate: i64,
    out: Vec<u8>,
}

impl CramBaiIndexer {
    pub fn new(sequence_lengths: Vec<i32>) -> Self {
        let mut out = b"BAI\x01".to_vec();
        out.extend_from_slice(&(sequence_lengths.len() as i32).to_le_bytes());
        Self {
            sequence_lengths,
            current: 0,
            bins: None,
            bins_seen: 0,
            linear: vec![0; MAX_LINEAR_INDEX_SIZE],
            largest_window: -1,
            meta: Meta::new(),
            no_coordinate: 0,
            out,
        }
    }

    /// `processContainer`: the container's entries, sorted.
    pub fn process_container(&mut self, mut entries: Vec<BaiEntry>) -> Result<(), String> {
        entries.sort_by(BaiEntry::compare);
        for entry in entries {
            self.process(entry)?;
        }
        Ok(())
    }

    /// `processBAIEntry`.
    fn process(&mut self, entry: BaiEntry) -> Result<(), String> {
        if entry.reference >= 0 {
            if entry.reference != self.current {
                self.advance(entry.reference);
            }
            if entry.reference != self.current {
                return Err(format!(
                    "htsjdk.samtools.SAMException: Unexpected reference {} when constructing index for reference {} for slice",
                    entry.reference, self.current
                ));
            }
        }
        // `recordMetaData(BAIEntry)`.
        self.meta.aligned += i64::from(entry.mapped);
        self.no_coordinate += i64::from(entry.unplaced);
        self.meta.unaligned += i64::from(entry.unmapped);
        let start = entry.slice_offset;
        if compare_unsigned(start, self.meta.first_offset) < 1 || self.meta.first_offset == -1 {
            self.meta.first_offset = start;
            self.meta.last_offset = start;
        }
        if entry.reference < 0 {
            return Ok(());
        }

        // The bin.
        let start0 = entry.alignment_start - 1;
        let mut end = entry.alignment_start + entry.alignment_span - 1;
        if end <= start0 {
            end = start0 + 1;
        }
        let bin = htsjdk_bam::bin::region_to_bin(start0, end) as usize;
        if self.bins.is_none() {
            let size = self
                .sequence_lengths
                .get(entry.reference as usize)
                .map_or(MAX_BINS, |&l| max_bin_for_length(l))
                + 1;
            self.bins = Some(vec![None; size]);
        }
        let bins = self.bins.as_mut().unwrap();
        let slot = bins.get_mut(bin).ok_or_else(|| {
            format!("java.lang.ArrayIndexOutOfBoundsException: Index {bin} out of bounds")
        })?;
        let pointer = ((entry.container_offset as u64) << 16) | entry.landmark_index as u64;
        let chunk = Chunk {
            start: pointer,
            end: pointer + 1,
        };
        match slot {
            None => {
                *slot = Some(vec![chunk]);
                self.bins_seen += 1;
            }
            Some(chunks) => {
                let last = chunks.last_mut().expect("a seen bin has a chunk");
                let (a, b) = (last.end >> 16, chunk.start >> 16);
                if a == b || a + 1 == b {
                    last.end = chunk.end;
                } else {
                    chunks.push(chunk);
                }
            }
        }

        // The linear index.
        let alignment_end = entry.alignment_start + entry.alignment_span;
        let (start_window, end_window) = if alignment_end == 0 {
            let w = linear_offset(entry.alignment_start - 1);
            (w, w)
        } else {
            (
                linear_offset(entry.alignment_start),
                linear_offset(alignment_end),
            )
        };
        if end_window as i64 > self.largest_window {
            self.largest_window = end_window as i64;
        }
        for window in start_window..=end_window {
            if self.linear[window] == 0 || pointer < self.linear[window] {
                self.linear[window] = pointer;
            }
        }
        Ok(())
    }

    /// `advanceToReference`: every reference before `next` written out.
    fn advance(&mut self, next: i32) {
        while self.current < next {
            self.write_current();
            self.current += 1;
            // `startNewReference`.
            self.bins = None;
            if self.bins_seen > 0 {
                self.linear.iter_mut().for_each(|w| *w = 0);
            }
            self.bins_seen = 0;
            self.largest_window = -1;
            self.meta = Meta::new();
        }
    }

    /// `processCurrentReference` and `BinaryBAMIndexWriter.writeReference`.
    fn write_current(&mut self) {
        if self.bins_seen == 0 {
            self.out.extend_from_slice(&0i64.to_le_bytes());
            return;
        }
        let bins = self.bins.as_ref().expect("bins were seen");
        let count = bins.iter().filter(|b| b.is_some()).count();
        self.out
            .extend_from_slice(&((count + 1) as i32).to_le_bytes());
        for (number, chunks) in bins.iter().enumerate() {
            let Some(chunks) = chunks else { continue };
            if number == MAX_BINS {
                continue;
            }
            self.out.extend_from_slice(&(number as i32).to_le_bytes());
            self.out
                .extend_from_slice(&(chunks.len() as i32).to_le_bytes());
            for c in chunks {
                self.out.extend_from_slice(&c.start.to_le_bytes());
                self.out.extend_from_slice(&c.end.to_le_bytes());
            }
        }
        self.out.extend_from_slice(&(MAX_BINS as i32).to_le_bytes());
        self.out.extend_from_slice(&2i32.to_le_bytes());
        self.out
            .extend_from_slice(&self.meta.first_offset.to_le_bytes());
        self.out
            .extend_from_slice(&self.meta.last_offset.to_le_bytes());
        self.out.extend_from_slice(&self.meta.aligned.to_le_bytes());
        self.out
            .extend_from_slice(&self.meta.unaligned.to_le_bytes());

        let length = (self.largest_window + 1).max(0) as usize;
        self.out.extend_from_slice(&(length as i32).to_le_bytes());
        let mut last = 0u64;
        for window in 0..length {
            if self.linear[window] == 0 {
                self.linear[window] = last;
            } else {
                last = self.linear[window];
            }
            self.out
                .extend_from_slice(&(self.linear[window] as i64).to_le_bytes());
        }
    }

    /// `finish`: every remaining reference, then the count of records with no coordinate.
    pub fn finish(mut self) -> Vec<u8> {
        let references = self.sequence_lengths.len() as i32;
        self.advance(references);
        self.out
            .extend_from_slice(&self.no_coordinate.to_le_bytes());
        self.out
    }
}

/// `LinearIndex.convertToLinearIndexOffset`.
fn linear_offset(position: i32) -> usize {
    let index = if position <= 0 { 0 } else { position - 1 };
    (index >> 14) as usize
}

/// `AbstractBAMFileIndex.getMaxBinNumberForSequenceLength`.
fn max_bin_for_length(length: i32) -> usize {
    (4681 + (length.max(0) >> 14)) as usize
}
