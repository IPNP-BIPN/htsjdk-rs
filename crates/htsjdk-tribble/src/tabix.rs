//! The `.tbi` beside a block-compressed feature file, which is not the `.idx` beside a plain one.
//!
//! [`crate::index_write`] builds the Tribble index a **plain** VCF gets. A BGZF one gets a tabix
//! index instead, and the two share neither their layout nor the machinery under them: the Tribble
//! index picks a bin width by estimating density, while tabix is the BAM's own binning scheme,
//! driven by the very same [`htsjdk_bam::index::BinningIndexBuilder`] a `.bai` is built with.
//!
//! # A feature is indexed one feature late
//!
//! `addFeature` is handed a feature and the position it starts at, and a chunk needs both ends. The
//! end of one feature is the START of the next, so nothing can be indexed until the next arrives
//! and the last one waits for `finalizeIndex`'s own final position. A port that closed each chunk
//! on the feature it was given would write chunks of length zero and an index that finds nothing.
//!
//! # The bin is the region's, not the feature's
//!
//! `TabixFeature.getIndexingBin` returns null, so the builder computes it, and it computes it from
//! a **zero-based half-open** region: `regionToBin(start - 1, end)`. A feature with no end (`end`
//! at or below zero) is treated as one base, `start - 1` to `start`.
//!
//! # The sequence order is the file's, and it is checked
//!
//! References are numbered in the order their first feature appears, not by any dictionary, and a
//! sequence that comes back after another has intervened is refused rather than merged. Within a
//! sequence, features must not go backwards by start.
//!
//! # What the header says, and what it does not
//!
//! The eight header integers are the magic, the sequence count, and the six fields of
//! [`TabixFormat`]. `numHeaderLinesToSkip` is written and, in htsjdk's own words, "does not appear
//! to be used"; it is carried here for the same reason, which is that it is bytes of the file.
//!
//! Ported from `htsjdk.tribble.index.tabix.TabixIndex`,
//! `htsjdk.tribble.index.tabix.TabixIndexCreator` and
//! `htsjdk.tribble.index.tabix.TabixFormat`.

use htsjdk_bam::bin::region_to_bin;
use htsjdk_bam::index::{BinningIndexBuilder, Chunk, IndexContent};
use htsjdk_bgzf::Deflater;

/// `TabixFormat`'s flag bits.
pub const ZERO_BASED: i32 = 0x10000;
pub const GENERIC_FLAGS: i32 = 0;
pub const SAM_FLAGS: i32 = 1;
pub const VCF_FLAGS: i32 = 2;
pub const UCSC_FLAGS: i32 = GENERIC_FLAGS | ZERO_BASED;

/// `TabixIndex.MAGIC`, which is `TBI\1` read back as a little-endian integer.
pub const MAGIC: [u8; 4] = *b"TBI\x01";

/// How the file being indexed is laid out, which is six integers of the index's header.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TabixFormat {
    pub flags: i32,
    /// One-based column holding the sequence name.
    pub sequence_column: i32,
    /// One-based column holding the start position.
    pub start_position_column: i32,
    /// One-based column holding the end position, or zero where there is none.
    pub end_position_column: i32,
    /// Lines beginning with this are ignored by a reader.
    pub meta_character: u8,
    /// Written, and htsjdk's own comment says it does not appear to be used.
    pub num_header_lines_to_skip: i32,
}

impl TabixFormat {
    pub const GFF: TabixFormat = TabixFormat::new(GENERIC_FLAGS, 1, 4, 5, b'#', 0);
    pub const BED: TabixFormat = TabixFormat::new(UCSC_FLAGS, 1, 2, 3, b'#', 0);
    pub const PSLTBL: TabixFormat = TabixFormat::new(UCSC_FLAGS, 15, 17, 18, b'#', 0);
    pub const SAM: TabixFormat = TabixFormat::new(SAM_FLAGS, 3, 4, 0, b'@', 0);
    pub const VCF: TabixFormat = TabixFormat::new(VCF_FLAGS, 1, 2, 0, b'#', 0);

    pub const fn new(
        flags: i32,
        sequence_column: i32,
        start_position_column: i32,
        end_position_column: i32,
        meta_character: u8,
        num_header_lines_to_skip: i32,
    ) -> Self {
        Self {
            flags,
            sequence_column,
            start_position_column,
            end_position_column,
            meta_character,
            num_header_lines_to_skip,
        }
    }
}

/// What a creator refuses, in the reference's own words.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TabixError {
    /// A sequence came back after another had intervened.
    SequenceOutOfOrder { feature: String },
    /// Two features in the same sequence went backwards by start.
    FeaturesOutOfOrder { previous: String, next: String },
    /// A chunk whose end is not after its start, which is a position that did not advance.
    EmptyChunk { start: i64, end: i64 },
}

impl TabixError {
    /// The exception class the reference throws, which is one class for all three.
    pub fn java_class(&self) -> &'static str {
        "java.lang.IllegalArgumentException"
    }

    pub fn message(&self) -> String {
        match self {
            TabixError::SequenceOutOfOrder { feature } => {
                format!("Sequence {feature} added out sequence of order")
            }
            TabixError::FeaturesOutOfOrder { previous, next } => {
                format!("Features added out of order: previous ({previous}) > next ({next})")
            }
            TabixError::EmptyChunk { start, end } => {
                format!("Feature start position {start} >= feature end position {end}")
            }
        }
    }
}

/// One feature on its way into the index.
///
/// `description` is the feature's own `toString`, which one refusal quotes verbatim. GATK hands
/// the creator a `VariantContext` and htsjdk prints whatever that class prints, so the text cannot
/// be composed here: the caller supplies it, and a caller that has nothing to say can pass the
/// contig.
///
/// `sequence_length` is what a sequence dictionary would have supplied, and `0` where there is
/// none. It is read only when a new sequence opens, and it decides how much the builder allocates,
/// never what it writes.
#[derive(Debug, Clone, Copy)]
pub struct FeatureRef<'a> {
    pub contig: &'a str,
    /// One-based inclusive, as `FeatureToBeIndexed` specifies.
    pub start: i32,
    /// One-based inclusive, or zero where the feature has no end.
    pub end: i32,
    pub description: &'a str,
    pub sequence_length: i32,
}

/// One feature waiting for the next one to close its chunk.
#[derive(Debug, Clone)]
struct PendingFeature {
    reference_index: i32,
    start: i32,
    end: i32,
    start_file_position: i64,
}

impl PendingFeature {
    /// `TabixFeature.toString`, which is what a refusal quotes.
    fn describe(&self, end_file_position: i64) -> String {
        format!(
            "TabixFeature{{referenceIndex={}, start={}, end={}, featureStartFilePosition={}, \
             featureEndFilePosition={}}}",
            self.reference_index, self.start, self.end, self.start_file_position, end_file_position
        )
    }
}

/// `TabixIndexCreator`: features and their file positions in, one index out.
pub struct TabixIndexCreator {
    format: TabixFormat,
    contents: Vec<Option<IndexContent>>,
    sequence_names: Vec<String>,
    current_sequence: Option<String>,
    builder: Option<BinningIndexBuilder>,
    previous: Option<PendingFeature>,
}

impl TabixIndexCreator {
    pub fn new(format: TabixFormat) -> Self {
        Self {
            format,
            contents: Vec::new(),
            sequence_names: Vec::new(),
            current_sequence: None,
            builder: None,
            previous: None,
        }
    }

    /// `addFeature`, whose file position is the BGZF virtual offset the feature starts at.
    pub fn add_feature(
        &mut self,
        feature: FeatureRef<'_>,
        file_position: i64,
    ) -> Result<(), TabixError> {
        let contig = feature.contig;
        let same = self.current_sequence.as_deref() == Some(contig);
        let reference_index = if same {
            self.sequence_names.len() as i32 - 1
        } else {
            if self.current_sequence.is_some() && self.sequence_names.iter().any(|n| n == contig) {
                return Err(TabixError::SequenceOutOfOrder {
                    feature: feature.description.to_string(),
                });
            }
            self.sequence_names.len() as i32
        };

        let this = PendingFeature {
            reference_index,
            start: feature.start,
            end: feature.end,
            start_file_position: file_position,
        };
        if let Some(previous) = self.previous.clone() {
            // `TabixFeature.compareTo`: the reference index first, then the start. The END is not
            // compared, so two features sharing a start are in order either way round.
            let backwards =
                (previous.reference_index, previous.start) > (this.reference_index, this.start);
            if backwards {
                return Err(TabixError::FeaturesOutOfOrder {
                    previous: previous.describe(-1),
                    next: this.describe(-1),
                });
            }
            self.finalize_feature(file_position)?;
        }
        self.previous = Some(this);
        if reference_index == self.sequence_names.len() as i32 {
            self.advance_to_reference(contig, feature.sequence_length);
        }
        Ok(())
    }

    /// `finalizeFeature`: the position the next feature starts at is where this one ended.
    fn finalize_feature(&mut self, end_file_position: i64) -> Result<(), TabixError> {
        let previous = self.previous.take().expect("a feature to finalize");
        if previous.start_file_position >= end_file_position {
            return Err(TabixError::EmptyChunk {
                start: previous.start_file_position,
                end: end_file_position,
            });
        }
        // `getIndexingBin` is null, so the builder computes it from a zero-based half-open region.
        let end = if previous.end <= 0 {
            previous.start
        } else {
            previous.end
        };
        let bin = region_to_bin(previous.start - 1, end);
        self.builder
            .as_mut()
            .expect("a builder for the current sequence")
            .process_feature(
                previous.start,
                previous.end,
                bin,
                Chunk {
                    start: previous.start_file_position as u64,
                    end: end_file_position as u64,
                },
            );
        Ok(())
    }

    /// `advanceToReference`: the previous sequence's content is generated and a builder opened.
    fn advance_to_reference(&mut self, contig: &str, sequence_length: i32) {
        if let Some(builder) = self.builder.take() {
            self.contents.push(builder.generate());
        }
        self.builder = Some(BinningIndexBuilder::new(
            self.sequence_names.len() as i32,
            sequence_length,
            true,
        ));
        self.sequence_names.push(contig.to_string());
        self.current_sequence = Some(contig.to_string());
    }

    /// `finalizeIndex`: the last feature is closed against the file's own end.
    pub fn finish(mut self, final_file_position: i64) -> Result<TabixIndex, TabixError> {
        if self.previous.is_some() {
            self.finalize_feature(final_file_position)?;
        }
        if let Some(builder) = self.builder.take() {
            self.contents.push(builder.generate());
        }
        Ok(TabixIndex {
            format: self.format,
            sequence_names: self.sequence_names,
            indices: self.contents,
        })
    }
}

/// `TabixIndex`, which is the header, the names, and one block per sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TabixIndex {
    pub format: TabixFormat,
    pub sequence_names: Vec<String>,
    /// One entry per sequence, `None` where the sequence carried no feature.
    pub indices: Vec<Option<IndexContent>>,
}

impl TabixIndex {
    /// `TabixIndex.write`: the little-endian body, before it is block compressed.
    pub fn write_body(&self) -> Vec<u8> {
        let mut out = Vec::new();
        out.extend_from_slice(&MAGIC);
        out.extend_from_slice(&(self.sequence_names.len() as i32).to_le_bytes());
        out.extend_from_slice(&self.format.flags.to_le_bytes());
        out.extend_from_slice(&self.format.sequence_column.to_le_bytes());
        out.extend_from_slice(&self.format.start_position_column.to_le_bytes());
        out.extend_from_slice(&self.format.end_position_column.to_le_bytes());
        out.extend_from_slice(&(self.format.meta_character as i32).to_le_bytes());
        out.extend_from_slice(&self.format.num_header_lines_to_skip.to_le_bytes());

        // The name block counts its null terminators, and the names are written raw: htsjdk's
        // `StringUtil.stringToBytes` takes the low byte of each char rather than encoding UTF-8.
        let name_block_size: usize = self.sequence_names.len()
            + self
                .sequence_names
                .iter()
                .map(|name| name.chars().count())
                .sum::<usize>();
        out.extend_from_slice(&(name_block_size as i32).to_le_bytes());
        for name in &self.sequence_names {
            out.extend(name.chars().map(|c| c as u32 as u8));
            out.push(0);
        }

        for index in &self.indices {
            match index {
                None => out.extend_from_slice(&0i32.to_le_bytes()),
                Some(content) => {
                    out.extend_from_slice(&(content.bins.len() as i32).to_le_bytes());
                    for bin in &content.bins {
                        out.extend_from_slice(&bin.bin_number.to_le_bytes());
                        out.extend_from_slice(&(bin.chunks.len() as i32).to_le_bytes());
                        for chunk in &bin.chunks {
                            out.extend_from_slice(&(chunk.start as i64).to_le_bytes());
                            out.extend_from_slice(&(chunk.end as i64).to_le_bytes());
                        }
                    }
                    out.extend_from_slice(&(content.linear_index.len() as i32).to_le_bytes());
                    for entry in &content.linear_index {
                        out.extend_from_slice(&entry.to_le_bytes());
                    }
                }
            }
        }
        out
    }

    /// `TabixIndex.write(Path)`: the body inside a BGZF stream, which is the `.tbi` itself.
    ///
    /// The reference wraps a `LittleEndianOutputStream` around a `BlockCompressedOutputStream`
    /// built with a null path, which takes the default compression level. The terminator block is
    /// part of the file: `close` writes it, and a reader that checks termination refuses a file
    /// without it.
    ///
    /// The deflater is htsjdk's own, which is the JDK's. A caller running under GATK writes the
    /// same body through [`Self::write_with`], because GATK replaces the static factory and every
    /// byte after the framing differs.
    pub fn write(&self) -> Vec<u8> {
        self.write_with(Deflater::Jdk)
    }

    /// The same file with the deflater named. See [`htsjdk_bgzf::Deflater`].
    pub fn write_with(&self, deflater: Deflater) -> Vec<u8> {
        let mut writer = htsjdk_bgzf::BgzfWriter::with_deflater(
            Vec::new(),
            htsjdk_bgzf::DEFAULT_COMPRESSION_LEVEL,
            deflater,
        );
        std::io::Write::write_all(&mut writer, &self.write_body()).expect("a vector never fails");
        writer.into_inner().expect("a vector never fails")
    }
}
