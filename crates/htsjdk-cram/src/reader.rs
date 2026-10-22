//! Reading a CRAM into SAM records: `CRAMFileReader`, `CRAMIterator` and the normalisation in
//! `Slice.normalizeCRAMRecords` and `CRAMCompressionRecord.toSAMRecord`, at htsjdk 4.2.0.
//!
//! [`crate::file`] walks the containers and [`crate::record_read`] decodes one record; this is
//! what turns a slice's raw records into the records a tool sees: mates linked, names generated,
//! bases restored against the reference, scores resolved, tags decoded.
//!
//! # The reference is fetched per slice, and lazily
//!
//! The header is read when the file is opened and a record is decoded only when its container is
//! reached, so a missing reference is not an error until the first mapped slice asks for bases.
//! That is why [`CramReader::next_container`] takes the reference and the constructor does not:
//! a tool that has already written output when the reference fails has written it in Picard too.
//!
//! # A reference region is a copy of the slice's span
//!
//! `CRAMReferenceRegion` holds the bases from the slice's alignment start for its span, uppercased,
//! and the MD5 the slice header carries is checked against exactly those bytes, truncated to what
//! the contig has.

use crate::file::{read_container, read_preamble, ContainerWalk, FileError};
use crate::read_features::{alignment_end, restore_bases_and_tags, ReadFeature};
use crate::record_flags::{MF_MATE_NEG_STRAND, MF_MATE_UNMAPPED};
use crate::record_read::{read_record, CramRecord, RecordReadError, RecordReaders, SliceContext};
use crate::substitution_matrix::SubstitutionMatrix;
use htsjdk_bam::cigar::Cigar;
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue, Tags};
use md5::{Digest, Md5};
use std::collections::HashMap;
use std::sync::Arc;

/// `CRAMCompressionRecord.DEFAULT_QUALITY_SCORE`, `'?' - '!'`.
const DEFAULT_QUALITY_SCORE: u8 = b'?' - b'!';
/// `CRAMCompressionRecord.MISSING_QUALITY_SCORE`, `-1`.
const MISSING_QUALITY_SCORE: u8 = 0xFF;
/// `CRAMCompressionRecord.NO_READGROUP_ID`.
const NO_READGROUP_ID: i32 = -1;

/// Where reference bases come from: `CRAMReferenceSource.getReferenceBases(record, false)`.
pub trait ReferenceSource {
    /// The whole sequence, uppercased, or `None` when the source has no such sequence.
    fn sequence_bases(&mut self, name: &str) -> Option<Arc<Vec<u8>>>;

    /// What the source throws instead of answering, for one that never has a sequence:
    /// `CRAMLazyReferenceSource`'s `IllegalArgumentException`. `None` for a source that answers
    /// `null` and leaves the caller to refuse.
    fn refusal(&self, _name: &str) -> Option<String> {
        None
    }
}

/// `CRAMLazyReferenceSource`, which `ReferenceSource.getDefaultCRAMReferenceSource` gives a
/// writer made with no reference: every request throws.
pub struct LazyReference;

impl ReferenceSource for LazyReference {
    fn sequence_bases(&mut self, _name: &str) -> Option<Arc<Vec<u8>>> {
        None
    }

    fn refusal(&self, name: &str) -> Option<String> {
        Some(format!(
            "A reference must be supplied that includes the reference sequence for {name}."
        ))
    }
}

/// `CRAMLazyReferenceSource`: what a reader without a reference gets, and it has nothing.
pub struct NoReference;

impl ReferenceSource for NoReference {
    fn sequence_bases(&mut self, _name: &str) -> Option<Arc<Vec<u8>>> {
        None
    }
}

/// `ReferenceSource(File)` over an indexed FASTA.
///
/// `ReferenceSequenceFileFactory` hands back an unindexed file when there is no `.fai`, and
/// `findBasesByName` answers nothing for one, so a FASTA without its index finds no sequence.
pub struct FastaReference {
    fasta: Option<htsjdk_bam::fasta_index::IndexedFasta<htsjdk_bam::fasta_index::PositionalFile>>,
    cache: HashMap<String, Arc<Vec<u8>>>,
}

impl FastaReference {
    pub fn open(path: &std::path::Path) -> Self {
        Self {
            fasta: htsjdk_bam::fasta_index::IndexedFasta::open(path).ok(),
            cache: HashMap::new(),
        }
    }
}

impl ReferenceSource for FastaReference {
    fn sequence_bases(&mut self, name: &str) -> Option<Arc<Vec<u8>>> {
        if let Some(bases) = self.cache.get(name) {
            return Some(Arc::clone(bases));
        }
        let fasta = self.fasta.as_mut()?;
        let size = fasta.index().get(name)?.size as i64;
        // `getSequence`, which throws for a contig the index lacks; `findBasesByName` swallows it.
        let mut bases = fasta.query(name, 1, size).ok()?;
        // `addToCache` uppercases in place, and every later fetch is of the uppercased copy.
        bases.make_ascii_uppercase();
        let bases = Arc::new(bases);
        self.cache.insert(name.to_string(), Arc::clone(&bases));
        Some(bases)
    }
}

/// What reading refuses, with the Java class that carries it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CramError {
    File(FileError),
    Record(RecordReadError),
    /// `IllegalArgumentException` from `CRAMReferenceRegion` and `ReferenceSource`.
    IllegalArgument(String),
    /// `CRAMException`.
    Cram(String),
    /// Anything else the reference throws on the way to a record, by class and message.
    Other {
        class: &'static str,
        message: String,
    },
}

impl CramError {
    /// The fully qualified class of what the reference throws.
    pub fn java_exception(&self) -> &'static str {
        match self {
            CramError::File(error) => error.java_exception(),
            CramError::Record(error) => error.java_exception(),
            CramError::IllegalArgument(_) => "java.lang.IllegalArgumentException",
            CramError::Cram(_) => "htsjdk.samtools.cram.CRAMException",
            CramError::Other { class, .. } => class,
        }
    }

    pub fn message(&self) -> String {
        match self {
            CramError::File(error) => error.message(),
            CramError::Record(error) => error.message(),
            CramError::IllegalArgument(message) | CramError::Cram(message) => message.clone(),
            CramError::Other { message, .. } => message.clone(),
        }
    }

    /// `Exception in thread "main"`'s line: the class, then the message.
    pub fn java_line(&self) -> String {
        format!("{}: {}", self.java_exception(), self.message())
    }
}

impl From<FileError> for CramError {
    fn from(error: FileError) -> Self {
        CramError::File(error)
    }
}

impl From<RecordReadError> for CramError {
    fn from(error: RecordReadError) -> Self {
        CramError::Record(error)
    }
}

/// `CRAMReferenceRegion`: the bases of one stretch of one contig.
#[derive(Default)]
struct ReferenceRegion {
    reference_index: Option<i32>,
    bases: Vec<u8>,
    /// Zero-based.
    start: i32,
    length: i32,
    valid: bool,
}

impl ReferenceRegion {
    /// `fetchReferenceBasesByRegion(referenceIndex, zeroBasedStart, length)`.
    fn fetch(
        &mut self,
        header: &SamHeader,
        source: &mut dyn ReferenceSource,
        reference_index: i32,
        zero_based_start: i32,
        length: i32,
    ) -> Result<(), CramError> {
        if zero_based_start < 0 {
            return Err(CramError::IllegalArgument("start must be >= 0".into()));
        }
        if self.valid
            && self.reference_index == Some(reference_index)
            && zero_based_start == self.start
            && length == self.length
        {
            return Ok(());
        }
        let sequence = usize::try_from(reference_index)
            .ok()
            .and_then(|index| header.sequences.get(index))
            .ok_or_else(|| {
                CramError::IllegalArgument(format!(
                    "Reference sequence index {reference_index} not found"
                ))
            })?;
        self.reference_index = Some(reference_index);
        if zero_based_start >= sequence.length {
            return Err(CramError::IllegalArgument(format!(
                "Requested start {zero_based_start} is beyond the sequence length {}",
                sequence.name
            )));
        }
        if let Some(message) = source.refusal(&sequence.name) {
            return Err(CramError::IllegalArgument(message));
        }
        let Some(all) = source.sequence_bases(&sequence.name) else {
            self.valid = false;
            return Err(CramError::IllegalArgument(format!(
                "Failure getting reference bases for sequence {}",
                sequence.name
            )));
        };
        // `ReferenceSource.getReferenceBasesByRegion`, against the bases rather than the
        // dictionary's length.
        if zero_based_start as usize >= all.len() {
            return Err(CramError::IllegalArgument(format!(
                "Requested start {zero_based_start} is beyond the sequence length {}",
                sequence.name
            )));
        }
        let end = all
            .len()
            .min(zero_based_start as usize + length.max(0) as usize);
        self.bases = all[zero_based_start as usize..end].to_vec();
        self.start = zero_based_start;
        self.length = self.bases.len() as i32;
        self.valid = true;
        Ok(())
    }
}

/// A CRAM being read: the header at once, the records a container at a time.
pub struct CramReader<'a> {
    cram: &'a [u8],
    major: u8,
    at: usize,
    finished: bool,
    /// The SAM header text as the file holds it.
    pub header_text: String,
    pub header: SamHeader,
    region: ReferenceRegion,
}

impl<'a> CramReader<'a> {
    /// `CRAMFileReader`'s constructor: the definition and the SAM header container.
    pub fn new(cram: &'a [u8]) -> Result<Self, CramError> {
        let (definition, block, at) = read_preamble(cram)?;
        let header_text = header_text(&block)?;
        let header = htsjdk_bam::reader::parse_header_text(&header_text);
        Ok(Self {
            cram,
            major: definition.major,
            at,
            finished: false,
            header_text,
            header,
            region: ReferenceRegion::default(),
        })
    }

    /// The next container's records, or `None` after the EOF container or the end of the bytes.
    pub fn next_container(
        &mut self,
        reference: &mut dyn ReferenceSource,
    ) -> Result<Option<Vec<BamRecord>>, CramError> {
        loop {
            if self.finished || self.at >= self.cram.len() {
                self.finished = true;
                return Ok(None);
            }
            let (container, end) = read_container(self.cram, self.at, self.major)?;
            self.at = end;
            if container.header.is_eof() && container.slices.is_empty() {
                self.finished = true;
                return Ok(None);
            }
            if container.slices.is_empty() {
                continue;
            }
            return self.container_records(&container, reference).map(Some);
        }
    }

    /// Every record, in file order: the whole-file read the tools do.
    pub fn read_all(
        &mut self,
        reference: &mut dyn ReferenceSource,
    ) -> Result<Vec<BamRecord>, CramError> {
        let mut records = Vec::new();
        while let Some(batch) = self.next_container(reference)? {
            records.extend(batch);
        }
        Ok(records)
    }

    /// `Container.getSAMRecords`.
    fn container_records(
        &mut self,
        container: &ContainerWalk,
        reference: &mut dyn ReferenceSource,
    ) -> Result<Vec<BamRecord>, CramError> {
        let compression = container
            .compression_header
            .as_ref()
            .expect("a data container has a compression header");
        let readers = RecordReaders::new(&compression.encodings)?;
        let matrix = SubstitutionMatrix::from_encoded(compression.preservation.substitution_matrix);
        let mut out = Vec::new();

        for slice in &container.slices {
            let context = SliceContext {
                reference_context: slice.header.reference_context_id,
                alignment_start: slice.header.alignment_start,
            };
            let blocks = slice.block_bytes();
            let mut streams = crate::external_codecs::SliceReadStreams::new(&blocks);

            // `deserializeCRAMRecords`.
            let mut records = Vec::with_capacity(slice.header.record_count.max(0) as usize);
            let mut previous = slice.header.alignment_start;
            for _ in 0..slice.header.record_count.max(0) {
                let record = read_record(&readers, compression, &context, &mut streams, previous)?;
                previous = record.alignment_start;
                records.push(record);
            }

            // `normalizeCRAMRecords`.
            let mut embedded = false;
            if compression.preservation.reference_required {
                if context.reference_context >= 0 {
                    self.region.fetch(
                        &self.header,
                        reference,
                        context.reference_context,
                        slice.header.alignment_start - 1,
                        slice.header.alignment_span,
                    )?;
                    validate_md5(
                        &self.region,
                        slice.header.alignment_span,
                        &slice.header.reference_md5,
                    )?;
                }
            } else if slice.header.embedded_reference_content_id >= 0 {
                if let Some(bases) = blocks
                    .external
                    .get(&slice.header.embedded_reference_content_id)
                {
                    embedded = true;
                    self.region.reference_index = Some(context.reference_context);
                    self.region.bases = bases.clone();
                    self.region.start = slice.header.alignment_start - 1;
                    self.region.length = bases.len() as i32;
                    self.region.valid = true;
                }
            }

            let mut normal = normalize(
                records,
                slice.header.global_record_counter,
                compression.preservation.reference_required
                    && context.is_multi_reference()
                    && !embedded,
                &matrix,
                &mut self.region,
                &self.header,
                reference,
            )?;
            for (record, cigar) in normal.drain(..) {
                out.push(to_bam_record(record, cigar, &self.header)?);
            }
        }
        Ok(out)
    }
}

/// The FILE_HEADER block's content: a little-endian length, then that many bytes of text.
fn header_text(block: &[u8]) -> Result<String, CramError> {
    let length = block
        .get(..4)
        .map(|b| i32::from_le_bytes([b[0], b[1], b[2], b[3]]))
        .ok_or(CramError::File(FileError::Truncated))?;
    let text = block
        .get(4..4 + length.max(0) as usize)
        .ok_or(CramError::File(FileError::Truncated))?;
    Ok(text.iter().map(|&b| b as char).collect())
}

/// `Slice.validateReferenceBases`, for a single-reference slice whose reference is required.
fn validate_md5(region: &ReferenceRegion, span: i32, expected: &[u8; 16]) -> Result<(), CramError> {
    let span = (span.max(0) as usize).min(region.bases.len());
    let actual: [u8; 16] = Md5::digest(&region.bases[..span]).into();
    if &actual != expected {
        let hex: String = expected.iter().map(|b| format!("{b:02x}")).collect();
        return Err(CramError::Cram(format!(
            "The MD5 for the reference failed to validate against the expected value {hex}. \
             This indicates that the supplied reference is not the one originally used to create the CRAM.."
        )));
    }
    Ok(())
}

/// A raw record and what normalisation adds to it.
struct Normal {
    record: CramRecord,
    /// The cigar the fused restore built, which `toSAMRecord` prefers to rebuilding one.
    cigar: Option<Cigar>,
    sequential_index: i64,
    next: Option<usize>,
    previous: Option<usize>,
    alignment_end: i32,
}

/// `Slice.normalizeCRAMRecords`, after the slice's reference region is in place.
#[allow(clippy::too_many_arguments)]
fn normalize(
    records: Vec<CramRecord>,
    global_record_counter: i64,
    per_record_reference: bool,
    matrix: &SubstitutionMatrix,
    region: &mut ReferenceRegion,
    header: &SamHeader,
    reference: &mut dyn ReferenceSource,
) -> Result<Vec<(CramRecord, Option<Cigar>)>, CramError> {
    let mut normal: Vec<Normal> = records
        .into_iter()
        .enumerate()
        .map(|(i, record)| {
            let end = if record.alignment_start != 0 {
                alignment_end(
                    &record.read_features,
                    record.alignment_start,
                    record.read_length,
                )
            } else {
                0
            };
            Normal {
                record,
                cigar: None,
                sequential_index: global_record_counter + i as i64,
                next: None,
                previous: None,
                alignment_end: end,
            }
        })
        .collect();

    // Mate pairing first.
    for i in 0..normal.len() {
        let flags = normal[i].record.flags;
        if flags.is_read_paired() && !flags.is_detached() && flags.has_mate_downstream() {
            let target = normal[i].sequential_index
                + i64::from(normal[i].record.records_to_next_fragment)
                + 1
                - global_record_counter;
            let target = usize::try_from(target)
                .ok()
                .filter(|&t| t < normal.len())
                .ok_or_else(|| CramError::Other {
                    class: "java.lang.IndexOutOfBoundsException",
                    message: format!("Index {target} out of bounds for length {}", normal.len()),
                })?;
            normal[i].next = Some(target);
            normal[target].previous = Some(i);
        }
    }
    for i in 0..normal.len() {
        if normal[i].previous.is_none() && normal[i].next.is_some() {
            restore_mate_info(&mut normal, i);
        }
    }

    // Generated names, propagated to the neighbours.
    for i in 0..normal.len() {
        if normal[i].record.read_name.is_none() {
            let name = normal[i].sequential_index.to_string();
            if let Some(next) = normal[i].next {
                normal[next].record.read_name = Some(name.clone());
            }
            if let Some(previous) = normal[i].previous {
                normal[previous].record.read_name = Some(name.clone());
            }
            normal[i].record.read_name = Some(name);
        }
    }

    // Bases, the cigar, and MD and NM where the record does not carry its own.
    for item in normal.iter_mut() {
        let record = &mut item.record;
        if record.flags.is_segment_unmapped() {
            continue;
        }
        if per_record_reference && !record.flags.is_unknown_bases() {
            region.fetch(
                header,
                reference,
                record.reference_index,
                record.alignment_start - 1,
                item.alignment_end - record.alignment_start + 1,
            )?;
        }
        // htslib's `cF` tag, written by `embed_ref=2`, says which of the two not to generate.
        let mut suppress_md = false;
        let mut suppress_nm = false;
        if let Some(i) = record
            .tags
            .iter()
            .rposition(|(id, _)| id[0] == b'c' && id[1] == b'F')
        {
            let (id, bytes) = record.tags.remove(i);
            let cf = match tag_value(id[2], &bytes) {
                Ok(TagValue::Int(v)) => v as i32,
                _ => 0,
            };
            suppress_md = cf & 1 != 0;
            suppress_nm = cf & 2 != 0;
        }
        let has = |name: &[u8; 2]| {
            record
                .tags
                .iter()
                .any(|(id, _)| id[0] == name[0] && id[1] == name[1])
        };
        let need_md = !has(b"MD") && !suppress_md;
        let need_nm = !has(b"NM") && !suppress_nm;
        let compute = (need_md || need_nm) && record.read_length > 0 && region.valid;

        let result = restore_bases_and_tags(
            &record.read_features,
            record.flags.is_unknown_bases(),
            record.alignment_start,
            record.read_length,
            region.valid.then_some(region.bases.as_slice()),
            region.start,
            matrix,
            compute,
        )
        .map_err(|error| CramError::Other {
            class: error.java_exception(),
            message: error.message(),
        })?;
        record.read_bases = result.read_bases;
        item.cigar = Some(result.cigar);
        if compute {
            if need_md {
                if let Some(md) = result.md {
                    let mut bytes = md.into_bytes();
                    bytes.push(0);
                    record.tags.push((*b"MDZ", bytes));
                }
            }
            if need_nm {
                if let Some(nm) = result.nm.filter(|&n| n >= 0) {
                    record.tags.push(integer_tag(*b"NM", nm));
                }
            }
        }
    }

    // Scores.
    for item in normal.iter_mut() {
        resolve_quality_scores(&mut item.record);
    }

    Ok(normal.into_iter().map(|n| (n.record, n.cigar)).collect())
}

/// `CRAMCompressionRecord.restoreMateInfo`, from the first segment of a chain.
fn restore_mate_info(normal: &mut [Normal], first: usize) {
    let mut current = first;
    while let Some(next) = normal[current].next {
        set_next_mate(normal, current, next);
        current = next;
    }
    let last = current;
    set_next_mate(normal, last, first);

    let template_length = compute_insert_size(&normal[first], &normal[last]);
    normal[first].record.template_size = template_length;
    normal[last].record.template_size = template_length.wrapping_neg();
}

/// `setNextMate`.
fn set_next_mate(normal: &mut [Normal], of: usize, next: usize) {
    let (start, unmapped, negative, reference) = {
        let n = &normal[next].record;
        (
            n.alignment_start,
            n.flags.is_segment_unmapped(),
            n.flags.is_negative_strand(),
            n.reference_index,
        )
    };
    let record = &mut normal[of].record;
    record.mate_alignment_start = start;
    set_bit(&mut record.flags.mate, MF_MATE_UNMAPPED, unmapped);
    set_bit(&mut record.flags.mate, MF_MATE_NEG_STRAND, negative);
    record.mate_reference_index = reference;
    if reference == -1 {
        record.mate_alignment_start = 0;
    }
}

fn set_bit(word: &mut i32, bit: i32, on: bool) {
    if on {
        *word |= bit;
    } else {
        *word &= !bit;
    }
}

/// `CRAMCompressionRecord.computeInsertSize` at htsjdk 5.0.0, which is htslib's: leftmost
/// mapped base to rightmost across both mates, positive for the leftmost read and ties broken by
/// first-of-pair.
fn compute_insert_size(first: &Normal, second: &Normal) -> i32 {
    let (a, b) = (&first.record, &second.record);
    if a.flags.is_segment_unmapped()
        || b.flags.is_segment_unmapped()
        || a.reference_index != b.reference_index
    {
        return 0;
    }
    let left = a.alignment_start.min(b.alignment_start);
    let right = first.alignment_end.max(second.alignment_end);
    let magnitude = right.wrapping_sub(left).wrapping_add(1);
    if a.alignment_start < b.alignment_start {
        return magnitude;
    }
    if a.alignment_start > b.alignment_start {
        return magnitude.wrapping_neg();
    }
    if a.flags.is_first_segment() {
        return magnitude;
    }
    magnitude.wrapping_neg()
}

/// A tag written as `ReadTag.deriveTypeFromValue` types an integer: the smallest type that holds
/// it, signed before unsigned at each width.
fn integer_tag(name: [u8; 2], value: i32) -> ([u8; 3], Vec<u8>) {
    let v = i64::from(value);
    let (kind, bytes) = if v > 0xFFFF {
        (b'i', value.to_le_bytes().to_vec())
    } else if v > i64::from(i16::MAX) {
        (b'S', (value as u16).to_le_bytes().to_vec())
    } else if v > 0xFF {
        (b's', (value as i16).to_le_bytes().to_vec())
    } else if v > i64::from(i8::MAX) {
        (b'C', vec![value as u8])
    } else if v >= i64::from(i8::MIN) {
        (b'c', vec![value as i8 as u8])
    } else if v >= i64::from(i16::MIN) {
        (b's', (value as i16).to_le_bytes().to_vec())
    } else {
        (b'i', value.to_le_bytes().to_vec())
    };
    ([name[0], name[1], kind], bytes)
}

/// `resolveQualityScores`.
fn resolve_quality_scores(record: &mut CramRecord) {
    let length = record.read_length.max(0) as usize;
    if !record.flags.is_force_preserve_quality_scores() {
        let mut missing = true;
        let mut scores = vec![DEFAULT_QUALITY_SCORE; length];
        for feature in &record.read_features {
            match feature {
                ReadFeature::BaseQualityScore { position, quality } => {
                    put(&mut scores, *position - 1, *quality as u8);
                    missing = false;
                }
                ReadFeature::Scores {
                    position,
                    scores: run,
                } => {
                    for (offset, score) in run.iter().enumerate() {
                        put(&mut scores, *position - 1 + offset as i32, *score);
                    }
                    missing = false;
                }
                ReadFeature::ReadBase {
                    position, quality, ..
                } => {
                    put(&mut scores, *position - 1, *quality as u8);
                    missing = false;
                }
                _ => {}
            }
        }
        record.quality_scores = if missing { Vec::new() } else { scores };
    } else {
        let mut missing = 0;
        for score in record.quality_scores.iter_mut() {
            if *score == MISSING_QUALITY_SCORE {
                *score = DEFAULT_QUALITY_SCORE;
                missing += 1;
            }
        }
        if missing == record.quality_scores.len() {
            record.quality_scores = Vec::new();
        }
    }
}

fn put(scores: &mut [u8], at: i32, value: u8) {
    if let Some(slot) = usize::try_from(at).ok().and_then(|i| scores.get_mut(i)) {
        *slot = value;
    }
}

/// `toSAMRecord`.
fn to_bam_record(
    record: CramRecord,
    cigar: Option<Cigar>,
    header: &SamHeader,
) -> Result<BamRecord, CramError> {
    let flags = record.flags;
    let mut bam_flags: u16 = 0;
    let mut set = |on: bool, bit: u16| {
        if on {
            bam_flags |= bit;
        }
    };
    let paired = flags.is_read_paired();
    set(paired, 0x1);
    set(flags.bam & 0x2 != 0, 0x2);
    set(flags.is_segment_unmapped(), 0x4);
    set(flags.is_negative_strand(), 0x10);
    set(flags.is_first_segment(), 0x40);
    set(flags.is_last_segment(), 0x80);
    set(flags.is_secondary_alignment(), 0x100);
    set(flags.bam & 0x200 != 0, 0x200);
    set(flags.bam & 0x400 != 0, 0x400);
    set(flags.bam & 0x800 != 0, 0x800);

    let mut out = BamRecord {
        read_name: record.read_name.unwrap_or_default(),
        ..BamRecord::default()
    };
    if record.reference_index == -1 {
        out.reference_index = -1;
        out.alignment_start = 0;
        // `SAMRecord.NO_MAPPING_QUALITY`, which is 0; 255 is the unknown one.
        out.mapping_quality = 0;
    } else {
        out.reference_index = record.reference_index;
        out.alignment_start = record.alignment_start;
        out.mapping_quality = record.mapping_quality as u8;
    }
    out.cigar = if flags.is_segment_unmapped() {
        Cigar::default()
    } else {
        cigar.unwrap_or_else(|| {
            crate::read_features::cigar_for_read_features(&record.read_features, record.read_length)
        })
    };
    if paired {
        out.mate_reference_index = record.mate_reference_index;
        out.mate_alignment_start = record.mate_alignment_start.max(0);
        if record.mate_alignment_start <= 0 {
            out.mate_alignment_start = 0;
        }
        if flags.mate & MF_MATE_NEG_STRAND != 0 {
            bam_flags |= 0x20;
        }
        if flags.mate & MF_MATE_UNMAPPED != 0 {
            bam_flags |= 0x8;
        }
    } else {
        out.mate_reference_index = -1;
        out.mate_alignment_start = 0;
    }
    out.flags = bam_flags;
    out.inferred_insert_size = record.template_size;
    out.read_bases = record.read_bases;
    out.base_qualities = record.quality_scores;

    let mut tags = Tags::new();
    for (id, bytes) in &record.tags {
        let value = tag_value(id[2], bytes).map_err(|message| CramError::Other {
            class: "htsjdk.samtools.SAMFormatException",
            message,
        })?;
        tags.insert(Tag::new(&[id[0], id[1]]), value);
    }
    if record.read_group_id != NO_READGROUP_ID {
        let group = usize::try_from(record.read_group_id)
            .ok()
            .and_then(|i| header.read_groups.get(i))
            .ok_or_else(|| CramError::Other {
                class: "java.lang.IndexOutOfBoundsException",
                message: format!(
                    "Index {} out of bounds for length {}",
                    record.read_group_id,
                    header.read_groups.len()
                ),
            })?;
        tags.insert(Tag::new(b"RG"), TagValue::Str(group.id.clone()));
    }
    out.tags = tags;
    Ok(out)
}

/// `ReadTag.readSingleValue`, then what `setAttribute` keeps of it: an array's unsigned flag is
/// lost, because `setAttribute(tag, byte[])` stores a signed array.
fn tag_value(kind: u8, bytes: &[u8]) -> Result<TagValue, String> {
    let int = |n: usize| -> Result<&[u8], String> {
        bytes
            .get(..n)
            .ok_or_else(|| "java.nio.BufferUnderflowException".to_string())
    };
    Ok(match kind {
        b'Z' => {
            let end = bytes.iter().position(|&b| b == 0).unwrap_or(bytes.len());
            TagValue::Str(bytes[..end].iter().map(|&b| b as char).collect())
        }
        b'A' => TagValue::Char(*int(1)?.first().unwrap()),
        b'I' => {
            let b = int(4)?;
            TagValue::Int(i64::from(u32::from_le_bytes([b[0], b[1], b[2], b[3]])))
        }
        b'i' => {
            let b = int(4)?;
            TagValue::Int(i64::from(i32::from_le_bytes([b[0], b[1], b[2], b[3]])))
        }
        b's' => {
            let b = int(2)?;
            TagValue::Int(i64::from(i16::from_le_bytes([b[0], b[1]])))
        }
        b'S' => {
            let b = int(2)?;
            TagValue::Int(i64::from(u16::from_le_bytes([b[0], b[1]])))
        }
        b'c' => TagValue::Int(i64::from(int(1)?[0] as i8)),
        b'C' => TagValue::Int(i64::from(int(1)?[0])),
        b'f' => {
            let b = int(4)?;
            TagValue::Float(f32::from_le_bytes([b[0], b[1], b[2], b[3]]))
        }
        b'H' => {
            let end = bytes.iter().position(|&b| b == 0).unwrap_or(bytes.len());
            let hex = &bytes[..end];
            let mut values = Vec::with_capacity(hex.len() / 2);
            for pair in hex.chunks(2) {
                let text = std::str::from_utf8(pair).map_err(|e| e.to_string())?;
                values.push(u8::from_str_radix(text, 16).map_err(|e| e.to_string())? as i8);
            }
            TagValue::ByteArray {
                values,
                unsigned: false,
            }
        }
        b'B' => {
            let kind = *int(1)?.first().unwrap();
            let count = {
                let b = bytes.get(1..5).ok_or("java.nio.BufferUnderflowException")?;
                i32::from_le_bytes([b[0], b[1], b[2], b[3]]).max(0) as usize
            };
            let body = &bytes[5..];
            let need = |width: usize| -> Result<&[u8], String> {
                body.get(..count * width)
                    .ok_or_else(|| "java.nio.BufferUnderflowException".to_string())
            };
            match kind.to_ascii_lowercase() {
                b'c' => TagValue::ByteArray {
                    values: need(1)?.iter().map(|&b| b as i8).collect(),
                    unsigned: false,
                },
                b's' => TagValue::ShortArray {
                    values: need(2)?
                        .chunks(2)
                        .map(|b| i16::from_le_bytes([b[0], b[1]]))
                        .collect(),
                    unsigned: false,
                },
                b'i' => TagValue::IntArray {
                    values: need(4)?
                        .chunks(4)
                        .map(|b| i32::from_le_bytes([b[0], b[1], b[2], b[3]]))
                        .collect(),
                    unsigned: false,
                },
                b'f' => TagValue::FloatArray(
                    need(4)?
                        .chunks(4)
                        .map(|b| f32::from_le_bytes([b[0], b[1], b[2], b[3]]))
                        .collect(),
                ),
                other => return Err(format!("Unrecognized tag array type: {}", other as char)),
            }
        }
        other => return Err(format!("Unrecognized tag type: {}", other as char)),
    })
}
