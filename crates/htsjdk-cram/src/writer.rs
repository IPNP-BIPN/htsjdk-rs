//! Writing SAM records as a CRAM 3.1 file, as htsjdk 5.0.0's `CRAMFileWriter` does with the
//! default `CRAMEncodingStrategy`, which is `CRAMCompressionProfile.NORMAL`.
//!
//! Ported from `htsjdk.samtools.CRAMContainerStreamWriter`, `ContainerFactory`, `SliceFactory`,
//! `CompressionHeaderFactory`, `CompressionHeaderEncodingMap`, `Slice`, `CramRecordWriter`,
//! `CRAMCompressionRecord(SAMRecord)` and `Container.write` at htsjdk 5.0.0.
//!
//! # What NORMAL means
//!
//! CRAM 3.1, 10,000 reads or 5,000,000 bases a slice, one slice a container. rANS Nx16 order 1 for
//! the small integer series, order 0 for `AP` and `RI`, fqzcomp for the scores, the name tokeniser
//! for the names, GZIP level 5 for everything else. A tag's block is GZIP or rANS order 0,
//! whichever a running trial of three blocks finds smaller, re-tried every seventy.
//!
//! # Three things htsjdk does that a reader of the specification would not
//!
//! * Every slice header's record counter is zero: `SliceFactory` computes the next counter and
//!   throws it away, so each slice's records count from zero again.
//! * A tag's external block, once created, stays in every later container's slices, empty: the
//!   encoding map that holds the tag compressors lives as long as the writer.
//! * `MD` and `NM` are dropped from a mapped record when they equal what a reader would regenerate,
//!   and kept verbatim when they do not; with no bases or no cigar they are dropped unchecked.

use crate::block::{CompressionMethod, ContentType};
use crate::compression_header::{crc32, raw_block, raw_block_with_id, CompressionHeader};
use crate::encoding_factory::{create_encoding, DataSeriesType as FactoryType, Encoding};
use crate::encoding_map::{DataSeries, EncodingDescriptor, EncodingId, EncodingMap};
use crate::external_codecs::SliceWriteStreams;
use crate::preservation_map::PreservationMap;
use crate::read_features::{alignment_end, create_read_features, Qualities, ReadFeature, NO_CODE};
use crate::reader::{CramError, ReferenceSource};
use crate::record_flags::{MF_MATE_NEG_STRAND, MF_MATE_UNMAPPED};
use crate::slice_header::SliceHeader;
use crate::substitution_matrix::{SubstitutionMatrix, SYMBOL_SPACE_SIZE};
use crate::tag_encoding_map::TagEncodingMap;
use crate::varint::{write_unsigned_itf8, write_unsigned_ltf8};
use htsjdk_bam::cigar::{Cigar, Op};
use htsjdk_bam::header::SamHeader;
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::TagValue;
use md5::{Digest, Md5};
use std::collections::{BTreeSet, HashMap};
use std::sync::Arc;

const MAJOR: u8 = 3;
const MINOR: u8 = 1;
const READS_PER_SLICE: usize = 10_000;
const BASES_PER_SLICE: i64 = READS_PER_SLICE as i64 * 500;
const SLICES_PER_CONTAINER: usize = 1;
const MINIMUM_SINGLE_REFERENCE_SLICE_SIZE: usize = 1000;
const GZIP_LEVEL: i32 = 5;

const UNINITIALIZED_REFERENCE_ID: i32 = -3;
const MULTIPLE_REFERENCE_ID: i32 = -2;
const UNMAPPED_UNPLACED_ID: i32 = -1;

const CF_QS_PRESERVED_AS_ARRAY: i32 = 0x1;
const CF_DETACHED: i32 = 0x2;
const CF_HAS_MATE_DOWNSTREAM: i32 = 0x4;
const CF_UNKNOWN_BASES: i32 = 0x8;
const MISSING_QUALITY_SCORE: u8 = 0xFF;
const NO_READGROUP_ID: i32 = -1;

/// `CRAMIO.ZERO_F_EOF_MARKER`: the version 3 EOF container, byte for byte.
const EOF_V3: [u8; 38] = [
    0x0f, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x0f, 0xe0, 0x45, 0x4f, 0x46, 0x00, 0x00, 0x00,
    0x00, 0x01, 0x00, 0x05, 0xbd, 0xd9, 0x4f, 0x00, 0x01, 0x00, 0x06, 0x06, 0x01, 0x00, 0x01, 0x00,
    0x01, 0x00, 0xee, 0x63, 0x01, 0x4b,
];

/// How a data series' block is compressed under NORMAL.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Compressor {
    Gzip,
    Rans0,
    Rans1,
    Fqzcomp,
    NameTokeniser,
}

/// The series NORMAL writes, in `DataSeries` declaration order, with htsjdk 5.0.0's content ids.
const WRITTEN: [(&str, FactoryType, i32, Compressor); 26] = [
    ("BF", FactoryType::Int, 15, Compressor::Rans1),
    ("CF", FactoryType::Int, 16, Compressor::Rans1),
    ("RI", FactoryType::Int, 33, Compressor::Rans0),
    ("RL", FactoryType::Int, 25, Compressor::Rans1),
    ("AP", FactoryType::Int, 17, Compressor::Rans0),
    ("RG", FactoryType::Int, 18, Compressor::Rans1),
    ("RN", FactoryType::ByteArray, 11, Compressor::NameTokeniser),
    ("NF", FactoryType::Int, 24, Compressor::Gzip),
    ("MF", FactoryType::Int, 21, Compressor::Rans1),
    ("NS", FactoryType::Int, 20, Compressor::Rans1),
    ("NP", FactoryType::Int, 23, Compressor::Gzip),
    ("TS", FactoryType::Int, 22, Compressor::Rans1),
    ("TL", FactoryType::Int, 32, Compressor::Rans1),
    ("MQ", FactoryType::Int, 19, Compressor::Rans1),
    ("FN", FactoryType::Int, 26, Compressor::Rans1),
    ("FP", FactoryType::Int, 28, Compressor::Gzip),
    ("FC", FactoryType::Byte, 27, Compressor::Rans1),
    ("BA", FactoryType::Byte, 30, Compressor::Rans1),
    ("QS", FactoryType::Byte, 12, Compressor::Fqzcomp),
    ("BS", FactoryType::Byte, 31, Compressor::Rans1),
    ("IN", FactoryType::ByteArray, 13, Compressor::Gzip),
    ("DL", FactoryType::Int, 29, Compressor::Gzip),
    ("RS", FactoryType::Int, 34, Compressor::Gzip),
    ("SC", FactoryType::ByteArray, 14, Compressor::Gzip),
    ("PD", FactoryType::Int, 35, Compressor::Gzip),
    ("HC", FactoryType::Int, 36, Compressor::Gzip),
];

fn itf8(value: i32) -> Vec<u8> {
    write_unsigned_itf8(value).0
}

/// The encoding descriptor `CompressionHeaderEncodingMap` gives a series.
fn descriptor(name: &str, content_id: i32) -> EncodingDescriptor {
    match name {
        "RN" | "IN" | "SC" => {
            let stop = if name == "RN" { 0u8 } else { b'\t' };
            let mut parameters = vec![stop];
            parameters.extend_from_slice(&itf8(content_id));
            EncodingDescriptor {
                id: EncodingId::ByteArrayStop,
                parameters,
            }
        }
        _ => EncodingDescriptor {
            id: EncodingId::External,
            parameters: itf8(content_id),
        },
    }
}

/// `TrialCompressor` over GZIP and rANS order 0, one per tag for the life of the writer.
struct TrialCompressor {
    sizes: [i64; 2],
    winner: Option<usize>,
    trials_remaining: i32,
    until_retrial: i32,
}

impl TrialCompressor {
    const TRIALS: i32 = 3;
    const SPAN: i32 = 70;

    fn new() -> Self {
        Self {
            sizes: [0; 2],
            winner: None,
            trials_remaining: Self::TRIALS,
            until_retrial: 0,
        }
    }

    fn candidate(index: usize, data: &[u8]) -> (Vec<u8>, CompressionMethod) {
        if index == 0 {
            (
                crate::gzip::compress(data, GZIP_LEVEL),
                CompressionMethod::Gzip,
            )
        } else {
            (
                crate::rans_nx16::compress(data, 0),
                CompressionMethod::RansNx16,
            )
        }
    }

    fn compress(&mut self, data: &[u8]) -> (Vec<u8>, CompressionMethod) {
        if data.is_empty() {
            return Self::candidate(self.winner.unwrap_or(0), data);
        }
        if let Some(winner) = self.winner {
            if self.until_retrial > 0 {
                self.until_retrial -= 1;
                return Self::candidate(winner, data);
            }
            self.trials_remaining = Self::TRIALS;
            self.winner = None;
            for size in self.sizes.iter_mut() {
                *size /= 2;
            }
        }
        let results = [Self::candidate(0, data), Self::candidate(1, data)];
        for (size, result) in self.sizes.iter_mut().zip(&results) {
            *size += result.0.len() as i64;
        }
        self.trials_remaining -= 1;
        if self.trials_remaining == 0 {
            let mut best = i64::MAX;
            let mut winner = 0;
            for (i, &size) in self.sizes.iter().enumerate() {
                if size < best {
                    best = size;
                    winner = i;
                }
            }
            self.winner = Some(winner);
            self.until_retrial = Self::SPAN;
            return results[winner].clone();
        }
        let best = if results[1].0.len() < results[0].0.len() {
            1
        } else {
            0
        };
        results[best].clone()
    }
}

/// One tag of a record as `ReadTag` holds it: the three-byte id and the value's bytes.
#[derive(Clone, Debug)]
struct Tag {
    id: [u8; 3],
    value: Vec<u8>,
    /// `getTagValueByteSize`, which differs from the value's length for nothing but is computed
    /// from the value rather than measured.
    size: usize,
}

impl Tag {
    fn key(&self) -> i32 {
        (i32::from(self.id[0]) << 16) | (i32::from(self.id[1]) << 8) | i32::from(self.id[2])
    }
}

/// A record on its way out: `CRAMCompressionRecord` built from a `SAMRecord`.
#[derive(Clone, Debug)]
struct WriteRecord {
    bam_flags: i32,
    cram_flags: i32,
    mate_flags: i32,
    read_name: Vec<u8>,
    reference_index: i32,
    alignment_start: i32,
    alignment_end: i32,
    read_length: i32,
    read_bases: Vec<u8>,
    features: Vec<ReadFeature>,
    template_size: i32,
    mapping_quality: i32,
    mate_alignment_start: i32,
    mate_reference_index: i32,
    quality_scores: Vec<u8>,
    read_group_id: i32,
    tags: Vec<Tag>,
    records_to_next_fragment: i32,
    tag_list: usize,
}

impl WriteRecord {
    fn is_placed(&self) -> bool {
        self.alignment_start != 0
    }

    fn is_unmapped(&self) -> bool {
        self.bam_flags & 0x4 != 0
    }
}

/// `CRAMReferenceRegion` as the writer uses it: whole contigs, by index.
struct WriterReference<'a> {
    source: &'a mut dyn ReferenceSource,
    current: Option<(i32, Arc<Vec<u8>>)>,
}

impl WriterReference<'_> {
    /// `fetchReferenceBases(referenceIndex)`, which tries `chr`-prefix variants of the name.
    fn contig(&mut self, header: &SamHeader, index: i32) -> Result<Arc<Vec<u8>>, CramError> {
        if let Some((current, bases)) = &self.current {
            if *current == index {
                return Ok(Arc::clone(bases));
            }
        }
        let sequence = usize::try_from(index)
            .ok()
            .and_then(|i| header.sequences.get(i))
            .ok_or_else(|| {
                CramError::IllegalArgument(format!("Reference sequence index {index} not found"))
            })?;
        let mut found = self.source.sequence_bases(&sequence.name);
        if found.is_none() {
            for variant in name_variants(&sequence.name) {
                found = self.source.sequence_bases(&variant);
                if found.is_some() {
                    break;
                }
            }
        }
        let Some(bases) = found else {
            let assembly = sequence
                .attributes
                .get("AS")
                .map_or("null".to_string(), str::to_string);
            let alternates = sequence.attributes.get("AN").map_or("[]".to_string(), |a| {
                format!("[{}]", a.split(',').collect::<Vec<_>>().join(", "))
            });
            return Err(CramError::IllegalArgument(format!(
                "A reference must be supplied (reference sequence SAMSequenceRecord(name={},length={},dict_index={},assembly={},alternate_names={}) not found).",
                sequence.name, sequence.length, index, assembly, alternates
            )));
        };
        self.current = Some((index, Arc::clone(&bases)));
        Ok(bases)
    }
}

/// `ReferenceSource.getVariants`.
fn name_variants(name: &str) -> Vec<String> {
    let mut variants = Vec::new();
    if name == "M" {
        variants.push("MT".to_string());
    }
    if name == "MT" {
        variants.push("M".to_string());
    }
    if name.len() >= 3 && name[..3].eq_ignore_ascii_case("chr") {
        variants.push(name[3..].to_string());
    } else {
        variants.push(format!("chr{name}"));
    }
    if name == "chrM" {
        variants.push("MT".to_string());
    }
    variants
}

/// `SequenceUtil.bases`: an IUPAC code's bit mask, and 0 for anything else.
fn iupac(byte: u8) -> u8 {
    let upper = byte.to_ascii_uppercase();
    match upper {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' => 8,
        b'M' => 3,
        b'R' => 5,
        b'W' => 9,
        b'S' => 6,
        b'Y' => 10,
        b'K' => 12,
        b'V' => 7,
        b'H' => 11,
        b'D' => 13,
        b'B' => 14,
        b'N' => 15,
        _ if byte == b'.' => 15,
        _ => 0,
    }
}

/// `SequenceUtil.calculateMdAndNm` over a whole contig.
fn calculate_md_nm(
    cigar: &Cigar,
    read: &[u8],
    reference: &[u8],
    alignment_start: i32,
) -> (String, i32) {
    let mut block_ref = alignment_start as i64 - 1;
    let mut block_read = 0usize;
    let mut matches = 0;
    let mut nm = 0;
    let mut md = String::new();
    for element in &cigar.elements {
        let length = element.length as i64;
        match element.op {
            Op::M | Op::Eq | Op::X => {
                let mut offset = 0;
                while offset < length {
                    let reference_index = block_ref + offset;
                    if reference_index >= reference.len() as i64 {
                        break;
                    }
                    let read_base = read.get(block_read + offset as usize).copied().unwrap_or(0);
                    let reference_base = reference[reference_index.max(0) as usize];
                    if iupac(read_base) == iupac(reference_base) || read_base == 0 {
                        matches += 1;
                    } else {
                        md.push_str(&matches.to_string());
                        md.push(reference_base as char);
                        matches = 0;
                        nm += 1;
                    }
                    offset += 1;
                }
                if offset < length {
                    break;
                }
                block_ref += length;
                block_read += length as usize;
            }
            Op::D => {
                md.push_str(&matches.to_string());
                md.push('^');
                let mut offset = 0;
                while offset < length {
                    let reference_index = block_ref + offset;
                    if reference_index >= reference.len() as i64
                        || reference[reference_index.max(0) as usize] == 0
                    {
                        break;
                    }
                    md.push(reference[reference_index as usize] as char);
                    offset += 1;
                }
                matches = 0;
                if offset < length {
                    break;
                }
                block_ref += length;
                nm += length as i32;
            }
            Op::I | Op::S => {
                block_read += length as usize;
                if element.op == Op::I {
                    nm += length as i32;
                }
            }
            Op::N => block_ref += length,
            _ => {}
        }
    }
    md.push_str(&matches.to_string());
    (md, nm)
}

/// `SequenceUtil.toBamReadBasesInPlace`.
fn to_bam_read_bases(bases: &[u8]) -> Vec<u8> {
    const SET: &[u8] = b"=ABCDGHKMNRSTVWY";
    bases
        .iter()
        .map(|&b| {
            for &base in SET {
                if b == base || b == base + 32 {
                    return base;
                }
            }
            b'N'
        })
        .collect()
}

/// `ReadTag.getIntegerType`.
fn integer_type(value: i64) -> Result<u8, CramError> {
    const MAX_INT: i64 = i32::MAX as i64;
    if value > 2 * MAX_INT + 1 {
        return Err(CramError::IllegalArgument(format!(
            "Integer attribute value too large: {value}"
        )));
    }
    Ok(if value > MAX_INT {
        b'I'
    } else if value > 65535 {
        b'i'
    } else if value > 32767 {
        b'S'
    } else if value > 255 {
        b's'
    } else if value > 127 {
        b'C'
    } else if value >= -128 {
        b'c'
    } else if value >= -32768 {
        b's'
    } else if value >= i64::from(i32::MIN) {
        b'i'
    } else {
        return Err(CramError::IllegalArgument(
            "Integer attribute value too negative to be encoded in BAM".into(),
        ));
    })
}

/// `ReadTag.deriveTypeFromValue` and `writeSingleValue`: the type letter and the bytes.
fn tag_of(name: [u8; 2], value: &TagValue) -> Result<Tag, CramError> {
    let (kind, bytes, size): (u8, Vec<u8>, usize) = match value {
        TagValue::Str(s) => {
            let mut bytes: Vec<u8> = s
                .chars()
                .map(|c| if c.is_ascii() { c as u8 } else { b'?' })
                .collect();
            let size = s.encode_utf16().count() + 1;
            bytes.push(0);
            (b'Z', bytes, size)
        }
        TagValue::Char(c) => (b'A', vec![*c], 1),
        TagValue::Float(f) => (b'f', f.to_le_bytes().to_vec(), 4),
        TagValue::Int(v) => {
            let kind = integer_type(*v)?;
            let bytes = match kind {
                b'I' | b'i' => (*v as i32).to_le_bytes().to_vec(),
                b's' | b'S' => (*v as i16).to_le_bytes().to_vec(),
                _ => vec![*v as u8],
            };
            let size = bytes.len();
            (kind, bytes, size)
        }
        TagValue::ByteArray { values, .. } => {
            let mut bytes = vec![b'c'];
            bytes.extend_from_slice(&(values.len() as i32).to_le_bytes());
            bytes.extend(values.iter().map(|&v| v as u8));
            (b'B', bytes, 5 + values.len())
        }
        TagValue::ShortArray { values, .. } => {
            let mut bytes = vec![b's'];
            bytes.extend_from_slice(&(values.len() as i32).to_le_bytes());
            for v in values {
                bytes.extend_from_slice(&v.to_le_bytes());
            }
            (b'B', bytes, 5 + 2 * values.len())
        }
        TagValue::IntArray { values, .. } => {
            let mut bytes = vec![b'i'];
            bytes.extend_from_slice(&(values.len() as i32).to_le_bytes());
            for v in values {
                bytes.extend_from_slice(&v.to_le_bytes());
            }
            (b'B', bytes, 5 + 4 * values.len())
        }
        TagValue::FloatArray(values) => {
            let mut bytes = vec![b'f'];
            bytes.extend_from_slice(&(values.len() as i32).to_le_bytes());
            for v in values {
                bytes.extend_from_slice(&v.to_le_bytes());
            }
            (b'B', bytes, 5 + 4 * values.len())
        }
    };
    Ok(Tag {
        id: [name[0], name[1], kind],
        value: bytes,
        size,
    })
}

/// `CRAMCompressionRecord(cramVersion, encodingStrategy, samRecord, referenceBases, ...)`.
fn to_write_record(
    record: &BamRecord,
    header: &SamHeader,
    reference: Option<&[u8]>,
    read_groups: &HashMap<String, i32>,
) -> Result<WriteRecord, CramError> {
    let flags = i32::from(record.flags);
    let paired = flags & 0x1 != 0;
    let unmapped = flags & 0x4 != 0;
    let mut bam = 0;
    let mut set = |on: bool, bit: i32| {
        if on {
            bam |= bit;
        }
    };
    set(paired, 0x1);
    set(paired && flags & 0x2 != 0, 0x2);
    set(unmapped, 0x4);
    set(paired && flags & 0x40 != 0, 0x40);
    set(paired && flags & 0x80 != 0, 0x80);
    set(flags & 0x10 != 0, 0x10);
    set(flags & 0x100 != 0, 0x100);
    set(flags & 0x800 != 0, 0x800);
    set(flags & 0x200 != 0, 0x200);
    set(flags & 0x400 != 0, 0x400);

    let mut cram = CF_DETACHED;
    let read_length = record.read_bases.len() as i32;
    let read_bases = to_bam_read_bases(&record.read_bases);

    let (features, end) = if unmapped {
        (Vec::new(), 0)
    } else {
        let qualities = if record.base_qualities.is_empty() {
            Qualities::Missing
        } else {
            Qualities::Present(&record.base_qualities)
        };
        let features = create_read_features(
            &record.cigar,
            record.alignment_start,
            &read_bases,
            qualities,
            reference.unwrap_or(&[]),
        )
        .map_err(|error| CramError::Other {
            class: error.java_exception(),
            message: error.message(),
        })?;
        let end = alignment_end(&features, record.alignment_start, read_length);
        (features, end)
    };

    let mut mate = 0;
    let (mate_start, mate_reference) = if paired {
        if flags & 0x8 != 0 {
            mate |= MF_MATE_UNMAPPED;
            bam |= 0x8;
        }
        if flags & 0x20 != 0 {
            mate |= MF_MATE_NEG_STRAND;
            bam |= 0x20;
        }
        (record.mate_alignment_start, record.mate_reference_index)
    } else {
        (0, -1)
    };
    if record.read_bases.is_empty() {
        cram |= CF_UNKNOWN_BASES;
    }
    let mut qualities = record.base_qualities.clone();
    if !qualities.is_empty() {
        cram |= CF_QS_PRESERVED_AS_ARRAY;
    } else if !features.is_empty() {
        qualities = vec![MISSING_QUALITY_SCORE; read_length.max(0) as usize];
        cram |= CF_QS_PRESERVED_AS_ARRAY;
    }

    let read_group = match record.tags.get(htsjdk_bam::tag::Tag::new(b"RG")) {
        Some(TagValue::Str(id)) => read_groups.get(id).copied().unwrap_or(NO_READGROUP_ID),
        _ => NO_READGROUP_ID,
    };

    // NM and MD: dropped when they are what a reader would regenerate.
    let mut strip_nm = !unmapped;
    let mut strip_md = !unmapped;
    if let Some(reference) = reference {
        if !unmapped && !record.cigar.is_empty() && !record.read_bases.is_empty() {
            let (md, nm) = calculate_md_nm(
                &record.cigar,
                &record.read_bases,
                reference,
                record.alignment_start,
            );
            if let Some(stored) = record.tags.get(htsjdk_bam::tag::Tag::new(b"NM")) {
                let stored = match stored {
                    TagValue::Int(v) => *v as i32,
                    TagValue::Float(f) => *f as i32,
                    _ => i32::MIN,
                };
                if stored != nm {
                    strip_nm = false;
                }
            }
            if let Some(stored) = record.tags.get(htsjdk_bam::tag::Tag::new(b"MD")) {
                if !matches!(stored, TagValue::Str(s) if *s == md) {
                    strip_md = false;
                }
            }
        }
    }
    let mut tags = Vec::new();
    for (tag, value) in record.tags.iter() {
        let name = tag.name();
        if &name == b"RG" || (strip_nm && &name == b"NM") || (strip_md && &name == b"MD") {
            continue;
        }
        tags.push(tag_of(name, value)?);
    }
    let _ = header;

    Ok(WriteRecord {
        bam_flags: bam,
        cram_flags: cram,
        mate_flags: mate,
        read_name: record.read_name.as_bytes().to_vec(),
        reference_index: record.reference_index,
        alignment_start: record.alignment_start,
        alignment_end: end,
        read_length,
        read_bases: if record.read_bases.is_empty() {
            Vec::new()
        } else {
            read_bases
        },
        features,
        template_size: record.inferred_insert_size,
        mapping_quality: i32::from(record.mapping_quality),
        mate_alignment_start: mate_start,
        mate_reference_index: mate_reference,
        quality_scores: qualities,
        read_group_id: read_group,
        tags,
        records_to_next_fragment: -1,
        tag_list: 0,
    })
}

/// `CRAMCompressionRecord.computeInsertSize` at 5.0.0, htslib's.
fn insert_size(first: &WriteRecord, second: &WriteRecord) -> i32 {
    if first.is_unmapped()
        || second.is_unmapped()
        || first.reference_index != second.reference_index
    {
        return 0;
    }
    let left = first.alignment_start.min(second.alignment_start);
    let right = first.alignment_end.max(second.alignment_end);
    let magnitude = right.wrapping_sub(left).wrapping_add(1);
    if first.alignment_start < second.alignment_start {
        return magnitude;
    }
    if first.alignment_start > second.alignment_start {
        return magnitude.wrapping_neg();
    }
    if first.bam_flags & 0x40 != 0 {
        return magnitude;
    }
    magnitude.wrapping_neg()
}

/// `Slice.linkMatesWithinSlice`: a pair whose stored template length is htslib's is attached.
fn link_mates(records: &mut [WriteRecord]) {
    let mut by_name: HashMap<Vec<u8>, usize> = HashMap::new();
    for i in 0..records.len() {
        let r = &records[i];
        if r.bam_flags & 0x1 == 0 || r.bam_flags & 0x100 != 0 || r.bam_flags & 0x800 != 0 {
            continue;
        }
        match by_name.get(&r.read_name).copied() {
            None => {
                by_name.insert(r.read_name.clone(), i);
            }
            Some(previous) => {
                let computed = insert_size(&records[previous], &records[i]);
                let name = records[i].read_name.clone();
                if records[previous].template_size != computed
                    || records[i].template_size != computed.wrapping_neg()
                {
                    by_name.remove(&name);
                    continue;
                }
                let p = &mut records[previous];
                p.cram_flags &= !CF_DETACHED;
                p.cram_flags |= CF_HAS_MATE_DOWNSTREAM;
                p.records_to_next_fragment = (i - previous - 1) as i32;
                let r = &mut records[i];
                r.cram_flags &= !CF_DETACHED;
                r.cram_flags &= !CF_HAS_MATE_DOWNSTREAM;
                by_name.remove(&name);
            }
        }
    }
}

/// A slice's bytes, its reference context, start and span, and its block count.
type SliceBytes = (Vec<u8>, (i32, i32, i32), i32);

/// The writer's state across containers.
pub struct CramWriter<'a> {
    header: SamHeader,
    reference: WriterReference<'a>,
    coordinate_sorted: bool,
    read_groups: HashMap<String, i32>,
    /// Every tag id any container has seen, which is every tag block every slice writes.
    tag_ids: BTreeSet<i32>,
    trials: HashMap<i32, TrialCompressor>,
    global_record_counter: i64,
    out: Vec<u8>,
    /// `CRAMBAIIndexer`, when the output is to be indexed.
    index: Option<crate::bai::CramBaiIndexer>,
}

impl<'a> CramWriter<'a> {
    /// `CRAMContainerStreamWriter.writeHeader`: the file definition and the SAM header container.
    ///
    /// `id` is what `SAMFileWriterFactory` passes, the output's `file:` URI, truncated to twenty
    /// bytes in the file definition.
    pub fn new(header: &SamHeader, id: &str, reference: &'a mut dyn ReferenceSource) -> Self {
        let header = stamp_sort_order(header);
        let mut out = Vec::new();
        out.extend_from_slice(b"CRAM");
        out.push(MAJOR);
        out.push(MINOR);
        let mut id_bytes = [0u8; 20];
        let bytes = id.as_bytes();
        let take = bytes.len().min(20);
        id_bytes[..take].copy_from_slice(&bytes[..take]);
        out.extend_from_slice(&id_bytes);

        // The SAM header container: one GZIP FILE_HEADER block of length-prefixed text.
        let text = header.encode_replacing_version();
        let mut content = (text.len() as i32).to_le_bytes().to_vec();
        content.extend_from_slice(text.as_bytes());
        let compressed = crate::gzip::compress(&content, GZIP_LEVEL);
        let block = block_bytes(
            CompressionMethod::Gzip,
            ContentType::FileHeader,
            0,
            &compressed,
            content.len(),
        );
        out.extend_from_slice(&container_header_bytes(
            block.len() as i32,
            -1,
            0,
            0,
            0,
            0,
            0,
            1,
            &[],
        ));
        out.extend_from_slice(&block);

        let read_groups = header
            .read_groups
            .iter()
            .enumerate()
            .map(|(i, g)| (g.id.clone(), i as i32))
            .collect();
        let coordinate_sorted = header.attributes.get("SO") == Some("coordinate");
        Self {
            header,
            reference: WriterReference {
                source: reference,
                current: None,
            },
            coordinate_sorted,
            read_groups,
            tag_ids: BTreeSet::new(),
            trials: HashMap::new(),
            global_record_counter: 0,
            out,
            index: None,
        }
    }

    /// Index the file as it is written, `SAMFileWriterFactory`'s `CREATE_INDEX` for a CRAM.
    pub fn with_index(mut self) -> Self {
        let lengths = self.header.sequences.iter().map(|s| s.length).collect();
        self.index = Some(crate::bai::CramBaiIndexer::new(lengths));
        self
    }

    /// Every record, then the EOF container: `writeAlignment` for each and `finish(true)`.
    pub fn write_all(self, records: &[BamRecord]) -> Result<Vec<u8>, CramError> {
        Ok(self.write_all_indexed(records)?.0)
    }

    /// [`Self::write_all`], with the BAI index when [`Self::with_index`] asked for one.
    pub fn write_all_indexed(
        mut self,
        records: &[BamRecord],
    ) -> Result<(Vec<u8>, Option<Vec<u8>>), CramError> {
        let mut current = UNINITIALIZED_REFERENCE_ID;
        let mut pending: Vec<&BamRecord> = Vec::new();
        let mut bases: i64 = 0;
        let mut slices: Vec<Vec<&BamRecord>> = Vec::new();
        for record in records {
            let next = record.reference_index;
            let updated =
                self.updated_reference_context(current, next, pending.len(), bases, slices.len())?;
            if updated == UNINITIALIZED_REFERENCE_ID {
                slices.push(std::mem::take(&mut pending));
                bases = 0;
                if slices.len() >= SLICES_PER_CONTAINER
                    || current == MULTIPLE_REFERENCE_ID
                    || current != next
                {
                    let staged = std::mem::take(&mut slices);
                    self.write_container(&staged)?;
                }
                current = next;
            } else {
                current = updated;
            }
            pending.push(record);
            bases += record.read_bases.len() as i64;
        }
        if !pending.is_empty() {
            slices.push(pending);
        }
        if !slices.is_empty() {
            self.write_container(&slices)?;
        }
        self.out.extend_from_slice(&EOF_V3);
        let index = self.index.map(crate::bai::CramBaiIndexer::finish);
        Ok((self.out, index))
    }

    /// `SliceFactory.getUpdatedReferenceContext`.
    fn updated_reference_context(
        &self,
        current: i32,
        next: i32,
        count: usize,
        bases: i64,
        staged_slices: usize,
    ) -> Result<i32, CramError> {
        let full = count >= READS_PER_SLICE || bases >= BASES_PER_SLICE;
        Ok(match current {
            UNINITIALIZED_REFERENCE_ID => {
                if count != 0 {
                    return Err(CramError::Cram(
                        "Reference context should have been initialized if records have previously been processed".into(),
                    ));
                }
                next
            }
            UNMAPPED_UNPLACED_ID => {
                if next == current {
                    if full {
                        UNINITIALIZED_REFERENCE_ID
                    } else {
                        UNMAPPED_UNPLACED_ID
                    }
                } else if self.coordinate_sorted {
                    return Err(CramError::Cram(
                        "Invalid coord-sorted input - unmapped records must be last".into(),
                    ));
                } else if full {
                    UNINITIALIZED_REFERENCE_ID
                } else {
                    MULTIPLE_REFERENCE_ID
                }
            }
            MULTIPLE_REFERENCE_ID => {
                if self.coordinate_sorted {
                    if count < MINIMUM_SINGLE_REFERENCE_SLICE_SIZE {
                        MULTIPLE_REFERENCE_ID
                    } else {
                        UNINITIALIZED_REFERENCE_ID
                    }
                } else if full {
                    UNINITIALIZED_REFERENCE_ID
                } else {
                    MULTIPLE_REFERENCE_ID
                }
            }
            _ => {
                if next == current {
                    if full {
                        UNINITIALIZED_REFERENCE_ID
                    } else {
                        next
                    }
                } else if count < MINIMUM_SINGLE_REFERENCE_SLICE_SIZE {
                    if staged_slices > 0 {
                        UNINITIALIZED_REFERENCE_ID
                    } else {
                        MULTIPLE_REFERENCE_ID
                    }
                } else {
                    UNINITIALIZED_REFERENCE_ID
                }
            }
        })
    }

    /// `ContainerFactory.makeContainer` and `Container.write`.
    fn write_container(&mut self, staged: &[Vec<&BamRecord>]) -> Result<(), CramError> {
        // `SliceFactory.convertToCRAMRecords`, every slice counting from zero.
        let mut slices: Vec<Vec<WriteRecord>> = Vec::with_capacity(staged.len());
        for slice in staged {
            let mut converted = Vec::with_capacity(slice.len());
            for record in slice {
                let reference = if record.reference_index != -1 {
                    Some(
                        self.reference
                            .contig(&self.header, record.reference_index)?,
                    )
                } else {
                    None
                };
                converted.push(to_write_record(
                    record,
                    &self.header,
                    reference.as_deref().map(Vec::as_slice),
                    &self.read_groups,
                )?);
            }
            slices.push(converted);
        }

        // `CompressionHeaderFactory.createCompressionHeader`.
        let dictionary = build_tag_dictionary(&mut slices);
        let mut tag_encodings = TagEncodingMap::default();
        let mut discovered: BTreeSet<i32> = BTreeSet::new();
        for record in slices.iter().flatten() {
            for tag in &record.tags {
                discovered.insert(tag.key());
            }
        }
        for &key in &discovered {
            if (key & 0xFFFF00) == ((b'S' as i32) << 16 | (b'A' as i32) << 8)
                || (key & 0xFFFF00) == ((b'X' as i32) << 16 | (b'A' as i32) << 8)
            {
                return Err(CramError::Other {
                    class: "java.lang.UnsupportedOperationException",
                    message: "an SA or XA tag in a CRAM output is compressed with bzip2, which is not ported".into(),
                });
            }
            tag_encodings.put(key, tag_encoding(key, slices.iter().flatten()));
            self.tag_ids.insert(key);
            self.trials.entry(key).or_insert_with(TrialCompressor::new);
        }
        let mut frequencies = vec![[0i64; SYMBOL_SPACE_SIZE]; SYMBOL_SPACE_SIZE];
        for record in slices.iter().flatten() {
            for feature in &record.features {
                if let ReadFeature::Substitution {
                    base,
                    reference_base,
                    ..
                } = feature
                {
                    if *reference_base == 0 || *base == 0 || *reference_base > 127 || *base > 127 {
                        return Err(CramError::IllegalArgument(format!(
                            "CRAM: Attempt to generate a substitution code for invalid reference base with value '{}'",
                            *reference_base as i8
                        )));
                    }
                    frequencies[*reference_base as usize][*base as usize] += 1;
                }
            }
        }
        let matrix = SubstitutionMatrix::from_frequencies(&frequencies);
        for record in slices.iter_mut().flatten() {
            for feature in record.features.iter_mut() {
                if let ReadFeature::Substitution {
                    base,
                    reference_base,
                    code,
                    ..
                } = feature
                {
                    if *code == NO_CODE {
                        *code = matrix
                            .code(*reference_base, *base)
                            .map_err(|e| CramError::IllegalArgument(e.message()))?
                            as i8;
                    }
                }
            }
        }
        let mut encodings = EncodingMap::default();
        for (name, _, content_id, _) in WRITTEN {
            let series = DataSeries::by_canonical_name(name.as_bytes().try_into().unwrap())
                .expect("a written series is a data series");
            encodings.put(series, descriptor(name, content_id));
        }
        let compression = CompressionHeader {
            preservation: PreservationMap {
                preserve_read_names: true,
                ap_delta: self.coordinate_sorted,
                reference_required: true,
                substitution_matrix: matrix.encoded(),
                tag_id_dictionary: dictionary,
            },
            encodings,
            tag_encodings,
        };

        let container_offset = self.out.len() as i64;
        let mut entries = Vec::new();
        let mut body = compression.write_block(MAJOR);
        let mut landmarks = Vec::new();
        let (mut record_count, mut block_count, mut base_count) = (0i32, 0i32, 0i64);
        let mut contexts: Vec<(i32, i32, i32)> = Vec::new();
        for records in slices.iter_mut() {
            landmarks.push(body.len() as i32);
            if self.index.is_some() {
                entries.extend(bai_entries(
                    records,
                    container_offset,
                    i64::from(*landmarks.last().unwrap()),
                    landmarks.len() as i32 - 1,
                ));
            }
            let (slice_bytes, context, blocks) = self.slice(records, &compression, &matrix)?;
            body.extend_from_slice(&slice_bytes);
            record_count += records.len() as i32;
            block_count += blocks;
            base_count += records
                .iter()
                .map(|r| i64::from(r.read_length))
                .sum::<i64>();
            contexts.push(context);
        }

        // `Container`: the slices' common reference context, or multiple.
        let distinct: BTreeSet<i32> = contexts.iter().map(|c| c.0).collect();
        let (reference_context, start, span) = if distinct.len() > 1 {
            (MULTIPLE_REFERENCE_ID, 0, 0)
        } else {
            let id = *distinct.iter().next().unwrap();
            if id >= 0 {
                let start = contexts.iter().map(|c| c.1).min().unwrap();
                let end_plus_one = contexts.iter().map(|c| c.1 + c.2).max().unwrap();
                (id, start, end_plus_one - start)
            } else {
                (id, 0, 0)
            }
        };
        self.out.extend_from_slice(&container_header_bytes(
            body.len() as i32,
            reference_context,
            start,
            span,
            record_count,
            self.global_record_counter,
            base_count,
            block_count,
            &landmarks,
        ));
        self.out.extend_from_slice(&body);
        self.global_record_counter += i64::from(record_count);
        if let Some(index) = self.index.as_mut() {
            index
                .process_container(entries)
                .map_err(|message| CramError::Other {
                    class: "htsjdk.samtools.SAMException",
                    message,
                })?;
        }
        Ok(())
    }

    /// `new Slice(records, compressionHeader, ...)`, `setReferenceMD5` and `Slice.write`: the
    /// slice's bytes, its reference context, start and span, and its block count.
    fn slice(
        &mut self,
        records: &mut [WriteRecord],
        compression: &CompressionHeader,
        _matrix: &SubstitutionMatrix,
    ) -> Result<SliceBytes, CramError> {
        let mut contexts: BTreeSet<i32> = BTreeSet::new();
        let mut start = i32::MAX;
        let mut end = 0;
        for record in records.iter() {
            if record.is_placed() {
                contexts.insert(record.reference_index);
                start = start.min(record.alignment_start);
                end = end.max(record.alignment_end);
            } else {
                contexts.insert(UNMAPPED_UNPLACED_ID);
            }
        }
        let (reference_context, alignment_start, alignment_span) = match contexts.len() {
            0 => (UNMAPPED_UNPLACED_ID, 0, 0),
            1 => {
                let id = *contexts.iter().next().unwrap();
                if id >= 0 {
                    (id, start, end - start + 1)
                } else {
                    (UNMAPPED_UNPLACED_ID, 0, 0)
                }
            }
            _ => (MULTIPLE_REFERENCE_ID, 0, 0),
        };

        // The codec context fqzcomp reads, before the mates are linked.
        let quality_lengths: Vec<i32> = records
            .iter()
            .map(|r| {
                if r.cram_flags & CF_QS_PRESERVED_AS_ARRAY != 0 {
                    r.read_length
                } else {
                    0
                }
            })
            .collect();
        let bam_flags: Vec<i32> = records.iter().map(|r| r.bam_flags).collect();
        link_mates(records);

        // `CramRecordWriter`.
        let mut streams = SliceWriteStreams::new();
        let mut external_ids: BTreeSet<i32> = WRITTEN.iter().map(|w| w.2).collect();
        external_ids.extend(self.tag_ids.iter().copied());
        for &id in &external_ids {
            streams.external(id);
        }
        let encodings = resolve_encodings(compression)?;
        let multi = reference_context == MULTIPLE_REFERENCE_ID;
        let mut previous = alignment_start;
        for record in records.iter() {
            write_record(
                &encodings,
                compression,
                multi,
                &mut streams,
                record,
                previous,
            )?;
            previous = record.alignment_start;
        }
        let blocks = streams.finish();

        let mut block_bytes_out = raw_block_with_id(ContentType::Core, 0, &blocks.core, MAJOR);
        for (&id, content) in &blocks.external {
            let (compressed, method) =
                self.compress_block(id, content, &quality_lengths, &bam_flags)?;
            block_bytes_out.extend_from_slice(&block_bytes(
                method,
                ContentType::External,
                id,
                &compressed,
                content.len(),
            ));
        }

        let mut md5 = [0u8; 16];
        if reference_context >= 0 {
            let contig = self.reference.contig(&self.header, reference_context)?;
            // `fetchReferenceBasesByRegion` for the slice, then its MD5 truncated to the contig.
            let from = (alignment_start - 1).max(0) as usize;
            if from >= contig.len() {
                return Err(CramError::IllegalArgument(format!(
                    "Requested start {} is beyond the sequence length {}",
                    alignment_start - 1,
                    self.header.sequences[reference_context as usize].name
                )));
            }
            let to = contig.len().min(from + alignment_span.max(0) as usize);
            md5 = Md5::digest(&contig[from..to]).into();
        }
        let header = SliceHeader {
            reference_context_id: reference_context,
            alignment_start,
            alignment_span,
            record_count: records.len() as i32,
            global_record_counter: 0,
            block_count: 1 + blocks.external.len() as i32,
            content_ids: blocks.external.keys().copied().collect(),
            embedded_reference_content_id: -1,
            reference_md5: md5,
            tags: Vec::new(),
        };
        let mut out = raw_block(ContentType::MappedSlice, &header.write(MAJOR), MAJOR);
        out.extend_from_slice(&block_bytes_out);
        Ok((
            out,
            (reference_context, alignment_start, alignment_span),
            1 + blocks.external.len() as i32,
        ))
    }

    /// `createCompressedBlockForWriter`: the compressor the encoding map holds for this id.
    fn compress_block(
        &mut self,
        id: i32,
        content: &[u8],
        quality_lengths: &[i32],
        bam_flags: &[i32],
    ) -> Result<(Vec<u8>, CompressionMethod), CramError> {
        if let Some(&(_, _, _, compressor)) = WRITTEN.iter().find(|w| w.2 == id) {
            return Ok(match compressor {
                Compressor::Gzip => (
                    crate::gzip::compress(content, GZIP_LEVEL),
                    CompressionMethod::Gzip,
                ),
                Compressor::Rans0 => (
                    crate::rans_nx16::compress(content, 0),
                    CompressionMethod::RansNx16,
                ),
                Compressor::Rans1 => (
                    crate::rans_nx16::compress(content, 1),
                    CompressionMethod::RansNx16,
                ),
                Compressor::Fqzcomp => {
                    if content.is_empty() {
                        (Vec::new(), CompressionMethod::Fqzcomp)
                    } else {
                        let sum: i64 = quality_lengths.iter().map(|&l| i64::from(l)).sum();
                        let (lengths, flags): (Vec<i32>, Vec<i32>) = if sum == content.len() as i64
                        {
                            quality_lengths
                                .iter()
                                .zip(bam_flags)
                                .filter(|(l, _)| **l != 0)
                                .map(|(l, f)| (*l, *f))
                                .unzip()
                        } else {
                            (vec![content.len() as i32], Vec::new())
                        };
                        // Without lengths that add up the flags are those of every record, and
                        // the one-record fallback reads the first of them.
                        let flags = if sum == content.len() as i64 {
                            flags
                        } else {
                            bam_flags.to_vec()
                        };
                        (
                            crate::fqzcomp::compress(content, &lengths, Some(&flags)),
                            CompressionMethod::Fqzcomp,
                        )
                    }
                }
                Compressor::NameTokeniser => (
                    crate::name_tokeniser::compress(content, crate::name_tokeniser::NAME_SEPARATOR)
                        .map_err(|e| CramError::Cram(e.to_string()))?,
                    CompressionMethod::NameTokeniser,
                ),
            });
        }
        let trial = self.trials.entry(id).or_insert_with(TrialCompressor::new);
        Ok(trial.compress(content))
    }
}

/// `Slice.getBAIEntries` for a slice being written, whose counts its constructor kept.
fn bai_entries(
    records: &[WriteRecord],
    container_offset: i64,
    slice_offset: i64,
    landmark_index: i32,
) -> Vec<crate::bai::BaiEntry> {
    use crate::bai::{BaiEntry, Span};
    let entry = |reference: i32, s: Span| BaiEntry {
        reference,
        alignment_start: s.start,
        alignment_span: s.span,
        mapped: s.mapped,
        unmapped: s.unmapped,
        unplaced: s.unplaced,
        container_offset,
        slice_offset,
        landmark_index,
    };
    let mut contexts: BTreeSet<i32> = BTreeSet::new();
    let (mut start, mut end) = (i32::MAX, 0);
    let (mut mapped, mut unmapped, mut unplaced) = (0, 0, 0);
    for r in records {
        if r.is_placed() {
            contexts.insert(r.reference_index);
            start = start.min(r.alignment_start);
            end = end.max(r.alignment_end);
            if r.is_unmapped() {
                unmapped += 1;
            } else {
                mapped += 1;
            }
        } else {
            contexts.insert(UNMAPPED_UNPLACED_ID);
        }
        if r.alignment_start == 0 {
            unplaced += 1;
        }
    }
    let counts = |start, span| Span {
        start,
        span,
        mapped,
        unmapped,
        unplaced,
    };
    if contexts.len() <= 1 {
        let id = contexts
            .iter()
            .next()
            .copied()
            .unwrap_or(UNMAPPED_UNPLACED_ID);
        if id < 0 {
            return vec![entry(UNMAPPED_UNPLACED_ID, counts(0, 0))];
        }
        return vec![entry(id, counts(start, end - start + 1))];
    }
    // `getMultiRefAlignmentSpans` and `mergeRecordSpan`.
    let mut spans: std::collections::BTreeMap<i32, Span> = std::collections::BTreeMap::new();
    for r in records {
        let (key, span) = if r.is_unmapped() {
            if r.alignment_start == 0 {
                (
                    UNMAPPED_UNPLACED_ID,
                    Span {
                        start: 0,
                        span: 0,
                        mapped: 0,
                        unmapped: 1,
                        unplaced: 1,
                    },
                )
            } else {
                (
                    r.reference_index,
                    Span {
                        start: r.alignment_start,
                        span: r.read_length,
                        mapped: 0,
                        unmapped: 1,
                        unplaced: 0,
                    },
                )
            }
        } else {
            (
                r.reference_index,
                Span {
                    start: r.alignment_start,
                    span: r.alignment_end - r.alignment_start,
                    mapped: 1,
                    unmapped: 0,
                    unplaced: 0,
                },
            )
        };
        spans
            .entry(key)
            .and_modify(|existing| *existing = Span::combine(*existing, span))
            .or_insert(span);
    }
    let mut out: Vec<BaiEntry> = spans
        .iter()
        .filter(|(k, _)| **k != UNMAPPED_UNPLACED_ID)
        .map(|(k, s)| entry(*k, *s))
        .collect();
    if let Some(s) = spans.get(&UNMAPPED_UNPLACED_ID) {
        out.push(entry(UNMAPPED_UNPLACED_ID, *s));
    }
    out
}

/// `SAMFileWriterImpl.setHeader`'s sort order stamp: a recognised `SO` kept, anything else
/// `unsorted`.
fn stamp_sort_order(header: &SamHeader) -> SamHeader {
    const RECOGNISED: [&str; 5] = [
        "unsorted",
        "queryname",
        "coordinate",
        "duplicate",
        "unknown",
    ];
    let mut stamped = header.clone();
    let order = match header.attributes.get("SO") {
        Some(value) if RECOGNISED.contains(&value) => value.to_string(),
        _ => "unsorted".to_string(),
    };
    stamped.set_sort_order(&order);
    stamped
}

/// `buildTagIdDictionary`: each record's tags sorted by id, the distinct id lists ordered by
/// length and then by their bytes as signed values, the empty list first.
fn build_tag_dictionary(slices: &mut [Vec<WriteRecord>]) -> Vec<Vec<[u8; 3]>> {
    fn compare(a: &[u8], b: &[u8]) -> std::cmp::Ordering {
        a.len().cmp(&b.len()).then_with(|| {
            for (x, y) in a.iter().zip(b) {
                if x != y {
                    return (*x as i8).cmp(&(*y as i8));
                }
            }
            std::cmp::Ordering::Equal
        })
    }
    let mut lists: Vec<Vec<u8>> = vec![Vec::new()];
    for record in slices.iter_mut().flatten() {
        record.tags.sort_by_key(Tag::key);
        let ids: Vec<u8> = record.tags.iter().flat_map(|t| t.id).collect();
        if !lists.contains(&ids) {
            lists.push(ids);
        }
    }
    lists.sort_by(|a, b| compare(a, b));
    for record in slices.iter_mut().flatten() {
        let ids: Vec<u8> = record.tags.iter().flat_map(|t| t.id).collect();
        record.tag_list = lists.iter().position(|l| *l == ids).unwrap();
    }
    lists
        .iter()
        .map(|l| l.chunks(3).map(|c| [c[0], c[1], c[2]]).collect())
        .collect()
}

/// `buildEncodingForTag`.
fn tag_encoding<'r>(
    key: i32,
    records: impl Iterator<Item = &'r WriteRecord>,
) -> EncodingDescriptor {
    let fixed = |size: i32| {
        // ByteArrayLen(CanonicalHuffmanInteger({size}, {0}), ExternalByteArray(key)).
        let mut huffman = itf8(1);
        huffman.extend_from_slice(&itf8(size));
        huffman.extend_from_slice(&itf8(1));
        huffman.extend_from_slice(&itf8(0));
        let external = itf8(key);
        let mut parameters = itf8(EncodingId::Huffman as i32);
        parameters.extend_from_slice(&itf8(huffman.len() as i32));
        parameters.extend_from_slice(&huffman);
        parameters.extend_from_slice(&itf8(EncodingId::External as i32));
        parameters.extend_from_slice(&itf8(external.len() as i32));
        parameters.extend_from_slice(&external);
        EncodingDescriptor {
            id: EncodingId::ByteArrayLen,
            parameters,
        }
    };
    match (key & 0xFF) as u8 {
        b'A' | b'c' | b'C' => fixed(1),
        b'I' | b'i' | b'f' => fixed(4),
        b's' | b'S' => fixed(2),
        b'Z' => {
            let mut parameters = vec![b'\t'];
            parameters.extend_from_slice(&itf8(key));
            EncodingDescriptor {
                id: EncodingId::ByteArrayStop,
                parameters,
            }
        }
        _ => {
            let mut min = usize::MAX;
            let mut max = 0usize;
            let mut values = Vec::new();
            for record in records {
                for tag in &record.tags {
                    if tag.key() == key {
                        min = min.min(tag.size);
                        max = max.max(tag.size);
                        values.extend_from_slice(&tag.value);
                    }
                }
            }
            if min == max {
                return fixed(min as i32);
            }
            if min > 100 {
                let mut used = [false; 256];
                for &b in &values {
                    used[b as usize] = true;
                }
                if let Some(unused) = used.iter().position(|u| !u) {
                    let mut parameters = vec![unused as u8];
                    parameters.extend_from_slice(&itf8(key));
                    return EncodingDescriptor {
                        id: EncodingId::ByteArrayStop,
                        parameters,
                    };
                }
            }
            // ByteArrayLen(ExternalInteger(key), ExternalByteArray(key)).
            let external = itf8(key);
            let mut parameters = itf8(EncodingId::External as i32);
            parameters.extend_from_slice(&itf8(external.len() as i32));
            parameters.extend_from_slice(&external);
            parameters.extend_from_slice(&itf8(EncodingId::External as i32));
            parameters.extend_from_slice(&itf8(external.len() as i32));
            parameters.extend_from_slice(&external);
            EncodingDescriptor {
                id: EncodingId::ByteArrayLen,
                parameters,
            }
        }
    }
}

/// The resolved writers: one encoding per series by name, the QS array, and one per tag.
struct Encodings {
    series: HashMap<&'static str, Encoding>,
    quality_array: Encoding,
    tags: HashMap<i32, Encoding>,
}

impl Encodings {
    fn get(&self, name: &str) -> Result<&Encoding, CramError> {
        self.series.get(name).ok_or_else(|| {
            CramError::IllegalArgument(format!(
                "Attempt to create data series writer for data series {name} for which no encoding can be found"
            ))
        })
    }
}

fn resolve_encodings(compression: &CompressionHeader) -> Result<Encodings, CramError> {
    let factory = |e: crate::encoding_factory::FactoryError| CramError::Cram(e.message());
    let mut series = HashMap::new();
    for (name, kind, _, _) in WRITTEN {
        let ds = DataSeries::by_canonical_name(name.as_bytes().try_into().unwrap()).unwrap();
        let descriptor = compression.encodings.get(ds).unwrap();
        series.insert(
            name,
            create_encoding(kind, descriptor.id, &descriptor.parameters).map_err(factory)?,
        );
    }
    let qs = compression
        .encodings
        .get(DataSeries::by_canonical_name(b"QS").unwrap())
        .unwrap();
    let quality_array =
        create_encoding(FactoryType::ByteArray, qs.id, &qs.parameters).map_err(factory)?;
    let mut tags = HashMap::new();
    for key in compression.tag_encodings.keys() {
        let descriptor = compression.tag_encodings.get(key).unwrap();
        tags.insert(
            key,
            create_encoding(
                FactoryType::ByteArray,
                descriptor.id,
                &descriptor.parameters,
            )
            .map_err(factory)?,
        );
    }
    Ok(Encodings {
        series,
        quality_array,
        tags,
    })
}

/// `CramRecordWriter.writeCRAMRecord` at htsjdk 5.0.0.
fn write_record(
    encodings: &Encodings,
    compression: &CompressionHeader,
    multi_reference: bool,
    streams: &mut SliceWriteStreams,
    r: &WriteRecord,
    previous_alignment_start: i32,
) -> Result<(), CramError> {
    use crate::codecs::{write_byte, write_byte_array, write_int};
    let codec = |e: crate::codecs::ReadError| CramError::Cram(e.message());
    write_int(encodings.get("BF")?, streams, r.bam_flags).map_err(codec)?;
    write_int(encodings.get("CF")?, streams, r.cram_flags & 0xFF).map_err(codec)?;
    if multi_reference {
        write_int(encodings.get("RI")?, streams, r.reference_index).map_err(codec)?;
    }
    write_int(encodings.get("RL")?, streams, r.read_length).map_err(codec)?;
    let ap = if compression.preservation.ap_delta {
        r.alignment_start.wrapping_sub(previous_alignment_start)
    } else {
        r.alignment_start
    };
    write_int(encodings.get("AP")?, streams, ap).map_err(codec)?;
    write_int(encodings.get("RG")?, streams, r.read_group_id).map_err(codec)?;
    write_byte_array(encodings.get("RN")?, streams, &r.read_name).map_err(codec)?;
    if r.cram_flags & CF_DETACHED != 0 {
        write_int(encodings.get("MF")?, streams, r.mate_flags & 0xFF).map_err(codec)?;
        write_int(encodings.get("NS")?, streams, r.mate_reference_index).map_err(codec)?;
        write_int(encodings.get("NP")?, streams, r.mate_alignment_start).map_err(codec)?;
        write_int(encodings.get("TS")?, streams, r.template_size).map_err(codec)?;
    } else if r.cram_flags & CF_HAS_MATE_DOWNSTREAM != 0 {
        write_int(encodings.get("NF")?, streams, r.records_to_next_fragment).map_err(codec)?;
    }
    write_int(encodings.get("TL")?, streams, r.tag_list as i32).map_err(codec)?;
    for tag in &r.tags {
        let encoding = encodings
            .tags
            .get(&tag.key())
            .ok_or_else(|| CramError::Cram(format!("no encoding for tag {}", tag.key())))?;
        write_byte_array(encoding, streams, &tag.value).map_err(codec)?;
    }
    if r.bam_flags & 0x4 == 0 {
        write_int(encodings.get("FN")?, streams, r.features.len() as i32).map_err(codec)?;
        let mut previous = 0;
        for f in &r.features {
            write_byte(encodings.get("FC")?, streams, f.operator() as i8).map_err(codec)?;
            write_int(encodings.get("FP")?, streams, f.position() - previous).map_err(codec)?;
            previous = f.position();
            match f {
                ReadFeature::ReadBase { base, quality, .. } => {
                    write_byte(encodings.get("BA")?, streams, *base as i8).map_err(codec)?;
                    write_byte(encodings.get("QS")?, streams, *quality).map_err(codec)?;
                }
                ReadFeature::Substitution { code, .. } => {
                    write_byte(encodings.get("BS")?, streams, *code).map_err(codec)?;
                }
                ReadFeature::Insertion { sequence, .. } => {
                    write_byte_array(encodings.get("IN")?, streams, sequence).map_err(codec)?;
                }
                ReadFeature::SoftClip { sequence, .. } => {
                    write_byte_array(encodings.get("SC")?, streams, sequence).map_err(codec)?;
                }
                ReadFeature::HardClip { length, .. } => {
                    write_int(encodings.get("HC")?, streams, *length).map_err(codec)?;
                }
                ReadFeature::Padding { length, .. } => {
                    write_int(encodings.get("PD")?, streams, *length).map_err(codec)?;
                }
                ReadFeature::Deletion { length, .. } => {
                    write_int(encodings.get("DL")?, streams, *length).map_err(codec)?;
                }
                ReadFeature::RefSkip { length, .. } => {
                    write_int(encodings.get("RS")?, streams, *length).map_err(codec)?;
                }
                ReadFeature::InsertBase { base, .. } => {
                    write_byte(encodings.get("BA")?, streams, *base as i8).map_err(codec)?;
                }
                ReadFeature::BaseQualityScore { quality, .. } => {
                    write_byte(encodings.get("QS")?, streams, *quality).map_err(codec)?;
                }
                other => {
                    return Err(CramError::Other {
                        class: "java.lang.RuntimeException",
                        message: format!(
                            "Unknown read feature operator: {}",
                            other.operator() as char
                        ),
                    })
                }
            }
        }
        write_int(encodings.get("MQ")?, streams, r.mapping_quality).map_err(codec)?;
        if r.cram_flags & CF_QS_PRESERVED_AS_ARRAY != 0 {
            write_byte_array(&encodings.quality_array, streams, &r.quality_scores)
                .map_err(codec)?;
        }
    } else {
        if r.cram_flags & CF_UNKNOWN_BASES == 0 {
            for &b in &r.read_bases {
                write_byte(encodings.get("BA")?, streams, b as i8).map_err(codec)?;
            }
        }
        if r.cram_flags & CF_QS_PRESERVED_AS_ARRAY != 0 {
            write_byte_array(&encodings.quality_array, streams, &r.quality_scores)
                .map_err(codec)?;
        }
    }
    Ok(())
}

/// `Block.write`: method, content type, content id, the two sizes, the data and its CRC.
fn block_bytes(
    method: CompressionMethod,
    content_type: ContentType,
    content_id: i32,
    compressed: &[u8],
    uncompressed: usize,
) -> Vec<u8> {
    let mut out = vec![method as u8, content_type as u8];
    out.extend_from_slice(&itf8(content_id));
    out.extend_from_slice(&itf8(compressed.len() as i32));
    out.extend_from_slice(&itf8(uncompressed as i32));
    out.extend_from_slice(compressed);
    out.extend_from_slice(&crc32(&out).to_le_bytes());
    out
}

/// `ContainerHeader.write` for version 3: the fields, then their CRC-32, little-endian.
#[allow(clippy::too_many_arguments)]
fn container_header_bytes(
    blocks_byte_size: i32,
    reference_context: i32,
    alignment_start: i32,
    alignment_span: i32,
    record_count: i32,
    global_record_counter: i64,
    base_count: i64,
    block_count: i32,
    landmarks: &[i32],
) -> Vec<u8> {
    let mut out = blocks_byte_size.to_le_bytes().to_vec();
    out.extend_from_slice(&itf8(reference_context));
    out.extend_from_slice(&itf8(alignment_start));
    out.extend_from_slice(&itf8(alignment_span));
    out.extend_from_slice(&itf8(record_count));
    out.extend_from_slice(&write_unsigned_ltf8(global_record_counter).0);
    out.extend_from_slice(&write_unsigned_ltf8(base_count).0);
    out.extend_from_slice(&itf8(block_count));
    out.extend_from_slice(&itf8(landmarks.len() as i32));
    for &landmark in landmarks {
        out.extend_from_slice(&itf8(landmark));
    }
    out.extend_from_slice(&crc32(&out).to_le_bytes());
    out
}

/// A whole file and its BAI index: `CREATE_INDEX` with a `.cram` output.
pub fn write_cram_indexed(
    header: &SamHeader,
    records: &[BamRecord],
    id: &str,
    reference: &mut dyn ReferenceSource,
) -> Result<(Vec<u8>, Vec<u8>), CramError> {
    let (cram, index) = CramWriter::new(header, id, reference)
        .with_index()
        .write_all_indexed(records)?;
    Ok((cram, index.expect("an indexing writer returns its index")))
}

/// A whole file: [`CramWriter::new`] then [`CramWriter::write_all`].
pub fn write_cram(
    header: &SamHeader,
    records: &[BamRecord],
    id: &str,
    reference: &mut dyn ReferenceSource,
) -> Result<Vec<u8>, CramError> {
    CramWriter::new(header, id, reference).write_all(records)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_eof_container_is_the_one_htsjdk_writes() {
        let header = ContainerHeader::read(&EOF_V3, 3).unwrap();
        assert!(header.is_eof());
        assert_eq!(
            header.byte_length + header.blocks_byte_size as usize,
            EOF_V3.len()
        );
    }

    use crate::container::ContainerHeader;
}
