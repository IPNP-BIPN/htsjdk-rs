//! Writing a CRAM record: the same data series as the read side, in the same order.
//!
//! Ported from `htsjdk.samtools.cram.encoding.writer.CramRecordWriter.writeCRAMRecord` at
//! htsjdk 4.2.0.
//!
//! [`crate::record_read`] is the half that consumes. This is the half that produces, and between
//! them they are a round trip: the order is prescribed and shared, so a difference in either shows
//! as bytes that do not match rather than as an error.
//!
//! # Three features can be read and not written
//!
//! `Bases` and `Scores` fall to the writer's default arm, which throws; the reader has branches for
//! both. htsjdk's own comment says it does not generate them. A `Substitution` whose code is
//! negative is resolved against the compression header's substitution matrix on the way out, so a
//! record carrying a base and a reference base rather than a code needs that matrix to be written
//! at all.
//!
//! # An unmapped record writes its bases one at a time
//!
//! Into the same series a mapped record's read feature writes a base into, and its quality scores
//! only if it also has bases: the branch is nested, exactly as it is on the read side.

use crate::codecs::{write_byte, write_byte_array, write_int, ReadError};
use crate::compression_header::CompressionHeader;
use crate::encoding_factory::{create_encoding, DataSeriesType, Encoding};
use crate::encoding_map::{DataSeries, EncodingMap};
use crate::external_codecs::SliceWriteStreams;
use crate::read_features::ReadFeature;
use crate::record_read::{CramRecord, RecordReadError, SliceContext};

/// What writing a record refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RecordWriteError {
    /// A data series the compression header does not name. The reference words this one
    /// differently from the reader's, and refuses where the reader returns null.
    NoEncoding {
        name: &'static str,
    },
    /// `Bases` or `Scores`, which the writer has no arm for.
    UnknownReadFeature {
        operator: u8,
    },
    /// A substitution carrying no code, which needs the compression header's matrix.
    SubstitutionWithoutCode,
    Read(ReadError),
    Record(RecordReadError),
}

impl RecordWriteError {
    pub fn message(&self) -> String {
        match self {
            RecordWriteError::NoEncoding { name } => format!(
                "Attempt to create data series writer for data series {name} for which no encoding \
                 can be found"
            ),
            RecordWriteError::UnknownReadFeature { operator } => {
                format!("Unknown read feature operator: {}", *operator as char)
            }
            RecordWriteError::SubstitutionWithoutCode => {
                "a substitution with no code needs the compression header's matrix".to_string()
            }
            RecordWriteError::Read(error) => error.message(),
            RecordWriteError::Record(error) => error.message(),
        }
    }

    pub fn java_exception(&self) -> &'static str {
        match self {
            RecordWriteError::NoEncoding { .. } => "IllegalArgumentException",
            RecordWriteError::UnknownReadFeature { .. }
            | RecordWriteError::SubstitutionWithoutCode => "RuntimeException",
            RecordWriteError::Read(error) => error.java_exception(),
            RecordWriteError::Record(error) => error.java_exception(),
        }
    }
}

impl From<ReadError> for RecordWriteError {
    fn from(error: ReadError) -> Self {
        RecordWriteError::Read(error)
    }
}

impl From<RecordReadError> for RecordWriteError {
    fn from(error: RecordReadError) -> Self {
        RecordWriteError::Record(error)
    }
}

/// The encodings a record is written through, resolved once from the compression header.
///
/// The reference refuses here rather than at the first use: a writer for a series with no encoding
/// cannot be built at all, where a reader for one is simply null and fails later.
pub struct RecordWriters {
    map: Vec<(&'static str, Encoding)>,
    quality_score_array: Option<Encoding>,
}

impl RecordWriters {
    pub fn new(encodings: &EncodingMap) -> Result<Self, RecordWriteError> {
        let mut map = Vec::new();
        for name in SERIES {
            let series = match DataSeries::by_canonical_name(name.as_bytes().try_into().unwrap()) {
                Some(series) => series,
                None => continue,
            };
            if let Some(descriptor) = encodings.get(series) {
                let encoding = create_encoding(
                    series_type(series.series_type()),
                    descriptor.id,
                    &descriptor.parameters,
                )
                .map_err(RecordReadError::from)?;
                map.push((name, encoding));
            }
        }

        let quality_score_array = DataSeries::by_canonical_name(b"QS")
            .and_then(|series| encodings.get(series))
            .map(|descriptor| {
                create_encoding(
                    DataSeriesType::ByteArray,
                    descriptor.id,
                    &descriptor.parameters,
                )
            })
            .transpose()
            .map_err(RecordReadError::from)?;

        Ok(Self {
            map,
            quality_score_array,
        })
    }

    fn get(&self, name: &'static str) -> Result<&Encoding, RecordWriteError> {
        self.map
            .iter()
            .find(|(series, _)| *series == name)
            .map(|(_, encoding)| encoding)
            .ok_or(RecordWriteError::NoEncoding { name })
    }
}

const SERIES: [&str; 29] = [
    "BF", "CF", "RI", "RL", "AP", "RG", "RN", "MF", "NS", "NP", "TS", "NF", "TL", "FN", "FC", "FP",
    "DL", "BB", "QQ", "BS", "IN", "SC", "HC", "PD", "RS", "BA", "QS", "MQ", "TM",
];

fn series_type(series_type: crate::encoding_map::DataSeriesType) -> DataSeriesType {
    match series_type {
        crate::encoding_map::DataSeriesType::Byte => DataSeriesType::Byte,
        crate::encoding_map::DataSeriesType::Int => DataSeriesType::Int,
        crate::encoding_map::DataSeriesType::Long => DataSeriesType::Long,
        crate::encoding_map::DataSeriesType::ByteArray => DataSeriesType::ByteArray,
    }
}

/// `writeCRAMRecord`, in the reference's order.
///
/// `tag_ids_index` is the record's index into the compression header's tag dictionary, which the
/// slice assigns rather than the record carrying it from anywhere.
pub fn write_record(
    writers: &RecordWriters,
    header: &CompressionHeader,
    slice: &SliceContext,
    streams: &mut SliceWriteStreams,
    record: &CramRecord,
    tag_ids_index: i32,
    previous_alignment_start: i32,
) -> Result<(), RecordWriteError> {
    write_int(writers.get("BF")?, streams, record.flags.bam)?;
    write_int(writers.get("CF")?, streams, record.flags.cram_flags())?;

    if slice.is_multi_reference() {
        write_int(writers.get("RI")?, streams, record.reference_index)?;
    }

    write_int(writers.get("RL")?, streams, record.read_length)?;

    if header.preservation.ap_delta {
        let delta = record
            .alignment_start
            .wrapping_sub(previous_alignment_start);
        write_int(writers.get("AP")?, streams, delta)?;
    } else {
        write_int(writers.get("AP")?, streams, record.alignment_start)?;
    }

    write_int(writers.get("RG")?, streams, record.read_group_id)?;

    if header.preservation.preserve_read_names {
        write_name(writers, streams, record)?;
    }

    if record.flags.is_detached() {
        write_int(writers.get("MF")?, streams, record.flags.mate_flags())?;
        // The same rule the reader has: a generated name goes after the mate flags.
        if !header.preservation.preserve_read_names {
            write_name(writers, streams, record)?;
        }
        write_int(writers.get("NS")?, streams, record.mate_reference_index)?;
        write_int(writers.get("NP")?, streams, record.mate_alignment_start)?;
        write_int(writers.get("TS")?, streams, record.template_size)?;
    } else if record.flags.has_mate_downstream() {
        write_int(writers.get("NF")?, streams, record.records_to_next_fragment)?;
    }

    write_int(writers.get("TL")?, streams, tag_ids_index)?;
    // A record's tag values would go here, one series per tag id, in the dictionary's order.

    if !record.flags.is_segment_unmapped() {
        write_read_features(writers, streams, &record.read_features)?;
        write_int(writers.get("MQ")?, streams, record.mapping_quality)?;
        if record.flags.is_force_preserve_quality_scores() {
            write_quality_score_array(writers, streams, &record.quality_scores)?;
        }
    } else if !record.flags.is_unknown_bases() {
        for base in &record.read_bases {
            write_byte(writers.get("BA")?, streams, *base as i8)?;
        }
        if record.flags.is_force_preserve_quality_scores() {
            write_quality_score_array(writers, streams, &record.quality_scores)?;
        }
    }

    Ok(())
}

fn write_name(
    writers: &RecordWriters,
    streams: &mut SliceWriteStreams,
    record: &CramRecord,
) -> Result<(), RecordWriteError> {
    let name = record.read_name.clone().unwrap_or_default();
    write_byte_array(writers.get("RN")?, streams, name.as_bytes())?;
    Ok(())
}

fn write_quality_score_array(
    writers: &RecordWriters,
    streams: &mut SliceWriteStreams,
    scores: &[u8],
) -> Result<(), RecordWriteError> {
    let encoding = writers
        .quality_score_array
        .as_ref()
        .ok_or(RecordWriteError::NoEncoding { name: "QS" })?;
    write_byte_array(encoding, streams, scores)?;
    Ok(())
}

/// The features, whose positions go out as deltas from the previous one.
fn write_read_features(
    writers: &RecordWriters,
    streams: &mut SliceWriteStreams,
    features: &[ReadFeature],
) -> Result<(), RecordWriteError> {
    write_int(writers.get("FN")?, streams, features.len() as i32)?;

    let mut previous_position = 0i32;
    for feature in features {
        let (operator, position) = operator_and_position(feature);
        write_byte(writers.get("FC")?, streams, operator as i8)?;
        write_int(
            writers.get("FP")?,
            streams,
            position.wrapping_sub(previous_position),
        )?;
        previous_position = position;

        match feature {
            ReadFeature::ReadBase { base, quality, .. } => {
                write_byte(writers.get("BA")?, streams, *base as i8)?;
                write_byte(writers.get("QS")?, streams, *quality)?;
            }
            ReadFeature::Substitution { code, .. } => {
                // A negative code means the matrix has not assigned one yet, which the reference
                // resolves here against the compression header's matrix.
                if *code < 0 {
                    return Err(RecordWriteError::SubstitutionWithoutCode);
                }
                write_byte(writers.get("BS")?, streams, *code)?;
            }
            ReadFeature::Insertion { sequence, .. } => {
                write_byte_array(writers.get("IN")?, streams, sequence)?
            }
            ReadFeature::SoftClip { sequence, .. } => {
                write_byte_array(writers.get("SC")?, streams, sequence)?
            }
            ReadFeature::HardClip { length, .. } => {
                write_int(writers.get("HC")?, streams, *length)?
            }
            ReadFeature::Padding { length, .. } => write_int(writers.get("PD")?, streams, *length)?,
            ReadFeature::Deletion { length, .. } => {
                write_int(writers.get("DL")?, streams, *length)?
            }
            ReadFeature::RefSkip { length, .. } => write_int(writers.get("RS")?, streams, *length)?,
            ReadFeature::InsertBase { base, .. } => {
                write_byte(writers.get("BA")?, streams, *base as i8)?
            }
            ReadFeature::BaseQualityScore { quality, .. } => {
                write_byte(writers.get("QS")?, streams, *quality)?
            }
            // The two the reader has branches for and the writer does not.
            ReadFeature::Bases { .. } | ReadFeature::Scores { .. } => {
                return Err(RecordWriteError::UnknownReadFeature { operator })
            }
        }
    }
    Ok(())
}

fn operator_and_position(feature: &ReadFeature) -> (u8, i32) {
    match feature {
        ReadFeature::ReadBase { position, .. } => (b'B', *position),
        ReadFeature::Substitution { position, .. } => (b'X', *position),
        ReadFeature::Insertion { position, .. } => (b'I', *position),
        ReadFeature::SoftClip { position, .. } => (b'S', *position),
        ReadFeature::HardClip { position, .. } => (b'H', *position),
        ReadFeature::Padding { position, .. } => (b'P', *position),
        ReadFeature::Deletion { position, .. } => (b'D', *position),
        ReadFeature::RefSkip { position, .. } => (b'N', *position),
        ReadFeature::InsertBase { position, .. } => (b'i', *position),
        ReadFeature::BaseQualityScore { position, .. } => (b'Q', *position),
        ReadFeature::Bases { position, .. } => (b'b', *position),
        ReadFeature::Scores { position, .. } => (b'q', *position),
    }
}
