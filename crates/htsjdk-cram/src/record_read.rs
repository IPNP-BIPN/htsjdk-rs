//! Reading a CRAM record: the data series, in the order the specification prescribes.
//!
//! Ported from `htsjdk.samtools.cram.encoding.reader.CramRecordReader.readCRAMRecord` at
//! htsjdk 4.2.0.
//!
//! Every field of a record comes from its own data series, and the series share streams. So the
//! order the reads happen in is not a detail of the implementation: read two series in the wrong
//! order and both come back wrong, with nothing to say so. That is why this module is a single
//! function in the reference's own order rather than a set of independent field readers.
//!
//! # The read name moves
//!
//! A preserved read name is read before the mate block. A generated one is read *inside* it, after
//! the mate flags, and only when the record is detached. The specification says so explicitly, and
//! the reference's comment repeats it.
//!
//! # The mate flags are propagated into the BAM flags
//!
//! Two bits, at different positions in the two words, because a writer is not required to have put
//! them in the BAM flags at all.
//!
//! # An unmapped record reads no read features
//!
//! Not zero of them: the series is not consulted. Its bases and scores come from elsewhere.

use crate::codecs::{read_byte, read_byte_array, read_int, ReadError};
use crate::compression_header::CompressionHeader;
use crate::encoding_factory::{create_encoding, DataSeriesType, Encoding, FactoryError};
use crate::encoding_map::{DataSeries, EncodingMap};
use crate::external_codecs::SliceReadStreams;
use crate::read_features::ReadFeature;
use crate::record_flags::{Flags, MATE_REVERSE_STRAND, MATE_UNMAPPED};

/// What reading a record refuses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RecordReadError {
    /// A data series the compression header does not name.
    MissingDataSeries {
        name: &'static str,
    },
    /// A tag list index past the end of the header's dictionary.
    TagListOutOfBounds {
        index: i32,
        length: usize,
    },
    /// A read feature operator the reader has no branch for.
    UnknownReadFeature {
        operator: u8,
    },
    Encoding(FactoryError),
    Read(ReadError),
}

impl RecordReadError {
    pub fn message(&self) -> String {
        match self {
            RecordReadError::MissingDataSeries { name } => {
                format!("Could not find Data Series Encoding for: {name}")
            }
            RecordReadError::TagListOutOfBounds { index, length } => {
                format!("Index {index} out of bounds for length {length}")
            }
            RecordReadError::UnknownReadFeature { operator } => {
                format!("Unknown read feature: {}", *operator as char)
            }
            RecordReadError::Encoding(error) => error.message(),
            RecordReadError::Read(error) => error.message(),
        }
    }

    pub fn java_exception(&self) -> &'static str {
        match self {
            RecordReadError::MissingDataSeries { .. } => "CRAMException",
            RecordReadError::TagListOutOfBounds { .. } => "ArrayIndexOutOfBoundsException",
            RecordReadError::UnknownReadFeature { .. } => "RuntimeException",
            RecordReadError::Encoding(error) => error.java_exception(),
            RecordReadError::Read(error) => error.java_exception(),
        }
    }
}

impl From<FactoryError> for RecordReadError {
    fn from(error: FactoryError) -> Self {
        RecordReadError::Encoding(error)
    }
}

impl From<ReadError> for RecordReadError {
    fn from(error: ReadError) -> Self {
        RecordReadError::Read(error)
    }
}

/// A record as the reader produces it, before normalisation.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct CramRecord {
    pub flags: Flags,
    pub read_name: Option<String>,
    pub read_length: i32,
    pub reference_index: i32,
    pub alignment_start: i32,
    pub read_group_id: i32,
    pub mate_reference_index: i32,
    pub mate_alignment_start: i32,
    pub template_size: i32,
    pub records_to_next_fragment: i32,
    pub mapping_quality: i32,
    pub read_features: Vec<ReadFeature>,
    /// The tag ids the dictionary gave, three bytes each, in its order.
    pub tags: Vec<[u8; 3]>,
}

/// What a slice knows about itself that the reader needs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SliceContext {
    /// The reference this slice is on, `-2` for multi-reference and `-1` for unmapped.
    pub reference_context: i32,
    pub alignment_start: i32,
}

impl SliceContext {
    pub const MULTIPLE_REFERENCE: i32 = -2;

    pub fn is_multi_reference(&self) -> bool {
        self.reference_context == Self::MULTIPLE_REFERENCE
    }
}

/// The encodings a record is read through, resolved once from the compression header.
///
/// Resolving them per record would be the same bytes and a great deal slower, but that is not why
/// they are held: the reference builds its readers once in the constructor, and a codec that holds
/// state across reads would behave differently if it were rebuilt.
pub struct RecordReaders {
    map: Vec<(&'static str, Encoding)>,
}

impl RecordReaders {
    /// Resolve every data series the reader can consult. A series the header does not name is left
    /// out, and asking for it later is the refusal the reference raises when it reaches one.
    pub fn new(encodings: &EncodingMap) -> Result<Self, RecordReadError> {
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
                )?;
                map.push((name, encoding));
            }
        }
        Ok(Self { map })
    }

    fn get(&self, name: &'static str) -> Result<&Encoding, RecordReadError> {
        self.map
            .iter()
            .find(|(series, _)| *series == name)
            .map(|(_, encoding)| encoding)
            .ok_or(RecordReadError::MissingDataSeries { name })
    }
}

/// The data series a record can be read from, by their two-character names.
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

/// `readCRAMRecord`, in the reference's order.
///
/// `previous_alignment_start` is the slice's own alignment start for the first record and the
/// previous record's for every one after it, which is what `Slice.getRecords` passes.
pub fn read_record(
    readers: &RecordReaders,
    header: &CompressionHeader,
    slice: &SliceContext,
    streams: &mut SliceReadStreams<'_>,
    previous_alignment_start: i32,
) -> Result<CramRecord, RecordReadError> {
    let mut record = CramRecord::default();

    record.flags.bam = read_int(readers.get("BF")?, streams)?;
    record.flags.cram = read_int(readers.get("CF")?, streams)?;

    // The reference index is only in the stream for a multi-reference slice; otherwise it is the
    // slice's own, whether that is a reference or the -1 of an unmapped slice.
    record.reference_index = if slice.is_multi_reference() {
        read_int(readers.get("RI")?, streams)?
    } else {
        slice.reference_context
    };

    record.read_length = read_int(readers.get("RL")?, streams)?;

    record.alignment_start = if header.preservation.ap_delta {
        // A negative delta is legal, and the corpus carries one.
        previous_alignment_start.wrapping_add(read_int(readers.get("AP")?, streams)?)
    } else {
        read_int(readers.get("AP")?, streams)?
    };

    record.read_group_id = read_int(readers.get("RG")?, streams)?;

    if header.preservation.preserve_read_names {
        record.read_name = Some(read_name(readers, streams)?);
    }

    record.mate_reference_index = -1;
    record.records_to_next_fragment = -1;

    if record.flags.is_detached() {
        record.flags.mate = read_int(readers.get("MF")?, streams)?;
        // A writer need not have put these two in the BAM flags, so they are propagated here. Note
        // the positions differ between the two words.
        if record.flags.mate & crate::record_flags::MF_MATE_NEG_STRAND != 0 {
            record.flags.bam |= MATE_REVERSE_STRAND;
        }
        if record.flags.mate & crate::record_flags::MF_MATE_UNMAPPED != 0 {
            record.flags.bam |= MATE_UNMAPPED;
        }
        // The specification prescribes that a generated read name is decoded AFTER the mate flags.
        if !header.preservation.preserve_read_names {
            record.read_name = Some(read_name(readers, streams)?);
        }
        record.mate_reference_index = read_int(readers.get("NS")?, streams)?;
        record.mate_alignment_start = read_int(readers.get("NP")?, streams)?;
        record.template_size = read_int(readers.get("TS")?, streams)?;
    } else if record.flags.has_mate_downstream() {
        record.records_to_next_fragment = read_int(readers.get("NF")?, streams)?;
    }

    // The tag list is an index into the header's dictionary, and the tags follow in its order.
    let tag_list = read_int(readers.get("TL")?, streams)?;
    let dictionary = &header.preservation.tag_id_dictionary;
    let ids =
        dictionary
            .get(tag_list.max(0) as usize)
            .ok_or(RecordReadError::TagListOutOfBounds {
                index: tag_list,
                length: dictionary.len(),
            })?;
    for id in ids {
        // One data series per tag id, named by the id itself rather than by a canonical name.
        let descriptor =
            header
                .tag_encodings
                .get(tag_id(id))
                .ok_or(RecordReadError::MissingDataSeries {
                    name: "a tag series",
                })?;
        let encoding = create_encoding(
            DataSeriesType::ByteArray,
            descriptor.id,
            &descriptor.parameters,
        )?;
        read_byte_array(&encoding, streams, None)?;
        record.tags.push(*id);
    }

    if !record.flags.is_segment_unmapped() {
        record.read_features = read_read_features(readers, streams)?;
        record.mapping_quality = read_int(readers.get("MQ")?, streams)?;
    }

    Ok(record)
}

/// `ReadTag.name3BytesToInt`: the three bytes of a tag id, packed low byte first.
fn tag_id(id: &[u8; 3]) -> i32 {
    i32::from(id[0]) | (i32::from(id[1]) << 8) | (i32::from(id[2]) << 16)
}

fn read_name(
    readers: &RecordReaders,
    streams: &mut SliceReadStreams<'_>,
) -> Result<String, RecordReadError> {
    let bytes = read_byte_array(readers.get("RN")?, streams, None)?;
    Ok(String::from_utf8_lossy(&bytes).to_string())
}

/// The read features, whose positions are deltas from the previous one.
fn read_read_features(
    readers: &RecordReaders,
    streams: &mut SliceReadStreams<'_>,
) -> Result<Vec<ReadFeature>, RecordReadError> {
    let count = read_int(readers.get("FN")?, streams)?;
    let mut features = Vec::new();
    let mut previous_position = 0i32;

    for _ in 0..count.max(0) {
        let operator = read_byte(readers.get("FC")?, streams)? as u8;
        let position = previous_position.wrapping_add(read_int(readers.get("FP")?, streams)?);
        previous_position = position;

        features.push(match operator {
            b'B' => ReadFeature::ReadBase {
                position,
                base: read_byte(readers.get("BA")?, streams)? as u8,
                quality: read_byte(readers.get("QS")?, streams)?,
            },
            b'X' => ReadFeature::Substitution {
                position,
                base: 0,
                reference_base: 0,
                code: read_byte(readers.get("BS")?, streams)?,
            },
            b'I' => ReadFeature::Insertion {
                position,
                sequence: read_byte_array(readers.get("IN")?, streams, None)?,
            },
            b'S' => ReadFeature::SoftClip {
                position,
                sequence: read_byte_array(readers.get("SC")?, streams, None)?,
            },
            b'H' => ReadFeature::HardClip {
                position,
                length: read_int(readers.get("HC")?, streams)?,
            },
            b'P' => ReadFeature::Padding {
                position,
                length: read_int(readers.get("PD")?, streams)?,
            },
            b'D' => ReadFeature::Deletion {
                position,
                length: read_int(readers.get("DL")?, streams)?,
            },
            b'N' => ReadFeature::RefSkip {
                position,
                length: read_int(readers.get("RS")?, streams)?,
            },
            b'i' => ReadFeature::InsertBase {
                position,
                base: read_byte(readers.get("BA")?, streams)? as u8,
            },
            b'Q' => ReadFeature::BaseQualityScore {
                position,
                quality: read_byte(readers.get("QS")?, streams)?,
            },
            b'b' => ReadFeature::Bases {
                position,
                bases: read_byte_array(readers.get("BB")?, streams, None)?,
            },
            b'q' => ReadFeature::Scores {
                position,
                scores: read_byte_array(readers.get("QQ")?, streams, None)?,
            },
            other => return Err(RecordReadError::UnknownReadFeature { operator: other }),
        });
    }

    Ok(features)
}
