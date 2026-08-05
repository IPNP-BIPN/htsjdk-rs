//! The CRAM compression header's data series encoding map: the second of its three maps, and the
//! first thing in CRAM that describes a record rather than a container.
//!
//! Ported from `htsjdk.samtools.cram.structure.CompressionHeaderEncodingMap`, `DataSeries`,
//! `DataSeriesType`, `EncodingID` and `EncodingDescriptor` at htsjdk 4.2.0.
//!
//! [`crate::preservation_map`] pinned the first map. This one says, for each of the record's data
//! series, which of the ten encodings carries it and with what parameters.
//!
//! # This map's size is a real count, where the preservation map's is the literal 5
//!
//! Two maps in the same header, one counted and one not. The count here excludes any data series
//! whose encoding is `NULL`.
//!
//! # The write order is the enum's ordinal order, not the order the constructor populates
//!
//! The map is a `TreeMap` keyed by the enum, so it sorts by ordinal; the constructor adds the
//! series alphabetically by canonical name. The two orders differ and only the first is in the
//! bytes. Measured: `BF, CF, RI, RL, AP, RG, RN, NF, MF, NS, NP, TS, TL, MQ, FN, FP, FC, BA, QS,
//! BS, IN, DL, RS, SC, PD, HC`.
//!
//! # htsjdk writes 26 of the 32 data series
//!
//! `BB` and `QQ` are never written by this implementation, `TC` and `TN` are obsolete, and `TM` and
//! `TV` exist only for tests. A port that writes all 32 writes a map no htsjdk-written CRAM
//! contains.
//!
//! # `TC` and `TN` are read and then dropped
//!
//! A CRAM from another writer that carries them gets a log warning and the entries never reach the
//! map, so **the reader's map can hold fewer entries than the count the file declared**.
//!
//! # The content ids are htsjdk's, not the specification's
//!
//! The spec does not prescribe them. This implementation numbers the data series 1 to 32 in enum
//! order, and a reader has to discover them from the map rather than assume them.
//!
//! # An unknown encoding id is an array index
//!
//! `EncodingID.values()[buffer.get()]` has no bounds check, so a tenth encoding is an
//! `ArrayIndexOutOfBoundsException` rather than a CRAM error. And `buffer.get()` is **signed**, so
//! an id byte of 255 is index -1. Measured: `Index 10 out of bounds for length 10` and
//! `Index -1 out of bounds for length 10`.

use crate::varint::RuntimeEof;

/// `DataSeriesType`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DataSeriesType {
    Byte,
    Int,
    Long,
    ByteArray,
}

impl DataSeriesType {
    /// The enum constant's name, which is what the golden records.
    pub fn name(&self) -> &'static str {
        match self {
            DataSeriesType::Byte => "BYTE",
            DataSeriesType::Int => "INT",
            DataSeriesType::Long => "LONG",
            DataSeriesType::ByteArray => "BYTE_ARRAY",
        }
    }
}

/// `DataSeries`, in declaration order, which is the order this map is written in.
///
/// The content id is this implementation's choice and not the specification's.
pub const DATA_SERIES: [(&str, DataSeriesType, i32); 32] = [
    ("BF", DataSeriesType::Int, 1),
    ("CF", DataSeriesType::Int, 2),
    ("RI", DataSeriesType::Int, 3),
    ("RL", DataSeriesType::Int, 4),
    ("AP", DataSeriesType::Int, 5),
    ("RG", DataSeriesType::Int, 6),
    ("RN", DataSeriesType::ByteArray, 7),
    ("NF", DataSeriesType::Int, 8),
    ("MF", DataSeriesType::Int, 9),
    ("NS", DataSeriesType::Int, 10),
    ("NP", DataSeriesType::Int, 11),
    ("TS", DataSeriesType::Int, 12),
    ("TL", DataSeriesType::Int, 13),
    ("TC", DataSeriesType::Int, 14),
    ("TN", DataSeriesType::Int, 15),
    ("MQ", DataSeriesType::Int, 16),
    ("FN", DataSeriesType::Int, 17),
    ("FP", DataSeriesType::Int, 18),
    ("FC", DataSeriesType::Byte, 19),
    ("BB", DataSeriesType::ByteArray, 20),
    ("QQ", DataSeriesType::ByteArray, 21),
    ("BA", DataSeriesType::Byte, 22),
    ("QS", DataSeriesType::Byte, 23),
    ("BS", DataSeriesType::Byte, 24),
    ("IN", DataSeriesType::ByteArray, 25),
    ("DL", DataSeriesType::Int, 26),
    ("RS", DataSeriesType::Int, 27),
    ("SC", DataSeriesType::ByteArray, 28),
    ("PD", DataSeriesType::Int, 29),
    ("HC", DataSeriesType::Int, 30),
    ("TM", DataSeriesType::Int, 31),
    ("TV", DataSeriesType::Int, 32),
];

/// `CompressionHeaderEncodingMap.DATASERIES_NOT_READ_BY_HTSJDK`: read, warned about, dropped.
pub const NOT_READ: [&str; 2] = ["TC", "TN"];

/// One data series, held as its ordinal so the map sorts the way the `TreeMap` does.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct DataSeries(pub usize);

impl DataSeries {
    /// `DataSeries.byCanonicalName`.
    pub fn by_canonical_name(name: &[u8; 2]) -> Option<Self> {
        DATA_SERIES
            .iter()
            .position(|(canonical, _, _)| canonical.as_bytes() == name)
            .map(DataSeries)
    }

    pub fn canonical_name(&self) -> &'static str {
        DATA_SERIES[self.0].0
    }

    pub fn series_type(&self) -> DataSeriesType {
        DATA_SERIES[self.0].1
    }

    /// The content id htsjdk assigns on write, which a reader must not assume.
    pub fn content_id(&self) -> i32 {
        DATA_SERIES[self.0].2
    }

    pub fn is_read_by_htsjdk(&self) -> bool {
        !NOT_READ.contains(&self.canonical_name())
    }
}

/// `EncodingID`, by the id the byte carries.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EncodingId {
    Null = 0,
    External = 1,
    Golomb = 2,
    Huffman = 3,
    ByteArrayLen = 4,
    ByteArrayStop = 5,
    Beta = 6,
    Subexponential = 7,
    GolombRice = 8,
    Gamma = 9,
}

/// How many there are, which is the bound the reader's array index is checked against by the JVM
/// and by nothing else.
pub const ENCODING_ID_COUNT: i32 = 10;

impl EncodingId {
    pub fn from_id(id: i32) -> Option<Self> {
        Some(match id {
            0 => EncodingId::Null,
            1 => EncodingId::External,
            2 => EncodingId::Golomb,
            3 => EncodingId::Huffman,
            4 => EncodingId::ByteArrayLen,
            5 => EncodingId::ByteArrayStop,
            6 => EncodingId::Beta,
            7 => EncodingId::Subexponential,
            8 => EncodingId::GolombRice,
            9 => EncodingId::Gamma,
            _ => return None,
        })
    }

    pub fn name(&self) -> &'static str {
        match self {
            EncodingId::Null => "NULL",
            EncodingId::External => "EXTERNAL",
            EncodingId::Golomb => "GOLOMB",
            EncodingId::Huffman => "HUFFMAN",
            EncodingId::ByteArrayLen => "BYTE_ARRAY_LEN",
            EncodingId::ByteArrayStop => "BYTE_ARRAY_STOP",
            EncodingId::Beta => "BETA",
            EncodingId::Subexponential => "SUBEXPONENTIAL",
            EncodingId::GolombRice => "GOLOMB_RICE",
            EncodingId::Gamma => "GAMMA",
        }
    }

    /// `EncodingID.isExternalEncoding`, which is true for `BYTE_ARRAY_LEN` even though that one can
    /// use a core sub-encoding for its length.
    pub fn is_external(&self) -> bool {
        matches!(
            self,
            EncodingId::External | EncodingId::ByteArrayLen | EncodingId::ByteArrayStop
        )
    }
}

/// `EncodingDescriptor`: an id and an untyped blob of parameters.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EncodingDescriptor {
    pub id: EncodingId,
    pub parameters: Vec<u8>,
}

/// What an encoding map is refused with.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EncodingMapError {
    /// `CRAMException` from `DataSeries.byCanonicalName`. The reader always takes exactly two bytes
    /// for a name, so the "exactly two characters" branch is unreachable from here and an
    /// unprintable byte is simply a name that does not exist.
    UnknownDataSeries([u8; 2]),
    /// The `ArrayIndexOutOfBoundsException` from `EncodingID.values()[buffer.get()]`. The index is
    /// what a **signed** byte read gives, so 255 arrives as -1.
    EncodingIdOutOfBounds(i32),
    /// The map ends inside one of its own entries.
    Truncated,
}

impl EncodingMapError {
    pub fn message(&self) -> String {
        match self {
            EncodingMapError::UnknownDataSeries(name) => format!(
                "Could not find Data Series Encoding for: {}",
                String::from_utf8_lossy(name)
            ),
            EncodingMapError::EncodingIdOutOfBounds(index) => {
                format!("Index {index} out of bounds for length {ENCODING_ID_COUNT}")
            }
            EncodingMapError::Truncated => {
                "the encoding map ends inside one of its entries".to_string()
            }
        }
    }
}

impl From<RuntimeEof> for EncodingMapError {
    fn from(_: RuntimeEof) -> Self {
        EncodingMapError::Truncated
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The content ids are this implementation's, and they run 1 to 32 in enum order.
    /// The content ids are this implementation's, and they run 1 to 32 in enum order.
    #[test]
    fn the_content_ids_are_the_ordinal_plus_one() {
        for (ordinal, (_, _, content_id)) in DATA_SERIES.iter().enumerate() {
            assert_eq!(*content_id, ordinal as i32 + 1);
        }
    }
}
