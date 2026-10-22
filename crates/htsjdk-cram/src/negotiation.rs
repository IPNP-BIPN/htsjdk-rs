//! Codec negotiation: which encoding a writer picks for each data series, and which compressor.
//!
//! Ported from `htsjdk.samtools.cram.build.CompressionHeaderFactory` and
//! `CompressionHeaderEncodingMap`'s default strategy at htsjdk 5.0.0.
//!
//! The reader takes what a file names. The writer chooses, and the choice is what makes one CRAM of
//! a set of records rather than another. Two things do the choosing: a fixed table for the data
//! series, and a measurement for the tag series, which are not known until the records are.
//!
//! # A tag's compressor is a trial
//!
//! Each tag gets a `TrialCompressor` over GZIP and rANS Nx16 order 0 for the life of the writer.
//! On the first blocks it returns the smaller output and adds both lengths to a running total;
//! after three it keeps the candidate with the smaller total for seventy blocks, then trials
//! again. Every comparison is the same one, [`trial_choice`]: rANS Nx16 only when strictly
//! smaller, so a tie goes to GZIP, the first candidate. (`SA` and `XA` add BZIP2, whose writer is
//! not ported, so the writer refuses them.)
//!
//! [`trial_choice`] takes lengths rather than data so the suite can feed it the lengths the
//! reference measured, and not depend on this crate's gzip matching the JDK's `Deflater`.
//!
//! # A tag's encoding comes from its type, and from the range of its values
//!
//! The fixed-width types get a length-prefixed array whose length is a single-symbol Huffman code
//! of zero bits, which writes no bits at all. The two variable-width types, `Z` and `B`, get the
//! **Every `Z` is a stop-byte array with a tab** for the stop byte, chosen rather than searched
//! for, so a `Z` whose text contains a tab is split by its own encoding; 4.2.0 gave a `Z` whose
//! values were all one size the fixed-size encoding instead. A `B` gets the fixed size when every
//! value is one size; otherwise one over a hundred bytes searches for a byte its data never uses,
//! and falls back to a length-prefixed array whose two halves are **both external and both on the
//! same content id**.

use crate::encoding_map::EncodingId;

/// A tag trial's candidates, in the order `getTagTrialCompressor` lists them, which is also their
/// index in the trial.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Compressor {
    Gzip,
    RansNx16,
}

impl Compressor {
    /// The `BlockCompressionMethod` constant's name.
    pub fn name(&self) -> &'static str {
        match self {
            Compressor::Gzip => "GZIP",
            Compressor::RansNx16 => "RANSNx16",
        }
    }
}

/// `TrialCompressor`'s comparison, between two outputs or two running totals: the first strictly
/// smallest candidate wins.
pub fn trial_choice<T: PartialOrd>(gzip: T, rans_nx16: T) -> Compressor {
    if rans_nx16 < gzip {
        Compressor::RansNx16
    } else {
        Compressor::Gzip
    }
}

/// `getUnusedByte`: the first byte value the data never uses, or `-1` if it uses all of them.
pub fn unused_byte(data: &[u8]) -> i32 {
    let mut used = [false; 256];
    for byte in data {
        used[*byte as usize] = true;
    }
    used.iter()
        .position(|seen| !seen)
        .map(|at| at as i32)
        .unwrap_or(-1)
}

/// `ALL_BYTES_USED`.
pub const ALL_BYTES_USED: i32 = -1;

/// The threshold above which a `B` tag looks for a stop byte at all.
pub const STOP_ENCODING_MIN_SIZE: usize = 100;

/// The stop byte a `Z` tag is given, chosen rather than searched for.
pub const Z_STOP_BYTE: u8 = b'\t';

/// An encoding as the factory produces it for a tag: an identifier and its parameters.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TagEncoding {
    pub id: EncodingId,
    pub parameters: Vec<u8>,
}

/// `buildTagEncodingForSize`: a byte-array-len whose length is a Huffman alphabet of one symbol at
/// zero bits, and whose bytes are external on the tag's own id.
///
/// A zero-bit code word writes nothing, so the length costs no bits at all: the size is constant
/// and the file says so once, in the encoding rather than in the data.
pub fn tag_encoding_for_size(size: i32, tag_id: i32) -> TagEncoding {
    let length = crate::huffman::serialize_integer_params(&[size], &[0]);
    let bytes = crate::external_codecs::serialize_external_params(tag_id);
    TagEncoding {
        id: EncodingId::ByteArrayLen,
        parameters: crate::external_codecs::serialize_byte_array_len_params(
            (EncodingId::Huffman, &length),
            (EncodingId::External, &bytes),
        ),
    }
}

/// A stop-byte array on the tag's own content id.
pub fn tag_encoding_stop(stop: u8, tag_id: i32) -> TagEncoding {
    TagEncoding {
        id: EncodingId::ByteArrayStop,
        parameters: crate::external_codecs::serialize_stop_params(stop, tag_id),
    }
}

/// A length-prefixed array whose two halves are both external and both on the tag's own id.
pub fn tag_encoding_both_external(tag_id: i32) -> TagEncoding {
    let id = crate::external_codecs::serialize_external_params(tag_id);
    TagEncoding {
        id: EncodingId::ByteArrayLen,
        parameters: crate::external_codecs::serialize_byte_array_len_params(
            (EncodingId::External, &id),
            (EncodingId::External, &id),
        ),
    }
}

/// `buildEncodingForTag`, without the compressor: the encoding a tag's type and value sizes give.
///
/// `sizes` is the byte size of every value of this tag in the container, and `data` all of their
/// bytes end to end, which is what the unused-byte search looks at.
pub fn tag_encoding(tag_type: u8, tag_id: i32, sizes: &[usize], data: &[u8]) -> TagEncoding {
    match tag_type {
        b'A' | b'c' | b'C' => tag_encoding_for_size(1, tag_id),
        b'I' | b'i' | b'f' => tag_encoding_for_size(4, tag_id),
        b's' | b'S' => tag_encoding_for_size(2, tag_id),
        b'Z' => tag_encoding_stop(Z_STOP_BYTE, tag_id),
        b'B' => {
            let min = sizes.iter().copied().min().unwrap_or(0);
            let max = sizes.iter().copied().max().unwrap_or(0);
            if min == max {
                return tag_encoding_for_size(min as i32, tag_id);
            }
            if min > STOP_ENCODING_MIN_SIZE {
                let unused = unused_byte(data);
                if unused > ALL_BYTES_USED {
                    return tag_encoding_stop(unused as u8, tag_id);
                }
            }
            tag_encoding_both_external(tag_id)
        }
        // The reference throws here; nothing in a valid file reaches it.
        _ => tag_encoding_both_external(tag_id),
    }
}

/// `buildTagIdDictionary`: the distinct sorted tag id lists a container's records use.
///
/// Two records whose tags are the same in a different order share one entry, because the ids are
/// sorted before they are compared. Group 0 is the empty list every container carries, and a
/// record with no tags takes it.
///
/// Returns the groups and the index each record was given, in the records' own order.
pub fn tag_id_dictionary(records: &[Vec<[u8; 3]>]) -> (Vec<Vec<[u8; 3]>>, Vec<i32>) {
    // A map that keeps insertion order, which is what the reference's LinkedHashMap does.
    let mut groups: Vec<Vec<[u8; 3]>> = vec![Vec::new()];
    let mut indexes = Vec::with_capacity(records.len());

    for tags in records {
        let mut sorted = tags.clone();
        sorted.sort_by_key(tag_id_as_int);
        let position = groups.iter().position(|group| *group == sorted);
        let index = match position {
            Some(index) => index,
            None => {
                groups.push(sorted);
                groups.len() - 1
            }
        };
        indexes.push(index as i32);
    }

    (groups, indexes)
}

/// `ReadTag.name3BytesToInt`: the three bytes of a tag id packed **high byte first**, so `XAc` is
/// 0x584163 and the type letter is the low byte, which is what [`tag_type`] reads back.
pub fn tag_id_as_int(id: &[u8; 3]) -> i32 {
    (i32::from(id[0]) << 16) | (i32::from(id[1]) << 8) | i32::from(id[2])
}

/// `getTagType`: the lowest byte of a packed tag id.
pub fn tag_type(tag_id: i32) -> u8 {
    (tag_id & 0xFF) as u8
}
