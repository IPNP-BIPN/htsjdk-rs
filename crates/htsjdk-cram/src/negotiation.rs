//! Codec negotiation: which encoding a writer picks for each data series, and which compressor.
//!
//! Ported from `htsjdk.samtools.cram.build.CompressionHeaderFactory` and
//! `CompressionHeaderEncodingMap`'s default strategy at htsjdk 4.2.0.
//!
//! The reader takes what a file names. The writer chooses, and the choice is what makes one CRAM of
//! a set of records rather than another. Two things do the choosing: a fixed table for the data
//! series, and a measurement for the tag series, which are not known until the records are.
//!
//! # The compressor is chosen by trying all three
//!
//! GZIP, rANS order 0 and rANS order 1 are each run over the data and the smallest output wins.
//! The tie-break is the order of the *comparisons*, rANS 0 then rANS 1 then GZIP, not the order of
//! the compressions: a thousand identical bytes compress to 29 under GZIP and 29 under rANS 0, and
//! rANS wins.
//!
//! The gzip length is not something this crate can produce: htsjdk's CRAM compressor is the JDK's
//! `Deflater`, whose output length is its zlib's business. So [`best_compressor`] takes the three
//! lengths rather than the data, and the suite feeds it the lengths the reference measured.
//!
//! # A tag's encoding comes from its type, and from the range of its values
//!
//! The fixed-width types get a length-prefixed array whose length is a single-symbol Huffman code
//! of zero bits, which writes no bits at all. The two variable-width types, `Z` and `B`, get the
//! same when every value is the same size, and otherwise: a `Z` becomes a stop-byte array with a
//! **tab** for the stop byte, chosen rather than searched for, so a `Z` whose text contains a tab
//! is split by its own encoding. A `B` over a hundred bytes searches for a byte its data never
//! uses, and falls back to a length-prefixed array whose two halves are **both external and both
//! on the same content id**, which the reference's own comment marks with three question marks.

use crate::encoding_map::EncodingId;

/// `BlockCompressionMethod`, as far as negotiation is concerned.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Compressor {
    Gzip,
    Rans,
}

impl Compressor {
    pub fn name(&self) -> &'static str {
        match self {
            Compressor::Gzip => "GZIP",
            Compressor::Rans => "RANS",
        }
    }
}

/// `getBestExternalCompressor`, given the three lengths it would have measured.
///
/// The comparison is `min(gzip, min(rans0, rans1))` and then three equality tests in the order
/// rANS 0, rANS 1, GZIP, so a tie goes to rANS however it arose.
pub fn best_compressor(gzip: usize, rans0: usize, rans1: usize) -> Compressor {
    let smallest = gzip.min(rans0.min(rans1));
    if smallest == rans0 || smallest == rans1 {
        Compressor::Rans
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
        b'Z' | b'B' => {
            let min = sizes.iter().copied().min().unwrap_or(0);
            let max = sizes.iter().copied().max().unwrap_or(0);
            if min == max {
                return tag_encoding_for_size(min as i32, tag_id);
            }
            if tag_type == b'Z' {
                return tag_encoding_stop(Z_STOP_BYTE, tag_id);
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
