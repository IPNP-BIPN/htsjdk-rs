//! The CRAM compression header's preservation map, the first of its three length-prefixed maps.
//!
//! Ported from `htsjdk.samtools.cram.structure.CompressionHeader`'s `internalRead`,
//! `internalWrite`, `parseDictionary` and `dictionaryToByteArray` at htsjdk 4.2.0.
//!
//! [`crate::block`] measured the compression header as the RAW block that follows the GZIP header
//! block in every CRAM. This is what is at the front of it: whether read names are kept, whether
//! alignment positions are deltas, whether a reference is required, the substitution matrix, and
//! which combinations of tag ids appear in the file.
//!
//! # The map size is a hardcoded 5, not a count
//!
//! `internalWrite` calls `ITF8.writeUnsignedITF8(5, mapBuffer)` and then writes exactly `RN`, `AP`,
//! `RR`, `SM` and `TD` in that order, whatever the header holds. The field is a constant wearing a
//! count's clothes, and the write order is htsjdk's rather than the specification's. The reader
//! accepts the keys in any order, so the order is only visible in the bytes, which is exactly what
//! a byte-identical port has to reproduce.
//!
//! # A boolean is `== 1`, not `!= 0`
//!
//! `preserveReadNames = buffer.get() == 1`. Measured over 0, 1, 2, 127 and 255: only 1 is true, and
//! nothing else raises anything. Three of the five keys go through it.
//!
//! # The substitution matrix and the tag dictionary are mandatory
//!
//! The check is after the loop rather than inside it, and it names both keys whichever one is
//! missing, so the message never says which. Measured: omitting either gives
//! `substitution matrix and tag ID dictionary must be present in the compression header`.
//!
//! # An unknown key is a plain `RuntimeException`
//!
//! Not a `CRAMException`, and it carries the two characters it did not recognise:
//! `Unknown preservation map key: ZZ`.
//!
//! # The first dictionary group is always empty
//!
//! Measured on a file whose every record carries three tags: the dictionary is
//! `00 4d445a 4e4d63 585866 00`, an empty group and then `MDZ, NMc, XXf`. Group 0 is the record
//! with no tags, and it is written even when no record is without tags.
//!
//! # The whole header does not depend on the reads
//!
//! Measured on four files differing in record count and read length: the compression header block's
//! raw content is the same 160 bytes and the same digest in all four. Only the tags moved it.

/// `SubstitutionMatrix.BASES_SIZE`, measured rather than assumed: five bases, five bytes.
pub const BASES_SIZE: usize = 5;

/// The count `internalWrite` writes, which is a constant and not a count.
pub const WRITTEN_MAP_SIZE: i32 = 5;

/// The five keys, in the order htsjdk writes them.
pub const RN: &[u8; 2] = b"RN";
pub const AP: &[u8; 2] = b"AP";
pub const RR: &[u8; 2] = b"RR";
pub const SM: &[u8; 2] = b"SM";
pub const TD: &[u8; 2] = b"TD";

use crate::varint::{read_unsigned_itf8, RuntimeEof};

/// What a preservation map is refused with.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PreservationMapError {
    /// `RuntimeException`, not a `CRAMException`, carrying the two bytes it did not recognise.
    UnknownKey([u8; 2]),
    /// The `CRAMException` thrown after the loop when either mandatory key is absent. It names both
    /// whichever one is missing.
    MissingMatrixOrDictionary,
    /// The map ends inside one of its own values.
    Truncated,
}

impl PreservationMapError {
    pub fn message(&self) -> String {
        match self {
            PreservationMapError::UnknownKey(key) => format!(
                "Unknown preservation map key: {}",
                String::from_utf8_lossy(key)
            ),
            PreservationMapError::MissingMatrixOrDictionary => {
                "substitution matrix and tag ID dictionary must be present in the compression header"
                    .to_string()
            }
            PreservationMapError::Truncated => {
                "the preservation map ends inside one of its values".to_string()
            }
        }
    }
}

impl From<RuntimeEof> for PreservationMapError {
    fn from(_: RuntimeEof) -> Self {
        PreservationMapError::Truncated
    }
}

/// The preservation map, as htsjdk holds it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PreservationMap {
    /// `RN`. Defaults to true, which is what an absent key leaves it as.
    pub preserve_read_names: bool,
    /// `AP`. Measured false on every unsorted file, because the writer writes what it used.
    pub ap_delta: bool,
    /// `RR`.
    pub reference_required: bool,
    /// `SM`, exactly [`BASES_SIZE`] bytes.
    pub substitution_matrix: [u8; BASES_SIZE],
    /// `TD`, parsed into groups. Group 0 is the empty one every file carries.
    pub tag_id_dictionary: Vec<Vec<[u8; 3]>>,
}

impl Default for PreservationMap {
    /// The three flags default to **true**, which is what the reader leaves them as when a key is
    /// absent. The two mandatory keys have no default: a map without them is refused.
    fn default() -> Self {
        Self {
            preserve_read_names: true,
            ap_delta: true,
            reference_required: true,
            substitution_matrix: [0u8; BASES_SIZE],
            tag_id_dictionary: Vec::new(),
        }
    }
}

impl PreservationMap {
    /// Read the map given its own bytes, without the length prefix that precedes it in the stream.
    pub fn read(map: &[u8]) -> Result<Self, PreservationMapError> {
        let mut at = 0usize;
        let itf8 = |at: &mut usize| -> Result<i32, PreservationMapError> {
            let (value, consumed) = read_unsigned_itf8(&map[(*at).min(map.len())..])?;
            *at += consumed;
            Ok(value)
        };
        let count = itf8(&mut at)?;

        let mut out = Self::default();
        let mut saw_matrix = false;
        let mut saw_dictionary = false;

        for _ in 0..count.max(0) {
            let key: [u8; 2] = map
                .get(at..at + 2)
                .ok_or(PreservationMapError::Truncated)?
                .try_into()
                .expect("two bytes");
            at += 2;
            match &key {
                b"RN" | b"AP" | b"RR" => {
                    // `buffer.get() == 1`, so everything else is false and nothing complains.
                    let value = *map.get(at).ok_or(PreservationMapError::Truncated)? == 1;
                    at += 1;
                    match &key {
                        b"RN" => out.preserve_read_names = value,
                        b"AP" => out.ap_delta = value,
                        _ => out.reference_required = value,
                    }
                }
                b"SM" => {
                    out.substitution_matrix = map
                        .get(at..at + BASES_SIZE)
                        .ok_or(PreservationMapError::Truncated)?
                        .try_into()
                        .expect("five bytes");
                    at += BASES_SIZE;
                    saw_matrix = true;
                }
                b"TD" => {
                    let size = itf8(&mut at)?.max(0) as usize;
                    let bytes = map
                        .get(at..at + size)
                        .ok_or(PreservationMapError::Truncated)?;
                    at += size;
                    out.tag_id_dictionary = parse_dictionary(bytes);
                    saw_dictionary = true;
                }
                _ => return Err(PreservationMapError::UnknownKey(key)),
            }
        }

        // After the loop, and it names both keys whichever one is missing.
        if !saw_matrix || !saw_dictionary {
            return Err(PreservationMapError::MissingMatrixOrDictionary);
        }
        Ok(out)
    }

    /// Read the map from the head of a compression header's content, past its length prefix.
    /// Returns the map and how many bytes it occupied, prefix included.
    pub fn read_prefixed(content: &[u8]) -> Result<(Self, usize), PreservationMapError> {
        let (size, consumed) = read_unsigned_itf8(content)?;
        let size = size.max(0) as usize;
        let bytes = content
            .get(consumed..consumed + size)
            .ok_or(PreservationMapError::Truncated)?;
        Ok((Self::read(bytes)?, consumed + size))
    }
}

/// `CompressionHeader.parseDictionary`: three bytes at a time until a zero, per group.
///
/// Nothing checks that the terminator falls on a three-byte boundary, and htsjdk computes a
/// `maxWidth` here that it never uses.
pub fn parse_dictionary(bytes: &[u8]) -> Vec<Vec<[u8; 3]>> {
    let mut groups = Vec::new();
    let mut at = 0usize;
    while at < bytes.len() {
        let mut group = Vec::new();
        while at < bytes.len() && bytes[at] != 0 {
            let end = (at + 3).min(bytes.len());
            let mut id = [0u8; 3];
            id[..end - at].copy_from_slice(&bytes[at..end]);
            group.push(id);
            at += 3;
        }
        at += 1;
        groups.push(group);
    }
    groups
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::varint::write_unsigned_itf8;

    /// `write_unsigned_itf8` returns the bytes and the bit count htsjdk's writer returns; the
    /// tests below build maps by hand and only want the bytes.
    fn push_itf8(value: i32, out: &mut Vec<u8>) {
        out.extend_from_slice(&write_unsigned_itf8(value).0);
    }

    /// Only 1 is true, and nothing else is an error.
    /// Only 1 is true, and nothing else is an error.
    #[test]
    fn a_boolean_is_one_rather_than_non_zero() {
        for (value, expected) in [
            (0u8, false),
            (1, true),
            (2, false),
            (127, false),
            (255, false),
        ] {
            let mut map = Vec::new();
            push_itf8(5, &mut map);
            map.extend_from_slice(RN);
            map.push(value);
            map.extend_from_slice(AP);
            map.push(1);
            map.extend_from_slice(RR);
            map.push(1);
            map.extend_from_slice(SM);
            map.extend_from_slice(&[0u8; BASES_SIZE]);
            map.extend_from_slice(TD);
            push_itf8(4, &mut map);
            map.extend_from_slice(b"NMi\x00");

            let read = PreservationMap::read(&map).expect("parses");
            assert_eq!(read.preserve_read_names, expected, "RN = {value}");
        }
    }

    /// The message names both keys whichever one is missing.
    #[test]
    fn the_two_mandatory_keys_share_one_message() {
        let mut map = Vec::new();
        push_itf8(1, &mut map);
        map.extend_from_slice(RN);
        map.push(1);
        assert_eq!(
            PreservationMap::read(&map),
            Err(PreservationMapError::MissingMatrixOrDictionary)
        );
        assert!(PreservationMapError::MissingMatrixOrDictionary
            .message()
            .contains("substitution matrix and tag ID dictionary"));
    }

    #[test]
    fn an_unknown_key_carries_the_two_characters_it_did_not_recognise() {
        let mut map = Vec::new();
        push_itf8(1, &mut map);
        map.extend_from_slice(b"ZZ");
        map.push(1);
        let error = PreservationMap::read(&map).expect_err("refused");
        assert_eq!(error, PreservationMapError::UnknownKey(*b"ZZ"));
        assert_eq!(error.message(), "Unknown preservation map key: ZZ");
    }

    /// Group 0 is the record with no tags, and it is there even when every record has tags.
    #[test]
    fn the_first_dictionary_group_is_empty() {
        let groups = parse_dictionary(b"\x00MDZNMcXXf\x00");
        assert_eq!(groups.len(), 2);
        assert!(groups[0].is_empty());
        assert_eq!(groups[1], vec![*b"MDZ", *b"NMc", *b"XXf"]);
    }
}
