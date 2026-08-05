//! The CRAM compression header's tag encoding map: the third of its maps, and the one that closes
//! the header.
//!
//! Ported from `htsjdk.samtools.cram.structure.CompressionHeader`'s tag encoding map and from
//! `htsjdk.samtools.cram.structure.ReadTag` at htsjdk 4.2.0.
//!
//! [`crate::encoding_map`] pinned the second map, whose keys are two-character data series names.
//! This one has the same entry shape and integer keys, and the integer **is the tag**: two
//! characters of name and one of type, packed into twenty-four bits.
//!
//! # The key is the tag, packed
//!
//! `nameType3BytesToInt` shifts `name[0]`, `name[1]` and the **type** into one int, so the key
//! carries the whole identity of the tag and the map needs no separate name field. Measured: `NMc`
//! is 5131619, `MDZ` is 5063770, and three spaces are 2105376.
//!
//! # The type is part of the key, so one name at two types is two entries
//!
//! Measured on a file carrying `XX` as a string on half its records and as an integer on the other
//! half: the map holds **two** entries, `XXZ` and `XXi`, with two external blocks.
//!
//! # The write order is numeric order of that key
//!
//! The map is a `TreeMap<Integer>`. Measured: a file whose records introduce `XX`, `NM`, `MD` in
//! that order produces **byte-identical** output to one that introduces `MD`, `NM`, `XX`, and both
//! write `MDZ, NMc, XXf`. The order the records arrived in is not in the file.
//!
//! # The collision guard cannot fire
//!
//! `putTagBlockCompression` refuses a tag id equal to a data series content id, and those are 1 to
//! 32, while the smallest printable tag packs to `0x202020` = 2105376. Measured: 1 and 32 are
//! refused, 33 and 2105376 are accepted. Real code that no input can reach.
//!
//! # The type is the value's, not the caller's
//!
//! Not a property of this map, but visible only through it: htsjdk narrows an integer attribute to
//! the smallest type that holds it. Measured, `NM` set to values 1 to 4 is written as **`NMc`** and
//! `NM` set to 100000 as **`NMi`**, from the same Java `Integer`. A port that writes the declared
//! type rather than the narrowed one produces a different key and therefore a different map.

use std::collections::BTreeMap;

use crate::encoding_map::{EncodingDescriptor, EncodingId, DATA_SERIES};
use crate::varint::{read_unsigned_itf8, write_unsigned_itf8, RuntimeEof};

/// `ReadTag.nameType3BytesToInt`: two name bytes and the type, big-endian into twenty-four bits.
pub fn name_type_to_int(name: [u8; 2], tag_type: u8) -> i32 {
    (i32::from(name[0]) << 16) | (i32::from(name[1]) << 8) | i32::from(tag_type)
}

/// `ReadTag.intToNameType3Bytes`.
pub fn int_to_name_type_3(value: i32) -> [u8; 3] {
    [
        ((value >> 16) & 0xFF) as u8,
        ((value >> 8) & 0xFF) as u8,
        (value & 0xFF) as u8,
    ]
}

/// `ReadTag.intToNameType4Bytes`, which inserts a colon before the type.
///
/// htsjdk carries a TODO asking for it to be merged with the three-byte form; they differ only by
/// that colon, and only this one is a display form.
pub fn int_to_name_type_4(value: i32) -> [u8; 4] {
    let three = int_to_name_type_3(value);
    [three[0], three[1], b':', three[2]]
}

/// `putTagBlockCompression`'s check: does this tag id collide with a data series content id?
///
/// It cannot, for any tag whose characters are printable. Kept because it is in the reference and
/// because a port that omits it would differ on a hand-built map.
pub fn overlaps_data_series_content_id(tag_id: i32) -> bool {
    DATA_SERIES
        .iter()
        .any(|(_, _, content_id)| *content_id == tag_id)
}

/// The message that guard raises.
pub fn overlap_message(tag_id: i32) -> String {
    format!("tagID {tag_id} overlaps with data series content ID")
}

/// What a tag encoding map is refused with.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TagEncodingMapError {
    /// The same unchecked, signed array index the data series map has.
    EncodingIdOutOfBounds(i32),
    /// The map ends inside one of its own entries.
    Truncated,
}

impl TagEncodingMapError {
    pub fn message(&self) -> String {
        match self {
            TagEncodingMapError::EncodingIdOutOfBounds(index) => {
                format!("Index {index} out of bounds for length 10")
            }
            TagEncodingMapError::Truncated => {
                "the tag encoding map ends inside one of its entries".to_string()
            }
        }
    }
}

impl From<RuntimeEof> for TagEncodingMapError {
    fn from(_: RuntimeEof) -> Self {
        TagEncodingMapError::Truncated
    }
}

/// The map, ordered by packed tag id, which is the order it is written in.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct TagEncodingMap {
    entries: BTreeMap<i32, EncodingDescriptor>,
}

impl TagEncodingMap {
    pub fn get(&self, tag_id: i32) -> Option<&EncodingDescriptor> {
        self.entries.get(&tag_id)
    }

    pub fn put(&mut self, tag_id: i32, descriptor: EncodingDescriptor) {
        self.entries.insert(tag_id, descriptor);
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// The tag ids it holds, in write order.
    pub fn tag_ids(&self) -> Vec<i32> {
        self.entries.keys().copied().collect()
    }

    /// Read the map given its own bytes, without the length prefix that precedes it.
    pub fn read(map: &[u8]) -> Result<Self, TagEncodingMapError> {
        let mut at = 0usize;
        let itf8 = |at: &mut usize| -> Result<i32, TagEncodingMapError> {
            let (value, consumed) = read_unsigned_itf8(&map[(*at).min(map.len())..])?;
            *at += consumed;
            Ok(value)
        };

        let count = itf8(&mut at)?;
        let mut out = Self::default();
        for _ in 0..count.max(0) {
            let key = itf8(&mut at)?;
            let raw = i32::from(*map.get(at).ok_or(TagEncodingMapError::Truncated)? as i8);
            at += 1;
            let id =
                EncodingId::from_id(raw).ok_or(TagEncodingMapError::EncodingIdOutOfBounds(raw))?;
            let length = itf8(&mut at)?.max(0) as usize;
            let parameters = map
                .get(at..at + length)
                .ok_or(TagEncodingMapError::Truncated)?
                .to_vec();
            at += length;
            out.entries
                .insert(key, EncodingDescriptor { id, parameters });
        }
        Ok(out)
    }

    /// The map's own bytes. Unlike the data series map, nothing is filtered out here: the count is
    /// `tagEncodingMap.size()` and every entry is written.
    pub fn write(&self) -> Vec<u8> {
        let mut out = Vec::new();
        push_itf8(self.entries.len() as i32, &mut out);
        for (tag_id, descriptor) in &self.entries {
            push_itf8(*tag_id, &mut out);
            out.push(descriptor.id as u8);
            push_itf8(descriptor.parameters.len() as i32, &mut out);
            out.extend_from_slice(&descriptor.parameters);
        }
        out
    }

    pub fn write_prefixed(&self) -> Vec<u8> {
        let map = self.write();
        let mut out = Vec::with_capacity(map.len() + 5);
        push_itf8(map.len() as i32, &mut out);
        out.extend_from_slice(&map);
        out
    }

    pub fn read_prefixed(content: &[u8]) -> Result<(Self, usize), TagEncodingMapError> {
        let (size, consumed) = read_unsigned_itf8(content)?;
        let size = size.max(0) as usize;
        let bytes = content
            .get(consumed..consumed + size)
            .ok_or(TagEncodingMapError::Truncated)?;
        Ok((Self::read(bytes)?, consumed + size))
    }
}

fn push_itf8(value: i32, out: &mut Vec<u8>) {
    out.extend_from_slice(&write_unsigned_itf8(value).0);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_key_is_the_tag_packed_into_twenty_four_bits() {
        assert_eq!(name_type_to_int(*b"NM", b'c'), 5131619);
        assert_eq!(name_type_to_int(*b"MD", b'Z'), 5063770);
        assert_eq!(name_type_to_int(*b"  ", b' '), 0x0020_2020);
        assert_eq!(int_to_name_type_3(5131619), *b"NMc");
        assert_eq!(int_to_name_type_4(5131619), *b"NM:c");
    }

    /// The type is part of the key, so the two forms of one name sort by their type character.
    #[test]
    fn one_name_at_two_types_is_two_keys_in_type_order() {
        let z = name_type_to_int(*b"XX", b'Z');
        let i = name_type_to_int(*b"XX", b'i');
        assert_eq!((z, i), (5789786, 5789801));
        assert!(z < i, "Z sorts before i because 0x5A is below 0x69");
    }

    /// A guard over a range no printable tag can enter.
    #[test]
    fn the_collision_guard_cannot_fire_for_a_printable_tag() {
        assert!(overlaps_data_series_content_id(1));
        assert!(overlaps_data_series_content_id(32));
        assert!(!overlaps_data_series_content_id(33));
        // The smallest printable tag, three spaces.
        assert!(!overlaps_data_series_content_id(name_type_to_int(
            *b"  ", b' '
        )));
        assert_eq!(
            overlap_message(1),
            "tagID 1 overlaps with data series content ID"
        );
    }

    /// The write order is the key's, so building the map in either order gives the same bytes.
    #[test]
    fn the_write_order_is_numeric_and_not_the_order_of_insertion() {
        let descriptor = EncodingDescriptor {
            id: EncodingId::ByteArrayLen,
            parameters: vec![1, 2],
        };
        let mut forwards = TagEncodingMap::default();
        let mut backwards = TagEncodingMap::default();
        let keys = [
            name_type_to_int(*b"XX", b'f'),
            name_type_to_int(*b"NM", b'c'),
            name_type_to_int(*b"MD", b'Z'),
        ];
        for key in keys {
            forwards.put(key, descriptor.clone());
        }
        for key in keys.iter().rev() {
            backwards.put(*key, descriptor.clone());
        }
        assert_eq!(forwards.write(), backwards.write());
        assert_eq!(
            forwards.tag_ids(),
            vec![keys[2], keys[1], keys[0]],
            "MDZ, NMc, XXf"
        );
    }

    #[test]
    fn an_empty_map_is_a_single_zero() {
        let map = TagEncodingMap::default();
        assert_eq!(map.write(), vec![0]);
        assert_eq!(TagEncodingMap::read(&[0]), Ok(map));
    }
}
