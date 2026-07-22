//! BAM auxiliary tags: the binary codec and, above all, the integer type promotion.
//!
//! Ported from `htsjdk.samtools.BinaryTagCodec` and `htsjdk.samtools.SAMBinaryTagAndValue`.
//!
//! This module is where a naive encoder produces a valid BAM that is not htsjdk's BAM. Three
//! independent traps live here, and all three yield files that `samtools` reads without
//! complaint:
//!
//! 1. **Integer width is chosen from the value, not the declared type.** `getIntegerType`
//!    picks the narrowest representation, and its ladder is not the obvious one: 300 is
//!    written as a *signed* short, while 200 is written as an *unsigned* byte.
//! 2. **Tag order is by the packed short, not by the tag string.** `makeBinaryTag` packs the
//!    *second* character into the high byte, so tags sort on their second letter first.
//! 3. **Strings are one byte per UTF-16 unit**, truncated, not encoded.

use std::fmt;

/// `BinaryTagCodec.FIXED_TAG_SIZE`: two bytes of name plus one of type.
pub const FIXED_TAG_SIZE: usize = 3;

/// `BinaryTagCodec.FIXED_BINARY_ARRAY_TAG_SIZE`: element type byte plus 4-byte count.
pub const FIXED_BINARY_ARRAY_TAG_SIZE: usize = 5;

const MAX_INT: i64 = i32::MAX as i64;
const MAX_UINT: i64 = MAX_INT * 2 + 1;
const MAX_SHORT: i64 = i16::MAX as i64;
const MAX_USHORT: i64 = MAX_SHORT * 2 + 1;
const MAX_BYTE: i64 = i8::MAX as i64;
const MAX_UBYTE: i64 = MAX_BYTE * 2 + 1;

/// A two-character tag name, held in htsjdk's packed form.
///
/// The packing is `(char[1] << 8) | char[0]`, so the on-disk little-endian bytes come out as
/// the two characters in reading order, but the **numeric** value used for ordering weights
/// the second character more heavily.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Tag(pub i16);

impl Tag {
    /// `SAMTag.makeBinaryTag`.
    pub fn new(name: &[u8; 2]) -> Self {
        Tag(((name[1] as i16) << 8) | name[0] as i16)
    }

    /// The two characters, in reading order, which is also their on-disk order.
    pub fn name(self) -> [u8; 2] {
        [(self.0 & 0xFF) as u8, ((self.0 >> 8) & 0xFF) as u8]
    }
}

impl fmt::Display for Tag {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let n = self.name();
        write!(f, "{}{}", n[0] as char, n[1] as char)
    }
}

/// The in-memory representation of a tag value.
///
/// Integers collapse to a single `Int` variant deliberately: htsjdk reaches every integral
/// Java box type through `((Number) value).longValue()`, so `Byte(100)`, `Short(100)` and
/// `Integer(100)` all encode to the same byte. The declared type has no influence. Arrays are
/// the opposite: their element width is taken from the array's class and never narrowed.
#[derive(Debug, Clone, PartialEq)]
pub enum TagValue {
    /// `'A'`, a single printable character.
    Char(u8),
    /// Any integral value. The on-disk type is derived by [`integer_type`].
    Int(i64),
    /// `'f'`.
    Float(f32),
    /// `'Z'`, stored as UTF-16 units because that is what htsjdk measures and truncates.
    Str(String),
    /// `'H'`. htsjdk **reads** this and never writes it, preferring `B` as more compact.
    Hex(Vec<u8>),
    /// `'B'` with element type `c`/`C`.
    ByteArray { values: Vec<i8>, unsigned: bool },
    /// `'B'` with element type `s`/`S`.
    ShortArray { values: Vec<i16>, unsigned: bool },
    /// `'B'` with element type `i`/`I`.
    IntArray { values: Vec<i32>, unsigned: bool },
    /// `'B'` with element type `f`. There is no unsigned float, and htsjdk ignores the flag.
    FloatArray(Vec<f32>),
}

/// A tag value could not be encoded.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TagError {
    /// `getIntegerType`: "Integer attribute value too large to be encoded in BAM".
    IntegerTooLarge(i64),
    /// `getIntegerType`: "Integer attribute value too negative to be encoded in BAM".
    IntegerTooNegative(i64),
    /// `'H'` is read but never written by htsjdk, so writing one is out of contract.
    HexIsNeverWritten,
}

/// `BinaryTagCodec.getIntegerType`.
///
/// Reproduced as the same ordered ladder rather than as a set of ranges, because the ladder's
/// order *is* the specification. Reading it as ranges:
///
/// | value | type | width |
/// |---|---|---|
/// | `-2^31 .. -32769` | `i` | 4 |
/// | `-32768 .. -129` | `s` | 2 |
/// | `-128 .. 127` | `c` | 1 |
/// | `128 .. 255` | `C` | 1, unsigned |
/// | `256 .. 32767` | `s` | 2 |
/// | `32768 .. 65535` | `S` | 2, unsigned |
/// | `65536 .. 2^31-1` | `i` | 4 |
/// | `2^31 .. 2^32-1` | `I` | 4, unsigned |
///
/// The non-obvious row is `256 .. 32767`, which takes the **signed** short even though the
/// unsigned one would hold it just as well. An encoder that reasoned "smallest type that
/// fits, preferring unsigned" would emit `S` there and be wrong by one byte of type code on
/// every such tag.
pub fn integer_type(val: i64) -> Result<u8, TagError> {
    if val > MAX_UINT {
        return Err(TagError::IntegerTooLarge(val));
    }
    if val > MAX_INT {
        return Ok(b'I');
    }
    if val > MAX_USHORT {
        return Ok(b'i');
    }
    if val > MAX_SHORT {
        return Ok(b'S');
    }
    if val > MAX_UBYTE {
        return Ok(b's');
    }
    if val > MAX_BYTE {
        return Ok(b'C');
    }
    if val >= i8::MIN as i64 {
        return Ok(b'c');
    }
    if val >= i16::MIN as i64 {
        return Ok(b's');
    }
    if val >= i32::MIN as i64 {
        return Ok(b'i');
    }
    Err(TagError::IntegerTooNegative(val))
}

/// `BinaryTagCodec.getTagValueType`.
pub fn tag_value_type(value: &TagValue) -> Result<u8, TagError> {
    Ok(match value {
        TagValue::Str(_) => b'Z',
        TagValue::Char(_) => b'A',
        TagValue::Float(_) => b'f',
        TagValue::Int(v) => integer_type(*v)?,
        TagValue::Hex(_) => return Err(TagError::HexIsNeverWritten),
        TagValue::ByteArray { .. }
        | TagValue::ShortArray { .. }
        | TagValue::IntArray { .. }
        | TagValue::FloatArray(_) => b'B',
    })
}

/// `StringUtil.stringToBytes`: one byte per UTF-16 code unit, truncated to the low 8 bits.
///
/// Not a UTF-8 encoding and not a lossy-ASCII conversion. `é` (U+00E9) becomes the single
/// byte `0xE9`, and a supplementary character becomes its two surrogate halves truncated to
/// two bytes. This is faithful to htsjdk, which is the only property that matters here.
fn string_to_bytes(s: &str) -> Vec<u8> {
    s.encode_utf16().map(|u| (u & 0xFF) as u8).collect()
}

/// `BinaryTagCodec.getBinaryValueSize`.
pub fn binary_value_size(value: &TagValue) -> Result<usize, TagError> {
    Ok(match value {
        // `String.length()` is UTF-16 units, matching `string_to_bytes` byte for byte.
        TagValue::Str(s) => s.encode_utf16().count() + 1,
        TagValue::Char(_) => 1,
        TagValue::Float(_) => 4,
        TagValue::Int(v) => match integer_type(*v)? {
            b'I' | b'i' => 4,
            b's' | b'S' => 2,
            b'c' | b'C' => 1,
            t => unreachable!("integer_type returned {}", t as char),
        },
        TagValue::Hex(b) => b.len() * 2 + 1,
        TagValue::ByteArray { values, .. } => values.len() + FIXED_BINARY_ARRAY_TAG_SIZE,
        TagValue::ShortArray { values, .. } => values.len() * 2 + FIXED_BINARY_ARRAY_TAG_SIZE,
        TagValue::IntArray { values, .. } => values.len() * 4 + FIXED_BINARY_ARRAY_TAG_SIZE,
        TagValue::FloatArray(values) => values.len() * 4 + FIXED_BINARY_ARRAY_TAG_SIZE,
    })
}

/// `BinaryTagCodec.getTagSize`.
pub fn tag_size(value: &TagValue) -> Result<usize, TagError> {
    Ok(FIXED_TAG_SIZE + binary_value_size(value)?)
}

/// `BinaryTagCodec.writeTag`.
pub fn write_tag(out: &mut Vec<u8>, tag: Tag, value: &TagValue) -> Result<(), TagError> {
    out.extend_from_slice(&tag.0.to_le_bytes());
    let ty = tag_value_type(value)?;
    out.push(ty);

    match value {
        TagValue::Str(s) => {
            out.extend_from_slice(&string_to_bytes(s));
            out.push(0);
        }
        TagValue::Char(c) => out.push(*c),
        TagValue::Float(f) => out.extend_from_slice(&f.to_le_bytes()),
        TagValue::Int(v) => match ty {
            // `writeUInt`: the low 32 bits, which for a value above i32::MAX is the
            // two's-complement pattern a signed write would also have produced.
            b'I' | b'i' => out.extend_from_slice(&(*v as i32).to_le_bytes()),
            b's' | b'S' => out.extend_from_slice(&(*v as i16).to_le_bytes()),
            b'c' | b'C' => out.push(*v as u8),
            _ => unreachable!(),
        },
        TagValue::Hex(_) => return Err(TagError::HexIsNeverWritten),
        TagValue::ByteArray { values, unsigned } => {
            out.push(if *unsigned { b'C' } else { b'c' });
            out.extend_from_slice(&(values.len() as i32).to_le_bytes());
            for v in values {
                out.push(*v as u8);
            }
        }
        TagValue::ShortArray { values, unsigned } => {
            out.push(if *unsigned { b'S' } else { b's' });
            out.extend_from_slice(&(values.len() as i32).to_le_bytes());
            for v in values {
                out.extend_from_slice(&v.to_le_bytes());
            }
        }
        TagValue::IntArray { values, unsigned } => {
            out.push(if *unsigned { b'I' } else { b'i' });
            out.extend_from_slice(&(values.len() as i32).to_le_bytes());
            for v in values {
                out.extend_from_slice(&v.to_le_bytes());
            }
        }
        TagValue::FloatArray(values) => {
            // No unsigned float: htsjdk writes 'f' unconditionally, ignoring the flag.
            out.push(b'f');
            out.extend_from_slice(&(values.len() as i32).to_le_bytes());
            for v in values {
                out.extend_from_slice(&v.to_le_bytes());
            }
        }
    }
    Ok(())
}

/// A tag list in the order htsjdk keeps it.
///
/// `SAMBinaryTagAndValue.insert` maintains a linked list sorted ascending by the packed short,
/// with an equal tag *replacing* rather than duplicating. `BAMRecordCodec.encode` then walks
/// that list, so the sort order is the write order.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Tags {
    entries: Vec<(Tag, TagValue)>,
}

impl Tags {
    pub fn new() -> Self {
        Self::default()
    }

    /// `SAMBinaryTagAndValue.insert`: sorted insert, replacing on an equal tag.
    pub fn insert(&mut self, tag: Tag, value: TagValue) {
        match self.entries.binary_search_by_key(&tag, |(t, _)| *t) {
            Ok(i) => self.entries[i] = (tag, value),
            Err(i) => self.entries.insert(i, (tag, value)),
        }
    }

    pub fn get(&self, tag: Tag) -> Option<&TagValue> {
        self.entries
            .binary_search_by_key(&tag, |(t, _)| *t)
            .ok()
            .map(|i| &self.entries[i].1)
    }

    /// `SAMRecord.setAttribute(tag, null)`: drop the tag if present.
    pub fn remove(&mut self, tag: Tag) {
        if let Ok(i) = self.entries.binary_search_by_key(&tag, |(t, _)| *t) {
            self.entries.remove(i);
        }
    }

    pub fn iter(&self) -> impl Iterator<Item = &(Tag, TagValue)> {
        self.entries.iter()
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Total on-disk size, as `BAMRecordCodec.encode` accumulates it.
    pub fn binary_size(&self) -> Result<usize, TagError> {
        self.entries.iter().map(|(_, v)| tag_size(v)).sum()
    }

    pub fn write(&self, out: &mut Vec<u8>) -> Result<(), TagError> {
        for (tag, value) in &self.entries {
            write_tag(out, *tag, value)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn t(s: &str) -> Tag {
        Tag::new(s.as_bytes().try_into().unwrap())
    }

    #[test]
    fn a_tag_round_trips_through_its_packed_form() {
        for name in ["NM", "MD", "RG", "AS", "ZZ", "aa"] {
            assert_eq!(t(name).name(), name.as_bytes());
            assert_eq!(t(name).to_string(), name);
        }
    }

    /// The packed short puts the *second* character in the high byte, so ordering is by the
    /// second letter first. Sorting by the tag string instead would reorder the tag block of
    /// every record that carries more than one tag.
    #[test]
    fn tags_order_by_their_second_character_first() {
        assert!(
            t("ZA") < t("AZ"),
            "ZA packs to {} and AZ to {}: the second character dominates",
            t("ZA").0,
            t("AZ").0
        );
        assert!(
            "AZ" < "ZA",
            "and the naive string order is the opposite, which is the whole trap"
        );
    }

    #[test]
    fn the_packed_bytes_are_the_characters_in_reading_order() {
        let mut out = Vec::new();
        write_tag(&mut out, t("NM"), &TagValue::Int(0)).unwrap();
        assert_eq!(&out[..2], b"NM", "on disk a tag reads left to right");
    }

    /// The full promotion ladder, boundary by boundary. Every row here is a value at which a
    /// plausible alternative rule would pick a different type.
    #[test]
    fn integer_promotion_follows_the_exact_ladder() {
        let cases: &[(i64, u8)] = &[
            (i32::MIN as i64, b'i'),
            (-32_769, b'i'),
            (-32_768, b's'),
            (-129, b's'),
            (-128, b'c'),
            (0, b'c'),
            (127, b'c'),
            (128, b'C'),
            (255, b'C'),
            (256, b's'),
            (32_767, b's'),
            (32_768, b'S'),
            (65_535, b'S'),
            (65_536, b'i'),
            (i32::MAX as i64, b'i'),
            (i32::MAX as i64 + 1, b'I'),
            (4_294_967_295, b'I'),
        ];
        for &(v, expect) in cases {
            assert_eq!(
                integer_type(v).unwrap(),
                expect,
                "value {v} must be type '{}'",
                expect as char
            );
        }
    }

    /// The single most counter-intuitive row, stated on its own so it cannot be lost in a
    /// table: 300 fits an unsigned short, and htsjdk still writes a signed one.
    #[test]
    fn a_value_of_three_hundred_takes_the_signed_short() {
        assert_eq!(integer_type(300).unwrap(), b's');
        assert_eq!(integer_type(200).unwrap(), b'C');
    }

    #[test]
    fn out_of_range_integers_are_refused() {
        assert_eq!(
            integer_type(4_294_967_296),
            Err(TagError::IntegerTooLarge(4_294_967_296))
        );
        let too_low = i32::MIN as i64 - 1;
        assert_eq!(
            integer_type(too_low),
            Err(TagError::IntegerTooNegative(too_low))
        );
    }

    #[test]
    fn declared_width_does_not_influence_the_encoding() {
        // Whatever Java box type held it, only the value matters.
        let mut a = Vec::new();
        let mut b = Vec::new();
        write_tag(&mut a, t("NM"), &TagValue::Int(100)).unwrap();
        write_tag(&mut b, t("NM"), &TagValue::Int(100i64)).unwrap();
        assert_eq!(a, b);
        assert_eq!(a, b"NM\x63\x64");
    }

    /// Arrays are the opposite rule: element width comes from the array type and is never
    /// narrowed, so a `[1, 2, 3]` int array stays 4 bytes per element.
    #[test]
    fn array_element_width_is_not_narrowed() {
        let mut out = Vec::new();
        write_tag(
            &mut out,
            t("ZI"),
            &TagValue::IntArray {
                values: vec![1, 2, 3],
                unsigned: false,
            },
        )
        .unwrap();
        assert_eq!(out.len(), FIXED_TAG_SIZE + FIXED_BINARY_ARRAY_TAG_SIZE + 12);
        assert_eq!(out[3], b'i');
    }

    #[test]
    fn unsigned_arrays_differ_only_in_the_type_letter() {
        let mk = |unsigned| {
            let mut o = Vec::new();
            write_tag(
                &mut o,
                t("ZB"),
                &TagValue::ByteArray {
                    values: vec![1, 2],
                    unsigned,
                },
            )
            .unwrap();
            o
        };
        let (s, u) = (mk(false), mk(true));
        assert_eq!(s[3], b'c');
        assert_eq!(u[3], b'C');
        assert_eq!(s[..3], u[..3]);
        assert_eq!(s[4..], u[4..]);
    }

    #[test]
    fn sizes_agree_with_what_is_actually_written() {
        let values = [
            TagValue::Char(b'X'),
            TagValue::Int(-5),
            TagValue::Int(300),
            TagValue::Int(70_000),
            TagValue::Int(3_000_000_000),
            TagValue::Float(1.5),
            TagValue::Str("hello".into()),
            TagValue::ByteArray {
                values: vec![1, 2, 3],
                unsigned: false,
            },
            TagValue::ShortArray {
                values: vec![1, 2],
                unsigned: true,
            },
            TagValue::IntArray {
                values: vec![7],
                unsigned: false,
            },
            TagValue::FloatArray(vec![1.0, 2.0]),
        ];
        for v in &values {
            let mut out = Vec::new();
            write_tag(&mut out, t("XX"), v).unwrap();
            assert_eq!(
                out.len(),
                tag_size(v).unwrap(),
                "declared size disagrees with written bytes for {v:?}"
            );
        }
    }

    /// The size accounting drives the record's `block_size` field. If it disagreed with the
    /// bytes actually written, every downstream record would be misparsed, so this is checked
    /// against reality rather than against a second copy of the same arithmetic.
    #[test]
    fn a_z_string_is_one_byte_per_utf16_unit() {
        let v = TagValue::Str("café".into());
        let mut out = Vec::new();
        write_tag(&mut out, t("XX"), &v).unwrap();
        // 4 UTF-16 units + NUL, not the 5 UTF-8 bytes + NUL.
        assert_eq!(binary_value_size(&v).unwrap(), 5);
        assert_eq!(out.len(), tag_size(&v).unwrap());
        assert_eq!(&out[3..], b"caf\xE9\x00");
    }

    #[test]
    fn hex_is_refused_because_htsjdk_never_writes_it() {
        let mut out = Vec::new();
        assert_eq!(
            write_tag(&mut out, t("XX"), &TagValue::Hex(vec![1, 2])),
            Err(TagError::HexIsNeverWritten)
        );
    }

    #[test]
    fn tags_are_kept_in_packed_short_order() {
        let mut tags = Tags::new();
        for name in ["NM", "AS", "ZA", "AZ", "MD"] {
            tags.insert(t(name), TagValue::Int(1));
        }
        let order: Vec<String> = tags.iter().map(|(k, _)| k.to_string()).collect();
        let mut expected = ["NM", "AS", "ZA", "AZ", "MD"];
        expected.sort_by_key(|n| t(n).0);
        assert_eq!(order, expected);
        // Concretely: ZA before AZ.
        let (za, az) = (
            order.iter().position(|s| s == "ZA").unwrap(),
            order.iter().position(|s| s == "AZ").unwrap(),
        );
        assert!(za < az);
    }

    #[test]
    fn inserting_the_same_tag_twice_replaces_it() {
        let mut tags = Tags::new();
        tags.insert(t("NM"), TagValue::Int(1));
        tags.insert(t("NM"), TagValue::Int(2));
        assert_eq!(tags.len(), 1);
        assert_eq!(tags.get(t("NM")), Some(&TagValue::Int(2)));
    }

    #[test]
    fn declared_binary_size_matches_the_written_block() {
        let mut tags = Tags::new();
        tags.insert(t("NM"), TagValue::Int(3));
        tags.insert(t("MD"), TagValue::Str("100".into()));
        tags.insert(
            t("BQ"),
            TagValue::ByteArray {
                values: vec![1, 2, 3],
                unsigned: true,
            },
        );
        let mut out = Vec::new();
        tags.write(&mut out).unwrap();
        assert_eq!(out.len(), tags.binary_size().unwrap());
    }
}
