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
//! 4. **There is no in-memory `H`.** It decodes into the same `byte[]` a `B` array does, so
//!    every branch that would write one back is dead, and an `H` tag read and rewritten comes
//!    out as a `B` array in both the binary and the text codec.

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
    /// `'B'` with element type `c`/`C`, and also what an `'H'` becomes.
    ///
    /// There is no in-memory `H`: htsjdk decodes one into a plain `byte[]`, which is the same
    /// thing a `B` array of signed bytes decodes into, and `getTagValueType` answers `'B'` for it.
    /// Both codecs carry a dead branch for writing `H` back, and `TextTagCodec`'s says so in a
    /// comment. Measured: `XX:H:48656C` re-encodes as `XX:B:c,72,101,108`.
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

    /// `BinaryTagCodec.readTags`: a whole tag block, in the order htsjdk gives it back.
    ///
    /// The block carries no count and no length of its own: it is read until the bytes run out,
    /// which is why the CRAM slice header's tag section can only be delimited by its block.
    pub fn read(bytes: &[u8]) -> Result<Self, TagReadError> {
        let mut cursor = Cursor { bytes, at: 0usize };
        let mut tags = Self::new();
        while cursor.at < cursor.bytes.len() {
            let tag = Tag(cursor.i16()?);
            let ty = cursor.u8()?;
            let value = if ty == b'B' {
                read_array(&mut cursor)?
            } else {
                read_single_value(ty, &mut cursor)?
            };
            // `insert` is htsjdk's own: sorted by the packed short, and an equal tag *replaces*.
            tags.insert(tag, value);
        }
        Ok(tags)
    }
}

/// What a tag block is refused with.
///
/// The three exception types are htsjdk's, and two of them are the JDK's rather than htsjdk's own:
/// a truncated block is a `BufferUnderflowException` with no message at all, which is what a
/// caller trying to tell "malformed" from "empty" apart has to work with.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TagReadError {
    /// `SAMFormatException: Unrecognized tag type: <c>`.
    UnrecognizedType(u8),
    /// `SAMFormatException: Unrecognized tag array type: <c>`.
    UnrecognizedArrayType(u8),
    /// `BufferUnderflowException`, which carries no message.
    BufferUnderflow,
    /// `NegativeArraySizeException`, whose message is the length itself.
    NegativeArraySize(i32),
    /// `NumberFormatException` out of `StringUtil.hexStringToBytes`, for an `H` value.
    OddHexLength(String),
    /// `NumberFormatException` out of `StringUtil.fromHexDigit`.
    NotAHexDigit(char),
}

impl TagReadError {
    /// The message htsjdk raises, or the empty string where the JDK raises none.
    pub fn message(&self) -> String {
        match self {
            // Both messages concatenate `(char) theByte`, so a type byte above 0x7F is reported
            // as the sign-extended character rather than as itself.
            TagReadError::UnrecognizedType(c) => {
                format!("Unrecognized tag type: {}", as_java_char(*c))
            }
            TagReadError::UnrecognizedArrayType(c) => {
                format!("Unrecognized tag array type: {}", as_java_char(*c))
            }
            TagReadError::BufferUnderflow => String::new(),
            TagReadError::NegativeArraySize(n) => n.to_string(),
            TagReadError::OddHexLength(s) => format!(
                "Hex representation of byte string does not have even number of hex chars: {s}"
            ),
            TagReadError::NotAHexDigit(c) => format!("Not a valid hex digit: {c}"),
        }
    }
}

impl TagReadError {
    /// The Java exception it arrives as, which is not always htsjdk's own.
    pub fn java_exception(&self) -> &'static str {
        match self {
            TagReadError::UnrecognizedType(_) | TagReadError::UnrecognizedArrayType(_) => {
                "SAMFormatException"
            }
            TagReadError::BufferUnderflow => "BufferUnderflowException",
            TagReadError::NegativeArraySize(_) => "NegativeArraySizeException",
            TagReadError::OddHexLength(_) | TagReadError::NotAHexDigit(_) => {
                "NumberFormatException"
            }
        }
    }
}

impl fmt::Display for TagReadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let message = self.message();
        if message.is_empty() {
            write!(f, "{}", self.java_exception())
        } else {
            write!(f, "{}: {message}", self.java_exception())
        }
    }
}

/// `(char) aSignedByte`, which sign-extends before it truncates.
///
/// The same cast the CRAM substitution matrix's messages go through. It matters here for `'A'`:
/// `readSingleValue` returns `(char) byteBuffer.get()`, so the byte `0xE9` becomes `U+FFE9` and
/// not `U+00E9`. The value written back is the low byte again, so the *bytes* survive a round trip
/// while the in-memory character never was the one in the file.
pub fn java_char(byte: u8) -> u16 {
    byte as i8 as i32 as u16
}

/// That character as Rust sees it, for the two messages that concatenate one.
fn as_java_char(byte: u8) -> char {
    char::from_u32(u32::from(java_char(byte))).expect("a UTF-16 unit below the surrogates")
}

/// A checked little-endian cursor. Every overrun is the one `BufferUnderflowException` htsjdk lets
/// through, so there is nothing to distinguish between the ways a block can be short.
struct Cursor<'a> {
    bytes: &'a [u8],
    at: usize,
}

impl Cursor<'_> {
    fn take(&mut self, n: usize) -> Result<&[u8], TagReadError> {
        let slice = self
            .bytes
            .get(self.at..self.at + n)
            .ok_or(TagReadError::BufferUnderflow)?;
        self.at += n;
        Ok(slice)
    }

    fn u8(&mut self) -> Result<u8, TagReadError> {
        Ok(self.take(1)?[0])
    }

    fn i16(&mut self) -> Result<i16, TagReadError> {
        Ok(i16::from_le_bytes(self.take(2)?.try_into().expect("two")))
    }

    fn i32(&mut self) -> Result<i32, TagReadError> {
        Ok(i32::from_le_bytes(self.take(4)?.try_into().expect("four")))
    }

    fn f32(&mut self) -> Result<f32, TagReadError> {
        Ok(f32::from_le_bytes(self.take(4)?.try_into().expect("four")))
    }

    /// `readNullTerminatedString`, whose scan for the terminator is what runs off the end.
    fn null_terminated(&mut self) -> Result<String, TagReadError> {
        let end = self.bytes[self.at..]
            .iter()
            .position(|b| *b == 0)
            .ok_or(TagReadError::BufferUnderflow)?;
        // `StringUtil.bytesToString` is the inverse of `stringToBytes`: one char per byte, so a
        // byte of 0xE9 is U+00E9 and not the two bytes UTF-8 would have wanted.
        let text = self.bytes[self.at..self.at + end]
            .iter()
            .map(|b| *b as char)
            .collect();
        self.at += end + 1;
        Ok(text)
    }
}

/// `StringUtil.hexStringToBytes`, which `H` goes through.
///
/// `Character.digit(c, 16)` also accepts non-ASCII digit forms, which cannot occur here: the text
/// came from bytes mapped one to one into `U+0000..U+00FF`, and Latin-1 holds no such form.
fn hex_string_to_bytes(text: &str) -> Result<Vec<u8>, TagReadError> {
    // `String.length()` counts UTF-16 units, and a byte of 0xE9 became one character that Rust
    // stores in two bytes, so counting `text.len()` here would refuse a string Java accepts.
    let digits: Vec<char> = text.chars().collect();
    if !digits.len().is_multiple_of(2) {
        return Err(TagReadError::OddHexLength(text.to_string()));
    }
    let mut out = Vec::with_capacity(digits.len() / 2);
    for pair in digits.chunks(2) {
        let hi = from_hex_digit(pair[0])?;
        let lo = from_hex_digit(pair[1])?;
        out.push(((hi << 4) | lo) as u8);
    }
    Ok(out)
}

fn from_hex_digit(c: char) -> Result<u32, TagReadError> {
    c.to_digit(16).ok_or(TagReadError::NotAHexDigit(c))
}

/// `BinaryTagCodec.readSingleValue`.
///
/// Every narrow integer widens to one type here. `c`, `C`, `s`, `S` and `i` all come back as a
/// Java `Integer`, so the width the file chose stops existing the moment it is read, and the type
/// a rewrite picks is derived from the *value* alone. Only `I` above `Integer.MAX_VALUE` comes
/// back as a `Long`, and [`TagValue::Int`] holds both because the rewrite cannot tell them apart.
fn read_single_value(ty: u8, cursor: &mut Cursor<'_>) -> Result<TagValue, TagReadError> {
    Ok(match ty {
        b'Z' => TagValue::Str(cursor.null_terminated()?),
        b'A' => TagValue::Char(cursor.u8()?),
        // `getInt() & 0xffffffffL` is already inside [0, 2^32-1], which is exactly the range
        // `isValidUnsignedIntegerAttribute` accepts, so the validation branch under this line
        // cannot fire for any input at all.
        b'I' => TagValue::Int(i64::from(cursor.i32()?) & 0xFFFF_FFFF),
        b'i' => TagValue::Int(i64::from(cursor.i32()?)),
        b's' => TagValue::Int(i64::from(cursor.i16()?)),
        b'S' => TagValue::Int(i64::from(cursor.i16()?) & 0xFFFF),
        b'c' => TagValue::Int(i64::from(cursor.u8()? as i8)),
        b'C' => TagValue::Int(i64::from(cursor.u8()?)),
        b'f' => TagValue::Float(cursor.f32()?),
        // A byte[], indistinguishable from what a signed B array gives back, which is why
        // rewriting an H tag produces a B one.
        b'H' => TagValue::ByteArray {
            values: hex_string_to_bytes(&cursor.null_terminated()?)?
                .into_iter()
                .map(|b| b as i8)
                .collect(),
            unsigned: false,
        },
        other => return Err(TagReadError::UnrecognizedType(other)),
    })
}

/// `BinaryTagCodec.readArray`.
///
/// The unsigned flag is the *case* of the element type letter, and it changes nothing but the
/// letter written back: the elements themselves stay signed. A `C` array holding `0xFF` comes back
/// as `-1`, and htsjdk's own `B-S` array of `65535` comes back as `-1` too.
///
/// The length is read before the element type is judged, so an unrecognized type is refused even
/// when the length that follows it is impossible.
///
/// The elements are grown rather than reserved: htsjdk allocates the whole array up front and
/// meets a hostile length with an `OutOfMemoryError`, which is not a behaviour worth reproducing
/// by aborting the process.
fn read_array(cursor: &mut Cursor<'_>) -> Result<TagValue, TagReadError> {
    let array_type = cursor.u8()?;
    let unsigned = array_type.is_ascii_uppercase();
    let length = cursor.i32()?;
    let lowered = array_type.to_ascii_lowercase();
    if !matches!(lowered, b'c' | b's' | b'i' | b'f') {
        return Err(TagReadError::UnrecognizedArrayType(array_type));
    }
    if length < 0 {
        return Err(TagReadError::NegativeArraySize(length));
    }
    let length = length as usize;
    Ok(match lowered {
        b'c' => TagValue::ByteArray {
            values: cursor.take(length)?.iter().map(|b| *b as i8).collect(),
            unsigned,
        },
        b's' => {
            let mut values = Vec::new();
            for _ in 0..length {
                values.push(cursor.i16()?);
            }
            TagValue::ShortArray { values, unsigned }
        }
        b'i' => {
            let mut values = Vec::new();
            for _ in 0..length {
                values.push(cursor.i32()?);
            }
            TagValue::IntArray { values, unsigned }
        }
        _ => {
            let mut values = Vec::new();
            for _ in 0..length {
                values.push(cursor.f32()?);
            }
            // There is no unsigned float: the flag is read and then dropped.
            TagValue::FloatArray(values)
        }
    })
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

    /// An `H` value has no type of its own, so it is written as the `B` array it decoded into.
    #[test]
    fn what_an_h_decoded_into_is_written_back_as_an_array() {
        let value = TagValue::ByteArray {
            values: vec![1, 2],
            unsigned: false,
        };
        assert_eq!(tag_value_type(&value).unwrap(), b'B');
        let mut out = Vec::new();
        write_tag(&mut out, t("XX"), &value).unwrap();
        assert_eq!(out, b"XXBc\x02\x00\x00\x00\x01\x02");
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

    /// A block with one tag in it, built the way a file has it.
    fn block(name: &str, ty: u8, value: &[u8]) -> Vec<u8> {
        let mut out = name.as_bytes().to_vec();
        out.push(ty);
        out.extend_from_slice(value);
        out
    }

    /// The width the file chose stops existing on the way in: five type letters, one value type.
    #[test]
    fn every_narrow_integer_widens_to_one_type() {
        let cases: &[(u8, &[u8], i64)] = &[
            (b'c', &[0xFB], -5),
            (b'C', &[0xC8], 200),
            (b's', &[0x2C, 0x01], 300),
            (b'S', &[0x40, 0x9C], 40_000),
            (b'i', &[0x70, 0x11, 0x01, 0x00], 70_000),
        ];
        for &(ty, bytes, expect) in cases {
            let tags = Tags::read(&block("NM", ty, bytes)).expect("reads");
            assert_eq!(
                tags.get(t("NM")),
                Some(&TagValue::Int(expect)),
                "type '{}'",
                ty as char
            );
        }
    }

    /// Reading and writing back is the identity for everything htsjdk itself writes, and for
    /// exactly two types that it never does.
    #[test]
    fn the_round_trip_breaks_only_where_htsjdk_never_writes() {
        // 'I' holding a small value: written back as 'c', which is one byte shorter.
        let small_i = block("IA", b'I', &[5, 0, 0, 0]);
        let mut out = Vec::new();
        Tags::read(&small_i).unwrap().write(&mut out).unwrap();
        assert_eq!(out, block("IA", b'c', &[5]));

        // 'H': read as bytes, written back as a 'B' array of signed bytes.
        let hex = block("Hh", b'H', b"48656C\0");
        let mut out = Vec::new();
        Tags::read(&hex).unwrap().write(&mut out).unwrap();
        assert_eq!(
            out,
            block("Hh", b'B', &[b'c', 3, 0, 0, 0, 0x48, 0x65, 0x6C]),
            "an H tag comes back as a B array, and grows by a byte doing it"
        );

        // And an 'I' that really needs its width keeps it.
        let big_i = block("IF", b'I', &[0xFF, 0xFF, 0xFF, 0xFF]);
        let mut out = Vec::new();
        Tags::read(&big_i).unwrap().write(&mut out).unwrap();
        assert_eq!(out, big_i);
    }

    /// `(char) aSignedByte` sign-extends, so the character is not the byte. The bytes still
    /// survive the round trip, because the write truncates it back.
    #[test]
    fn a_high_char_tag_is_sign_extended_in_memory_only() {
        assert_eq!(java_char(0xE9), 0xFFE9);
        assert_eq!(java_char(b'Q'), 0x0051);
        let bytes = block("CB", b'A', &[0xE9]);
        let tags = Tags::read(&bytes).expect("reads");
        assert_eq!(tags.get(t("CB")), Some(&TagValue::Char(0xE9)));
        let mut out = Vec::new();
        tags.write(&mut out).unwrap();
        assert_eq!(out, bytes);
    }

    /// The unsigned flag is the case of the type letter and nothing else: the elements stay
    /// signed, so a `C` array holding 0xFF comes back as -1.
    #[test]
    fn an_unsigned_array_keeps_signed_elements() {
        let tags =
            Tags::read(&block("BC", b'B', &[b'C', 3, 0, 0, 0, 0xFF, 0x00, 0x7F])).expect("reads");
        assert_eq!(
            tags.get(t("BC")),
            Some(&TagValue::ByteArray {
                values: vec![-1, 0, 127],
                unsigned: true,
            })
        );
    }

    /// Sorted by the packed short, so by the second character first, whatever order the file had.
    #[test]
    fn a_block_comes_back_sorted_by_the_packed_short() {
        let mut bytes = Vec::new();
        for (name, value) in [("ZA", 1u8), ("AZ", 2), ("NM", 3), ("MD", 4)] {
            bytes.extend_from_slice(&block(name, b'c', &[value]));
        }
        let tags = Tags::read(&bytes).expect("reads");
        let order: Vec<String> = tags.iter().map(|(tag, _)| tag.to_string()).collect();
        assert_eq!(order, ["ZA", "MD", "NM", "AZ"]);
    }

    /// A repeated tag replaces rather than duplicates, and the last one in the file wins.
    #[test]
    fn a_repeated_tag_leaves_one_entry_and_the_last_value() {
        let mut bytes = block("NM", b'c', &[1]);
        bytes.extend_from_slice(&block("NM", b'c', &[2]));
        bytes.extend_from_slice(&block("NM", b'c', &[3]));
        let tags = Tags::read(&bytes).expect("reads");
        assert_eq!(tags.len(), 1);
        assert_eq!(tags.get(t("NM")), Some(&TagValue::Int(3)));
    }

    #[test]
    fn an_empty_block_is_an_empty_list_and_not_an_error() {
        assert_eq!(Tags::read(&[]), Ok(Tags::new()));
    }

    /// Every way a block can be malformed, with the message the reference raises.
    #[test]
    fn the_read_errors_are_the_reference_errors() {
        let cases: &[(Vec<u8>, TagReadError, &str)] = &[
            (
                block("XX", b'q', &[1]),
                TagReadError::UnrecognizedType(b'q'),
                "Unrecognized tag type: q",
            ),
            (
                block("XX", b'B', &[b'q', 1, 0, 0, 0, 1]),
                TagReadError::UnrecognizedArrayType(b'q'),
                "Unrecognized tag array type: q",
            ),
            (block("XX", b'Z', b"ab"), TagReadError::BufferUnderflow, ""),
            (b"XX".to_vec(), TagReadError::BufferUnderflow, ""),
            (
                block("XX", b'B', &[b'c', 0xFF, 0xFF, 0xFF, 0xFF]),
                TagReadError::NegativeArraySize(-1),
                "-1",
            ),
            (
                block("XX", b'H', b"486\0"),
                TagReadError::OddHexLength("486".into()),
                "Hex representation of byte string does not have even number of hex chars: 486",
            ),
            (
                block("XX", b'H', b"4G\0"),
                TagReadError::NotAHexDigit('G'),
                "Not a valid hex digit: G",
            ),
        ];
        for (bytes, expect, message) in cases {
            let error = Tags::read(bytes).expect_err("refused");
            assert_eq!(&error, expect);
            assert_eq!(error.message(), *message);
        }
    }
}
