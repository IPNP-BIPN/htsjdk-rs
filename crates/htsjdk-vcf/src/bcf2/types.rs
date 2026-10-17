//! `BCF2Type` and the integer-sizing helpers of `BCF2Utils`.

use super::Bcf2Error;

/// `BCF2Type`, in declaration order: `determineIntegerType(int[])` compares by ordinal.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Bcf2Type {
    Missing,
    Int8,
    Int16,
    Int32,
    Float,
    Char,
}

/// `BCF2Utils.MAX_ALLELES_IN_GENOTYPES`.
pub const MAX_ALLELES_IN_GENOTYPES: usize = 127;
/// `BCF2Utils.OVERFLOW_ELEMENT_MARKER`.
pub const OVERFLOW_ELEMENT_MARKER: i32 = 15;
/// `BCF2Utils.MAX_INLINE_ELEMENTS`.
pub const MAX_INLINE_ELEMENTS: i32 = 14;

impl Bcf2Type {
    /// The type ID in the low nibble of a descriptor.
    pub fn id(self) -> i32 {
        match self {
            Bcf2Type::Missing => 0,
            Bcf2Type::Int8 => 1,
            Bcf2Type::Int16 => 2,
            Bcf2Type::Int32 => 3,
            Bcf2Type::Float => 5,
            Bcf2Type::Char => 7,
        }
    }

    /// `BCF2Utils.ID_TO_ENUM`: an array of eight with holes at 4 and 6, which are `null` and
    /// only fail when the type is used.
    pub fn from_id(id: i32) -> Option<Bcf2Type> {
        match id {
            0 => Some(Bcf2Type::Missing),
            1 => Some(Bcf2Type::Int8),
            2 => Some(Bcf2Type::Int16),
            3 => Some(Bcf2Type::Int32),
            5 => Some(Bcf2Type::Float),
            7 => Some(Bcf2Type::Char),
            _ => None,
        }
    }

    pub fn size_in_bytes(self) -> usize {
        match self {
            Bcf2Type::Missing => 0,
            Bcf2Type::Int8 | Bcf2Type::Char => 1,
            Bcf2Type::Int16 => 2,
            Bcf2Type::Int32 | Bcf2Type::Float => 4,
        }
    }

    /// `getMissingBytes()`, as the `int` the readers compare against: the INT8 and INT16 values
    /// are sign-extended, which is what makes `read() == getMissingBytes()` work for them.
    pub fn missing_bytes(self) -> i32 {
        match self {
            Bcf2Type::Missing | Bcf2Type::Char => 0,
            Bcf2Type::Int8 => 0xFFFF_FF80u32 as i32,
            Bcf2Type::Int16 => 0xFFFF_8000u32 as i32,
            Bcf2Type::Int32 => i32::MIN,
            Bcf2Type::Float => 0x7F80_0001,
        }
    }

    /// `withinRange`: the minimum is one above the missing value, so `-128` is not an INT8.
    /// The non-integer types have a range of `[0, 0]`.
    pub fn within_range(self, v: i64) -> bool {
        let (min, max) = match self {
            Bcf2Type::Int8 => (-127, 127),
            Bcf2Type::Int16 => (-32767, 32767),
            Bcf2Type::Int32 => (-2_147_483_647, 2_147_483_647),
            _ => (0, 0),
        };
        v >= min && v <= max
    }

    pub fn is_integer_type(self) -> bool {
        matches!(self, Bcf2Type::Int8 | Bcf2Type::Int16 | Bcf2Type::Int32)
    }

    /// `write(value, out)`: little-endian, the value's low bytes only, so an INT8 written with a
    /// value that does not fit is truncated rather than refused.
    pub fn write(self, value: i32, out: &mut Vec<u8>) -> Result<(), Bcf2Error> {
        match self {
            Bcf2Type::Missing => Err(Bcf2Error::illegal_argument("Cannot write MISSING type")),
            Bcf2Type::Int8 | Bcf2Type::Char => {
                out.push(value as u8);
                Ok(())
            }
            Bcf2Type::Int16 => {
                out.extend_from_slice(&(value as u16).to_le_bytes());
                Ok(())
            }
            Bcf2Type::Int32 | Bcf2Type::Float => {
                out.extend_from_slice(&value.to_le_bytes());
                Ok(())
            }
        }
    }

    /// `read(in)`, sign-extending INT8 and INT16. A read past the end of a
    /// `ByteArrayInputStream` returns -1, which `readByte` turns into the byte `0xFF`.
    pub fn read(self, cursor: &mut Cursor) -> Result<i32, Bcf2Error> {
        match self {
            Bcf2Type::Missing => Err(Bcf2Error::illegal_argument("Cannot read MISSING type")),
            Bcf2Type::Int8 | Bcf2Type::Char => Ok(cursor.byte() as i8 as i32),
            Bcf2Type::Int16 => {
                let b2 = cursor.byte() as i32;
                let b1 = cursor.byte() as i32;
                Ok(((b1 << 8) | b2) as i16 as i32)
            }
            Bcf2Type::Int32 | Bcf2Type::Float => {
                let b4 = cursor.byte() as u32;
                let b3 = cursor.byte() as u32;
                let b2 = cursor.byte() as u32;
                let b1 = cursor.byte() as u32;
                Ok((b1 << 24 | b2 << 16 | b3 << 8 | b4) as i32)
            }
        }
    }
}

impl std::fmt::Display for Bcf2Type {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Bcf2Type::Missing => "MISSING",
            Bcf2Type::Int8 => "INT8",
            Bcf2Type::Int16 => "INT16",
            Bcf2Type::Int32 => "INT32",
            Bcf2Type::Float => "FLOAT",
            Bcf2Type::Char => "CHAR",
        })
    }
}

/// A `ByteArrayInputStream` over one block.
#[derive(Debug, Clone)]
pub struct Cursor<'a> {
    bytes: &'a [u8],
    at: usize,
}

impl<'a> Cursor<'a> {
    pub fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, at: 0 }
    }

    /// `(byte) (read() & 0xFF)`: `0xFF` once the block is exhausted.
    pub fn byte(&mut self) -> u8 {
        match self.bytes.get(self.at) {
            Some(&b) => {
                self.at += 1;
                b
            }
            None => 0xFF,
        }
    }

    /// `read(byte[])`: as many as are left, and zeroes after them in the destination, which a
    /// fresh Java array holds.
    pub fn bytes(&mut self, n: usize) -> Vec<u8> {
        let end = (self.at + n).min(self.bytes.len());
        let mut out = self.bytes[self.at..end].to_vec();
        self.at = end;
        out.resize(n, 0);
        out
    }

    pub fn remaining(&self) -> usize {
        self.bytes.len() - self.at
    }
}

/// `BCF2Utils.encodeTypeDescriptor`, with the size masked to four bits: a negative size gives
/// the overflow nibble with no count after it.
pub fn encode_type_descriptor(n_elements: i32, t: Bcf2Type) -> u8 {
    (((0x0F & n_elements) << 4) | (t.id() & 0x0F)) as u8
}

pub fn decode_size(descriptor: u8) -> i32 {
    ((descriptor & 0xF0) >> 4) as i32
}

pub fn decode_type_id(descriptor: u8) -> i32 {
    (descriptor & 0x0F) as i32
}

/// `determineIntegerType(int)`.
pub fn determine_integer_type(value: i32) -> Result<Bcf2Type, Bcf2Error> {
    for t in [Bcf2Type::Int8, Bcf2Type::Int16, Bcf2Type::Int32] {
        if t.within_range(value as i64) {
            return Ok(t);
        }
    }
    Err(Bcf2Error::tribble(format!(
        "Integer cannot be encoded in allowable range of even INT32: {value}"
    )))
}

/// `determineIntegerType(int[])`: from the extremes, both starting at zero.
pub fn determine_integer_type_of_array(values: &[i32]) -> Result<Bcf2Type, Bcf2Error> {
    let max = values.iter().copied().fold(0, i32::max);
    let min = values.iter().copied().fold(0, i32::min);
    let max_type = determine_integer_type(max)?;
    let min_type = determine_integer_type(min)?;
    Ok(max_type.max(min_type))
}

/// `determineIntegerType(List<Integer>)`, which stops at the first value needing 32 bits.
pub fn determine_integer_type_of_list(values: &[i32]) -> Result<Bcf2Type, Bcf2Error> {
    let mut max_type = Bcf2Type::Int8;
    for &value in values {
        match determine_integer_type(value)? {
            Bcf2Type::Int8 => {}
            Bcf2Type::Int16 => max_type = Bcf2Type::Int16,
            _ => return Ok(Bcf2Type::Int32),
        }
    }
    Ok(max_type)
}

/// `maxIntegerType`.
pub fn max_integer_type(t1: Bcf2Type, t2: Bcf2Type) -> Result<Bcf2Type, Bcf2Error> {
    match t1 {
        Bcf2Type::Int8 => Ok(t2),
        Bcf2Type::Int16 => Ok(if t2 == Bcf2Type::Int32 { t2 } else { t1 }),
        Bcf2Type::Int32 => Ok(t1),
        other => Err(Bcf2Error::tribble(format!(
            "BUG: unexpected BCF2Type {other}"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_missing_value_is_outside_the_range() {
        assert_eq!(determine_integer_type(-127).unwrap(), Bcf2Type::Int8);
        assert_eq!(determine_integer_type(-128).unwrap(), Bcf2Type::Int16);
        assert_eq!(determine_integer_type(32768).unwrap(), Bcf2Type::Int32);
        assert!(determine_integer_type(i32::MIN).is_err());
    }

    #[test]
    fn reads_sign_extend() {
        let bytes = [0x80, 0x00, 0x80, 0x01, 0x00, 0x00, 0x80];
        let mut c = Cursor::new(&bytes);
        assert_eq!(
            Bcf2Type::Int8.read(&mut c).unwrap(),
            Bcf2Type::Int8.missing_bytes()
        );
        assert_eq!(
            Bcf2Type::Int16.read(&mut c).unwrap(),
            Bcf2Type::Int16.missing_bytes()
        );
        assert_eq!(Bcf2Type::Int32.read(&mut c).unwrap(), i32::MIN + 1);
        // Past the end every byte is 0xFF.
        assert_eq!(Bcf2Type::Int16.read(&mut c).unwrap(), -1);
    }

    #[test]
    fn a_negative_size_is_the_overflow_nibble() {
        assert_eq!(encode_type_descriptor(-1, Bcf2Type::Int8), 0xF1);
        assert_eq!(encode_type_descriptor(3, Bcf2Type::Char), 0x37);
    }
}
