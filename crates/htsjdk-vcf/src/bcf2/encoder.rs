//! `BCF2Encoder`: the byte buffer one block is built in.

use super::types::{
    determine_integer_type, encode_type_descriptor, Bcf2Type, MAX_INLINE_ELEMENTS,
    OVERFLOW_ELEMENT_MARKER,
};
use super::Bcf2Error;

/// `BCF2Encoder`. `getRecordBytes` hands the buffer over and starts a new one.
#[derive(Debug, Default)]
pub struct Bcf2Encoder {
    buffer: Vec<u8>,
}

/// `Float.floatToIntBits((float) value)`: the narrowing rounds to nearest, and every NaN
/// collapses onto the canonical `0x7fc00000`.
pub fn float_to_int_bits(value: f64) -> i32 {
    let f = value as f32;
    if f.is_nan() {
        0x7fc0_0000
    } else {
        f.to_bits() as i32
    }
}

impl Bcf2Encoder {
    pub fn new() -> Self {
        Self::default()
    }

    /// `getRecordBytes`.
    pub fn record_bytes(&mut self) -> Vec<u8> {
        std::mem::take(&mut self.buffer)
    }

    /// `encodeTypedMissing`: a descriptor of size zero.
    pub fn encode_typed_missing(&mut self, t: Bcf2Type) -> Result<(), Bcf2Error> {
        self.encode_type(0, t)
    }

    /// `encodeTypedInt(v)`: the smallest type that holds it.
    pub fn encode_typed_int(&mut self, v: i32) -> Result<(), Bcf2Error> {
        let t = determine_integer_type(v)?;
        self.encode_typed_int_as(v, t)
    }

    pub fn encode_typed_int_as(&mut self, v: i32, t: Bcf2Type) -> Result<(), Bcf2Error> {
        self.encode_type(1, t)?;
        self.encode_raw_int(v, t)
    }

    /// `encodeTypedString(byte[])`.
    pub fn encode_typed_string(&mut self, s: &[u8]) -> Result<(), Bcf2Error> {
        self.encode_type(s.len() as i32, Bcf2Type::Char)?;
        self.buffer.extend_from_slice(s);
        Ok(())
    }

    /// `encodeTyped(List<Integer>, type)`, the only list form the writer uses with an integer
    /// type.
    pub fn encode_typed_ints(&mut self, v: &[i32], t: Bcf2Type) -> Result<(), Bcf2Error> {
        self.encode_type(v.len() as i32, t)?;
        for &i in v {
            self.encode_raw_int(i, t)?;
        }
        Ok(())
    }

    pub fn encode_raw_missing_value(&mut self, t: Bcf2Type) -> Result<(), Bcf2Error> {
        t.write(t.missing_bytes(), &mut self.buffer)
    }

    pub fn encode_raw_char(&mut self, c: u8) {
        self.buffer.push(c);
    }

    pub fn encode_raw_float(&mut self, value: f64) -> Result<(), Bcf2Error> {
        Bcf2Type::Float.write(float_to_int_bits(value), &mut self.buffer)
    }

    /// `encodeType`: inline up to fourteen elements, otherwise the overflow nibble and a typed
    /// count.
    pub fn encode_type(&mut self, size: i32, t: Bcf2Type) -> Result<(), Bcf2Error> {
        if size <= MAX_INLINE_ELEMENTS {
            self.buffer.push(encode_type_descriptor(size, t));
            Ok(())
        } else {
            self.buffer
                .push(encode_type_descriptor(OVERFLOW_ELEMENT_MARKER, t));
            self.encode_typed_int(size)
        }
    }

    /// `encodeRawInt` and `encodeRawBytes`, which are the same call.
    pub fn encode_raw_int(&mut self, value: i32, t: Bcf2Type) -> Result<(), Bcf2Error> {
        t.write(value, &mut self.buffer)
    }

    /// `encodeRawString(s, sizeToWrite)`: the UTF-8 bytes, cut or padded with `\0` to the size.
    pub fn encode_raw_string(&mut self, bytes: &[u8], size_to_write: usize) {
        for i in 0..size_to_write {
            self.buffer.push(bytes.get(i).copied().unwrap_or(0));
        }
    }
}
