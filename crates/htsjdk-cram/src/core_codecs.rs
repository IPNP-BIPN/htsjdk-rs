//! The three integer codecs written on the bit stream: Beta, Gamma and Subexponential.
//!
//! Ported from `htsjdk.samtools.cram.encoding.core.BetaIntegerCodec`, `GammaIntegerCodec` and
//! `SubexponentialIntegerCodec` at htsjdk 4.2.0.
//!
//! [`crate::bit_stream`] is the floor. These are the three codecs the encoding map can name that
//! are written on it rather than on an external byte block, and each is a different bargain
//! between a fixed width and a variable one.
//!
//! # Every codec carries an offset
//!
//! Added before storage and subtracted after. It is how a range that starts below zero is stored
//! in bits that cannot hold a sign: `offset = -MIN` puts the whole range at or above zero.
//!
//! # Beta is the only one with an upper bound
//!
//! It refuses what does not fit, with two messages naming the value, the offset and the limit.
//! Gamma and Subexponential grow instead.
//!
//! # Gamma and Subexponential compute a bit length with a floating-point log
//!
//! `1 + (int)(Math.log(v) / Math.log(2))` decides how many bits a value takes, so the bytes
//! written depend on a `double` division landing on the right side of an integer at every power of
//! two. Measured across the powers of two up to `2^31 - 1`, it does; the port uses the same
//! arithmetic rather than an integer bit count, because "it agrees on the values measured" is not
//! the same claim as "it is the same function".
//!
//! # Subexponential has two regimes, split at `2^k`
//!
//! Below it the value goes in `k` bits behind a single `0`. At or above it, `b = floor(log2(v))`
//! bits go behind `b - k + 1` ones and a `0`, and the top bit is implied rather than written.

use crate::bit_stream::{BitError, BitInputStream, BitOutputStream};

/// What a core codec refuses.
///
/// All four are `IllegalArgumentException`, and the wording is the reference's: two of the
/// messages have quirks a port has to keep, a double space in Gamma's and "less then" in
/// Subexponential's.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CodecError {
    /// Beta, when the offset value is below zero.
    BetaNotPositive { value: i32, offset: i32 },
    /// Beta, when the offset value does not fit the declared width.
    BetaAboveLimit { value: i32, offset: i32, limit: i64 },
    /// Gamma, whose length prefix cannot encode a value with no bits.
    GammaNotPositive { value: i32, offset: i32 },
    /// Subexponential, whose message names neither the offset nor the sum.
    SubexponentialBelowOffset { value: i32 },
    /// The bit stream underneath refused.
    Bits(BitError),
}

impl CodecError {
    pub fn message(&self) -> String {
        match self {
            CodecError::BetaNotPositive { value, offset } => {
                format!("Value {value} plus offset {offset} must be positive")
            }
            CodecError::BetaAboveLimit {
                value,
                offset,
                limit,
            } => format!(
                "Value {value} plus offset {offset} is greater than or equal to limit {limit}"
            ),
            // The two spaces after the full stop are the reference's.
            CodecError::GammaNotPositive { value, offset } => format!(
                "Gamma codec handles only positive values.  Value {value} + Offset {offset} <= 0"
            ),
            CodecError::SubexponentialBelowOffset { value } => {
                format!("Value is less then offset: {value}")
            }
            CodecError::Bits(error) => error.message(),
        }
    }

    pub fn java_exception(&self) -> &'static str {
        match self {
            CodecError::Bits(error) => error.java_exception(),
            _ => "IllegalArgumentException",
        }
    }
}

impl From<BitError> for CodecError {
    fn from(error: BitError) -> Self {
        CodecError::Bits(error)
    }
}

/// `1 + (int)(Math.log(value) / Math.log(2))`, the bit length Gamma writes.
///
/// Kept as the floating-point computation rather than replaced by a leading-zero count: the two
/// agree on every value measured, and agreeing is not being the same function.
fn beta_code_length(value: i64) -> i32 {
    1 + ((value as f64).ln() / std::f64::consts::LN_2) as i32
}

/// `BetaIntegerCodec`: an offset and a fixed width.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BetaIntegerCodec {
    offset: i32,
    bits_per_value: i32,
    /// `1 << bitsPerValue`, held as a long because 32 bits does not fit an int.
    value_limit: i64,
}

impl BetaIntegerCodec {
    pub fn new(offset: i32, bits_per_value: i32) -> Self {
        Self {
            offset,
            bits_per_value,
            value_limit: 1i64 << bits_per_value,
        }
    }

    /// `getAndCheckOffsetValue`, which is the only bound any of the three enforces.
    fn checked(&self, value: i32) -> Result<i32, CodecError> {
        let shifted = value.wrapping_add(self.offset);
        if shifted < 0 {
            return Err(CodecError::BetaNotPositive {
                value,
                offset: self.offset,
            });
        }
        if i64::from(shifted) >= self.value_limit {
            return Err(CodecError::BetaAboveLimit {
                value,
                offset: self.offset,
                limit: self.value_limit,
            });
        }
        Ok(shifted)
    }

    pub fn write(&self, out: &mut BitOutputStream, value: i32) -> Result<(), CodecError> {
        let shifted = self.checked(value)?;
        out.write_int_bits(shifted, self.bits_per_value)?;
        Ok(())
    }

    pub fn read(&self, input: &mut BitInputStream<'_>) -> Result<i32, CodecError> {
        Ok(input.read_bits(self.bits_per_value)? - self.offset)
    }
}

/// `GammaIntegerCodec`: Elias gamma, with an offset.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GammaIntegerCodec {
    offset: i32,
}

impl GammaIntegerCodec {
    pub fn new(offset: i32) -> Self {
        Self { offset }
    }

    /// `write`: `length - 1` zeros, then the value in `length` bits with its top bit included.
    pub fn write(&self, out: &mut BitOutputStream, value: i32) -> Result<(), CodecError> {
        if value.wrapping_add(self.offset) < 1 {
            return Err(CodecError::GammaNotPositive {
                value,
                offset: self.offset,
            });
        }
        let shifted = i64::from(value) + i64::from(self.offset);
        let length = beta_code_length(shifted);
        if length > 1 {
            out.write_long_bits(0, length - 1)?;
        }
        out.write_long_bits(shifted, length)?;
        Ok(())
    }

    /// `read`: count the leading zeros, then take that many more bits and put the top one back.
    pub fn read(&self, input: &mut BitInputStream<'_>) -> Result<i32, CodecError> {
        let mut length = 1;
        while !input.read_bit()? {
            length += 1;
        }
        let read = input.read_bits(length - 1)?;
        // The leading one was consumed as the terminator, so it is restored rather than read.
        let value = read | (1 << (length - 1));
        Ok(value - self.offset)
    }
}

/// `SubexponentialIntegerCodec`: two regimes, split at `2^k`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SubexponentialIntegerCodec {
    offset: i32,
    k: i32,
}

impl SubexponentialIntegerCodec {
    pub fn new(offset: i32, k: i32) -> Self {
        Self { offset, k }
    }

    pub fn write(&self, out: &mut BitOutputStream, value: i32) -> Result<(), CodecError> {
        if value.wrapping_add(self.offset) < 0 {
            return Err(CodecError::SubexponentialBelowOffset { value });
        }
        let shifted = i64::from(value) + i64::from(self.offset);
        let (b, u) = if shifted < (1i64 << self.k) {
            (self.k, 0)
        } else {
            // Note the missing `1 +`: this is floor(log2(v)), one less than Gamma's length.
            let b = beta_code_length(shifted) - 1;
            (b, b - self.k + 1)
        };
        out.write_bits(true, u as u64)?;
        out.write_bit(false)?;
        // Only the low `b` bits go out; above the split, the top bit is implied by `u`.
        out.write_long_bits(shifted, b)?;
        Ok(())
    }

    pub fn read(&self, input: &mut BitInputStream<'_>) -> Result<i32, CodecError> {
        let mut u = 0;
        while input.read_bit()? {
            u += 1;
        }
        let (b, n) = if u == 0 {
            let b = self.k;
            (b, input.read_bits(b)?)
        } else {
            let b = u + self.k - 1;
            (b, (1 << b) | input.read_bits(b)?)
        };
        let _ = b;
        Ok(n - self.offset)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn written(f: impl FnOnce(&mut BitOutputStream) -> Result<(), CodecError>) -> Vec<u8> {
        let mut out = BitOutputStream::new();
        f(&mut out).expect("written");
        out.into_bytes()
    }

    /// Beta is a fixed width, left-aligned like everything else on this stream.
    #[test]
    fn beta_writes_a_fixed_width() {
        let codec = BetaIntegerCodec::new(0, 4);
        assert_eq!(written(|out| codec.write(out, 0)), [0x00]);
        assert_eq!(written(|out| codec.write(out, 1)), [0x10]);
        assert_eq!(written(|out| codec.write(out, 15)), [0xF0]);
        assert_eq!(
            written(|out| BetaIntegerCodec::new(0, 32).write(out, i32::MAX)),
            [0x7F, 0xFF, 0xFF, 0xFF]
        );
    }

    /// The offset is what lets a negative range be stored at all.
    #[test]
    fn the_offset_moves_the_range() {
        let codec = BetaIntegerCodec::new(10, 4);
        assert_eq!(written(|out| codec.write(out, -10)), [0x00]);
        assert_eq!(written(|out| codec.write(out, 5)), [0xF0]);
        // And it can move a range down as well as up.
        assert_eq!(
            written(|out| BetaIntegerCodec::new(-5, 4).write(out, 5)),
            [0x00]
        );
    }

    /// Beta is the only one of the three that refuses anything for being too large.
    #[test]
    fn beta_refuses_what_does_not_fit() {
        let codec = BetaIntegerCodec::new(0, 4);
        let mut out = BitOutputStream::new();
        assert_eq!(
            codec.write(&mut out, 16).unwrap_err().message(),
            "Value 16 plus offset 0 is greater than or equal to limit 16"
        );
        assert_eq!(
            codec.write(&mut out, -1).unwrap_err().message(),
            "Value -1 plus offset 0 must be positive"
        );
    }

    /// Gamma: `length - 1` zeros, then the value with its top bit.
    #[test]
    fn gamma_writes_a_length_prefix_of_zeros() {
        let codec = GammaIntegerCodec::new(0);
        assert_eq!(written(|out| codec.write(out, 1)), [0x80]);
        assert_eq!(written(|out| codec.write(out, 2)), [0x40]);
        assert_eq!(written(|out| codec.write(out, 3)), [0x60]);
        assert_eq!(written(|out| codec.write(out, 4)), [0x20]);
        assert_eq!(written(|out| codec.write(out, 16)), [0x08, 0x00]);
        assert_eq!(written(|out| codec.write(out, 31)), [0x0F, 0x80]);
    }

    /// Gamma cannot encode a value with no bits, so zero and below are refused.
    #[test]
    fn gamma_refuses_zero_and_below() {
        let codec = GammaIntegerCodec::new(0);
        let mut out = BitOutputStream::new();
        assert_eq!(
            codec.write(&mut out, 0).unwrap_err().message(),
            "Gamma codec handles only positive values.  Value 0 + Offset 0 <= 0",
            "and the two spaces are the reference's"
        );
        assert!(GammaIntegerCodec::new(1).write(&mut out, 0).is_ok());
    }

    /// Subexponential: below `2^k` the value goes in `k` bits behind a single zero.
    #[test]
    fn subexponential_has_two_regimes() {
        let codec = SubexponentialIntegerCodec::new(0, 2);
        assert_eq!(written(|out| codec.write(out, 0)), [0x00]);
        assert_eq!(written(|out| codec.write(out, 3)), [0x60]);
        // At the split the prefix appears, and the top bit stops being written.
        assert_eq!(written(|out| codec.write(out, 4)), [0x80]);
        assert_eq!(written(|out| codec.write(out, 7)), [0xB0]);
        assert_eq!(written(|out| codec.write(out, 8)), [0xC0]);
        assert_eq!(written(|out| codec.write(out, 32)), [0xF0, 0x00]);
    }

    #[test]
    fn subexponential_refuses_below_its_offset() {
        let codec = SubexponentialIntegerCodec::new(0, 2);
        let mut out = BitOutputStream::new();
        assert_eq!(
            codec.write(&mut out, -1).unwrap_err().message(),
            "Value is less then offset: -1",
            "and the typo is the reference's"
        );
    }

    /// Values in a row pack against each other with no alignment between them.
    #[test]
    fn values_pack_against_each_other() {
        let beta = BetaIntegerCodec::new(0, 3);
        assert_eq!(
            written(|out| {
                for value in 1..=4 {
                    beta.write(out, value)?;
                }
                Ok(())
            }),
            [0x29, 0xC0]
        );
        let gamma = GammaIntegerCodec::new(0);
        assert_eq!(
            written(|out| {
                for value in 1..=4 {
                    gamma.write(out, value)?;
                }
                Ok(())
            }),
            [0xA6, 0x40]
        );
    }

    /// Everything written comes back, which is the only property a codec owes on its own.
    #[test]
    fn every_codec_round_trips() {
        let values: &[i32] = &[0, 1, 2, 3, 4, 7, 8, 15, 16, 31, 32, 1023, 1024, i32::MAX];
        for value in values {
            for k in 0..4 {
                let codec = SubexponentialIntegerCodec::new(0, k);
                let bytes = written(|out| codec.write(out, *value));
                let mut input = BitInputStream::new(&bytes);
                assert_eq!(
                    codec.read(&mut input).expect("back"),
                    *value,
                    "subexp k={k}"
                );
            }
            if *value > 0 {
                let codec = GammaIntegerCodec::new(0);
                let bytes = written(|out| codec.write(out, *value));
                let mut input = BitInputStream::new(&bytes);
                assert_eq!(codec.read(&mut input).expect("back"), *value, "gamma");
            }
        }
        for bits in 1..=32 {
            let codec = BetaIntegerCodec::new(0, bits);
            // A 32-bit Beta still takes an int, so the widest value it can hold is i32::MAX.
            let value = (((1i64 << bits) - 1).min(i64::from(i32::MAX))) as i32;
            let bytes = written(|out| codec.write(out, value));
            let mut input = BitInputStream::new(&bytes);
            assert_eq!(codec.read(&mut input).expect("back"), value, "beta {bits}");
        }
    }
}
