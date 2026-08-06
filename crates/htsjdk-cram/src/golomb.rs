//! Golomb, Golomb-Rice and Golomb-Long: the three codecs the CRAM specification is removing.
//!
//! Ported from `htsjdk.samtools.cram.encoding.core.experimental` at htsjdk 4.2.0.
//!
//! They are reachable all the same. The encoding factory dispatches to them by identifier, so a
//! file that names one has to be read, and a port that skips them cannot claim to read every legal
//! file. htsjdk marks them experimental and logs a warning when one is built; nothing else guards
//! them, and neither does this.
//!
//! # The quotient is unary
//!
//! A value costs one bit per multiple of `m` it contains, so the corpus behind this module keeps
//! every value small next to its `m`. One long of `2^32` with `m` of 10 writes four hundred million
//! ones.
//!
//! # Golomb does not round-trip a negative
//!
//! With offset 0 and `m` 4, `-1` is written as `60` and read back as `3`. Java's division truncates
//! towards zero and its remainder takes the sign of the dividend, so a negative quotient writes no
//! unary bits at all and a negative remainder is written as its low bits. Nothing reports it. The
//! offset exists to keep values at or above zero, and this is what happens when it does not.
//!
//! # Golomb-Rice's parameter is not `m`
//!
//! `GolombRiceIntegerEncoding` calls it `m` and hands it to the codec as `log2m`, so an encoding
//! built with 8 divides by 256. It also accepts what Golomb refuses: `m < 2` is an exception there
//! and a different encoding here.
//!
//! # The remainder is written at one of two widths
//!
//! `ceiling - 1` or `ceiling`, chosen by comparing the remainder against `2^ceiling - m`, which is
//! what lets an `m` that is not a power of two avoid wasting a bit. `ceiling` comes from
//! `(int)(Math.log(m) / Math.log(2) + 1)`, and the comparisons are against `Math.pow`, so an
//! integer remainder is promoted to a double before being compared and again before being
//! subtracted.

use crate::bit_stream::{BitError, BitInputStream, BitOutputStream};

/// What the three experimental codecs refuse.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GolombError {
    /// Golomb and Golomb-Long, on a divisor below two. Golomb-Rice has no such check.
    MTooSmall,
    /// `read(length)` on Golomb and Golomb-Long.
    MultiValueRead,
    /// `read(length)` on Golomb-Rice, which words the same refusal differently.
    NotImplemented,
    /// The bit stream underneath.
    Bits(BitError),
}

impl GolombError {
    pub fn message(&self) -> String {
        match self {
            GolombError::MTooSmall => "M parameter must be at least 2.".to_string(),
            GolombError::MultiValueRead => "Multi-value read method not defined.".to_string(),
            GolombError::NotImplemented => "Not implemented.".to_string(),
            GolombError::Bits(error) => error.message(),
        }
    }

    pub fn java_exception(&self) -> &'static str {
        match self {
            GolombError::MTooSmall => "IllegalArgumentException",
            GolombError::MultiValueRead | GolombError::NotImplemented => "RuntimeException",
            GolombError::Bits(error) => error.java_exception(),
        }
    }
}

impl From<BitError> for GolombError {
    fn from(error: BitError) -> Self {
        GolombError::Bits(error)
    }
}

/// `(int)(Math.log(m) / Math.log(2) + 1)`.
///
/// Not the same expression as the core codecs' `1 + (int)(Math.log(v) / Math.log(2))`: the addition
/// is inside the truncation here and outside it there. They agree on every `m` measured, and
/// agreeing is not being the same function, so both are kept as written.
fn ceiling(m: i32) -> i32 {
    ((f64::from(m).ln() / std::f64::consts::LN_2) + 1.0) as i32
}

/// `Math.pow(2, ceiling) - m`, the threshold both widths are chosen by. Held as a double because
/// that is what the reference compares against.
fn threshold(ceiling: i32, m: i32) -> f64 {
    2f64.powi(ceiling) - f64::from(m)
}

/// `GolombIntegerCodec`: a unary quotient and a remainder at one of two widths.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GolombIntegerCodec {
    offset: i32,
    m: i32,
}

impl GolombIntegerCodec {
    /// The divisor is checked in the constructor, which is where the reference checks it.
    pub fn new(offset: i32, m: i32) -> Result<Self, GolombError> {
        if m < 2 {
            return Err(GolombError::MTooSmall);
        }
        Ok(Self { offset, m })
    }

    pub fn write(&self, out: &mut BitOutputStream, value: i32) -> Result<(), GolombError> {
        let shifted = value.wrapping_add(self.offset);
        let quotient = shifted / self.m;
        let remainder = shifted % self.m;
        let ceiling = ceiling(self.m);

        // `write(bit, quotient)` counts up to the quotient, so a negative one writes nothing at
        // all. That is the first half of why a negative value does not survive the round trip.
        if quotient > 0 {
            out.write_bits(true, quotient as u64)?;
        }
        out.write_bit(false)?;

        if f64::from(remainder) < threshold(ceiling, self.m) {
            out.write_int_bits(remainder, ceiling - 1)?;
        } else {
            out.write_int_bits(
                (f64::from(remainder) + threshold(ceiling, self.m)) as i32,
                ceiling,
            )?;
        }
        Ok(())
    }

    pub fn read(&self, input: &mut BitInputStream<'_>) -> Result<i32, GolombError> {
        let mut quotient = 0i32;
        while input.read_bit()? {
            quotient += 1;
        }

        let ceiling = ceiling(self.m);
        let mut remainder = input.read_bits(ceiling - 1)?;
        if f64::from(remainder) >= threshold(ceiling, self.m) {
            remainder <<= 1;
            remainder |= input.read_bits(1)?;
            // A compound assignment against a double truncates back to int, which is the
            // reference's arithmetic and not a convenience of the port.
            remainder = (f64::from(remainder) - threshold(ceiling, self.m)) as i32;
        }

        Ok(quotient
            .wrapping_mul(self.m)
            .wrapping_add(remainder)
            .wrapping_sub(self.offset))
    }

    /// `read(length)`, which the reference defines only to refuse.
    pub fn read_with_length(&self, _length: i32) -> Result<i32, GolombError> {
        Err(GolombError::MultiValueRead)
    }
}

/// `GolombRiceIntegerCodec`: the same shape with a power of two for a divisor, and no check on it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GolombRiceIntegerCodec {
    offset: i32,
    log2m: i32,
    m: i32,
    mask: i64,
}

impl GolombRiceIntegerCodec {
    /// `log2m` is what `GolombRiceIntegerEncoding` calls `m`. An encoding built with 8 divides by
    /// 256, and nothing in the reference reconciles the two names.
    pub fn new(offset: i32, log2m: i32) -> Self {
        Self {
            offset,
            log2m,
            m: 1i32.wrapping_shl(log2m as u32),
            // `~(~0 << log2m)` on an int, then widened to a long, so a count of 0 or 32 gives 0.
            mask: i64::from(!(!0i32).wrapping_shl(log2m as u32)),
        }
    }

    pub fn write(&self, out: &mut BitOutputStream, value: i32) -> Result<(), GolombError> {
        // The addition happens in int and is only then widened, so it overflows where an int does.
        let shifted = i64::from(value.wrapping_add(self.offset));
        let quotient = ((shifted as u64) >> (self.log2m as u32 & 63)) as i64;

        if quotient > 0 {
            out.write_bits(true, quotient as u64)?;
        }
        out.write_bit(false)?;

        let remainder = shifted & self.mask;
        for bit in (0..self.log2m).rev() {
            out.write_bit(remainder & (1i64 << bit) != 0)?;
        }
        Ok(())
    }

    pub fn read(&self, input: &mut BitInputStream<'_>) -> Result<i32, GolombError> {
        let mut unary = 0i32;
        while input.read_bit()? {
            unary += 1;
        }
        let remainder = input.read_bits(self.log2m)?;
        Ok(unary
            .wrapping_mul(self.m)
            .wrapping_add(remainder)
            .wrapping_sub(self.offset))
    }

    /// `read(length)`, refused in the reference's other wording.
    pub fn read_with_length(&self, _length: i32) -> Result<i32, GolombError> {
        Err(GolombError::NotImplemented)
    }
}

/// `GolombLongCodec`: the same arithmetic on a long, and one cast that did not follow it there.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GolombLongCodec {
    offset: i64,
    m: i32,
}

impl GolombLongCodec {
    pub fn new(offset: i64, m: i32) -> Result<Self, GolombError> {
        if m < 2 {
            return Err(GolombError::MTooSmall);
        }
        Ok(Self { offset, m })
    }

    pub fn write(&self, out: &mut BitOutputStream, value: i64) -> Result<(), GolombError> {
        let shifted = value.wrapping_add(self.offset);
        let quotient = shifted / i64::from(self.m);
        let remainder = shifted % i64::from(self.m);
        let ceiling = ceiling(self.m);

        if quotient > 0 {
            out.write_bits(true, quotient as u64)?;
        }
        out.write_bit(false)?;

        if (remainder as f64) < threshold(ceiling, self.m) {
            // The narrow branch writes a long and the wide one casts to an int. The asymmetry is
            // the reference's.
            out.write_long_bits(remainder, ceiling - 1)?;
        } else {
            out.write_int_bits(
                (remainder as f64 + threshold(ceiling, self.m)) as i32,
                ceiling,
            )?;
        }
        Ok(())
    }

    pub fn read(&self, input: &mut BitInputStream<'_>) -> Result<i64, GolombError> {
        let mut quotient = 0i64;
        while input.read_bit()? {
            quotient += 1;
        }

        let ceiling = ceiling(self.m);
        let mut remainder = i64::from(input.read_bits(ceiling - 1)?);
        if (remainder as f64) >= threshold(ceiling, self.m) {
            remainder <<= 1;
            remainder |= i64::from(input.read_bits(1)?);
            remainder = (remainder as f64 - threshold(ceiling, self.m)) as i64;
        }

        Ok(quotient
            .wrapping_mul(i64::from(self.m))
            .wrapping_add(remainder)
            .wrapping_sub(self.offset))
    }

    pub fn read_with_length(&self, _length: i32) -> Result<i64, GolombError> {
        Err(GolombError::MultiValueRead)
    }
}

/// The encoding parameters of all three: the offset and the divisor, both ITF8.
pub fn serialize_params(offset: i32, m: i32) -> Vec<u8> {
    let mut out = crate::varint::write_unsigned_itf8(offset).0;
    out.extend_from_slice(&crate::varint::write_unsigned_itf8(m).0);
    out
}

/// The offset and the divisor, read back.
pub fn parse_params(bytes: &[u8]) -> Result<(i32, i32), crate::varint::RuntimeEof> {
    let (offset, used) = crate::varint::read_unsigned_itf8(bytes)?;
    let (m, _) = crate::varint::read_unsigned_itf8(bytes.get(used..).unwrap_or(&[]))?;
    Ok((offset, m))
}
