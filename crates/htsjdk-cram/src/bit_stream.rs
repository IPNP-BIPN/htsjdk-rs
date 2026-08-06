//! The bit stream the core codecs are written on: the floor under the encodings.
//!
//! Ported from `htsjdk.samtools.cram.io.DefaultBitOutputStream` and `DefaultBitInputStream` at
//! htsjdk 4.2.0.
//!
//! [`crate::read_features`] pins the record model in both directions. What carries it is a set of
//! codecs the encoding map names per data series, and three of them are written on a **bit** stream
//! rather than a byte one. That stream is this module, the way [`crate::varint`] was the floor
//! under the frames.
//!
//! # Bits go in most significant first
//!
//! A value of n bits is left-aligned into a one-byte buffer held back until it fills, so the first
//! bit of the stream is the top bit of the first byte. Writing `1` in one bit produces `0x80`, not
//! `0x01`.
//!
//! # Flush pads with zeros on the right, and the padding is data
//!
//! The buffer was left-aligned to begin with, so what is left of it goes out as it stands. Nothing
//! records how many bits were real: only the count of values a reader expects says where the
//! stream ends. Measured, a stream of one `true` bit and a stream of `0x80` in eight bits are the
//! same byte.
//!
//! # A multi-byte write splits at the top
//!
//! `write(long, n)` writes whole bytes from the **most significant** end while at least eight bits
//! remain, then the remainder. The leftover bits of a twelve-bit write are the low four, written
//! last: `0xABC` in twelve bits is `ab c0`.
//!
//! # The bounds are checked with three different wordings, and one case is not checked at all
//!
//! Each overload has its own message. And `write(byte, 0)` against a **non-empty** buffer indexes
//! a mask table of eight entries with 8, so writing nothing is a no-op or an exception depending
//! on what was written before it.

/// What a bit stream refuses, with the exception htsjdk raises.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BitError {
    /// `RuntimeIOException`, from the `long` overload, whose message quotes the value too.
    ExpectingOneToSixtyFourBits { value: i64, bits: i32 },
    /// `RuntimeIOException`, from the `int` overload.
    ExpectingOneToThirtyTwoBits,
    /// `RuntimeIOException`, from the `byte` overload.
    ExpectingZeroToEightBits,
    /// `ArrayIndexOutOfBoundsException` from `bitMasks[8 - nofBits]`, which only a zero-bit write
    /// against a non-empty buffer reaches.
    MaskIndexOutOfBounds,
    /// `RuntimeEOFException`, which carries no more than this.
    EndOfStream,
    /// `RuntimeException` from `readLongBits`, whose message has a typo the port keeps.
    MoreThanSixtyFourBitsRequested,
}

impl BitError {
    pub fn message(&self) -> String {
        match self {
            BitError::ExpectingOneToSixtyFourBits { value, bits } => {
                format!("Expecting 1 to 64 bits, got: value={value}, nofBits={bits}")
            }
            BitError::ExpectingOneToThirtyTwoBits => "Expecting 1 to 32 bits.".to_string(),
            BitError::ExpectingZeroToEightBits => "Expecting 0 to 8 bits.".to_string(),
            BitError::MaskIndexOutOfBounds => "Index 8 out of bounds for length 8".to_string(),
            BitError::EndOfStream => "End of stream.".to_string(),
            BitError::MoreThanSixtyFourBitsRequested => {
                "More then 64 bits are requested in one read from bit stream.".to_string()
            }
        }
    }

    /// The Java exception it arrives as. Three of the six are the JDK's rather than htsjdk's.
    pub fn java_exception(&self) -> &'static str {
        match self {
            BitError::ExpectingOneToSixtyFourBits { .. }
            | BitError::ExpectingOneToThirtyTwoBits
            | BitError::ExpectingZeroToEightBits => "RuntimeIOException",
            BitError::MaskIndexOutOfBounds => "ArrayIndexOutOfBoundsException",
            BitError::EndOfStream => "RuntimeEOFException",
            BitError::MoreThanSixtyFourBitsRequested => "RuntimeException",
        }
    }
}

/// `DefaultBitOutputStream`.
#[derive(Debug, Default)]
pub struct BitOutputStream {
    out: Vec<u8>,
    buffer_byte: u8,
    buffered_bits: i32,
}

impl BitOutputStream {
    pub fn new() -> Self {
        Self::default()
    }

    /// How many bits are held back, which is what [`Self::flush`] will pad out.
    pub fn buffered_bits(&self) -> i32 {
        self.buffered_bits
    }

    /// `write(byte, nofBits)`, the one every other overload ends in.
    ///
    /// The value is masked to the requested width, so `0xFF` in three bits is `0b111` and the bits
    /// above it are dropped rather than refused.
    pub fn write_byte_bits(&mut self, value: u8, bits: i32) -> Result<(), BitError> {
        if !(0..=8).contains(&bits) {
            return Err(BitError::ExpectingZeroToEightBits);
        }
        if bits == 8 {
            self.write_whole_byte(value);
            return Ok(());
        }
        if self.buffered_bits == 0 {
            // A zero-bit write lands here and does nothing at all.
            self.buffer_byte = ((value as u32) << (8 - bits)) as u8;
            self.buffered_bits = bits;
            return Ok(());
        }
        if bits == 0 {
            // `bitMasks[8 - 0]` is one past the end of a table of eight.
            return Err(BitError::MaskIndexOutOfBounds);
        }
        // `bitContainer & ~bitMasks[8 - nofBits]`: keep the low `bits` bits.
        let container = value & !mask(8 - bits);
        let remaining = 8 - self.buffered_bits - bits;
        if remaining < 0 {
            let shift = -remaining;
            self.buffer_byte |= container >> shift;
            self.out.push(self.buffer_byte);
            self.buffer_byte = ((container as u32) << (8 - shift)) as u8;
            self.buffered_bits = shift;
        } else if remaining == 0 {
            self.buffer_byte |= container;
            self.out.push(self.buffer_byte);
            self.buffered_bits = 0;
        } else {
            self.buffer_byte |= container << remaining;
            self.buffered_bits = 8 - remaining;
        }
        Ok(())
    }

    /// `writeByte`, which is where a whole byte straddles the buffer.
    fn write_whole_byte(&mut self, value: u8) {
        if self.buffered_bits == 0 {
            self.out.push(value);
        } else {
            self.buffer_byte |= value >> self.buffered_bits;
            self.out.push(self.buffer_byte);
            self.buffer_byte = ((value as u32) << (8 - self.buffered_bits)) as u8;
        }
    }

    /// `write(long, nofBits)`.
    pub fn write_long_bits(&mut self, value: i64, bits: i32) -> Result<(), BitError> {
        if bits == 0 {
            return Ok(());
        }
        if !(1..=64).contains(&bits) {
            return Err(BitError::ExpectingOneToSixtyFourBits { value, bits });
        }
        if bits <= 8 {
            return self.write_byte_bits(value as u8, bits);
        }
        // Whole bytes from the top while at least eight bits remain, then the leftover low bits.
        let mut i = bits - 8;
        while i >= 0 {
            self.write_whole_byte(((value as u64) >> i) as u8);
            i -= 8;
        }
        if bits % 8 != 0 {
            self.write_byte_bits(value as u8, bits % 8)?;
        }
        Ok(())
    }

    /// `write(int, nofBits)`, which is `write_int_LSB_0`.
    pub fn write_int_bits(&mut self, value: i32, bits: i32) -> Result<(), BitError> {
        if bits == 0 {
            return Ok(());
        }
        if !(1..=32).contains(&bits) {
            return Err(BitError::ExpectingOneToThirtyTwoBits);
        }
        if bits <= 8 {
            return self.write_byte_bits(value as u8, bits);
        }
        let mut i = bits - 8;
        while i >= 0 {
            self.write_whole_byte(((value as u32) >> i) as u8);
            i -= 8;
        }
        if bits % 8 != 0 {
            self.write_byte_bits(value as u8, bits % 8)?;
        }
        Ok(())
    }

    /// `write(boolean)`.
    pub fn write_bit(&mut self, bit: bool) -> Result<(), BitError> {
        self.write_byte_bits(u8::from(bit), 1)
    }

    /// `write(boolean, repeat)`, which is that one in a loop.
    pub fn write_bits(&mut self, bit: bool, repeat: u64) -> Result<(), BitError> {
        for _ in 0..repeat {
            self.write_bit(bit)?;
        }
        Ok(())
    }

    /// `flush`: the partial byte goes out as it stands, zero-padded on the right.
    pub fn flush(&mut self) {
        if self.buffered_bits > 0 {
            self.out.push(self.buffer_byte);
        }
        self.buffered_bits = 0;
        self.buffer_byte = 0;
    }

    /// Flush and take the bytes.
    pub fn into_bytes(mut self) -> Vec<u8> {
        self.flush();
        self.out
    }
}

/// `bitMasks[i]`, which is `~(0xFF >>> i)`: the top `i` bits set.
fn mask(index: i32) -> u8 {
    !((0xFFu32 >> index) as u8)
}

/// `DefaultBitInputStream`.
pub struct BitInputStream<'a> {
    bytes: &'a [u8],
    at: usize,
    byte_buffer: i32,
    buffered_bits: i32,
}

impl<'a> BitInputStream<'a> {
    pub fn new(bytes: &'a [u8]) -> Self {
        Self {
            bytes,
            at: 0,
            byte_buffer: 0,
            buffered_bits: 0,
        }
    }

    /// `in.read()`: the next byte, or the end of the stream.
    fn next_byte(&mut self) -> Result<i32, BitError> {
        match self.bytes.get(self.at) {
            Some(byte) => {
                self.at += 1;
                Ok(i32::from(*byte))
            }
            None => Err(BitError::EndOfStream),
        }
    }

    /// `readBit`.
    pub fn read_bit(&mut self) -> Result<bool, BitError> {
        self.buffered_bits -= 1;
        if self.buffered_bits >= 0 {
            return Ok((self.byte_buffer >> self.buffered_bits) & 1 == 1);
        }
        // The count is set before the read, so a stream that ends here has already been counted
        // down; the exception is what stops it mattering.
        self.buffered_bits = 7;
        self.byte_buffer = self.next_byte()?;
        Ok((self.byte_buffer >> 7) & 1 == 1)
    }

    /// `readBits`.
    pub fn read_bits(&mut self, mut n: i32) -> Result<i32, BitError> {
        if n == 0 {
            return Ok(0);
        }
        let mut x = 0i32;
        while n > self.buffered_bits {
            n -= self.buffered_bits;
            // Java's `<<` masks the count to five bits, and a count of 32 does reach here.
            x |= right_bits(self.buffered_bits, self.byte_buffer).wrapping_shl(n as u32);
            self.byte_buffer = self.next_byte()?;
            self.buffered_bits = 8;
        }
        self.buffered_bits -= n;
        Ok(x | right_bits(n, self.byte_buffer >> self.buffered_bits))
    }

    /// `readLongBits`, which reads its first byte before it looks at what it has.
    pub fn read_long_bits(&mut self, mut n: i32) -> Result<i64, BitError> {
        if n > 64 {
            return Err(BitError::MoreThanSixtyFourBitsRequested);
        }
        if n == 0 {
            return Ok(0);
        }
        let mut x = 0i64;
        let mut byte_buffer = i64::from(self.byte_buffer);
        if self.buffered_bits == 0 {
            byte_buffer = i64::from(self.next_byte()?);
            self.buffered_bits = 8;
        }
        byte_buffer &= long_mask(self.buffered_bits);
        while n > self.buffered_bits {
            n -= self.buffered_bits;
            x |= byte_buffer.wrapping_shl(n as u32);
            byte_buffer = i64::from(self.next_byte()?);
            self.buffered_bits = 8;
        }
        self.buffered_bits -= n;
        self.byte_buffer = (byte_buffer & long_mask(self.buffered_bits)) as i32;
        Ok(x | (byte_buffer >> self.buffered_bits))
    }
}

/// `rightBits`, whose mask is `(1 << n) - 1` with Java's shift semantics: a count of 32 masks to
/// 0, so the mask is 0 and the answer is 0 rather than the whole word.
fn right_bits(n: i32, x: i32) -> i32 {
    x & (1i32.wrapping_shl(n as u32) - 1)
}

/// `masks[n]`, which is `(1 << n) - 1` for zero to eight.
fn long_mask(n: i32) -> i64 {
    if n == 0 {
        0
    } else {
        (1i64 << n) - 1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn written(f: impl FnOnce(&mut BitOutputStream) -> Result<(), BitError>) -> Vec<u8> {
        let mut out = BitOutputStream::new();
        f(&mut out).expect("written");
        out.into_bytes()
    }

    /// The first bit of the stream is the top bit of the first byte.
    #[test]
    fn a_one_bit_value_is_the_high_bit_of_the_byte() {
        assert_eq!(written(|out| out.write_byte_bits(1, 1)), [0x80]);
        assert_eq!(written(|out| out.write_byte_bits(1, 8)), [0x01]);
        assert_eq!(written(|out| out.write_byte_bits(1, 4)), [0x10]);
    }

    /// The value is masked to the width, not refused for exceeding it.
    #[test]
    fn a_value_too_wide_for_its_width_is_masked() {
        assert_eq!(written(|out| out.write_byte_bits(0xFF, 3)), [0xE0]);
        assert_eq!(written(|out| out.write_byte_bits(0x07, 3)), [0xE0]);
    }

    /// Two writes that cross the buffered byte, which is where the shifts matter.
    #[test]
    fn a_write_straddles_the_buffered_byte() {
        assert_eq!(
            written(|out| {
                out.write_byte_bits(0x0A, 4)?;
                out.write_byte_bits(0x05, 4)
            }),
            [0xA5]
        );
        assert_eq!(
            written(|out| {
                out.write_byte_bits(0x15, 5)?;
                out.write_byte_bits(0x0A, 5)
            }),
            [0xAA, 0x80]
        );
        // A whole byte written against a one-bit buffer splits across two bytes.
        assert_eq!(
            written(|out| {
                out.write_bit(true)?;
                out.write_byte_bits(0xFF, 8)
            }),
            [0xFF, 0x80]
        );
    }

    /// A multi-byte write splits at the top: the leftover bits are the low ones, written last.
    #[test]
    fn a_multi_byte_write_leaves_the_low_bits_for_last() {
        assert_eq!(written(|out| out.write_long_bits(0xABC, 12)), [0xAB, 0xC0]);
        assert_eq!(written(|out| out.write_long_bits(0x1FF, 9)), [0xFF, 0x80]);
        assert_eq!(
            written(|out| out.write_long_bits(0x0123_4567_89AB_CDEF, 64)),
            [0x01, 0x23, 0x45, 0x67, 0x89, 0xAB, 0xCD, 0xEF]
        );
        assert_eq!(
            written(|out| out.write_int_bits(0xDEAD_BEEFu32 as i32, 32)),
            [0xDE, 0xAD, 0xBE, 0xEF]
        );
    }

    /// The padding a flush leaves cannot be told from data.
    #[test]
    fn a_flush_pads_with_zeros_that_read_back_as_data() {
        let padded = written(|out| out.write_bit(true));
        assert_eq!(padded, [0x80]);
        let mut input = BitInputStream::new(&padded);
        assert!(input.read_bit().expect("the real bit"));
        assert_eq!(
            input.read_bits(7).expect("the padding"),
            0,
            "and nothing says these seven bits were never written"
        );
    }

    /// Nothing written is nothing at all, not a zero byte.
    #[test]
    fn an_empty_stream_is_empty() {
        assert_eq!(written(|_| Ok(())), Vec::<u8>::new());
        assert_eq!(written(|out| out.write_byte_bits(1, 0)), Vec::<u8>::new());
    }

    /// A zero-bit write is a no-op or an exception depending on what came before it.
    #[test]
    fn a_zero_bit_write_against_a_partial_buffer_is_an_index_out_of_a_table_of_eight() {
        let mut out = BitOutputStream::new();
        out.write_byte_bits(0x0A, 4).expect("four bits");
        assert_eq!(
            out.write_byte_bits(1, 0),
            Err(BitError::MaskIndexOutOfBounds)
        );
        assert_eq!(
            BitError::MaskIndexOutOfBounds.message(),
            "Index 8 out of bounds for length 8"
        );
    }

    #[test]
    fn the_bounds_are_checked_with_their_own_wording() {
        let mut out = BitOutputStream::new();
        assert_eq!(
            out.write_long_bits(1, 65),
            Err(BitError::ExpectingOneToSixtyFourBits { value: 1, bits: 65 })
        );
        assert_eq!(
            out.write_long_bits(1, 65).unwrap_err().message(),
            "Expecting 1 to 64 bits, got: value=1, nofBits=65"
        );
        assert_eq!(
            out.write_int_bits(1, 33),
            Err(BitError::ExpectingOneToThirtyTwoBits)
        );
        assert_eq!(
            out.write_byte_bits(1, 9),
            Err(BitError::ExpectingZeroToEightBits)
        );
    }

    /// What comes back out, byte boundaries and all.
    #[test]
    fn the_reader_assembles_across_byte_boundaries() {
        let mut input = BitInputStream::new(&[0xAB, 0xCD]);
        assert_eq!(input.read_bits(5).expect("five"), 21);
        assert_eq!(input.read_bits(5).expect("five"), 15);
        assert_eq!(input.read_bits(6).expect("six"), 13);

        let mut wide = BitInputStream::new(&[0xAB, 0xCD]);
        assert_eq!(wide.read_bits(12).expect("twelve"), 2748);

        let mut long = BitInputStream::new(&[0xAB, 0xCD]);
        assert_eq!(long.read_long_bits(5).expect("five"), 21);
        assert_eq!(long.read_long_bits(11).expect("eleven"), 973);
    }

    /// End of stream is an exception, not a zero, and there is no way to ask first.
    #[test]
    fn reading_past_the_end_is_an_exception() {
        let mut input = BitInputStream::new(&[0xAB]);
        assert_eq!(input.read_bits(8).expect("the byte"), 0xAB);
        assert_eq!(input.read_bit(), Err(BitError::EndOfStream));

        let mut empty = BitInputStream::new(&[]);
        assert_eq!(empty.read_bit(), Err(BitError::EndOfStream));
        assert_eq!(BitError::EndOfStream.message(), "End of stream.");

        let mut wide = BitInputStream::new(&[0xAB]);
        assert_eq!(
            wide.read_long_bits(65),
            Err(BitError::MoreThanSixtyFourBitsRequested)
        );
        assert_eq!(
            BitError::MoreThanSixtyFourBitsRequested.message(),
            "More then 64 bits are requested in one read from bit stream.",
            "the typo is the reference's"
        );
    }

    /// A round trip over every width, which is what a codec above this will rely on.
    #[test]
    fn every_width_round_trips() {
        for bits in 1..=32 {
            let value = ((1i64 << bits) - 1) as i32;
            let bytes = written(|out| out.write_int_bits(value, bits));
            let mut input = BitInputStream::new(&bytes);
            assert_eq!(input.read_bits(bits).expect("back"), value, "{bits} bits");
        }
    }
}
