//! ITF8 and LTF8: the variable-length integers every other CRAM structure is measured in.
//!
//! Ported from `htsjdk.samtools.cram.io.ITF8` and `htsjdk.samtools.cram.io.LTF8` at htsjdk 4.2.0.
//!
//! This is the floor of the CRAM port. A container header is a run of ITF8s, a slice header is a
//! run of ITF8s, and a compression header's encoding parameters are ITF8s inside a byte array.
//! Nothing above them can be checked until they are, and they have two properties a reading of the
//! specification does not give you.
//!
//! # The five-byte ITF8 stores four bits twice, and the reader believes one copy
//!
//! The writer puts `(value >> 4) & 0xFF` in byte four and `value & 0xFF` in byte five, so bits 4 to
//! 7 appear in **both**. The reader takes byte four whole and masks byte five to its low nibble:
//!
//! ```text
//! ((b1 & 15) << 28) | b2 << 20 | b3 << 12 | b4 << 4 | (15 & b5)
//! ```
//!
//! Measured, `f0 00 00 01 12` and `f0 00 00 01 f2` both read **18**: the high nibble of byte five
//! is discarded. A stream whose two copies disagree is not an error, it resolves silently to byte
//! four's, and a port that reads the fifth byte whole answers differently on exactly those streams.
//!
//! # A truncated stream is a wrong number, not a refusal
//!
//! `InputStream.read()` returns `-1` at end of stream and nothing checks it after the first byte,
//! so the arithmetic proceeds with it. Measured, the two-byte form `80` reads **-1** and the
//! five-byte form `f0 00` reads **-1**. Only a stream that is empty at the *first* byte throws, and
//! it throws `RuntimeEOFException` with a null message. So a truncated CRAM does not fail here; it
//! produces a number, and the failure surfaces somewhere else or not at all.
//!
//! That is why [`read_unsigned_itf8`] returns the number rather than an error for a short slice:
//! reproducing the refusal a port would naturally add would be a divergence, and a silently wrong
//! number is what the reference produces.
//!
//! # Negative values go through unharmed and are not "unsigned"
//!
//! `writeUnsignedITF8` takes an `int`; a negative one has its high bits set, so it always takes the
//! five-byte form and round-trips through a reader that returns an `int`. Measured, `-1` writes
//! `ff ff ff ff ff` and reads back `-1`. The name says unsigned and the type does not.

/// `ITF8.MAX_BYTES`.
pub const ITF8_MAX_BYTES: usize = 5;

/// What a read at a genuinely empty position throws.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RuntimeEof;

impl RuntimeEof {
    pub fn class(&self) -> &'static str {
        "htsjdk.samtools.util.RuntimeEOFException"
    }

    /// The message is `null`: the exception is constructed with no argument.
    pub fn message(&self) -> &'static str {
        "null"
    }
}

/// A cursor whose `read` returns `-1` past the end, which is the behaviour the arithmetic below
/// depends on.
struct Stream<'a> {
    bytes: &'a [u8],
    at: usize,
}

impl Stream<'_> {
    /// `InputStream.read()`: the byte as 0..255, or -1 at end of stream.
    fn read(&mut self) -> i32 {
        match self.bytes.get(self.at) {
            Some(byte) => {
                self.at += 1;
                i32::from(*byte)
            }
            None => {
                self.at += 1;
                -1
            }
        }
    }
}

/// `ITF8.readUnsignedITF8(InputStream)`.
///
/// Returns the number of bytes consumed alongside the value, which the stream form does not need
/// and a caller decoding a header does.
pub fn read_unsigned_itf8(bytes: &[u8]) -> Result<(i32, usize), RuntimeEof> {
    let mut stream = Stream { bytes, at: 0 };
    let b1 = stream.read();
    if b1 == -1 {
        return Err(RuntimeEof);
    }
    let value = if (b1 & 128) == 0 {
        b1
    } else if (b1 & 64) == 0 {
        ((b1 & 127) << 8) | stream.read()
    } else if (b1 & 32) == 0 {
        let b2 = stream.read();
        let b3 = stream.read();
        ((b1 & 63) << 16) | (b2 << 8) | b3
    } else if (b1 & 16) == 0 {
        ((b1 & 31) << 24) | (stream.read() << 16) | (stream.read() << 8) | stream.read()
    } else {
        // The fifth byte contributes its low nibble only. Its high nibble is written by the
        // encoder and thrown away here; see the module doc.
        ((b1 & 15) << 28)
            | (stream.read() << 20)
            | (stream.read() << 12)
            | (stream.read() << 4)
            | (15 & stream.read())
    };
    Ok((value, stream.at))
}

/// `ITF8.writeUnsignedITF8`, returning the bytes and the **bit count** the Java returns.
pub fn write_unsigned_itf8(value: i32) -> (Vec<u8>, i32) {
    let unsigned = value as u32;
    if (unsigned >> 7) == 0 {
        (vec![value as u8], 8)
    } else if (unsigned >> 14) == 0 {
        (vec![((value >> 8) | 0x80) as u8, (value & 0xFF) as u8], 16)
    } else if (unsigned >> 21) == 0 {
        (
            vec![
                ((value >> 16) | 0xC0) as u8,
                ((value >> 8) & 0xFF) as u8,
                (value & 0xFF) as u8,
            ],
            24,
        )
    } else if (unsigned >> 28) == 0 {
        (
            vec![
                ((value >> 24) | 0xE0) as u8,
                ((value >> 16) & 0xFF) as u8,
                ((value >> 8) & 0xFF) as u8,
                (value & 0xFF) as u8,
            ],
            32,
        )
    } else {
        (
            vec![
                ((value >> 28) | 0xF0) as u8,
                ((value >> 20) & 0xFF) as u8,
                ((value >> 12) & 0xFF) as u8,
                ((value >> 4) & 0xFF) as u8,
                (value & 0xFF) as u8,
            ],
            40,
        )
    }
}

/// `LTF8.readUnsignedLTF8(InputStream)`.
///
/// Nine arms rather than five, and the last two are the ones with no length bits left: `0xFE`
/// means eight bytes follow and `0xFF` means eight bytes follow as well, the difference being
/// which bits of the first byte are part of the value. In the `0xFF` arm none are.
pub fn read_unsigned_ltf8(bytes: &[u8]) -> Result<(i64, usize), RuntimeEof> {
    let mut stream = Stream { bytes, at: 0 };
    let b1 = stream.read();
    if b1 == -1 {
        return Err(RuntimeEof);
    }
    let next = |stream: &mut Stream| i64::from(stream.read());
    let value = if (b1 & 128) == 0 {
        i64::from(b1)
    } else if (b1 & 64) == 0 {
        (i64::from(b1 & 127) << 8) | next(&mut stream)
    } else if (b1 & 32) == 0 {
        (i64::from(b1 & 63) << 16) | (next(&mut stream) << 8) | next(&mut stream)
    } else if (b1 & 16) == 0 {
        (i64::from(b1 & 31) << 24)
            | (next(&mut stream) << 16)
            | (next(&mut stream) << 8)
            | next(&mut stream)
    } else if (b1 & 8) == 0 {
        (i64::from(b1 & 15) << 32)
            | (next(&mut stream) << 24)
            | (next(&mut stream) << 16)
            | (next(&mut stream) << 8)
            | next(&mut stream)
    } else if (b1 & 4) == 0 {
        (i64::from(b1 & 7) << 40)
            | (next(&mut stream) << 32)
            | (next(&mut stream) << 24)
            | (next(&mut stream) << 16)
            | (next(&mut stream) << 8)
            | next(&mut stream)
    } else if (b1 & 2) == 0 {
        (i64::from(b1 & 3) << 48)
            | (next(&mut stream) << 40)
            | (next(&mut stream) << 32)
            | (next(&mut stream) << 24)
            | (next(&mut stream) << 16)
            | (next(&mut stream) << 8)
            | next(&mut stream)
    } else if (b1 & 1) == 0 {
        (next(&mut stream) << 48)
            | (next(&mut stream) << 40)
            | (next(&mut stream) << 32)
            | (next(&mut stream) << 24)
            | (next(&mut stream) << 16)
            | (next(&mut stream) << 8)
            | next(&mut stream)
    } else {
        (next(&mut stream) << 56)
            | (next(&mut stream) << 48)
            | (next(&mut stream) << 40)
            | (next(&mut stream) << 32)
            | (next(&mut stream) << 24)
            | (next(&mut stream) << 16)
            | (next(&mut stream) << 8)
            | next(&mut stream)
    };
    Ok((value, stream.at))
}

/// `LTF8.writeUnsignedLTF8`, returning the bytes and the bit count.
pub fn write_unsigned_ltf8(value: i64) -> (Vec<u8>, i32) {
    let unsigned = value as u64;
    let byte = |shift: i32| ((value >> shift) & 0xFF) as u8;
    if (unsigned >> 7) == 0 {
        (vec![value as u8], 8)
    } else if (unsigned >> 14) == 0 {
        (vec![((value >> 8) | 0x80) as u8, byte(0)], 16)
    } else if (unsigned >> 21) == 0 {
        (vec![((value >> 16) | 0xC0) as u8, byte(8), byte(0)], 24)
    } else if (unsigned >> 28) == 0 {
        (
            vec![((value >> 24) | 0xE0) as u8, byte(16), byte(8), byte(0)],
            32,
        )
    } else if (unsigned >> 35) == 0 {
        (
            vec![
                ((value >> 32) | 0xF0) as u8,
                byte(24),
                byte(16),
                byte(8),
                byte(0),
            ],
            40,
        )
    } else if (unsigned >> 42) == 0 {
        (
            vec![
                ((value >> 40) | 0xF8) as u8,
                byte(32),
                byte(24),
                byte(16),
                byte(8),
                byte(0),
            ],
            48,
        )
    } else if (unsigned >> 49) == 0 {
        (
            vec![
                ((value >> 48) | 0xFC) as u8,
                byte(40),
                byte(32),
                byte(24),
                byte(16),
                byte(8),
                byte(0),
            ],
            56,
        )
    } else if (unsigned >> 56) == 0 {
        (
            vec![
                0xFE,
                byte(48),
                byte(40),
                byte(32),
                byte(24),
                byte(16),
                byte(8),
                byte(0),
            ],
            64,
        )
    } else {
        (
            vec![
                0xFF,
                byte(56),
                byte(48),
                byte(40),
                byte(32),
                byte(24),
                byte(16),
                byte(8),
                byte(0),
            ],
            72,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The redundancy: two streams that differ in one byte and read identically.
    #[test]
    fn the_high_nibble_of_the_fifth_byte_is_discarded() {
        let agreeing = read_unsigned_itf8(&[0xF0, 0x00, 0x00, 0x01, 0x12]).unwrap();
        let disagreeing = read_unsigned_itf8(&[0xF0, 0x00, 0x00, 0x01, 0xF2]).unwrap();
        assert_eq!(agreeing, (18, 5));
        assert_eq!(disagreeing, (18, 5));
    }

    /// A truncated stream is a number and not a refusal.
    #[test]
    fn truncation_reads_minus_one_rather_than_failing() {
        assert_eq!(read_unsigned_itf8(&[0x80]).unwrap().0, -1);
        assert_eq!(read_unsigned_itf8(&[0xF0, 0x00]).unwrap().0, -1);
        assert_eq!(read_unsigned_itf8(&[]), Err(RuntimeEof));
    }

    /// Every value round-trips, negative ones included, which is what makes "unsigned" a misnomer.
    #[test]
    fn every_int_round_trips_including_the_negative_ones() {
        for value in [0, 1, 127, 128, 16383, 16384, i32::MAX, -1, -2, i32::MIN] {
            let (bytes, _) = write_unsigned_itf8(value);
            assert_eq!(read_unsigned_itf8(&bytes).unwrap().0, value, "{value}");
        }
    }

    #[test]
    fn every_long_round_trips_too() {
        for value in [0, 1, 127, 128, 268435456, i64::MAX, -1, i64::MIN] {
            let (bytes, _) = write_unsigned_ltf8(value);
            assert_eq!(read_unsigned_ltf8(&bytes).unwrap().0, value, "{value}");
        }
    }
}
