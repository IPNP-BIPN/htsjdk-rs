/*
 * The bit stream the core codecs are written on: the floor under the encodings.
 *
 * The record model is pinned in both directions. What carries it is a set of codecs the encoding
 * map names per data series, and three of them are written on a BIT stream rather than a byte one.
 * That stream is this file's subject, the way ITF8 was the floor under the frames.
 *
 * Six things here are decisions rather than layout.
 *
 *   - BITS ARE WRITTEN MOST SIGNIFICANT FIRST, into a one-byte buffer held back until it fills. A
 *     value of n bits is left-aligned into that buffer, so the first bit of the stream is the top
 *     bit of the first byte;
 *   - FLUSH PADS THE PARTIAL BYTE WITH ZEROS ON THE RIGHT, because the buffer was left-aligned to
 *     begin with. The padding is therefore indistinguishable from data, and only the count of
 *     values expected says where the stream really ends;
 *   - A MULTI-BYTE WRITE SPLITS AT THE TOP, not at the bottom: `write(long, n)` writes whole bytes
 *     from the most significant end while at least eight bits remain, then the remainder. So the
 *     leftover bits of a 12-bit write are the LOW four, written last;
 *   - THE READER BUFFERS A WHOLE BYTE AND COUNTS DOWN, and `readBits` assembles across byte
 *     boundaries by shifting what it has left and reading another byte;
 *   - END OF STREAM IS AN EXCEPTION, NOT A ZERO: RuntimeEOFException, and there is no way to ask
 *     the stream how many bits are left;
 *   - THE BOUNDS ARE CHECKED WITH DIFFERENT WORDING PER OVERLOAD, and one of them is not checked
 *     at all: `write(byte, 0)` against a non-empty buffer indexes a mask table of eight entries
 *     with 8.
 *
 * Output:
 *
 *     write\t<label>\t<what was written>\t<bytes after flush, hex>\t<bits buffered before flush>
 *     read\t<label>\t<input hex>\t<what was read>\t<values>
 *     err\t<label>\t<what it was given>\t<class>\t<message>
 *
 * Usage: CramBitStreamDump
 */

import htsjdk.samtools.cram.io.DefaultBitInputStream;
import htsjdk.samtools.cram.io.DefaultBitOutputStream;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.util.StringJoiner;

public class CramBitStreamDump {

    public static void main(final String[] args) {
        System.out.println("# CramBitStreamDump: the bit stream the core codecs are written on");

        // One value of n bits, left-aligned into the first byte.
        for (int bits = 1; bits <= 8; bits++) {
            writeByte("one-value-" + bits + "-bits", (byte) 1, bits);
        }
        writeByte("all-ones-3-bits", (byte) 0x07, 3);
        writeByte("high-bit-of-a-byte", (byte) 0x80, 8);
        writeByte("a-byte-in-5-bits", (byte) 0x1F, 5);
        // The value is masked to the requested width, so the bits above it are dropped.
        writeByte("a-byte-too-wide-for-its-width", (byte) 0xFF, 3);

        // Two writes that straddle the buffered byte.
        writeTwo("straddle-4-and-4", (byte) 0x0A, 4, (byte) 0x05, 4);
        writeTwo("straddle-5-and-5", (byte) 0x15, 5, (byte) 0x0A, 5);
        writeTwo("straddle-7-and-7", (byte) 0x55, 7, (byte) 0x2A, 7);
        writeTwo("straddle-1-and-8", (byte) 1, 1, (byte) 0xFF, 8);
        writeTwo("straddle-8-and-1", (byte) 0xFF, 8, (byte) 1, 1);

        // The multi-byte writes, which split at the top.
        writeLong("long-12-bits", 0xABCL, 12);
        writeLong("long-16-bits", 0xABCDL, 16);
        writeLong("long-24-bits", 0xABCDEFL, 24);
        writeLong("long-32-bits", 0xDEADBEEFL, 32);
        writeLong("long-64-bits", 0x0123456789ABCDEFL, 64);
        writeLong("long-9-bits", 0x1FFL, 9);
        writeLong("long-1-bit", 1L, 1);
        writeInt("int-12-bits", 0xABC, 12);
        writeInt("int-32-bits", 0xDEADBEEF, 32);

        // Booleans, which are one-bit writes.
        writeBits("seven-true", true, 7);
        writeBits("eight-true", true, 8);
        writeBits("nine-true", true, 9);
        writeBits("eight-false", false, 8);

        // Nothing at all: no bits buffered, no bytes out.
        writeNothing("nothing");

        // Reading back, from bytes the dump lays out itself.
        readBits("read-one-bit-at-a-time", new byte[] {(byte) 0xB5}, new int[] {1, 1, 1, 1, 1, 1,
                1, 1});
        readBits("read-4-and-4", new byte[] {(byte) 0xAB}, new int[] {4, 4});
        readBits("read-across-a-byte", new byte[] {(byte) 0xAB, (byte) 0xCD}, new int[] {5, 5, 6});
        readBits("read-12-bits", new byte[] {(byte) 0xAB, (byte) 0xCD}, new int[] {12});
        readBits("read-16-bits", new byte[] {(byte) 0xAB, (byte) 0xCD}, new int[] {16});
        readBits("read-zero-bits", new byte[] {(byte) 0xAB}, new int[] {0, 8});
        readLongBits("read-long-32", new byte[] {(byte) 0xDE, (byte) 0xAD, (byte) 0xBE,
                (byte) 0xEF}, new int[] {32});
        readLongBits("read-long-64", new byte[] {0x01, 0x23, 0x45, 0x67, (byte) 0x89, (byte) 0xAB,
                (byte) 0xCD, (byte) 0xEF}, new int[] {64});
        readLongBits("read-long-across", new byte[] {(byte) 0xAB, (byte) 0xCD}, new int[] {5, 11});
        readLongBits("read-long-zero", new byte[] {(byte) 0xAB}, new int[] {0, 8});

        // The padding a flush leaves is indistinguishable from data.
        readBits("read-the-padding", new byte[] {(byte) 0x80}, new int[] {1, 7});

        // End of stream, which is an exception rather than a zero.
        readBits("read-past-the-end", new byte[] {(byte) 0xAB}, new int[] {8, 1});
        readBits("read-past-the-end-wide", new byte[] {(byte) 0xAB}, new int[] {9});
        readLongBits("read-long-past-the-end", new byte[] {(byte) 0xAB}, new int[] {9});
        readBits("read-from-nothing", new byte[0], new int[] {1});

        // What each overload refuses, and with which wording.
        errLong("long-65-bits", 1L, 65);
        errLong("long-negative-bits", 1L, -1);
        errInt("int-33-bits", 1, 33);
        errByte("byte-9-bits", (byte) 1, 9);
        errByte("byte-negative-bits", (byte) 1, -1);
        errLongRead("read-long-65-bits", new byte[] {(byte) 0xAB}, 65);

        // A zero-bit write is a no-op against an empty buffer and an index out of a table of
        // eight against a full one.
        writeByte("zero-bits-into-an-empty-buffer", (byte) 1, 0);
        errZeroBitsAfter("zero-bits-into-a-partial-buffer");
    }

    /** One `write(byte, nofBits)`. */
    static void writeByte(final String label, final byte value, final int bits) {
        final ByteArrayOutputStream sink = new ByteArrayOutputStream();
        try (final DefaultBitOutputStream out = new DefaultBitOutputStream(sink)) {
            out.write(value, bits);
        } catch (final Throwable t) {
            err(label, String.format("byte 0x%02X in %d bits", value & 0xFF, bits), t);
            return;
        }
        System.out.printf("write\t%s\tbyte 0x%02X in %d bits\t%s\t%d%n", label, value & 0xFF, bits,
                hex(sink.toByteArray()), bits % 8);
    }

    /** Two writes in a row, which is where the buffer straddles. */
    static void writeTwo(final String label, final byte first, final int firstBits,
            final byte second, final int secondBits) {
        final ByteArrayOutputStream sink = new ByteArrayOutputStream();
        try (final DefaultBitOutputStream out = new DefaultBitOutputStream(sink)) {
            out.write(first, firstBits);
            out.write(second, secondBits);
        } catch (final Throwable t) {
            err(label, String.format("byte 0x%02X in %d then 0x%02X in %d", first & 0xFF,
                    firstBits, second & 0xFF, secondBits), t);
            return;
        }
        System.out.printf("write\t%s\tbyte 0x%02X in %d then 0x%02X in %d\t%s\t%d%n", label,
                first & 0xFF, firstBits, second & 0xFF, secondBits, hex(sink.toByteArray()),
                (firstBits + secondBits) % 8);
    }

    static void writeLong(final String label, final long value, final int bits) {
        final ByteArrayOutputStream sink = new ByteArrayOutputStream();
        try (final DefaultBitOutputStream out = new DefaultBitOutputStream(sink)) {
            out.write(value, bits);
        } catch (final Throwable t) {
            err(label, String.format("long 0x%X in %d bits", value, bits), t);
            return;
        }
        System.out.printf("write\t%s\tlong 0x%X in %d bits\t%s\t%d%n", label, value, bits,
                hex(sink.toByteArray()), bits % 8);
    }

    static void writeInt(final String label, final int value, final int bits) {
        final ByteArrayOutputStream sink = new ByteArrayOutputStream();
        try (final DefaultBitOutputStream out = new DefaultBitOutputStream(sink)) {
            out.write(value, bits);
        } catch (final Throwable t) {
            err(label, String.format("int 0x%X in %d bits", value, bits), t);
            return;
        }
        System.out.printf("write\t%s\tint 0x%X in %d bits\t%s\t%d%n", label, value, bits,
                hex(sink.toByteArray()), bits % 8);
    }

    static void writeBits(final String label, final boolean bit, final int repeat) {
        final ByteArrayOutputStream sink = new ByteArrayOutputStream();
        try (final DefaultBitOutputStream out = new DefaultBitOutputStream(sink)) {
            out.write(bit, repeat);
        } catch (final Throwable t) {
            err(label, String.format("%b repeated %d", bit, repeat), t);
            return;
        }
        System.out.printf("write\t%s\t%b repeated %d\t%s\t%d%n", label, bit, repeat,
                hex(sink.toByteArray()), repeat % 8);
    }

    static void writeNothing(final String label) {
        final ByteArrayOutputStream sink = new ByteArrayOutputStream();
        try (final DefaultBitOutputStream out = new DefaultBitOutputStream(sink)) {
            // nothing at all
            out.flush();
        }
        System.out.printf("write\t%s\tnothing\t%s\t0%n", label, hex(sink.toByteArray()));
    }

    /** A run of `readBits` over bytes the dump supplies. */
    static void readBits(final String label, final byte[] input, final int[] widths) {
        final StringJoiner joiner = new StringJoiner(",");
        final StringJoiner asked = new StringJoiner(",");
        for (final int width : widths) {
            asked.add(Integer.toString(width));
        }
        try (final DefaultBitInputStream in =
                new DefaultBitInputStream(new ByteArrayInputStream(input))) {
            for (final int width : widths) {
                joiner.add(width == 1 ? Boolean.toString(in.readBit())
                        : Integer.toString(in.readBits(width)));
            }
        } catch (final Throwable t) {
            err(label, String.format("%s readBits %s", hex(input), asked.toString()), t);
            return;
        }
        System.out.printf("read\t%s\t%s\treadBits %s\t%s%n", label, hex(input), asked.toString(),
                joiner.toString());
    }

    static void readLongBits(final String label, final byte[] input, final int[] widths) {
        final StringJoiner joiner = new StringJoiner(",");
        final StringJoiner asked = new StringJoiner(",");
        for (final int width : widths) {
            asked.add(Integer.toString(width));
        }
        try (final DefaultBitInputStream in =
                new DefaultBitInputStream(new ByteArrayInputStream(input))) {
            for (final int width : widths) {
                joiner.add(Long.toString(in.readLongBits(width)));
            }
        } catch (final Throwable t) {
            err(label, String.format("%s readLongBits %s", hex(input), asked.toString()), t);
            return;
        }
        System.out.printf("read\t%s\t%s\treadLongBits %s\t%s%n", label, hex(input),
                asked.toString(), joiner.toString());
    }

    static void errLong(final String label, final long value, final int bits) {
        final ByteArrayOutputStream sink = new ByteArrayOutputStream();
        try (final DefaultBitOutputStream out = new DefaultBitOutputStream(sink)) {
            out.write(value, bits);
            System.out.printf("write\t%s\tlong 0x%X in %d bits\t%s\t%d%n", label, value, bits,
                    hex(sink.toByteArray()), 0);
        } catch (final Throwable t) {
            err(label, String.format("long 0x%X in %d bits", value, bits), t);
        }
    }

    static void errInt(final String label, final int value, final int bits) {
        final ByteArrayOutputStream sink = new ByteArrayOutputStream();
        try (final DefaultBitOutputStream out = new DefaultBitOutputStream(sink)) {
            out.write(value, bits);
            System.out.printf("write\t%s\tint 0x%X in %d bits\t%s\t%d%n", label, value, bits,
                    hex(sink.toByteArray()), 0);
        } catch (final Throwable t) {
            err(label, String.format("int 0x%X in %d bits", value, bits), t);
        }
    }

    static void errByte(final String label, final byte value, final int bits) {
        final ByteArrayOutputStream sink = new ByteArrayOutputStream();
        try (final DefaultBitOutputStream out = new DefaultBitOutputStream(sink)) {
            out.write(value, bits);
            System.out.printf("write\t%s\tbyte 0x%02X in %d bits\t%s\t%d%n", label, value & 0xFF,
                    bits, hex(sink.toByteArray()), 0);
        } catch (final Throwable t) {
            err(label, String.format("byte 0x%02X in %d bits", value & 0xFF, bits), t);
        }
    }

    static void errLongRead(final String label, final byte[] input, final int bits) {
        try (final DefaultBitInputStream in =
                new DefaultBitInputStream(new ByteArrayInputStream(input))) {
            final long value = in.readLongBits(bits);
            System.out.printf("read\t%s\t%s\treadLongBits %d\t%d%n", label, hex(input), bits,
                    value);
        } catch (final Throwable t) {
            err(label, String.format("%s readLongBits %d", hex(input), bits), t);
        }
    }

    /** Four bits, then a write of zero bits against the buffer they left behind. */
    static void errZeroBitsAfter(final String label) {
        final ByteArrayOutputStream sink = new ByteArrayOutputStream();
        try (final DefaultBitOutputStream out = new DefaultBitOutputStream(sink)) {
            out.write((byte) 0x0A, 4);
            out.write((byte) 1, 0);
            System.out.printf("write\t%s\t4 bits then 0 bits\t%s\t4%n", label,
                    hex(sink.toByteArray()));
        } catch (final Throwable t) {
            err(label, "byte 0x0A in 4 bits then byte 0x01 in 0 bits", t);
        }
    }

    /** Every refusal carries what it was given, so no test has to rebuild an input from a label. */
    static void err(final String label, final String input, final Throwable t) {
        System.out.printf("err\t%s\t%s\t%s\t%s%n", label, input, t.getClass().getSimpleName(),
                String.valueOf(t.getMessage()));
    }

    static String hex(final byte[] bytes) {
        final StringBuilder b = new StringBuilder(bytes.length * 2);
        for (final byte value : bytes) {
            b.append(String.format("%02x", value));
        }
        return b.length() == 0 ? "-" : b.toString();
    }
}
