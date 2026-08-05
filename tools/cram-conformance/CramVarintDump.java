/*
 * ITF8 and LTF8: the variable-length integers every other CRAM structure is measured in.
 *
 * This is the floor of H.3. A container header is a run of ITF8s, a slice header is a run of ITF8s,
 * and the compression header's encoding parameters are ITF8s inside a byte array. Nothing above
 * them can be checked until they are, and they have two properties a reading of the specification
 * does not give you.
 *
 *   - THE FIVE-BYTE ITF8 STORES FOUR BITS TWICE. The writer puts (value >> 4) & 0xFF in byte four
 *     and value & 0xFF in byte five, so bits 4 to 7 appear in both. The reader takes byte four
 *     whole and masks byte five to its low nibble: `| inputStream.read() << 4 | (15 & read())`. So
 *     a hand-made stream whose two copies disagree is not an error, it resolves silently to byte
 *     four's copy, and a port that reads the fifth byte whole answers differently on exactly those
 *     streams;
 *   - THE ENCODING IS NOT SIGN-AWARE AND IS USED FOR SIGNED VALUES ANYWAY. `writeUnsignedITF8`
 *     takes an int, and a negative one has its high bits set, so it always takes the five-byte
 *     form and round-trips through a reader that returns an int. What it does NOT do is round-trip
 *     through the eight-byte forms of LTF8, where the same bit pattern is a large positive long.
 *
 * The dump reports the bytes, so a divergence names a position rather than a number, and it reports
 * what the reference reads back from streams it did not write.
 *
 * Output:
 *
 *     itf8\t<value>\t<bytes written, hex>\t<bits returned by the writer>\t<value read back>
 *     ltf8\t<value>\t<bytes written, hex>\t<bits returned by the writer>\t<value read back>
 *     itf8read\t<label>\t<bytes, hex>\t<value read or E:class:message>
 *     ltf8read\t<label>\t<bytes, hex>\t<value read or E:class:message>
 *
 * Usage: CramVarintDump
 */

import htsjdk.samtools.cram.io.ITF8;
import htsjdk.samtools.cram.io.LTF8;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;

public class CramVarintDump {

    public static void main(final String[] args) throws Exception {
        System.out.println("# CramVarintDump: the variable-length integers CRAM is measured in");

        // Every boundary of the five-way branch, from both sides, plus the values a container
        // header actually carries.
        final int[] ints = {
                0, 1, 63, 127, 128, 129, 255, 256,
                16383, 16384, 16385,
                2097151, 2097152, 2097153,
                268435455, 268435456, 268435457,
                Integer.MAX_VALUE - 1, Integer.MAX_VALUE,
                -1, -2, -128, -129, Integer.MIN_VALUE, Integer.MIN_VALUE + 1,
        };
        for (final int value : ints) {
            itf8(value);
        }

        // LTF8's branch has nine arms rather than five, so the boundaries are further apart.
        final long[] longs = {
                0L, 1L, 127L, 128L, 16383L, 16384L, 2097151L, 2097152L,
                268435455L, 268435456L,
                34359738367L, 34359738368L,
                4398046511103L, 4398046511104L,
                562949953421311L, 562949953421312L,
                72057594037927935L, 72057594037927936L,
                Long.MAX_VALUE - 1, Long.MAX_VALUE,
                -1L, -2L, Long.MIN_VALUE,
        };
        for (final long value : longs) {
            ltf8(value);
        }

        // Streams the writer would never produce. The first pair is the redundant nibble: two
        // five-byte ITF8s whose byte four and byte five disagree about bits 4 to 7.
        itf8Read("five-byte-agreeing", 0xF0, 0x00, 0x00, 0x01, 0x12);
        itf8Read("five-byte-nibble-disagrees", 0xF0, 0x00, 0x00, 0x01, 0xF2);
        itf8Read("five-byte-high-nibble-only", 0xF0, 0x00, 0x00, 0x00, 0xF0);
        // A one-byte form with the continuation bit clear but a value above 127 is impossible;
        // these are the shortest legal forms of values that also have longer spellings.
        itf8Read("one-byte-zero", 0x00);
        itf8Read("two-byte-zero", 0x80, 0x00);
        itf8Read("three-byte-zero", 0xC0, 0x00, 0x00);
        itf8Read("four-byte-zero", 0xE0, 0x00, 0x00, 0x00);
        itf8Read("five-byte-zero", 0xF0, 0x00, 0x00, 0x00, 0x00);
        // Truncated: fewer bytes than the first one promises.
        itf8Read("truncated-two", 0x80);
        itf8Read("truncated-five", 0xF0, 0x00);
        itf8Read("empty", new int[0]);

        ltf8Read("ltf8-eight-byte-zero", 0xFE, 0, 0, 0, 0, 0, 0, 0);
        ltf8Read("ltf8-nine-byte-zero", 0xFF, 0, 0, 0, 0, 0, 0, 0, 0);
        ltf8Read("ltf8-truncated", 0xFF, 0);
        ltf8Read("ltf8-empty", new int[0]);
    }

    static void itf8(final int value) throws Exception {
        final ByteArrayOutputStream out = new ByteArrayOutputStream();
        final int bits = ITF8.writeUnsignedITF8(value, out);
        final byte[] bytes = out.toByteArray();
        final int read = ITF8.readUnsignedITF8(new ByteArrayInputStream(bytes));
        System.out.printf("itf8\t%d\t%s\t%d\t%d%n", value, hex(bytes), bits, read);
    }

    static void ltf8(final long value) throws Exception {
        final ByteArrayOutputStream out = new ByteArrayOutputStream();
        final int bits = LTF8.writeUnsignedLTF8(value, out);
        final byte[] bytes = out.toByteArray();
        final long read = LTF8.readUnsignedLTF8(new ByteArrayInputStream(bytes));
        System.out.printf("ltf8\t%d\t%s\t%d\t%d%n", value, hex(bytes), bits, read);
    }

    static void itf8Read(final String label, final int... values) {
        final byte[] bytes = toBytes(values);
        String outcome;
        try {
            outcome = Integer.toString(ITF8.readUnsignedITF8(new ByteArrayInputStream(bytes)));
        } catch (final Throwable t) {
            outcome = "E:" + t.getClass().getName() + ":" + t.getMessage();
        }
        System.out.printf("itf8read\t%s\t%s\t%s%n", label, hex(bytes), outcome);
    }

    static void ltf8Read(final String label, final int... values) {
        final byte[] bytes = toBytes(values);
        String outcome;
        try {
            outcome = Long.toString(LTF8.readUnsignedLTF8(new ByteArrayInputStream(bytes)));
        } catch (final Throwable t) {
            outcome = "E:" + t.getClass().getName() + ":" + t.getMessage();
        }
        System.out.printf("ltf8read\t%s\t%s\t%s%n", label, hex(bytes), outcome);
    }

    static byte[] toBytes(final int[] values) {
        final byte[] bytes = new byte[values.length];
        for (int i = 0; i < values.length; i++) {
            bytes[i] = (byte) values[i];
        }
        return bytes;
    }

    static String hex(final byte[] bytes) {
        final StringBuilder b = new StringBuilder(bytes.length * 2);
        for (final byte value : bytes) {
            b.append(String.format("%02x", value));
        }
        return b.length() == 0 ? "-" : b.toString();
    }
}
