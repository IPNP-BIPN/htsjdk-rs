/*
 * The read side of the BAM tag codec: what comes back out, and what does not come back the same.
 *
 * The write side is pinned (BinaryTagCodec.writeTag, the integer promotion ladder, the packed tag
 * order). `readTags` is its inverse only in the loose sense: it is the function the CRAM slice
 * header's tag section is decoded with, and it is where a value's on-disk TYPE stops existing.
 *
 * Six things are decisions rather than layout.
 *
 *   - EVERY NARROW INTEGER WIDENS TO ONE JAVA TYPE. c, C, s, S and i all come back as Integer, so
 *     the width the file chose is gone the moment it is read. Only I above Integer.MAX_VALUE comes
 *     back as a Long. A port that keeps the on-disk width in its value type has more information
 *     than htsjdk has, and will not re-narrow the same way;
 *   - THE ROUND TRIP IS NOT THE IDENTITY, AND EXACTLY TWO TYPES BREAK IT. An 'I' holding a small
 *     value reads as an Integer and rewrites as 'c'; an 'H' reads as byte[] and rewrites as a 'B'
 *     array. htsjdk never writes either form, so only a foreign file can produce them, and reading
 *     one and writing it back changes the bytes;
 *   - 'A' IS A SIGNED BYTE CAST TO A CHAR. `(char) byteBuffer.get()` sign-extends, so the byte 0xE9
 *     becomes U+FFE9 rather than U+00E9. The bytes survive the rewrite because the write truncates
 *     back, but the in-memory character is not the one in the file;
 *   - THE 'I' RANGE CHECK CANNOT FIRE. `getInt() & 0xffffffffL` is already inside [0, 2^32-1], which
 *     is exactly what isValidUnsignedIntegerAttribute accepts, so the validation branch under it is
 *     unreachable for every possible input;
 *   - THE LIST COMES BACK SORTED, AND A DUPLICATE TAG REPLACES RATHER THAN DUPLICATES. Sorting is by
 *     the packed short, which weights the SECOND character, so ZA precedes AZ. Of two entries with
 *     the same tag, the LAST one in the file is the one that survives;
 *   - THE UNSIGNED FLAG OF AN ARRAY IS THE CASE OF ITS TYPE LETTER, and it is carried by the class
 *     of the returned node rather than by its value.
 *
 * Output:
 *
 *     sizes\t<fixed tag size>\t<fixed binary array tag size>
 *     value\t<label>\t<hex in>\t<tag>\t<class>\t<unsigned array>\t<value>\t<hex out>\t<stable>
 *     order\t<label>\t<hex in>\t<tags in the order they come back>
 *     dup\t<label>\t<hex in>\t<entries>\t<value that survived>
 *     readerr\t<label>\t<hex in>\t<class>\t<message>
 *     strtag\t<name>\t<packed short>\t<what makeStringTag gives back>
 *
 * Usage: BinaryTagCodecReadDump
 */

import htsjdk.samtools.BinaryTagCodec;
import htsjdk.samtools.SAMBinaryTagAndValue;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.TextTagCodec;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.BinaryCodec;

import java.io.ByteArrayOutputStream;
import java.util.Map;
import java.util.StringJoiner;

public class BinaryTagCodecReadDump {

    public static void main(final String[] args) {
        System.out.println("# BinaryTagCodecReadDump: what comes back out of a BAM tag block");

        System.out.printf("sizes\t%d\t%d%n", 3, 5);

        // Every scalar type, including the two htsjdk reads and never writes.
        value("char", new Buf().tag("CA").u8('A').u8('Q'));
        value("char-high", new Buf().tag("CB").u8('A').u8(0xE9));
        value("c", new Buf().tag("Ic").u8('c').u8(-5));
        value("c-min", new Buf().tag("Id").u8('c').u8(-128));
        value("C", new Buf().tag("IC").u8('C').u8(200));
        value("C-max", new Buf().tag("ID").u8('C').u8(255));
        value("s", new Buf().tag("Is").u8('s').i16(300));
        value("s-neg", new Buf().tag("It").u8('s').i16(-300));
        value("S", new Buf().tag("IS").u8('S').i16(40000));
        value("S-max", new Buf().tag("IT").u8('S').i16(65535));
        value("i", new Buf().tag("Ii").u8('i').i32(70000));
        value("i-neg", new Buf().tag("Ij").u8('i').i32(-70000));
        value("i-min", new Buf().tag("Ik").u8('i').i32(Integer.MIN_VALUE));
        value("I-small", new Buf().tag("IA").u8('I').i32(5));
        value("I-int-max", new Buf().tag("IB").u8('I').i32(Integer.MAX_VALUE));
        value("I-boundary", new Buf().tag("IE").u8('I').i32(Integer.MIN_VALUE));
        value("I-max", new Buf().tag("IF").u8('I').i32(-1));
        value("f", new Buf().tag("Ff").u8('f').f32(1.5f));
        value("f-negzero", new Buf().tag("Fg").u8('f').f32(-0.0f));
        value("Z", new Buf().tag("Zz").u8('Z').str("hello"));
        value("Z-empty", new Buf().tag("Zy").u8('Z').str(""));
        value("Z-high", new Buf().tag("Zx").u8('Z').u8(0xE9).u8(0));
        value("H", new Buf().tag("Hh").u8('H').str("48656C"));
        value("H-empty", new Buf().tag("Hi").u8('H').str(""));
        value("H-odd", new Buf().tag("Hj").u8('H').str("486"));
        value("H-not-hex", new Buf().tag("Hk").u8('H').str("4G"));
        value("H-lowercase", new Buf().tag("Hl").u8('H').str("ff0a"));

        // Arrays: the element type letter, its case, and the length that precedes the elements.
        value("B-c", new Buf().tag("Bc").u8('B').u8('c').i32(3).u8(-1).u8(0).u8(127));
        value("B-C", new Buf().tag("BC").u8('B').u8('C').i32(3).u8(255).u8(0).u8(127));
        value("B-s", new Buf().tag("Bs").u8('B').u8('s').i32(2).i16(-300).i16(300));
        value("B-S", new Buf().tag("BS").u8('B').u8('S').i32(2).i16(65535).i16(1));
        value("B-i", new Buf().tag("Bi").u8('B').u8('i').i32(2).i32(-70000).i32(70000));
        value("B-I", new Buf().tag("BI").u8('B').u8('I').i32(1).i32(-1));
        value("B-f", new Buf().tag("Bf").u8('B').u8('f').i32(2).f32(1.5f).f32(-2.5f));
        value("B-empty", new Buf().tag("Be").u8('B').u8('c').i32(0));

        // The order the list comes back in, from bytes deliberately not in that order.
        order("second-character-first", new Buf()
                .tag("ZA").u8('c').u8(1)
                .tag("AZ").u8('c').u8(2)
                .tag("NM").u8('c').u8(3)
                .tag("MD").u8('c').u8(4));
        order("already-sorted", new Buf()
                .tag("AZ").u8('c').u8(1)
                .tag("MD").u8('c').u8(2)
                .tag("NM").u8('c').u8(3)
                .tag("ZA").u8('c').u8(4));

        // Two entries with the same tag: how many come back, and which value.
        dup("duplicate-tag", new Buf()
                .tag("NM").u8('c').u8(1)
                .tag("NM").u8('c').u8(2));
        dup("duplicate-tag-three", new Buf()
                .tag("NM").u8('c').u8(1)
                .tag("NM").u8('c').u8(2)
                .tag("NM").u8('c').u8(3));
        dup("duplicate-around-another", new Buf()
                .tag("NM").u8('c').u8(1)
                .tag("AZ").u8('c').u8(9)
                .tag("NM").u8('c').u8(2));

        // What a malformed tag block is refused with.
        readErr("unknown-type", new Buf().tag("XX").u8('q').u8(1));
        readErr("unknown-array-type", new Buf().tag("XX").u8('B').u8('q').i32(1).u8(1));
        readErr("unterminated-string", new Buf().tag("XX").u8('Z').u8('a').u8('b'));
        readErr("truncated-after-name", new Buf().tag("XX"));
        readErr("truncated-value", new Buf().tag("XX").u8('i').u8(1).u8(2));
        readErr("negative-array-length", new Buf().tag("XX").u8('B').u8('c').i32(-1));
        readErr("array-longer-than-block", new Buf().tag("XX").u8('B').u8('c').i32(99).u8(1));
        readErr("empty-block", new Buf());

        // There is no in-memory H, so the text codec cannot write one either: the branch that
        // would is dead, and its own comment says so. A text tag read and written back proves it
        // without reference to CRAM or to BAM.
        textRoundTrip("XX:H:48656C");
        textRoundTrip("XX:H:");
        textRoundTrip("XX:B:c,1,2");
        textRoundTrip("XX:Z:hello");
        textRoundTrip("XX:i:5");

        // makeStringTag is a lookup into an array indexed by the packed short.
        strTag("NM");
        strTag("AZ");
        strTag("ZA");
        strTag("\u0000\u0001");
        strTag("A\u0080");
        strTag("A\u00FF");
        strTag("N");
        strTag("NMD");
    }

    /** One tag, decoded, then written back with the codec that wrote the file in the first place. */
    static void value(final String label, final Buf buf) {
        final byte[] in = buf.done();
        final SAMBinaryTagAndValue head;
        try {
            head = BinaryTagCodec.readTags(in, 0, in.length, ValidationStringency.SILENT);
        } catch (final Throwable t) {
            System.out.printf("readerr\t%s\t%s\t%s\t%s%n", label, hex(in),
                    t.getClass().getSimpleName(), escape(String.valueOf(t.getMessage())));
            return;
        }
        if (head == null) {
            System.out.printf("readerr\t%s\t%s\tnull\tno tag came back%n", label, hex(in));
            return;
        }

        String out;
        try {
            final ByteArrayOutputStream sink = new ByteArrayOutputStream();
            final BinaryCodec codec = new BinaryCodec(sink);
            new BinaryTagCodec(codec).writeTag(head.tag, head.value, head.isUnsignedArray());
            codec.close();
            out = hex(sink.toByteArray());
        } catch (final Throwable t) {
            out = t.getClass().getSimpleName();
        }

        System.out.printf("value\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%n", label, hex(in),
                SAMTag.makeStringTag(head.tag), head.value.getClass().getSimpleName(),
                head.isUnsignedArray(), describe(head.value), out,
                out.equals(hex(in)) ? "same" : "changed");
    }

    /** The order a whole block comes back in. */
    static void order(final String label, final Buf buf) {
        final byte[] in = buf.done();
        final StringJoiner tags = new StringJoiner(",");
        for (SAMBinaryTagAndValue tag = BinaryTagCodec.readTags(in, 0, in.length,
                ValidationStringency.SILENT); tag != null; tag = tag.getNext()) {
            tags.add(SAMTag.makeStringTag(tag.tag) + "=" + describe(tag.value));
        }
        System.out.printf("order\t%s\t%s\t%s%n", label, hex(in), tags.toString());
    }

    /** How many entries a repeated tag leaves behind, and which value survives. */
    static void dup(final String label, final Buf buf) {
        final byte[] in = buf.done();
        int entries = 0;
        String survivor = "-";
        for (SAMBinaryTagAndValue tag = BinaryTagCodec.readTags(in, 0, in.length,
                ValidationStringency.SILENT); tag != null; tag = tag.getNext()) {
            entries++;
            if ("NM".equals(SAMTag.makeStringTag(tag.tag))) {
                survivor = describe(tag.value);
            }
        }
        System.out.printf("dup\t%s\t%s\t%d\t%s%n", label, hex(in), entries, survivor);
    }

    /** What a malformed block is refused with, or what it silently gives back. */
    static void readErr(final String label, final Buf buf) {
        final byte[] in = buf.done();
        try {
            final SAMBinaryTagAndValue head = BinaryTagCodec.readTags(in, 0, in.length,
                    ValidationStringency.SILENT);
            int entries = 0;
            for (SAMBinaryTagAndValue tag = head; tag != null; tag = tag.getNext()) {
                entries++;
            }
            System.out.printf("readerr\t%s\t%s\tnone\t%d entries, first %s%n", label, hex(in),
                    entries, head == null ? "-" : describe(head.value));
        } catch (final Throwable t) {
            System.out.printf("readerr\t%s\t%s\t%s\t%s%n", label, hex(in),
                    t.getClass().getSimpleName(), escape(String.valueOf(t.getMessage())));
        }
    }

    /** One SAM text tag, decoded and encoded again by the codec that wrote it. */
    static void textRoundTrip(final String text) {
        try {
            final Map.Entry<String, Object> entry = new TextTagCodec().decode(text);
            final String back = new TextTagCodec().encode(entry.getKey(), entry.getValue());
            System.out.printf("textrt\t%s\t%s\t%s\t%s%n", text, entry.getValue().getClass()
                    .getSimpleName(), back, text.equals(back) ? "same" : "changed");
        } catch (final Throwable t) {
            System.out.printf("textrt\t%s\t%s\t%s\t-%n", text, t.getClass().getSimpleName(),
                    escape(String.valueOf(t.getMessage())));
        }
    }

    /** The packed short, and what makeStringTag does with it. */
    static void strTag(final String name) {
        final short packed;
        try {
            packed = SAMTag.makeBinaryTag(name);
        } catch (final Throwable t) {
            System.out.printf("strtag\t%s\t-\t%s: %s%n", escape(name), t.getClass().getSimpleName(),
                    escape(String.valueOf(t.getMessage())));
            return;
        }
        String back;
        try {
            back = escape(SAMTag.makeStringTag(packed));
        } catch (final Throwable t) {
            back = t.getClass().getSimpleName();
        }
        System.out.printf("strtag\t%s\t%d\t%s%n", escape(name), packed, back);
    }

    /** A value as text: arrays element by element, a character as its code point. */
    static String describe(final Object value) {
        if (value instanceof byte[]) {
            final StringJoiner j = new StringJoiner(",");
            for (final byte b : (byte[]) value) {
                j.add(Byte.toString(b));
            }
            return "[" + j + "]";
        }
        if (value instanceof short[]) {
            final StringJoiner j = new StringJoiner(",");
            for (final short s : (short[]) value) {
                j.add(Short.toString(s));
            }
            return "[" + j + "]";
        }
        if (value instanceof int[]) {
            final StringJoiner j = new StringJoiner(",");
            for (final int i : (int[]) value) {
                j.add(Integer.toString(i));
            }
            return "[" + j + "]";
        }
        if (value instanceof float[]) {
            final StringJoiner j = new StringJoiner(",");
            for (final float f : (float[]) value) {
                j.add(Float.toString(f));
            }
            return "[" + j + "]";
        }
        if (value instanceof Character) {
            return String.format("U+%04X", (int) (Character) value);
        }
        return escape(String.valueOf(value));
    }

    /** Non-printable characters as \\uXXXX, so a golden stays a text file. */
    static String escape(final String text) {
        final StringBuilder b = new StringBuilder();
        for (int i = 0; i < text.length(); i++) {
            final char c = text.charAt(i);
            if (c < 0x20 || c > 0x7E) {
                b.append(String.format("\\u%04X", (int) c));
            } else {
                b.append(c);
            }
        }
        return b.length() == 0 ? "-" : b.toString();
    }

    static String hex(final byte[] bytes) {
        final StringBuilder b = new StringBuilder(bytes.length * 2);
        for (final byte value : bytes) {
            b.append(String.format("%02x", value));
        }
        return b.length() == 0 ? "-" : b.toString();
    }

    /** A little-endian byte builder, so the dump controls exactly what the reader is handed. */
    static class Buf {
        private final ByteArrayOutputStream out = new ByteArrayOutputStream();

        Buf tag(final String name) {
            out.write(name.charAt(0));
            out.write(name.charAt(1));
            return this;
        }

        Buf u8(final int value) {
            out.write(value & 0xFF);
            return this;
        }

        Buf i16(final int value) {
            return u8(value).u8(value >> 8);
        }

        Buf i32(final int value) {
            return i16(value).i16(value >> 16);
        }

        Buf f32(final float value) {
            return i32(Float.floatToRawIntBits(value));
        }

        Buf str(final String text) {
            for (int i = 0; i < text.length(); i++) {
                u8(text.charAt(i));
            }
            return u8(0);
        }

        byte[] done() {
            return out.toByteArray();
        }
    }
}
