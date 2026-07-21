/*
 * Dumps SAMRecord.getSAMString() for records chosen to hit the text writer's choices.
 *
 * Output: <case name> <TAB> <the SAM line, newline stripped>
 *
 * The interesting property under test is that TextTagCodec collapses every integer tag width
 * to `i`, discarding exactly what BinaryTagCodec's promotion ladder computes.
 */

import htsjdk.samtools.*;

public class SamTextDump {

    static SAMFileHeader header;

    static void emit(final String name, final SAMRecord r) {
        System.out.println(name + "\t" + r.getSAMString().replace("\n", ""));
    }

    static SAMRecord base() {
        final SAMRecord r = new SAMRecord(header);
        r.setReadName("read1");
        r.setReferenceIndex(0);
        r.setAlignmentStart(100);
        r.setMappingQuality(60);
        r.setCigarString("4M");
        r.setMateReferenceIndex(0);
        r.setMateAlignmentStart(300);
        r.setInferredInsertSize(250);
        r.setReadBases("ACGT".getBytes());
        r.setBaseQualities(new byte[]{30, 31, 32, 33});
        r.setFlags(99);
        return r;
    }

    public static void main(final String[] args) {
        header = new SAMFileHeader();
        final SAMSequenceDictionary d = new SAMSequenceDictionary();
        d.addSequence(new SAMSequenceRecord("chr1", 250000000));
        d.addSequence(new SAMSequenceRecord("chr2", 200000000));
        header.setSequenceDictionary(d);

        emit("plain", base());

        final SAMRecord otherRef = base();
        otherRef.setMateReferenceIndex(1);
        emit("mate_other_reference", otherRef);

        final SAMRecord unplaced = base();
        unplaced.setReferenceIndex(-1);
        unplaced.setAlignmentStart(0);
        unplaced.setMateReferenceIndex(-1);
        unplaced.setMateAlignmentStart(0);
        unplaced.setCigarString("*");
        unplaced.setReadUnmappedFlag(true);
        emit("unplaced", unplaced);

        final SAMRecord noSeq = base();
        noSeq.setReadBases(SAMRecord.NULL_SEQUENCE);
        noSeq.setBaseQualities(SAMRecord.NULL_QUALS);
        emit("no_sequence", noSeq);

        final SAMRecord noQual = base();
        noQual.setBaseQualities(SAMRecord.NULL_QUALS);
        emit("no_quals", noQual);

        // Every integer width, which the text form collapses to `i`.
        final long[] values = {Integer.MIN_VALUE, -32769, -32768, -129, -128, -1, 0, 1, 127,
                               128, 200, 255, 256, 300, 32767, 32768, 65535, 65536,
                               2147483647L, 2147483648L, 4294967295L};
        for (final long v : values) {
            final SAMRecord r = base();
            if (v >= Integer.MIN_VALUE && v <= Integer.MAX_VALUE) r.setAttribute("XI", (int) v);
            else r.setAttribute("XI", v);
            emit("int_" + v, r);
        }

        final SAMRecord charTag = base();
        charTag.setAttribute("XA", 'Z');
        emit("tag_char", charTag);

        for (final float f : new float[]{0f, -0f, 1f, -1f, 0.5f, 3.14159f, 1e10f, 1e-10f,
                                          Float.MIN_VALUE, Float.MAX_VALUE,
                                          Float.NaN, Float.POSITIVE_INFINITY, Float.NEGATIVE_INFINITY}) {
            final SAMRecord r = base();
            r.setAttribute("XF", f);
            emit("float_" + Float.floatToRawIntBits(f), r);
        }

        for (final String s : new String[]{"", "hello", "with space", "100"}) {
            final SAMRecord r = base();
            r.setAttribute("XS", s);
            emit("str_" + s.length() + "_" + s.hashCode(), r);
        }

        final SAMRecord byteArr = base();
        byteArr.setAttribute("XB", new byte[]{-1, 0, 1, 127});
        emit("arr_byte_signed", byteArr);
        final SAMRecord ubyteArr = base();
        ubyteArr.setUnsignedArrayAttribute("XB", new byte[]{-1, 0, 1, 127});
        emit("arr_byte_unsigned", ubyteArr);
        final SAMRecord shortArr = base();
        shortArr.setAttribute("XB", new short[]{-1, 0, 300, 32767});
        emit("arr_short_signed", shortArr);
        final SAMRecord ushortArr = base();
        ushortArr.setUnsignedArrayAttribute("XB", new short[]{-1, 0, 300, 32767});
        emit("arr_short_unsigned", ushortArr);
        final SAMRecord intArr = base();
        intArr.setAttribute("XB", new int[]{-1, 0, 100000, Integer.MAX_VALUE});
        emit("arr_int_signed", intArr);
        final SAMRecord uintArr = base();
        uintArr.setUnsignedArrayAttribute("XB", new int[]{-1, 0, 100000, Integer.MAX_VALUE});
        emit("arr_int_unsigned", uintArr);
        final SAMRecord floatArr = base();
        floatArr.setAttribute("XB", new float[]{1.0f, -2.5f, 1e10f});
        emit("arr_float", floatArr);

        final SAMRecord ordered = base();
        for (final String t : new String[]{"ZA", "AZ", "NM", "MD", "AS"}) ordered.setAttribute(t, 1);
        emit("tag_order", ordered);

        for (final String c : new String[]{"4M", "2M2I", "1S2M1S", "2H4M2H", "4=", "4X",
                                           "1M1I1D1N1S1H1P1=1X"}) {
            final SAMRecord r = base();
            r.setCigarString(c);
            final int len = TextCigarCodec.decode(c).getReadLength();
            final byte[] bases = new byte[len];
            java.util.Arrays.fill(bases, (byte) 'A');
            r.setReadBases(bases);
            final byte[] q = new byte[len];
            java.util.Arrays.fill(q, (byte) 30);
            r.setBaseQualities(q);
            emit("cigar_" + c, r);
        }

        for (final int flags : new int[]{0, 4, 16, 99, 147, 1024, 2048, 4095}) {
            final SAMRecord r = base();
            r.setFlags(flags);
            emit("flags_" + flags, r);
        }
    }
}
