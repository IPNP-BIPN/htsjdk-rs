/*
 * Dumps the exact bytes htsjdk's BAMRecordCodec produces for a set of records chosen to hit
 * the places where the format leaves a choice open.
 *
 * Run inside the pinned oracle container. Output is one line per case:
 *
 *     <case name> <TAB> <hex of the encoded record, including its block_size field>
 *
 * The Rust side builds the same records and must reproduce the hex exactly. Nothing here is
 * asserted by this program: it is a recorder, and the comparison happens on the Rust side so
 * that a bug in the harness cannot quietly define the expected answer.
 */

import htsjdk.samtools.*;
import java.io.ByteArrayOutputStream;
import java.util.ArrayList;
import java.util.List;

public class BamCodecDump {

    static SAMFileHeader header;
    static BAMRecordCodec codec;
    static ByteArrayOutputStream sink;

    static void reset() {
        sink = new ByteArrayOutputStream();
        codec = new BAMRecordCodec(header);
        codec.setOutputStream(sink);
    }

    static void emit(final String name, final SAMRecord rec) {
        reset();
        codec.encode(rec);
        final byte[] bytes = sink.toByteArray();
        final StringBuilder sb = new StringBuilder();
        for (final byte b : bytes) sb.append(String.format("%02x", b));
        System.out.println(name + "\t" + sb);
    }

    static SAMRecord base() {
        final SAMRecord r = new SAMRecord(header);
        r.setReadName("read1");
        r.setReferenceIndex(0);
        r.setAlignmentStart(100);
        r.setMappingQuality(60);
        r.setCigarString("4M");
        r.setMateReferenceIndex(-1);
        r.setMateAlignmentStart(0);
        r.setInferredInsertSize(0);
        r.setReadBases("ACGT".getBytes());
        r.setBaseQualities(new byte[]{30, 30, 30, 30});
        r.setFlags(0);
        return r;
    }

    public static void main(final String[] args) {
        header = new SAMFileHeader();
        final SAMSequenceDictionary dict = new SAMSequenceDictionary();
        dict.addSequence(new SAMSequenceRecord("chr1", 250000000));
        dict.addSequence(new SAMSequenceRecord("chr2", 200000000));
        header.setSequenceDictionary(dict);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);

        emit("plain", base());

        // ---- alignment positions, which drive the indexing bin --------------------------
        final int[] starts = {1, 2, 16384, 16385, 16386, 131072, 131073, 1048576,
                              8388608, 67108864, 100000000, 250000000};
        for (final int start : starts) {
            final SAMRecord r = base();
            r.setAlignmentStart(start);
            emit("start_" + start, r);
        }
        // Alignments that straddle bin boundaries at every level. The reference span is made
        // with a deletion rather than with matches, so a 64 Mb span costs four read bases
        // instead of 64 million. Same bin arithmetic, four orders of magnitude less golden.
        final int[] spans = {1, 100, 16384, 131072, 1048576, 8388608, 67108864};
        for (final int span : spans) {
            final SAMRecord r = base();
            r.setAlignmentStart(16300);
            r.setCigarString("4M" + span + "D");
            emit("span_" + span, r);
        }

        // ---- unmapped and unplaced -------------------------------------------------------
        final SAMRecord unmapped = base();
        unmapped.setReadUnmappedFlag(true);
        unmapped.setReferenceIndex(-1);
        unmapped.setAlignmentStart(0);
        unmapped.setCigarString("*");
        emit("unmapped", unmapped);

        final SAMRecord placedUnmapped = base();
        placedUnmapped.setReadUnmappedFlag(true);
        emit("placed_unmapped", placedUnmapped);

        // ---- read shapes ------------------------------------------------------------------
        final String[] seqs = {"A", "AC", "ACG", "ACGT", "ACGTN", "=ACMGRSVTWYHKDBN",
                               "acgt", "ACGTacgtNn.="};
        for (final String seq : seqs) {
            final SAMRecord r = base();
            r.setReadBases(seq.getBytes());
            r.setCigarString(seq.length() + "M");
            final byte[] quals = new byte[seq.length()];
            for (int i = 0; i < quals.length; i++) quals[i] = (byte) (i % 60);
            r.setBaseQualities(quals);
            emit("seq_" + seq, r);
        }

        // No qualities at all: htsjdk writes 0xFF repeated.
        final SAMRecord noQual = base();
        noQual.setBaseQualities(SAMRecord.NULL_QUALS);
        emit("no_quals", noQual);

        // ---- read names -------------------------------------------------------------------
        final String[] names = {"a", "read/1", "x".repeat(100), "x".repeat(254)};
        for (int i = 0; i < names.length; i++) {
            final SAMRecord r = base();
            r.setReadName(names[i]);
            emit("name_" + i + "_len" + names[i].length(), r);
        }

        // ---- cigars ------------------------------------------------------------------------
        final String[] cigars = {"4M", "2M2I", "2M2D2M", "1S2M1S", "2H4M2H", "4=", "4X",
                                 "1M1I1D1N1S1H1P1=1X"};
        for (final String c : cigars) {
            final SAMRecord r = base();
            r.setCigarString(c);
            final int readLen = new Cigar(TextCigarCodec.decode(c).getCigarElements()).getReadLength();
            final byte[] bases = new byte[readLen];
            java.util.Arrays.fill(bases, (byte) 'A');
            r.setReadBases(bases);
            final byte[] quals = new byte[readLen];
            java.util.Arrays.fill(quals, (byte) 30);
            r.setBaseQualities(quals);
            emit("cigar_" + c, r);
        }

        // ---- the integer tag promotion ladder -----------------------------------------------
        final long[] values = {
                Integer.MIN_VALUE, -32769L, -32768L, -32767L, -129L, -128L, -127L, -1L, 0L,
                1L, 126L, 127L, 128L, 129L, 200L, 254L, 255L, 256L, 257L, 300L,
                32766L, 32767L, 32768L, 32769L, 65534L, 65535L, 65536L, 65537L,
                2147483646L, 2147483647L, 2147483648L, 4294967294L, 4294967295L};
        for (final long v : values) {
            final SAMRecord r = base();
            // Feed the narrowest Java box type that holds it, to confirm the declared type has
            // no influence on the encoding.
            if (v >= Integer.MIN_VALUE && v <= Integer.MAX_VALUE) {
                r.setAttribute("XI", (int) v);
            } else {
                r.setAttribute("XI", v);
            }
            emit("int_" + v, r);
        }
        // The same value through different Java box types must encode identically.
        for (final int v : new int[]{100, 30000}) {
            final SAMRecord byteBox = base();
            byteBox.setAttribute("XI", (byte) v);
            emit("box_byte_" + v, byteBox);
            final SAMRecord shortBox = base();
            shortBox.setAttribute("XI", (short) v);
            emit("box_short_" + v, shortBox);
            final SAMRecord intBox = base();
            intBox.setAttribute("XI", v);
            emit("box_int_" + v, intBox);
            final SAMRecord longBox = base();
            longBox.setAttribute("XI", (long) v);
            emit("box_long_" + v, longBox);
        }

        // ---- other tag types -------------------------------------------------------------
        final SAMRecord charTag = base();
        charTag.setAttribute("XA", 'Z');
        emit("tag_char", charTag);

        final SAMRecord floatTag = base();
        floatTag.setAttribute("XF", 1.5f);
        emit("tag_float", floatTag);

        for (final float f : new float[]{0.0f, -0.0f, 1.0f, -1.0f, 3.14159f, Float.MIN_VALUE,
                                         Float.MAX_VALUE, Float.NaN, Float.POSITIVE_INFINITY}) {
            final SAMRecord r = base();
            r.setAttribute("XF", f);
            emit("float_" + Float.floatToRawIntBits(f), r);
        }

        for (final String s : new String[]{"", "a", "hello world", "100", "é"}) {
            final SAMRecord r = base();
            r.setAttribute("XS", s);
            emit("tag_str_" + s.length() + "_" + s.hashCode(), r);
        }

        // Arrays, signed and unsigned, at each element width.
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
        floatArr.setAttribute("XB", new float[]{1.0f, -2.5f});
        emit("arr_float", floatArr);
        final SAMRecord emptyArr = base();
        emptyArr.setAttribute("XB", new int[]{});
        emit("arr_empty", emptyArr);

        // ---- tag ordering ------------------------------------------------------------------
        // Inserted in an order that is neither sorted nor reverse-sorted, by tag names whose
        // packed-short order differs from their string order.
        final SAMRecord ordered = base();
        for (final String t : new String[]{"ZA", "AZ", "NM", "MD", "AS", "XS", "SA", "aa", "Aa"}) {
            ordered.setAttribute(t, 1);
        }
        emit("tag_order", ordered);

        // A mixed bag of realistic tags on one record.
        final SAMRecord realistic = base();
        realistic.setAttribute("NM", 2);
        realistic.setAttribute("MD", "2A1");
        realistic.setAttribute("AS", 100);
        realistic.setAttribute("XS", -30);
        realistic.setAttribute("RG", "rg1");
        realistic.setAttribute("PG", "MarkDuplicates");
        emit("realistic", realistic);

        // ---- the long-cigar displacement ------------------------------------------------
        for (final int n : new int[]{65535, 65536, 65537}) {
            final List<CigarElement> elements = new ArrayList<>(n);
            for (int i = 0; i < n; i++) {
                elements.add(new CigarElement(1, i % 2 == 0 ? CigarOperator.M : CigarOperator.I));
            }
            final Cigar c = new Cigar(elements);
            final SAMRecord r = base();
            r.setCigar(c);
            final byte[] bases = new byte[c.getReadLength()];
            java.util.Arrays.fill(bases, (byte) 'A');
            r.setReadBases(bases);
            final byte[] quals = new byte[c.getReadLength()];
            java.util.Arrays.fill(quals, (byte) 30);
            r.setBaseQualities(quals);
            emit("longcigar_" + n, r);
        }

        // The same, with tags around it, to pin where CG lands in the sorted list.
        final List<CigarElement> elements = new ArrayList<>();
        for (int i = 0; i < 65536; i++) {
            elements.add(new CigarElement(1, i % 2 == 0 ? CigarOperator.M : CigarOperator.I));
        }
        final Cigar c = new Cigar(elements);
        final SAMRecord withTags = base();
        withTags.setCigar(c);
        final byte[] bases = new byte[c.getReadLength()];
        java.util.Arrays.fill(bases, (byte) 'A');
        withTags.setReadBases(bases);
        final byte[] quals = new byte[c.getReadLength()];
        java.util.Arrays.fill(quals, (byte) 30);
        withTags.setBaseQualities(quals);
        withTags.setAttribute("AG", 1);
        withTags.setAttribute("CH", 1);
        withTags.setAttribute("NM", 1);
        emit("longcigar_with_tags", withTags);

        // ---- flags, mates, insert sizes -------------------------------------------------
        for (final int flags : new int[]{0, 1, 4, 16, 99, 147, 1024, 2048, 4095}) {
            final SAMRecord r = base();
            r.setFlags(flags);
            emit("flags_" + flags, r);
        }
        for (final int isize : new int[]{Integer.MIN_VALUE, -1000, -1, 0, 1, 1000, Integer.MAX_VALUE}) {
            final SAMRecord r = base();
            r.setInferredInsertSize(isize);
            emit("isize_" + isize, r);
        }
        final SAMRecord mated = base();
        mated.setMateReferenceIndex(1);
        mated.setMateAlignmentStart(5000);
        mated.setInferredInsertSize(300);
        emit("mated", mated);

        for (final int mq : new int[]{0, 1, 60, 254, 255}) {
            final SAMRecord r = base();
            r.setMappingQuality(mq);
            emit("mapq_" + mq, r);
        }
    }
}
