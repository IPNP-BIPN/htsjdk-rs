/*
 * The CRAM slice header: the last frame before CRAM becomes reads.
 *
 * The compression header's three maps are pinned. A container's blocks then hold one or more
 * slices, and each slice begins with a MAPPED_SLICE block whose raw content is this header: an
 * alignment context, a record count, a block count, the external content ids the slice uses, an
 * embedded reference id, an MD5, and, from version 3, a run of BAM-encoded tags.
 *
 * Six things are decisions rather than layout.
 *
 *   - THE BLOCK COUNT DOES NOT COUNT THE HEADER BLOCK. `getNumberOfBlocks` returns
 *     `1 + numberOfExternalBlocks`: the core block plus the externals, and measured, it equals
 *     exactly the number of blocks that FOLLOW the header. A reader that counts the header block
 *     among them reads one block too few and stops before the last one;
 *   - THE SLICE HEADER CARRIES SIX TAGS, AND FOUR OF THEM DIGEST NOTHING. B1 and S1 are a SHA-1,
 *     B5 and S5 a SHA-512, and on an unmapped slice all four are the digest of the EMPTY string,
 *     byte for byte identical in every file. That is 168 bytes of constant per slice. Only BD and
 *     SD, four bytes apiece, vary with the reads, and they do not move when only the tags change;
 *   - THE TAG SECTION HAS NO LENGTH. It is read with readFully to the end of the block, so the
 *     slice header block's own length is the only thing that delimits it, and a header with no tags
 *     is indistinguishable from one whose tags are zero bytes long;
 *   - THE MD5 IS SIXTEEN ZEROES WHEN THERE IS NONE, not an absent field. `createSliceHeaderBlock`
 *     writes `new byte[16]` when the reference MD5 is null, so the field is always present and its
 *     emptiness is a value;
 *   - AN ABSENT EMBEDDED REFERENCE IS -1, written as an ITF8, which is the five-byte form. The
 *     commonest value of this field is also its longest encoding;
 *   - THE ALIGNMENT CONTEXT CARRIES THREE MAGIC NUMBERS: a reference id of -1 is unmapped-unplaced
 *     and -2 is multiple-reference, and both force the start and span to 0;
 *   - TAGS ARE VERSION-GATED ON BOTH SIDES. Below major version 3 they are neither written nor
 *     read, so the same slice is shorter in a 2.1 file by however many bytes its tags occupied.
 *
 * Output:
 *
 *     sizes\t<md5 byte size>\t<embedded reference absent id>\t<no alignment start>\t<no alignment span>
 *     slice\t<label>\t<index>\t<refContextId>\t<start>\t<span>\t<records>\t<counter>\t<blocks>\t<contentIds>\t<embeddedRefId>\t<md5 hex>\t<tag bytes hex>
 *     hdrbytes\t<label>\t<index>\t<the slice header block's raw content, hex>
 *     blockcount\t<label>\t<index>\t<declared block count>\t<blocks actually present after the header>
 *     counts\t<label>\t<slices>
 *     err\t<label>\t<class>\t<message>
 *
 * Usage: CramSliceHeaderDump
 */

import htsjdk.samtools.*;
import htsjdk.samtools.cram.io.CramInt;
import htsjdk.samtools.cram.io.ITF8;
import htsjdk.samtools.cram.io.LTF8;
import htsjdk.samtools.BinaryTagCodec;
import htsjdk.samtools.SAMBinaryTagAndValue;
import htsjdk.samtools.SAMTagUtil;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.cram.structure.CramHeader;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.StringJoiner;

public class CramSliceHeaderDump {

    private static final int MD5_BYTE_SIZE = 16;

    public static void main(final String[] args) throws Exception {
        System.out.println("# CramSliceHeaderDump: the last frame before CRAM becomes reads");

        System.out.printf("sizes\t%d\t%d\t%d\t%d%n", MD5_BYTE_SIZE, -1, 0, 0);

        emit("four-unmapped", build(4, 8, false));
        emit("one-unmapped", build(1, 8, false));
        emit("long-reads", build(2, 400, false));
        emit("many-reads", build(40, 20, false));
        emit("tagged", build(4, 8, true));
    }

    static byte[] build(final int count, final int length, final boolean tagged) {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        header.addSequence(new SAMSequenceRecord("chr1", 100000));

        final ByteArrayOutputStream out = new ByteArrayOutputStream();
        try (final SAMFileWriter writer =
                new SAMFileWriterFactory().makeCRAMWriter(header, out, (java.io.File) null)) {
            for (int i = 0; i < count; i++) {
                final SAMRecord record = new SAMRecord(header);
                record.setReadName("read" + i);
                record.setReadUnmappedFlag(true);
                record.setReferenceIndex(-1);
                record.setAlignmentStart(0);
                final byte[] bases = new byte[length];
                Arrays.fill(bases, (byte) "ACGT".charAt(i % 4));
                record.setReadBases(bases);
                final byte[] quals = new byte[length];
                Arrays.fill(quals, (byte) (30 + i % 5));
                record.setBaseQualities(quals);
                if (tagged) {
                    record.setAttribute("NM", i);
                    record.setAttribute("MD", "8");
                }
                writer.addAlignment(record);
            }
        }
        return out.toByteArray();
    }

    /** Walk the data container's blocks and parse every MAPPED_SLICE header among them. */
    static void emit(final String label, final byte[] cram) throws Exception {
        final byte[] blocks = dataContainerBlocks(cram);
        int at = 0;
        int index = 0;
        int slices = 0;

        while (at < blocks.length) {
            final int[] header = blockHeaderAt(blocks, at);
            final int method = header[0];
            final int contentType = header[1];
            final int headerLength = header[3];
            final int compressedSize = header[4];

            if (contentType == 2) { // MAPPED_SLICE
                if (method != 0) {
                    System.out.printf("err\t%s\tnot-raw\tthe slice header block is method %d%n",
                            label, method);
                    return;
                }
                final byte[] content = Arrays.copyOfRange(blocks, at + headerLength,
                        at + headerLength + compressedSize);
                System.out.printf("hdrbytes\t%s\t%d\t%s%n", label, index, hex(content));
                final int declaredBlocks = parse(label, index, content);

                // How many blocks actually follow this header before the next slice or the end.
                int following = 0;
                int cursor = at + headerLength + compressedSize + 4;
                while (cursor < blocks.length) {
                    final int[] next = blockHeaderAt(blocks, cursor);
                    if (next[1] == 2) {
                        break;
                    }
                    following++;
                    cursor += next[3] + next[4] + 4;
                }
                System.out.printf("blockcount\t%s\t%d\t%d\t%d%n", label, index, declaredBlocks,
                        following);
                index++;
                slices++;
            }
            at += headerLength + compressedSize + 4;
        }
        System.out.printf("counts\t%s\t%d%n", label, slices);
    }

    /** The slice header's own fields, in the order the parser reads them. */
    static int parse(final String label, final int index, final byte[] content) throws Exception {
        final ByteArrayInputStream in = new ByteArrayInputStream(content);
        final int refContextId = ITF8.readUnsignedITF8(in);
        final int alignmentStart = ITF8.readUnsignedITF8(in);
        final int alignmentSpan = ITF8.readUnsignedITF8(in);
        final int records = ITF8.readUnsignedITF8(in);
        final long counter = LTF8.readUnsignedLTF8(in);
        final int blocks = ITF8.readUnsignedITF8(in);

        final int contentIdCount = ITF8.readUnsignedITF8(in);
        final List<Integer> contentIds = new ArrayList<>();
        for (int i = 0; i < contentIdCount; i++) {
            contentIds.add(ITF8.readUnsignedITF8(in));
        }
        final int embeddedReferenceId = ITF8.readUnsignedITF8(in);

        final byte[] md5 = new byte[MD5_BYTE_SIZE];
        if (in.read(md5) != MD5_BYTE_SIZE) {
            System.out.printf("err\t%s\tshort\tthe MD5 field%n", label);
            return blocks;
        }
        // Whatever remains is the tag section, which carries no length of its own.
        final byte[] tags = new byte[in.available()];
        if (tags.length > 0 && in.read(tags) != tags.length) {
            System.out.printf("err\t%s\tshort\tthe tag section%n", label);
            return blocks;
        }

        final StringJoiner ids = new StringJoiner(",");
        for (final int id : contentIds) {
            ids.add(Integer.toString(id));
        }
        System.out.printf("slice\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%s\t%s%n", label, index,
                refContextId, alignmentStart, alignmentSpan, records, counter, blocks,
                contentIds.isEmpty() ? "-" : ids.toString(), embeddedReferenceId, hex(md5),
                hex(tags));

        // The tag section carries no length of its own, so it is whatever is left. Decoded here
        // with the same codec the reader uses, rather than guessed at from the bytes.
        if (tags.length > 0) {
            SAMBinaryTagAndValue tag = BinaryTagCodec.readTags(tags, 0, tags.length,
                    ValidationStringency.DEFAULT_STRINGENCY);
            while (tag != null) {
                final Object value = tag.value;
                System.out.printf("slicetag\t%s\t%d\t%s\t%s\t%s%n", label, index,
                        SAMTagUtil.getSingleton().makeStringTag((short) tag.tag),
                        value == null ? "<null>" : value.getClass().getSimpleName(),
                        describe(value));
                tag = tag.getNext();
            }
        }
        return blocks;
    }

    /** A tag value as text, with a byte array shown in hex. */
    static String describe(final Object value) {
        if (value instanceof byte[]) {
            return hex((byte[]) value);
        }
        return String.valueOf(value);
    }

    /** method, content type, content id, header length, compressed size. */
    static int[] blockHeaderAt(final byte[] bytes, final int at) {
        final ByteArrayInputStream in = new ByteArrayInputStream(bytes, at, bytes.length - at);
        final int method = in.read();
        final int contentType = in.read();
        final int contentId = ITF8.readUnsignedITF8(in);
        final int compressedSize = ITF8.readUnsignedITF8(in);
        ITF8.readUnsignedITF8(in);
        final int headerLength = (bytes.length - at) - in.available();
        return new int[] {method, contentType, contentId, headerLength, compressedSize};
    }

    /** The blocks of the second container, which is the first one holding records. */
    static byte[] dataContainerBlocks(final byte[] cram) {
        int at = CramHeader.CRAM_HEADER_LENGTH;
        for (int container = 0; container < 2; container++) {
            final ByteArrayInputStream in = new ByteArrayInputStream(cram, at, cram.length - at);
            final int blocksByteSize = CramInt.readInt32(in);
            ITF8.readUnsignedITF8(in);
            ITF8.readUnsignedITF8(in);
            ITF8.readUnsignedITF8(in);
            ITF8.readUnsignedITF8(in);
            LTF8.readUnsignedLTF8(in);
            LTF8.readUnsignedLTF8(in);
            ITF8.readUnsignedITF8(in);
            final int landmarkCount = ITF8.readUnsignedITF8(in);
            for (int i = 0; i < landmarkCount; i++) {
                ITF8.readUnsignedITF8(in);
            }
            CramInt.readInt32(in);
            final int headerLength = (cram.length - at) - in.available();
            if (container == 1) {
                return Arrays.copyOfRange(cram, at + headerLength,
                        at + headerLength + blocksByteSize);
            }
            at += headerLength + blocksByteSize;
        }
        throw new IllegalStateException("no data container");
    }

    static String hex(final byte[] bytes) {
        final StringBuilder b = new StringBuilder(bytes.length * 2);
        for (final byte value : bytes) {
            b.append(String.format("%02x", value));
        }
        return b.length() == 0 ? "-" : b.toString();
    }
}
