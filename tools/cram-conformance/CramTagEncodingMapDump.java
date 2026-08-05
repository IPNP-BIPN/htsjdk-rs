/*
 * The CRAM compression header's tag encoding map: the third and last of its maps, and the one
 * that closes the header.
 *
 * cram-encoding-map pinned the second map, whose keys are two-character data series names. This one
 * has the same entry shape but its keys are integers, and the integer is the tag itself: the two
 * characters of its name and its type character, packed into twenty-four bits.
 *
 * Five things are decisions rather than layout.
 *
 *   - THE KEY IS THE TAG, PACKED. `nameType3BytesToInt` shifts name[0], name[1] and the TYPE into
 *     one int, so the map's key carries the whole identity of the tag and needs no separate name
 *     field. NMc is 0x4E4D63, which is 5131107;
 *   - THE TYPE IS PART OF THE KEY, so the same tag name at two types is two entries with two
 *     external blocks, not one entry that carries a type;
 *   - THE WRITE ORDER IS NUMERIC ORDER OF THAT KEY, because the map is a TreeMap<Integer>. For
 *     printable tag names that is lexicographic order by name then type, which is not the order the
 *     records introduced the tags in;
 *   - THE COLLISION GUARD CANNOT FIRE. `putTagBlockCompression` refuses a tag id that equals a data
 *     series content id, and those are 1 to 32, while the smallest printable tag packs to 0x202020
 *     = 2105376. The check is real code that no input can reach;
 *   - THE SIZE IS A REAL COUNT, like the encoding map's and unlike the preservation map's literal 5.
 *     Three maps in one header and two counting conventions.
 *
 * Output:
 *
 *     id\t<name>\t<type>\t<packed key>\t<intToNameType3Bytes>\t<intToNameType4Bytes>
 *     tmap\t<label>\t<declared byte size>\t<declared map size>\t<the whole map, hex>
 *     tentry\t<label>\t<index>\t<key>\t<name and type>\t<encoding id>\t<param length>\t<params, hex>
 *     order\t<label>\t<the tags in write order>
 *     guard\t<tag id>\t<accepted or the message>
 *     err\t<label>\t<class>\t<message>
 *
 * Usage: CramTagEncodingMapDump
 */

import htsjdk.samtools.*;
import htsjdk.samtools.cram.io.CramInt;
import htsjdk.samtools.cram.io.ITF8;
import htsjdk.samtools.cram.io.LTF8;
import htsjdk.samtools.cram.structure.CompressionHeaderEncodingMap;
import htsjdk.samtools.cram.structure.CramHeader;
import htsjdk.samtools.cram.structure.ReadTag;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.util.Arrays;
import java.util.StringJoiner;

public class CramTagEncodingMapDump {

    public static void main(final String[] args) throws Exception {
        System.out.println("# CramTagEncodingMapDump: the compression header's third map");

        // The key is the tag. Every one of these is a name and a type packed into 24 bits.
        for (final String[] tag : new String[][] {
                {"NM", "c"}, {"NM", "i"}, {"MD", "Z"}, {"XX", "f"}, {"AS", "i"},
                {"aa", "A"}, {"  ", " "}, {"~~", "~"}, {"ZZ", "B"}}) {
            final int key = ReadTag.nameType3BytesToInt(tag[0], tag[1].charAt(0));
            System.out.printf("id\t%s\t%s\t%d\t%s\t%s%n", tag[0], tag[1], key,
                    escape(ReadTag.intToNameType3Bytes(key)),
                    escape(ReadTag.intToNameType4Bytes(key)));
        }

        // A tag map from a real file, and one from a file whose tags arrive out of order.
        emit("tagged", build(new String[][] {{"NM", "i"}, {"MD", "Z"}, {"XX", "f"}}));
        emit("reverse-order", build(new String[][] {{"XX", "f"}, {"NM", "i"}, {"MD", "Z"}}));
        // The same name at two types, one per record, because a second setAttribute on one record
        // replaces the first rather than adding to it.
        emit("same-name-two-types", buildAlternating("XX"));
        // A value too large for a signed byte, which forces the type htsjdk would otherwise narrow.
        emit("large-integer", buildWithValue("NM", 100000));
        emit("untagged", build(new String[][] {}));

        // The guard that cannot fire from any input: a data series content id is 1 to 32, and the
        // smallest printable tag packs to 2105376.
        for (final int tagId : new int[] {1, 32, 33, 2105376}) {
            guard(tagId);
        }
    }

    /** The same tag name at type i on even records and type Z on odd ones. */
    static byte[] buildAlternating(final String name) {
        final SAMFileHeader header = newHeader();
        final ByteArrayOutputStream out = new ByteArrayOutputStream();
        try (final SAMFileWriter writer =
                new SAMFileWriterFactory().makeCRAMWriter(header, out, (java.io.File) null)) {
            for (int i = 0; i < 4; i++) {
                final SAMRecord record = newRecord(header, i);
                if (i % 2 == 0) {
                    record.setAttribute(name, 100000);
                } else {
                    record.setAttribute(name, "value");
                }
                writer.addAlignment(record);
            }
        }
        return out.toByteArray();
    }

    /** One tag carrying one integer value, so the type htsjdk chooses for it is visible. */
    static byte[] buildWithValue(final String name, final int value) {
        final SAMFileHeader header = newHeader();
        final ByteArrayOutputStream out = new ByteArrayOutputStream();
        try (final SAMFileWriter writer =
                new SAMFileWriterFactory().makeCRAMWriter(header, out, (java.io.File) null)) {
            for (int i = 0; i < 4; i++) {
                final SAMRecord record = newRecord(header, i);
                record.setAttribute(name, value);
                writer.addAlignment(record);
            }
        }
        return out.toByteArray();
    }

    static SAMFileHeader newHeader() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        header.addSequence(new SAMSequenceRecord("chr1", 100000));
        return header;
    }

    static SAMRecord newRecord(final SAMFileHeader header, final int i) {
        final SAMRecord record = new SAMRecord(header);
        record.setReadName("read" + i);
        record.setReadUnmappedFlag(true);
        record.setReferenceIndex(-1);
        record.setAlignmentStart(0);
        final byte[] bases = new byte[8];
        Arrays.fill(bases, (byte) "ACGT".charAt(i % 4));
        record.setReadBases(bases);
        final byte[] quals = new byte[8];
        Arrays.fill(quals, (byte) 30);
        record.setBaseQualities(quals);
        return record;
    }

    static byte[] build(final String[][] tags) {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        header.addSequence(new SAMSequenceRecord("chr1", 100000));

        final ByteArrayOutputStream out = new ByteArrayOutputStream();
        try (final SAMFileWriter writer =
                new SAMFileWriterFactory().makeCRAMWriter(header, out, (java.io.File) null)) {
            for (int i = 0; i < 4; i++) {
                final SAMRecord record = new SAMRecord(header);
                record.setReadName("read" + i);
                record.setReadUnmappedFlag(true);
                record.setReferenceIndex(-1);
                record.setAlignmentStart(0);
                final byte[] bases = new byte[8];
                Arrays.fill(bases, (byte) "ACGT".charAt(i % 4));
                record.setReadBases(bases);
                final byte[] quals = new byte[8];
                Arrays.fill(quals, (byte) 30);
                record.setBaseQualities(quals);
                for (final String[] tag : tags) {
                    switch (tag[1]) {
                        case "i":
                            record.setAttribute(tag[0], i + 1);
                            break;
                        case "Z":
                            record.setAttribute(tag[0], "value");
                            break;
                        case "f":
                            record.setAttribute(tag[0], (float) i);
                            break;
                        default:
                            record.setAttribute(tag[0], (byte) i);
                            break;
                    }
                }
                writer.addAlignment(record);
            }
        }
        return out.toByteArray();
    }

    static void emit(final String label, final byte[] cram) throws Exception {
        final byte[] prefixed = tagMapBytes(cram);
        final ByteArrayInputStream in = new ByteArrayInputStream(prefixed);
        final int byteSize = ITF8.readUnsignedITF8(in);
        final byte[] map = new byte[byteSize];
        if (byteSize > 0 && in.read(map) != byteSize) {
            System.out.printf("err\t%s\tshort\tthe tag map is shorter than declared%n", label);
            return;
        }

        final ByteArrayInputStream mapIn = new ByteArrayInputStream(map);
        final int mapSize = ITF8.readUnsignedITF8(mapIn);
        System.out.printf("tmap\t%s\t%d\t%d\t%s%n", label, byteSize, mapSize, hex(map));

        final StringJoiner order = new StringJoiner(",");
        for (int i = 0; i < mapSize; i++) {
            final int key = ITF8.readUnsignedITF8(mapIn);
            final int encodingId = mapIn.read();
            final int paramLen = ITF8.readUnsignedITF8(mapIn);
            final byte[] params = new byte[paramLen];
            if (paramLen > 0 && mapIn.read(params) != paramLen) {
                System.out.printf("err\t%s\tshort\tparameters for %d%n", label, key);
                return;
            }
            final String name = ReadTag.intToNameType3Bytes(key);
            order.add(name);
            System.out.printf("tentry\t%s\t%d\t%d\t%s\t%d\t%d\t%s%n", label, i, key, escape(name),
                    encodingId, paramLen, hex(params));
        }
        System.out.printf("order\t%s\t%s%n", label, order.length() == 0 ? "-" : order.toString());
    }

    /** `putTagBlockCompression`'s check, which no printable tag can reach. */
    static void guard(final int tagId) {
        try {
            new CompressionHeaderEncodingMap(new htsjdk.samtools.cram.structure.CRAMEncodingStrategy())
                    .putTagBlockCompression(tagId, null);
            System.out.printf("guard\t%d\taccepted%n", tagId);
        } catch (final Exception e) {
            System.out.printf("guard\t%d\t%s%n", tagId, escape(e.getMessage()));
        }
    }

    /**
     * The tag encoding map, with its own length prefix: the third length-prefixed run inside the
     * compression header block's raw content.
     */
    static byte[] tagMapBytes(final byte[] cram) {
        final byte[] content = compressionHeaderContent(cram);
        final ByteArrayInputStream in = new ByteArrayInputStream(content);
        for (int map = 0; map < 2; map++) {
            final int size = ITF8.readUnsignedITF8(in);
            in.skip(size);
        }
        final int consumed = content.length - in.available();
        return Arrays.copyOfRange(content, consumed, content.length);
    }

    static byte[] compressionHeaderContent(final byte[] cram) {
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
                final byte[] blockBytes = Arrays.copyOfRange(cram, at + headerLength,
                        at + headerLength + blocksByteSize);
                final ByteArrayInputStream block = new ByteArrayInputStream(blockBytes);
                block.read();
                block.read();
                ITF8.readUnsignedITF8(block);
                final int compressedSize = ITF8.readUnsignedITF8(block);
                ITF8.readUnsignedITF8(block);
                final int blockHeaderLength = blockBytes.length - block.available();
                return Arrays.copyOfRange(blockBytes, blockHeaderLength,
                        blockHeaderLength + compressedSize);
            }
            at += headerLength + blocksByteSize;
        }
        throw new IllegalStateException("no data container");
    }

    /** Every character outside printable ASCII as \\uXXXX, so a golden stays a text file. */
    static String escape(final String text) {
        if (text == null) {
            return "<null>";
        }
        final StringBuilder out = new StringBuilder(text.length());
        for (int i = 0; i < text.length(); i++) {
            final char c = text.charAt(i);
            if (c >= 0x20 && c <= 0x7e) {
                out.append(c);
            } else {
                out.append(String.format("\\u%04x", (int) c));
            }
        }
        return out.toString();
    }

    static String hex(final byte[] bytes) {
        final StringBuilder b = new StringBuilder(bytes.length * 2);
        for (final byte value : bytes) {
            b.append(String.format("%02x", value));
        }
        return b.length() == 0 ? "-" : b.toString();
    }
}
