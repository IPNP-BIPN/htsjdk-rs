/*
 * The CRAM compression header's data series encoding map: the second of its three maps, and the
 * first thing in CRAM that describes a record rather than a container.
 *
 * cram-preservation-map pinned the first map. This is the one that says, for each of the record's
 * data series, which of the ten encodings carries it and with what parameters.
 *
 * Six things are decisions rather than layout.
 *
 *   - THIS MAP'S SIZE IS A REAL COUNT, WHERE THE PRESERVATION MAP'S IS THE LITERAL 5. Two maps in
 *     the same header, one counted and one not, and the count here excludes any data series whose
 *     encoding is NULL;
 *   - THE WRITE ORDER IS THE DataSeries ENUM'S ORDINAL ORDER, not the order the constructor
 *     populates them in. The map is a TreeMap keyed by the enum, so it sorts by ordinal, while the
 *     constructor adds them alphabetically by canonical name. The two orders are different and only
 *     the first is in the bytes;
 *   - htsjdk WRITES 26 OF THE 32 DATA SERIES. BB and QQ are never written by this implementation,
 *     TC and TN are obsolete, and TM and TV exist only for tests. A port that writes all 32 writes
 *     a map no htsjdk-written CRAM contains;
 *   - TC AND TN ARE READ AND THEN DROPPED. A CRAM from another writer that carries them gets a log
 *     warning and the entries never reach the map, so the reader's map can be smaller than the
 *     count the file declared;
 *   - THE CONTENT IDS ARE htsjdk'S, NOT THE SPECIFICATION'S. The spec does not prescribe them; this
 *     implementation numbers the data series 1 to 32 in enum order, and a reader must discover them
 *     from the map rather than assume them;
 *   - AN UNKNOWN ENCODING ID IS AN ARRAY INDEX. `EncodingID.values()[buffer.get()]` has no bounds
 *     check, so a tenth encoding is an ArrayIndexOutOfBoundsException with a null message rather
 *     than a CRAM error naming the id.
 *
 * Output:
 *
 *     sizes\t<data series count>\t<encoding id count>
 *     series\t<ordinal>\t<canonical name>\t<type>\t<content id>
 *     encid\t<id>\t<name>\t<is external>
 *     map\t<label>\t<declared byte size>\t<declared map size>\t<the whole map, hex>
 *     entry\t<label>\t<index>\t<series>\t<encoding id>\t<param length>\t<params, hex>
 *     written\t<the series htsjdk writes, in write order>
 *     ignored\t<the series read and dropped>
 *     refuse\t<case>\t<class>\t<message>
 *     err\t<label>\t<class>\t<message>
 *
 * Usage: CramEncodingMapDump
 */

import htsjdk.samtools.*;
import htsjdk.samtools.cram.io.CramInt;
import htsjdk.samtools.cram.io.ITF8;
import htsjdk.samtools.cram.io.LTF8;
import htsjdk.samtools.cram.structure.CompressionHeaderEncodingMap;
import htsjdk.samtools.cram.structure.CramHeader;
import htsjdk.samtools.cram.structure.DataSeries;
import htsjdk.samtools.cram.structure.EncodingID;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.util.Arrays;
import java.util.StringJoiner;

public class CramEncodingMapDump {

    public static void main(final String[] args) throws Exception {
        System.out.println("# CramEncodingMapDump: the data series encoding map");

        System.out.printf("sizes\t%d\t%d%n", DataSeries.values().length, EncodingID.values().length);

        // The data series, in the ordinal order the TreeMap writes them in, with the content ids
        // this implementation assigns and the specification does not.
        for (final DataSeries series : DataSeries.values()) {
            System.out.printf("series\t%d\t%s\t%s\t%d%n", series.ordinal(),
                    series.getCanonicalName(), series.getType(),
                    series.getExternalBlockContentId());
        }
        for (final EncodingID id : EncodingID.values()) {
            System.out.printf("encid\t%d\t%s\t%s%n", id.getId(), id.name(),
                    id.isExternalEncoding());
        }

        emit("four-unmapped", build(4, 8, false));
        emit("tagged", build(4, 8, true));

        // The series htsjdk actually writes, in the order it writes them, taken from a real file.
        final byte[] map = encodingMapBytes(build(4, 8, false));
        final ByteArrayInputStream in = new ByteArrayInputStream(map);
        ITF8.readUnsignedITF8(in);
        final StringJoiner written = new StringJoiner(",");
        final int count = ITF8.readUnsignedITF8(in);
        for (int i = 0; i < count; i++) {
            written.add("" + (char) in.read() + (char) in.read());
            in.read();
            final int paramLen = ITF8.readUnsignedITF8(in);
            in.skip(paramLen);
        }
        System.out.printf("written\t%s%n", written);

        final StringJoiner ignored = new StringJoiner(",");
        for (final DataSeries series : CompressionHeaderEncodingMap.DATASERIES_NOT_READ_BY_HTSJDK) {
            ignored.add(series.getCanonicalName());
        }
        System.out.printf("ignored\t%s%n", ignored);

        // What the reader refuses, and with what.
        refuse("unknown-canonical-name", handBuilt("ZZ", 1, new byte[] {1}));
        // The reader always takes exactly two bytes for a name, so the "exactly two characters"
        // branch of byCanonicalName is unreachable from here: a control byte is just a name that
        // does not exist.
        refuse("name-with-a-control-byte", handBuiltRaw(new byte[] {1, 'Z', 1, 1, 1}));
        refuse("encoding-id-past-the-end", handBuilt("BF", 10, new byte[] {1}));
        refuse("encoding-id-negative", handBuilt("BF", 255, new byte[] {1}));
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

    static void emit(final String label, final byte[] cram) throws Exception {
        final byte[] prefixed = encodingMapBytes(cram);
        final ByteArrayInputStream in = new ByteArrayInputStream(prefixed);
        final int byteSize = ITF8.readUnsignedITF8(in);
        final byte[] map = new byte[byteSize];
        if (in.read(map) != byteSize) {
            System.out.printf("err\t%s\tshort\tthe encoding map is shorter than declared%n", label);
            return;
        }

        final ByteArrayInputStream mapIn = new ByteArrayInputStream(map);
        final int mapSize = ITF8.readUnsignedITF8(mapIn);
        System.out.printf("map\t%s\t%d\t%d\t%s%n", label, byteSize, mapSize, hex(map));

        for (int i = 0; i < mapSize; i++) {
            final String name = "" + (char) mapIn.read() + (char) mapIn.read();
            final int encodingId = mapIn.read();
            final int paramLen = ITF8.readUnsignedITF8(mapIn);
            final byte[] params = new byte[paramLen];
            if (mapIn.read(params) != paramLen && paramLen != 0) {
                System.out.printf("err\t%s\tshort\tparameters for %s%n", label, name);
                return;
            }
            System.out.printf("entry\t%s\t%d\t%s\t%d\t%d\t%s%n", label, i, name, encodingId,
                    paramLen, hex(params));
        }
    }

    static void refuse(final String label, final byte[] prefixedMap) {
        try {
            new CompressionHeaderEncodingMap(new ByteArrayInputStream(prefixedMap));
            System.out.printf("refuse\t%s\tnone\taccepted%n", label);
        } catch (final Throwable t) {
            System.out.printf("refuse\t%s\t%s\t%s%n", label, t.getClass().getName(),
                    escape(t.getMessage()));
        }
    }

    /**
     * A message with every character outside printable ASCII written as \\uXXXX. The reader always
     * takes exactly two bytes for a name, so a name can contain a control byte and the message that
     * quotes it back would put that byte in the golden.
     */
    static String escape(final String message) {
        if (message == null) {
            return "<null>";
        }
        final StringBuilder out = new StringBuilder(message.length());
        for (int i = 0; i < message.length(); i++) {
            final char c = message.charAt(i);
            if (c >= 0x20 && c <= 0x7e) {
                out.append(c);
            } else {
                out.append(String.format("\\u%04x", (int) c));
            }
        }
        return out.toString();
    }

    /** One entry, with its length prefix, as the reading constructor expects it. */
    static byte[] handBuilt(final String name, final int encodingId, final byte[] params)
            throws Exception {
        final ByteArrayOutputStream map = new ByteArrayOutputStream();
        ITF8.writeUnsignedITF8(1, map);
        map.write(name.getBytes());
        map.write(encodingId);
        ITF8.writeUnsignedITF8(params.length, map);
        map.write(params);
        return handBuiltRaw(map.toByteArray());
    }

    static byte[] handBuiltRaw(final byte[] map) throws Exception {
        final ByteArrayOutputStream out = new ByteArrayOutputStream();
        ITF8.writeUnsignedITF8(map.length, out);
        out.write(map);
        return out.toByteArray();
    }

    /**
     * The encoding map, with its own length prefix, from a real CRAM: it is what follows the
     * preservation map inside the compression header block's raw content.
     */
    static byte[] encodingMapBytes(final byte[] cram) {
        final byte[] content = compressionHeaderContent(cram);
        final ByteArrayInputStream in = new ByteArrayInputStream(content);
        final int preservationSize = ITF8.readUnsignedITF8(in);
        in.skip(preservationSize);
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

    static String hex(final byte[] bytes) {
        final StringBuilder b = new StringBuilder(bytes.length * 2);
        for (final byte value : bytes) {
            b.append(String.format("%02x", value));
        }
        return b.length() == 0 ? "-" : b.toString();
    }
}
