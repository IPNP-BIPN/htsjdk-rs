/*
 * The CRAM file definition and the container header: the first structures built out of the ITF8s
 * the cram-varint suite pinned.
 *
 * A CRAM file is a 26-byte definition, then containers. The definition is fixed-width and the
 * container header is not: it mixes little-endian int32s with ITF8s and LTF8s, and its last field
 * exists only from major version 3.
 *
 * Four things a reading of the specification does not give you.
 *
 *   - THE FILE ID IS ZERO-PADDED TO EXACTLY 20 BYTES AND SILENTLY TRUNCATED. `CramHeader` fills the
 *     array with zeros and copies `min(id.length(), 20)`, so a longer id loses its tail with no
 *     error and a shorter one is padded rather than terminated;
 *   - THE CHECKSUM COVERS THE HEADER AND NOT THE CONTAINER. The CRC32 is computed over the bytes of
 *     the container header itself, up to but not including the checksum, and it is written LITTLE
 *     ENDIAN, which is the opposite of the CRC in a BGZF block's gzip trailer;
 *   - THE CHECKSUM IS ABSENT BELOW VERSION 3, so the same container is four bytes shorter in a 2.1
 *     file and a reader that always consumes them is off by four from the first container on;
 *   - `writeContainerHeader` RETURNS A BYTE COUNT COMPUTED FROM BIT COUNTS: every writer returns
 *     bits and the caller does `(bits + 7) / 8`. For ITF8 and LTF8 those are always multiples of
 *     eight, so the rounding never rounds; the expression is there and does nothing.
 *
 * The SAM header container is a shape of its own: `makeSAMFileHeaderContainer` fixes the alignment
 * context to unmapped-unplaced, one block, no records, no landmarks and no checksum, so the first
 * container of every CRAM has the same header but for its size.
 *
 * Output:
 *
 *     def\t<label>\t<magic>\t<major>.<minor>\t<id bytes, hex>
 *     hdr\t<label>\t<index>\t<blocksByteSize>\t<refContextId>\t<start>\t<span>\t<records>\t<counter>\t<bases>\t<blocks>\t<landmarks>\t<checksum>
 *     bytes\t<label>\t<index>\t<the container header's own bytes, hex>
 *     file\t<label>\t<total bytes>\t<container count>
 *     err\t<label>\t<class>\t<message>
 *
 * Usage: CramContainerDump
 */

import htsjdk.samtools.*;
import htsjdk.samtools.cram.common.CRAMVersion;
import htsjdk.samtools.cram.io.CramInt;
import htsjdk.samtools.cram.io.ITF8;
import htsjdk.samtools.cram.io.LTF8;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.cram.structure.CramHeader;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.StringJoiner;

public class CramContainerDump {

    public static void main(final String[] args) throws Exception {
        System.out.println("# CramContainerDump: the file definition and the container header");

        // Unmapped reads need no reference, which keeps this suite about the container rather than
        // about reference compression.
        emit("four-unmapped", build(4, 8));
        emit("one-unmapped", build(1, 8));
        emit("long-reads", build(2, 400));
        emit("no-reads", build(0, 0));

        // The file definition on its own, for ids the writer has to pad or truncate.
        for (final String id : new String[] {"", "short", "exactly-twenty-byte!", "far-too-long-to-fit-in-twenty"}) {
            final CramHeader header = new CramHeader(new CRAMVersion(3, 0), id);
            System.out.printf("id\t%s\t%s%n", escape(id), hex(header.getId()));
        }
    }

    /** A CRAM of `count` unmapped reads whose bases are `length` long. */
    static byte[] build(final int count, final int length) {
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
                Arrays.fill(quals, (byte) 30);
                record.setBaseQualities(quals);
                writer.addAlignment(record);
            }
        }
        return out.toByteArray();
    }

    static void emit(final String label, final byte[] cram) throws Exception {
        // The file definition: 4 magic, 2 version, 20 id.
        System.out.printf("def\t%s\t%s\t%d.%d\t%s%n", label,
                new String(Arrays.copyOfRange(cram, 0, 4)),
                cram[4] & 0xFF, cram[5] & 0xFF,
                hex(Arrays.copyOfRange(cram, 6, 26)));

        final CRAMVersion version = new CRAMVersion(cram[4] & 0xFF, cram[5] & 0xFF);
        int at = CramHeader.CRAM_HEADER_LENGTH;
        int index = 0;

        // Walk the containers by hand rather than through the reader, because the layout is what
        // this suite is about and the reader would hide it.
        while (at < cram.length) {
            final int start = at;
            final ByteArrayInputStream in =
                    new ByteArrayInputStream(cram, at, cram.length - at);

            final int blocksByteSize = CramInt.readInt32(in);
            final int refContextId = ITF8.readUnsignedITF8(in);
            final int alignmentStart = ITF8.readUnsignedITF8(in);
            final int alignmentSpan = ITF8.readUnsignedITF8(in);
            final int recordCount = ITF8.readUnsignedITF8(in);
            final long globalRecordCounter = LTF8.readUnsignedLTF8(in);
            final long baseCount = LTF8.readUnsignedLTF8(in);
            final int blockCount = ITF8.readUnsignedITF8(in);
            final int landmarkCount = ITF8.readUnsignedITF8(in);
            final List<Integer> landmarks = new ArrayList<>();
            for (int i = 0; i < landmarkCount; i++) {
                landmarks.add(ITF8.readUnsignedITF8(in));
            }
            final int checksum = version.getMajor() >= 3 ? CramInt.readInt32(in) : 0;

            final int consumed = (cram.length - at) - in.available();
            final StringJoiner joined = new StringJoiner(",");
            for (final int landmark : landmarks) {
                joined.add(Integer.toString(landmark));
            }

            System.out.printf("hdr\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%d%n",
                    label, index, blocksByteSize, refContextId, alignmentStart, alignmentSpan,
                    recordCount, globalRecordCounter, baseCount, blockCount,
                    landmarks.isEmpty() ? "-" : joined.toString(), checksum);
            System.out.printf("bytes\t%s\t%d\t%s%n", label, index,
                    hex(Arrays.copyOfRange(cram, start, start + consumed)));

            at = start + consumed + blocksByteSize;
            index++;
            if (blocksByteSize == 0 && recordCount == 0 && blockCount == 1) {
                // The EOF container has a fixed shape and no blocks worth walking past; it is the
                // last thing in the file either way.
                if (at >= cram.length) {
                    break;
                }
            }
        }
        System.out.printf("file\t%s\t%d\t%d%n", label, cram.length, index);
    }

    static String hex(final byte[] bytes) {
        final StringBuilder b = new StringBuilder(bytes.length * 2);
        for (final byte value : bytes) {
            b.append(String.format("%02x", value));
        }
        return b.length() == 0 ? "-" : b.toString();
    }

    static String escape(final String s) {
        return s.isEmpty() ? "<empty>" : s;
    }
}
