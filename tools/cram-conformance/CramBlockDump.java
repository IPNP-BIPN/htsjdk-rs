/*
 * The CRAM block: what a container's bytes are actually made of.
 *
 * The container header pinned by cram-container gives a size and a count; the blocks are what that
 * size counts. Each is a five-field header, then the compressed content, then a CRC32 from version
 * 3 up.
 *
 * Four things worth measuring rather than reading.
 *
 *   - THE CRC COVERS THE HEADER AND THE CONTENT TOGETHER, not the content alone, and it is a
 *     CRC32InputStream wrapped around the whole read. So a block cannot be verified without
 *     re-reading its own header, and the four checksum bytes are outside the compressedSize the
 *     header declares;
 *   - IT IS ABSENT BELOW VERSION 3, exactly as the container header's is, so a 2.1 block is four
 *     bytes shorter and the two version-dependent lengths compound;
 *   - A CONTENT ID IS ONLY LEGAL ON AN EXTERNAL BLOCK. Block's constructor throws "Cannot set a
 *     Content ID for non-external blocks" otherwise, so the field is present in every block and
 *     meaningful in one kind;
 *   - THE SAM HEADER BLOCK IS GZIP AND THE COMPRESSION HEADER BLOCK IS RAW. Neither is a choice
 *     the writer makes per file: createGZIPFileHeaderBlock and createRawCompressionHeaderBlock fix
 *     them, so the first two blocks of every CRAM have known methods.
 *
 * The compression methods are RAW 0, GZIP 1, BZIP2 2, LZMA 3, RANS 4, RANGE 5, and the content
 * types are FILE_HEADER 0, COMPRESSION_HEADER 1, MAPPED_SLICE 2, RESERVED 3, EXTERNAL 4, CORE 5.
 * Both are dumped by id rather than by name so a port compares numbers.
 *
 * Output:
 *
 *     blk\t<label>\t<container>\t<index>\t<method>\t<contentType>\t<contentId>\t<compressedSize>\t<uncompressedSize>\t<crc32>
 *     hdrbytes\t<label>\t<container>\t<index>\t<the block header's own bytes, hex>
 *     content\t<label>\t<container>\t<index>\t<sha256 of the compressed content>
 *     counts\t<label>\t<blocks seen>\t<methods used>
 *     err\t<label>\t<class>\t<message>
 *
 * Usage: CramBlockDump
 */

import htsjdk.samtools.*;
import htsjdk.samtools.cram.common.CRAMVersion;
import htsjdk.samtools.cram.io.CramInt;
import htsjdk.samtools.cram.io.ITF8;
import htsjdk.samtools.cram.structure.CramHeader;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.StringJoiner;
import java.util.TreeSet;

public class CramBlockDump {

    public static void main(final String[] args) throws Exception {
        System.out.println("# CramBlockDump: the blocks a container's byte size counts");

        emit("four-unmapped", build(4, 8));
        emit("one-unmapped", build(1, 8));
        emit("long-reads", build(2, 400));
        emit("no-reads", build(0, 0));
        emit("many-reads", build(40, 20));
    }

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
                Arrays.fill(quals, (byte) (30 + i % 5));
                record.setBaseQualities(quals);
                writer.addAlignment(record);
            }
        }
        return out.toByteArray();
    }

    static void emit(final String label, final byte[] cram) throws Exception {
        final CRAMVersion version = new CRAMVersion(cram[4] & 0xFF, cram[5] & 0xFF);
        final boolean v3 = version.getMajor() >= 3;
        int at = CramHeader.CRAM_HEADER_LENGTH;
        int container = 0;
        int blocksSeen = 0;
        final TreeSet<Integer> methods = new TreeSet<>();

        while (at < cram.length) {
            final ByteArrayInputStream in = new ByteArrayInputStream(cram, at, cram.length - at);
            final int blocksByteSize = CramInt.readInt32(in);
            ITF8.readUnsignedITF8(in); // reference context
            final int alignmentStart = ITF8.readUnsignedITF8(in);
            ITF8.readUnsignedITF8(in); // span
            final int recordCount = ITF8.readUnsignedITF8(in);
            htsjdk.samtools.cram.io.LTF8.readUnsignedLTF8(in);
            htsjdk.samtools.cram.io.LTF8.readUnsignedLTF8(in);
            final int blockCount = ITF8.readUnsignedITF8(in);
            final int landmarkCount = ITF8.readUnsignedITF8(in);
            for (int i = 0; i < landmarkCount; i++) {
                ITF8.readUnsignedITF8(in);
            }
            if (v3) {
                CramInt.readInt32(in);
            }
            final int headerLength = (cram.length - at) - in.available();

            // The EOF container's single block is not worth walking into: it is a fixed sequence
            // and the file ends there.
            final boolean eof = alignmentStart == 0x454F46 && recordCount == 0;

            int blockAt = at + headerLength;
            final int blocksEnd = blockAt + blocksByteSize;
            for (int index = 0; index < blockCount && blockAt < blocksEnd; index++) {
                final int start = blockAt;
                final ByteArrayInputStream bin =
                        new ByteArrayInputStream(cram, blockAt, cram.length - blockAt);
                final int method = bin.read();
                final int contentType = bin.read();
                final int contentId = ITF8.readUnsignedITF8(bin);
                final int compressedSize = ITF8.readUnsignedITF8(bin);
                final int uncompressedSize = ITF8.readUnsignedITF8(bin);
                final int blockHeaderLength = (cram.length - blockAt) - bin.available();

                final byte[] content = Arrays.copyOfRange(
                        cram, start + blockHeaderLength, start + blockHeaderLength + compressedSize);

                int crc = 0;
                if (v3) {
                    crc = CramInt.readInt32(new ByteArrayInputStream(
                            cram, start + blockHeaderLength + compressedSize, 4));
                }

                System.out.printf("blk\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d%n", label, container,
                        index, method, contentType, contentId, compressedSize, uncompressedSize,
                        crc);
                System.out.printf("hdrbytes\t%s\t%d\t%d\t%s%n", label, container, index,
                        hex(Arrays.copyOfRange(cram, start, start + blockHeaderLength)));
                System.out.printf("content\t%s\t%d\t%d\t%s%n", label, container, index,
                        sha256(content));

                methods.add(method);
                blocksSeen++;
                blockAt = start + blockHeaderLength + compressedSize + (v3 ? 4 : 0);
            }

            at += headerLength + blocksByteSize;
            container++;
            if (eof) {
                break;
            }
        }

        final StringJoiner joined = new StringJoiner(",");
        for (final int method : methods) {
            joined.add(Integer.toString(method));
        }
        System.out.printf("counts\t%s\t%d\t%s%n", label, blocksSeen, joined);
    }

    static String hex(final byte[] bytes) {
        final StringBuilder b = new StringBuilder(bytes.length * 2);
        for (final byte value : bytes) {
            b.append(String.format("%02x", value));
        }
        return b.length() == 0 ? "-" : b.toString();
    }

    /** A digest rather than the bytes: the content is what the compressors produce and this suite
     *  is about the framing around it. The digest still fails if the framing put the boundary in
     *  the wrong place. */
    static String sha256(final byte[] bytes) throws Exception {
        return hex(MessageDigest.getInstance("SHA-256").digest(bytes));
    }
}
