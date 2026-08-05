/*
 * The CRAM compression header's preservation map: the first thing inside the RAW block that
 * cram-block measured as the second block of every file.
 *
 * The compression header is three length-prefixed maps in a row: the preservation map, the data
 * series encoding map, and the tag encoding map. This suite is the first of them, and it is where
 * a CRAM says whether read names are kept, whether alignment positions are deltas, whether a
 * reference is required, what its substitution matrix is and which tag id combinations appear.
 *
 * Six things are decisions rather than layout.
 *
 *   - THE MAP SIZE IS A HARDCODED 5, not a count. `internalWrite` calls
 *     `ITF8.writeUnsignedITF8(5, mapBuffer)` and then writes exactly RN, AP, RR, SM and TD in that
 *     order, whatever the header holds. The field is a constant wearing a count's clothes;
 *   - THE WRITE ORDER IS NOT THE SPECIFICATION'S ORDER. htsjdk writes RN, AP, RR, SM, TD. The
 *     reader accepts them in any order, so this is only visible on the bytes, which is exactly what
 *     a byte-identical port has to reproduce;
 *   - A BOOLEAN IS `== 1`, NOT `!= 0`. `preserveReadNames = buffer.get() == 1`, so a 2 reads as
 *     FALSE and no error is raised. Three of the five keys go through it;
 *   - SM AND TD ARE MANDATORY, and the check is after the loop rather than in it: a header that
 *     omits either throws a CRAMException naming both, whatever the map size said;
 *   - AN UNKNOWN KEY IS A PLAIN RuntimeException, not a CRAMException, and it carries the two
 *     characters it did not recognise;
 *   - THE TAG ID DICTIONARY IS A RUN OF THREE-BYTE IDS TERMINATED BY A ZERO, per group, with the
 *     groups concatenated. `parseDictionary` computes a `maxWidth` it never uses, and reads three
 *     bytes at a time with no check that the terminator falls on a boundary.
 *
 * Output:
 *
 *     block\t<label>\t<compression header block's raw content length>\t<sha256>
 *     pmap\t<label>\t<declared byte size>\t<declared map size>\t<the whole map, hex>
 *     key\t<label>\t<index>\t<two-character key>\t<payload, hex>
 *     flags\t<label>\t<preserveReadNames>\t<apDelta>\t<referenceRequired>
 *     sm\t<label>\t<the substitution matrix bytes, hex>
 *     td\t<label>\t<dictionary length>\t<dictionary, hex>
 *     tdgroup\t<label>\t<group index>\t<the three-byte ids in it, comma separated>
 *     boolean\t<byte written>\t<what the reader made of it>
 *     refuse\t<case>\t<class>\t<message>
 *     sizes\t<basesSize>
 *     err\t<label>\t<class>\t<message>
 *
 * Usage: CramPreservationMapDump
 */

import htsjdk.samtools.*;
import htsjdk.samtools.cram.common.CRAMVersion;
import htsjdk.samtools.cram.io.CramInt;
import htsjdk.samtools.cram.io.ITF8;
import htsjdk.samtools.cram.io.LTF8;
import htsjdk.samtools.cram.structure.CompressionHeader;
import htsjdk.samtools.cram.structure.CramHeader;
import htsjdk.samtools.cram.structure.SubstitutionMatrix;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.StringJoiner;

public class CramPreservationMapDump {

    private static final CRAMVersion VERSION_3_0 = new CRAMVersion(3, 0);

    public static void main(final String[] args) throws Exception {
        System.out.println("# CramPreservationMapDump: the compression header's first map");

        System.out.printf("sizes\t%d%n", SubstitutionMatrix.BASES_SIZE);

        emit("four-unmapped", build(4, 8, false));
        emit("one-unmapped", build(1, 8, false));
        emit("long-reads", build(2, 400, false));
        emit("many-reads", build(40, 20, false));
        emit("tagged", build(4, 8, true));

        // A boolean is `== 1`, so everything else is false and nothing complains.
        for (final int value : new int[] {0, 1, 2, 127, 255}) {
            try {
                final CompressionHeader header = readHeader(handBuilt(value, "RN", true, true));
                System.out.printf("boolean\t%d\t%s%n", value, header.isPreserveReadNames());
            } catch (final Exception e) {
                System.out.printf("boolean\t%d\tERR %s%n", value, e.getClass().getSimpleName());
            }
        }

        // The refusals, each on a header that is well formed but for the one thing being tested.
        refuse("unknown-key", unknownKey());
        refuse("no-substitution-matrix", handBuilt(1, "RN", false, true));
        refuse("no-tag-dictionary", handBuilt(1, "RN", true, false));
    }

    /** A CRAM of `count` unmapped reads whose bases are `length` long, optionally with tags. */
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
                    record.setAttribute("XX", (float) i);
                }
                writer.addAlignment(record);
            }
        }
        return out.toByteArray();
    }

    static void emit(final String label, final byte[] cram) throws Exception {
        final byte[] content = compressionHeaderContent(cram);
        System.out.printf("block\t%s\t%d\t%s%n", label, content.length, sha256(content));

        // The preservation map is the first length-prefixed run of the compression header.
        final ByteArrayInputStream in = new ByteArrayInputStream(content);
        final int byteSize = ITF8.readUnsignedITF8(in);
        final byte[] map = new byte[byteSize];
        if (in.read(map) != byteSize) {
            System.out.printf("err\t%s\tshort\tthe preservation map is shorter than declared%n", label);
            return;
        }
        final ByteArrayInputStream mapIn = new ByteArrayInputStream(map);
        final int mapSize = ITF8.readUnsignedITF8(mapIn);
        System.out.printf("pmap\t%s\t%d\t%d\t%s%n", label, byteSize, mapSize, hex(map));

        for (int i = 0; i < mapSize; i++) {
            final String key = "" + (char) mapIn.read() + (char) mapIn.read();
            final byte[] payload;
            switch (key) {
                case "RN":
                case "AP":
                case "RR":
                    payload = new byte[] {(byte) mapIn.read()};
                    break;
                case "SM":
                    payload = new byte[SubstitutionMatrix.BASES_SIZE];
                    if (mapIn.read(payload) != payload.length) {
                        System.out.printf("err\t%s\tshort\tsubstitution matrix%n", label);
                        return;
                    }
                    System.out.printf("sm\t%s\t%s%n", label, hex(payload));
                    break;
                case "TD": {
                    final int size = ITF8.readUnsignedITF8(mapIn);
                    payload = new byte[size];
                    if (mapIn.read(payload) != size) {
                        System.out.printf("err\t%s\tshort\ttag dictionary%n", label);
                        return;
                    }
                    System.out.printf("td\t%s\t%d\t%s%n", label, size, hex(payload));
                    emitDictionaryGroups(label, payload);
                    break;
                }
                default:
                    System.out.printf("err\t%s\tunknown-key\t%s%n", label, key);
                    return;
            }
            System.out.printf("key\t%s\t%d\t%s\t%s%n", label, i, key, hex(payload));
        }

        final CompressionHeader header = new CompressionHeader(VERSION_3_0,
                new ByteArrayInputStream(compressionHeaderBlockBytes(cram)));
        System.out.printf("flags\t%s\t%s\t%s\t%s%n", label, header.isPreserveReadNames(),
                header.isAPDelta(), header.isReferenceRequired());
    }

    /** `parseDictionary`: three bytes at a time until a zero, per group. */
    static void emitDictionaryGroups(final String label, final byte[] bytes) {
        int at = 0;
        int group = 0;
        while (at < bytes.length) {
            final StringJoiner joined = new StringJoiner(",");
            while (at < bytes.length && bytes[at] != 0) {
                joined.add(new String(Arrays.copyOfRange(bytes, at, at + 3)));
                at += 3;
            }
            at++;
            System.out.printf("tdgroup\t%s\t%d\t%s%n", label, group,
                    joined.length() == 0 ? "-" : joined.toString());
            group++;
        }
    }

    static void refuse(final String label, final byte[] headerContent) {
        try {
            readHeader(headerContent);
            System.out.printf("refuse\t%s\tnone\taccepted%n", label);
        } catch (final Exception e) {
            System.out.printf("refuse\t%s\t%s\t%s%n", label, e.getClass().getName(),
                    e.getMessage() == null ? "<null>" : e.getMessage().replace('\t', ' '));
        }
    }

    /** Wrap a compression header's content in a RAW block and read it back as htsjdk does. */
    static CompressionHeader readHeader(final byte[] content) throws Exception {
        final ByteArrayOutputStream out = new ByteArrayOutputStream();
        htsjdk.samtools.cram.structure.block.Block.createRawCompressionHeaderBlock(content)
                .write(VERSION_3_0, out);
        return new CompressionHeader(VERSION_3_0, new ByteArrayInputStream(out.toByteArray()));
    }

    /**
     * A compression header built by hand: a preservation map with the boolean under test, then two
     * empty maps. `withMatrix` and `withDictionary` decide whether the two mandatory keys appear.
     */
    static byte[] handBuilt(final int booleanValue, final String booleanKey,
            final boolean withMatrix, final boolean withDictionary) throws Exception {
        final ByteArrayOutputStream map = new ByteArrayOutputStream();
        int keys = 3;
        if (withMatrix) {
            keys++;
        }
        if (withDictionary) {
            keys++;
        }
        ITF8.writeUnsignedITF8(keys, map);
        for (final String key : new String[] {"RN", "AP", "RR"}) {
            map.write(key.getBytes());
            map.write(key.equals(booleanKey) ? booleanValue : 1);
        }
        if (withMatrix) {
            map.write("SM".getBytes());
            map.write(new byte[SubstitutionMatrix.BASES_SIZE]);
        }
        if (withDictionary) {
            map.write("TD".getBytes());
            // One group holding one three-byte id, then its terminator.
            final byte[] dictionary = new byte[] {'N', 'M', 'i', 0};
            ITF8.writeUnsignedITF8(dictionary.length, map);
            map.write(dictionary);
        }
        return wrapMaps(map.toByteArray());
    }

    /** The same, but with a key htsjdk does not know. */
    static byte[] unknownKey() throws Exception {
        final ByteArrayOutputStream map = new ByteArrayOutputStream();
        ITF8.writeUnsignedITF8(1, map);
        map.write("ZZ".getBytes());
        map.write(1);
        return wrapMaps(map.toByteArray());
    }

    /** A preservation map followed by an empty encoding map and an empty tag encoding map. */
    static byte[] wrapMaps(final byte[] preservationMap) throws Exception {
        final ByteArrayOutputStream out = new ByteArrayOutputStream();
        ITF8.writeUnsignedITF8(preservationMap.length, out);
        out.write(preservationMap);
        for (int i = 0; i < 2; i++) {
            final ByteArrayOutputStream empty = new ByteArrayOutputStream();
            ITF8.writeUnsignedITF8(0, empty);
            ITF8.writeUnsignedITF8(empty.size(), out);
            out.write(empty.toByteArray());
        }
        return out.toByteArray();
    }

    /** The raw content of the compression header block, which cram-block measured as RAW. */
    static byte[] compressionHeaderContent(final byte[] cram) {
        final byte[] blockBytes = compressionHeaderBlockBytes(cram);
        final ByteArrayInputStream in = new ByteArrayInputStream(blockBytes);
        in.read(); // method
        in.read(); // content type
        ITF8.readUnsignedITF8(in); // content id
        final int compressedSize = ITF8.readUnsignedITF8(in);
        ITF8.readUnsignedITF8(in); // uncompressed size
        final int headerLength = blockBytes.length - in.available();
        return Arrays.copyOfRange(blockBytes, headerLength, headerLength + compressedSize);
    }

    /** The bytes of the compression header block, starting at its own header. */
    static byte[] compressionHeaderBlockBytes(final byte[] cram) {
        int at = CramHeader.CRAM_HEADER_LENGTH;
        // The first container is the SAM header; the second begins with the compression header
        // block, which cram-block measured as RAW in every file.
        for (int container = 0; container < 2; container++) {
            final ByteArrayInputStream in =
                    new ByteArrayInputStream(cram, at, cram.length - at);
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

    static String sha256(final byte[] bytes) throws Exception {
        return hex(MessageDigest.getInstance("SHA-256").digest(bytes));
    }
}
