/*
 * The compression header as a whole: three length-prefixed maps inside one raw block.
 *
 * Each of the three has been measured on its own. This is what joins them, and it is the last
 * structure between a container's bytes and the codecs that read its slices.
 *
 * Six things here are decisions rather than layout.
 *
 *   - THE ORDER IS THE WRITER'S, not the reader's. The write emits RN, AP, RR, SM and TD in that
 *     order behind a hardcoded count of 5; the read accepts them in any order, so only the bytes
 *     record which order was used;
 *   - THE BLOCK IS ALWAYS RAW, whatever the CRAM version, and it is read back through
 *     getRawContent rather than through a decompressor;
 *   - THE CONTENT TYPE IS CHECKED AND THE MESSAGE NAMES WHAT WAS FOUND, so a block of the wrong
 *     kind is refused rather than parsed;
 *   - TWO KEYS ARE REQUIRED AND THE OTHER THREE ARE NOT. A header without SM or without TD is
 *     refused after the whole map has been read, with one message covering both;
 *   - AN UNKNOWN KEY IS FATAL, and its message carries the two characters that were not
 *     recognised;
 *   - THE VERSION CHANGES THE BLOCK AND NOTHING INSIDE IT. A 3.0 block carries a four-byte
 *     checksum that a 2.1 block does not, so the same header is 178 bytes in one and 174 in the
 *     other, identical up to those four.
 *
 * Output:
 *
 *     header\t<version>\t<rn>\t<ap>\t<rr>\t<tags>\t<block hex>
 *     back\t<version>\t<rn>\t<ap>\t<rr>\t<tags>\t<rewritten: same, or its hex>\t<tag encoding map>
 *     section\t<name>\t<hex>
 *     err\t<what>\t<detail>\t<class>\t<message>
 *
 * Usage: CramCompressionHeaderDump
 */

import htsjdk.samtools.cram.common.CRAMVersion;
import htsjdk.samtools.cram.common.CramVersions;
import htsjdk.samtools.cram.structure.CompressionHeader;
import htsjdk.samtools.cram.structure.CompressionHeaderEncodingMap;
import htsjdk.samtools.cram.structure.CRAMEncodingStrategy;
import htsjdk.samtools.cram.structure.EncodingDescriptor;
import htsjdk.samtools.cram.structure.SubstitutionMatrix;
import htsjdk.samtools.cram.structure.block.Block;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.util.Map;
import java.util.StringJoiner;

public class CramCompressionHeaderDump {

    /** A substitution matrix of five bytes, one per base, distinct so a swap would show. */
    static final byte[] MATRIX = new byte[] {0x1b, 0x2d, 0x39, 0x4e, 0x63};

    public static void main(final String[] args) {
        System.out.println("# CramCompressionHeaderDump: three maps inside one raw block");

        for (final CRAMVersion version : new CRAMVersion[] {CramVersions.CRAM_v2_1,
                CramVersions.CRAM_v3}) {
            header(version, true, true, true, dictionary("OQZ", "XAZ"));
            header(version, false, false, false, dictionary("OQZ", "XAZ"));
        }
        header(CramVersions.CRAM_v3, true, false, true, dictionary());
        header(CramVersions.CRAM_v3, false, true, false, dictionary("MDZ"));

        // The three sections on their own, taken out of the block's content so a port can check
        // each against the suite that already covers it.
        sections(CramVersions.CRAM_v3, true, true, true, dictionary("OQZ", "XAZ"));

        // What it refuses.
        errWrongBlockType();
        errNoSubstitutionMatrix();
        errNoTagDictionary();
        errUnknownKey();
    }

    static byte[][][] dictionary(final String... tagSets) {
        final byte[][][] dictionary = new byte[tagSets.length][][];
        for (int i = 0; i < tagSets.length; i++) {
            final byte[] bytes = tagSets[i].getBytes();
            final byte[][] set = new byte[bytes.length / 3][];
            for (int j = 0; j < set.length; j++) {
                set[j] = new byte[] {bytes[j * 3], bytes[j * 3 + 1], bytes[j * 3 + 2]};
            }
            dictionary[i] = set;
        }
        return dictionary;
    }

    static CompressionHeader build(final boolean readNames, final boolean apDelta,
            final boolean referenceRequired, final byte[][][] tagDictionary) {
        final CompressionHeader header = new CompressionHeader(
                new CompressionHeaderEncodingMap(new CRAMEncodingStrategy()), apDelta, readNames,
                referenceRequired);
        header.setSubstitutionMatrix(new SubstitutionMatrix(MATRIX));
        header.setTagIdDictionary(tagDictionary);
        return header;
    }

    static void header(final CRAMVersion version, final boolean readNames, final boolean apDelta,
            final boolean referenceRequired, final byte[][][] tagDictionary) {
        final CompressionHeader header = build(readNames, apDelta, referenceRequired,
                tagDictionary);
        final ByteArrayOutputStream out = new ByteArrayOutputStream();
        header.write(version, out);
        final byte[] block = out.toByteArray();
        System.out.printf("header\t%s\t%b\t%b\t%b\t%s\t%s%n", version, readNames, apDelta,
                referenceRequired, tags(tagDictionary), hex(block));

        final CompressionHeader back = new CompressionHeader(version,
                new ByteArrayInputStream(block));
        // Written again from what was read, which is the property a byte-identical port needs and
        // the one place a lost field would show without being looked for.
        final ByteArrayOutputStream again = new ByteArrayOutputStream();
        back.write(version, again);
        System.out.printf("back\t%s\t%b\t%b\t%b\t%s\t%s\t%s%n", version,
                back.isPreserveReadNames(), back.isAPDelta(), back.isReferenceRequired(),
                tags(back.getTagIDDictionary()),
                java.util.Arrays.equals(block, again.toByteArray()) ? "same" : hex(again.toByteArray()),
                tagEncodings(back));
    }

    /** The block's content, split at the boundaries the three maps' own length prefixes give. */
    static void sections(final CRAMVersion version, final boolean readNames, final boolean apDelta,
            final boolean referenceRequired, final byte[][][] tagDictionary) {
        final ByteArrayOutputStream out = new ByteArrayOutputStream();
        build(readNames, apDelta, referenceRequired, tagDictionary).write(version, out);
        final Block block = Block.read(version, new ByteArrayInputStream(out.toByteArray()));
        final byte[] content = block.getRawContent();
        System.out.printf("section\tcontent\t%s%n", hex(content));
        System.out.printf("section\tcontent-type\t%s%n", block.getContentType().name());
        System.out.printf("section\tcompression\t%s%n", block.getCompressionMethod().name());
    }

    static void errWrongBlockType() {
        try {
            final ByteArrayOutputStream out = new ByteArrayOutputStream();
            Block.createRawSliceHeaderBlock(new byte[] {1, 2, 3}).write(CramVersions.CRAM_v3, out);
            new CompressionHeader(CramVersions.CRAM_v3, new ByteArrayInputStream(
                    out.toByteArray()));
            System.out.printf("err\twrong-block-type\tMAPPED_SLICE_HEADER\t-\t-%n");
        } catch (final Throwable t) {
            System.out.printf("err\twrong-block-type\tMAPPED_SLICE_HEADER\t%s\t%s%n",
                    t.getClass().getSimpleName(), String.valueOf(t.getMessage()));
        }
    }

    /** A preservation map holding everything but SM, wrapped and read back. */
    static void errNoSubstitutionMatrix() {
        readRaw("no-substitution-matrix", withKeys(true, true, true, false, true));
    }

    static void errNoTagDictionary() {
        readRaw("no-tag-dictionary", withKeys(true, true, true, true, false));
    }

    static void errUnknownKey() {
        readRaw("unknown-key", withUnknownKey());
    }

    /** Wrap a compression header's content in a block and read it back. */
    static void readRaw(final String what, final byte[] content) {
        try {
            final ByteArrayOutputStream out = new ByteArrayOutputStream();
            Block.createRawCompressionHeaderBlock(content).write(CramVersions.CRAM_v3, out);
            new CompressionHeader(CramVersions.CRAM_v3, new ByteArrayInputStream(
                    out.toByteArray()));
            System.out.printf("err\t%s\t%s\t-\t-%n", what, hex(content));
        } catch (final Throwable t) {
            System.out.printf("err\t%s\t%s\t%s\t%s%n", what, hex(content),
                    t.getClass().getSimpleName(), String.valueOf(t.getMessage()));
        }
    }

    /** A compression header's content with the chosen preservation keys and nothing after it. */
    static byte[] withKeys(final boolean rn, final boolean ap, final boolean rr, final boolean sm,
            final boolean td) {
        final ByteArrayOutputStream map = new ByteArrayOutputStream();
        int count = 0;
        if (rn) {
            map.writeBytes("RN".getBytes());
            map.write(1);
            count++;
        }
        if (ap) {
            map.writeBytes("AP".getBytes());
            map.write(1);
            count++;
        }
        if (rr) {
            map.writeBytes("RR".getBytes());
            map.write(1);
            count++;
        }
        if (sm) {
            map.writeBytes("SM".getBytes());
            map.writeBytes(MATRIX);
            count++;
        }
        if (td) {
            map.writeBytes("TD".getBytes());
            map.write(4);
            map.writeBytes(new byte[] {'O', 'Q', 'Z', 0});
            count++;
        }
        return prefixed(count, map.toByteArray());
    }

    static byte[] withUnknownKey() {
        final ByteArrayOutputStream map = new ByteArrayOutputStream();
        map.writeBytes("ZZ".getBytes());
        map.write(1);
        return prefixed(1, map.toByteArray());
    }

    /** The preservation map, behind its entry count and its byte length. */
    static byte[] prefixed(final int count, final byte[] entries) {
        final ByteArrayOutputStream withCount = new ByteArrayOutputStream();
        htsjdk.samtools.cram.io.ITF8.writeUnsignedITF8(count, withCount);
        withCount.writeBytes(entries);
        final byte[] map = withCount.toByteArray();

        final ByteArrayOutputStream out = new ByteArrayOutputStream();
        htsjdk.samtools.cram.io.ITF8.writeUnsignedITF8(map.length, out);
        out.writeBytes(map);
        return out.toByteArray();
    }

    static String tags(final byte[][][] dictionary) {
        if (dictionary == null) {
            return "null";
        }
        final StringJoiner outer = new StringJoiner(";");
        for (final byte[][] set : dictionary) {
            final StringBuilder builder = new StringBuilder();
            for (final byte[] tag : set) {
                builder.append(new String(tag));
            }
            outer.add(builder.length() == 0 ? "." : builder.toString());
        }
        return outer.length() == 0 ? "-" : outer.toString();
    }

    static String tagEncodings(final CompressionHeader header) {
        final StringJoiner joiner = new StringJoiner(";");
        for (final Map.Entry<Integer, EncodingDescriptor> entry
                : header.getTagEncodingMap().entrySet()) {
            joiner.add(entry.getKey() + "=" + entry.getValue().getEncodingID().name());
        }
        return joiner.length() == 0 ? "-" : joiner.toString();
    }

    static String hex(final byte[] bytes) {
        final StringBuilder builder = new StringBuilder(bytes.length * 2);
        for (final byte value : bytes) {
            builder.append(String.format("%02x", value));
        }
        return builder.length() == 0 ? "-" : builder.toString();
    }
}
