/*
 * SliceBlocks: the core block and the external blocks of one slice, in the order they are written.
 *
 * A slice's data is one core block of bits and any number of external blocks of bytes, each named
 * by a content id. This is what holds them, and it is the join between the codecs and the file.
 *
 * Five things here are decisions rather than layout.
 *
 *   - THE ORDER IS BY CONTENT ID AND NOT BY INSERTION. The externals live in a TreeMap, so writing
 *     is core first and then ascending content id, whatever order they were added in;
 *   - THE READER TAKES A BLOCK COUNT AND NOT AN ORDER, so it accepts them in any order and a
 *     stream that puts its core block last is read exactly as one that puts it first;
 *   - A DUPLICATE CONTENT ID IS FATAL, and its message names the id and the type of both the new
 *     block and the one already there;
 *   - A BLOCK OF ANY OTHER TYPE IS FATAL, with a message naming the type it found;
 *   - A STREAM WITHOUT A CORE BLOCK IS FATAL, and that is checked after all the blocks have been
 *     read rather than while reading them.
 *
 * Output:
 *
 *     write\t<version>\t<added order>\t<written order>\t<hex>
 *     read\t<version>\t<blocks>\t<core size>\t<external ids>
 *     err\t<what>\t<detail>\t<class>\t<message>
 *
 * Usage: CramSliceBlocksDump
 */

import htsjdk.samtools.cram.common.CRAMVersion;
import htsjdk.samtools.cram.common.CramVersions;
import htsjdk.samtools.cram.structure.CompressorCache;
import htsjdk.samtools.cram.structure.SliceBlocks;
import htsjdk.samtools.cram.structure.block.Block;
import htsjdk.samtools.cram.structure.block.BlockCompressionMethod;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.StringJoiner;

public class CramSliceBlocksDump {

    public static void main(final String[] args) {
        System.out.println("# CramSliceBlocksDump: the core block and the externals of one slice");

        for (final CRAMVersion version : new CRAMVersion[] {CramVersions.CRAM_v2_1,
                CramVersions.CRAM_v3}) {
            write(version, new int[] {1, 2, 3});
            write(version, new int[] {3, 2, 1});
        }
        write(CramVersions.CRAM_v3, new int[] {7});
        write(CramVersions.CRAM_v3, new int[] {300, 2, 128});
        write(CramVersions.CRAM_v3, new int[] {0, 1});

        // Read a stream back, including one whose core block comes last.
        readBack(CramVersions.CRAM_v3, false);
        readBack(CramVersions.CRAM_v3, true);
        readBack(CramVersions.CRAM_v2_1, false);

        // What it refuses.
        errDuplicate();
        errNoCore();
        errWrongType();
    }

    /** A core block and one external block per content id, added in the order given. */
    static void write(final CRAMVersion version, final int[] contentIds) {
        final Block core = Block.createRawCoreDataBlock(new byte[] {(byte) 0xAA, (byte) 0x55});
        final List<Block> externals = new ArrayList<>();
        for (final int contentId : contentIds) {
            externals.add(external(contentId));
        }
        final SliceBlocks blocks = new SliceBlocks(core, externals);

        final ByteArrayOutputStream out = new ByteArrayOutputStream();
        blocks.writeBlocks(version, out);
        System.out.printf("write\t%s\t%s\t%s\t%s%n", version, ids(contentIds),
                idList(blocks.getExternalContentIDs()), hex(out.toByteArray()));
    }

    /** One external block whose content is its own id, so a swap between two would show. */
    static Block external(final int contentId) {
        final byte[] content = new byte[] {(byte) contentId, (byte) (contentId + 1)};
        return Block.createExternalBlock(BlockCompressionMethod.RAW, contentId, content,
                content.length);
    }

    static void readBack(final CRAMVersion version, final boolean coreLast) {
        final Block core = Block.createRawCoreDataBlock(new byte[] {1, 2, 3, 4});
        final List<Block> externals = Arrays.asList(external(2), external(1));

        final ByteArrayOutputStream out = new ByteArrayOutputStream();
        if (!coreLast) {
            core.write(version, out);
        }
        for (final Block block : externals) {
            block.write(version, out);
        }
        if (coreLast) {
            core.write(version, out);
        }

        final SliceBlocks blocks = new SliceBlocks(version, externals.size() + 1,
                new ByteArrayInputStream(out.toByteArray()));
        System.out.printf("read\t%s\t%s\t%d\t%s%n", version, coreLast ? "core-last" : "core-first",
                blocks.getCoreBlock().getUncompressedContentSize(),
                idList(blocks.getExternalContentIDs()));
    }

    /** Two external blocks with the same content id. */
    static void errDuplicate() {
        try {
            new SliceBlocks(Block.createRawCoreDataBlock(new byte[] {1}),
                    Arrays.asList(external(4), external(4)));
            System.out.printf("err\tduplicate\tid=4\t-\t-%n");
        } catch (final Throwable t) {
            System.out.printf("err\tduplicate\tid=4\t%s\t%s%n", t.getClass().getSimpleName(),
                    String.valueOf(t.getMessage()));
        }
    }

    /** A stream of external blocks only, which is refused after all of them have been read. */
    static void errNoCore() {
        try {
            final ByteArrayOutputStream out = new ByteArrayOutputStream();
            external(1).write(CramVersions.CRAM_v3, out);
            external(2).write(CramVersions.CRAM_v3, out);
            new SliceBlocks(CramVersions.CRAM_v3, 2, new ByteArrayInputStream(out.toByteArray()));
            System.out.printf("err\tno-core\t2 external blocks\t-\t-%n");
        } catch (final Throwable t) {
            System.out.printf("err\tno-core\t2 external blocks\t%s\t%s%n",
                    t.getClass().getSimpleName(), String.valueOf(t.getMessage()));
        }
    }

    /** A compression header block in a slice's stream. */
    static void errWrongType() {
        try {
            final ByteArrayOutputStream out = new ByteArrayOutputStream();
            Block.createRawCompressionHeaderBlock(new byte[] {1}).write(CramVersions.CRAM_v3, out);
            new SliceBlocks(CramVersions.CRAM_v3, 1, new ByteArrayInputStream(out.toByteArray()));
            System.out.printf("err\twrong-type\tCOMPRESSION_HEADER\t-\t-%n");
        } catch (final Throwable t) {
            System.out.printf("err\twrong-type\tCOMPRESSION_HEADER\t%s\t%s%n",
                    t.getClass().getSimpleName(), String.valueOf(t.getMessage()));
        }
    }

    static String ids(final int[] contentIds) {
        final StringJoiner joiner = new StringJoiner(",");
        for (final int contentId : contentIds) {
            joiner.add(Integer.toString(contentId));
        }
        return joiner.toString();
    }

    static String idList(final List<Integer> contentIds) {
        final StringJoiner joiner = new StringJoiner(",");
        for (final Integer contentId : contentIds) {
            joiner.add(contentId.toString());
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
