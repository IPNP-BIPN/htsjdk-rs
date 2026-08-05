/*
 * htsjdk.tribble.index: the .idx file, and the blocks a query resolves to.
 *
 * The last named consumer of H.2, and the item gatk-rs G1.3 named when it closed: the Tribble
 * index is what turns a Feature file into a random-access source rather than a linear read.
 *
 * The dump comes before the port on purpose. The index type identifiers are defined through a
 * circular reference — `LinearIndex.INDEX_TYPE` reads a field of the `IndexType` enum whose own
 * constructor is handed `LinearIndex.INDEX_TYPE` — so the values cannot be read reliably out of
 * the source. They are measured here instead.
 *
 * Three behaviours of the query that a port would otherwise have to guess at:
 *
 *   - THE LONGEST FEATURE IS SUBTRACTED FROM THE START. `adjustedPosition = max(start -
 *     longestFeature, 0)`, so a query gets blocks holding features that begin before the interval
 *     and reach into it. That is why the index records the longest feature per contig at all;
 *   - THE ANSWER IS ALWAYS ONE BLOCK OR NONE, never a list. Linear-index blocks are adjacent by
 *     definition, so the query merges from the first bin to the last into a single block;
 *   - AN EMPTY ANSWER HAS TWO CAUSES. Off the end of the bin list, and a merged size of zero.
 *
 * The bytes travel base64 because the on-disk layout is the thing being reproduced: the header's
 * strings are NUL-terminated rather than length-prefixed, and a chromosome's N blocks are written
 * as N+1 longs whose differences are the sizes.
 *
 * WITH ONE FIELD MASKED. `indexedFileTS` is the source file's modification time, so the raw bytes
 * differ on every run — measured, by running this dump twice and diffing. A golden built on them
 * would fail intermittently, or worse be "fixed" by regenerating it. The eight bytes are zeroed
 * and the offset is reported, so the rest of the layout stays under test and the one field that
 * cannot be stable is visibly absent rather than quietly poisoning the file.
 *
 * Output:
 *
 *     idx\t<label>\t<timestamp offset>\t<base64 of the .idx, timestamp zeroed>
 *     header\t<label>\t<magic>\t<type>\t<version>\t<flags>\t<properties>\t<path shape>\t<size shape>\t<md5 shape>
 *     chr\t<label>\t<name>\t<binWidth>\t<nBins>\t<longestFeature>\t<unused>\t<nFeatures>\t<block positions>   (linear)
 *     chr\t<label>\t<name>\t<count>\t<start,end,pos,size;...>                                     (interval tree)
 *     query\t<label>\t<contig>:<start>-<end>\t<start,size|...  or `none`>
 *
 * Usage: TribbleIndexDump
 */

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.index.Block;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.index.linear.LinearIndex;
import htsjdk.tribble.util.LittleEndianInputStream;

import java.io.BufferedInputStream;
import java.io.ByteArrayInputStream;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Base64;
import java.util.List;

public class TribbleIndexDump {

    public static void main(final String[] args) throws Exception {
        System.out.println("# TribbleIndexDump: the .idx bytes, and what a query resolves to");

        final Path dir = Path.of("tribble-index-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        // A small BED with features of two very different lengths on one contig, so the
        // longest-feature back-off is visible, plus a second contig.
        final Path bed = dir.resolve("features.bed");
        try (final PrintWriter out = new PrintWriter(Files.newBufferedWriter(bed))) {
            out.println("chr1\t100\t110\tshort1");
            out.println("chr1\t200\t210\tshort2");
            // A long one, so longestFeature is not the common case.
            out.println("chr1\t300\t900\tlong1");
            out.println("chr1\t20000\t20010\tfar");
            out.println("chr2\t50\t60\tother");
        }

        index("linear-bed", bed, IndexFactory.IndexType.LINEAR);
        index("interval-bed", bed, IndexFactory.IndexType.INTERVAL_TREE);

        // A second file whose features all sit inside one bin, so the merged block is one bin wide.
        final Path dense = dir.resolve("dense.bed");
        try (final PrintWriter out = new PrintWriter(Files.newBufferedWriter(dense))) {
            for (int i = 0; i < 20; i++) {
                out.printf("chr1\t%d\t%d\td%d%n", 100 + i, 110 + i, i);
            }
        }
        index("linear-dense", dense, IndexFactory.IndexType.LINEAR);
    }

    static void index(final String label, final Path source, final IndexFactory.IndexType type)
            throws Exception {
        final FeatureCodec<? extends Feature, ?> codec = new BEDCodec();
        final Index index = IndexFactory.createIndex(source, codec, type);
        final Path idx = source.resolveSibling(source.getFileName() + "." + label + ".idx");
        index.write(idx);

        final byte[] bytes = Files.readAllBytes(idx);
        // magic(4) + type(4) + version(4) + NUL-terminated path + size(8), then the timestamp.
        int cursor = 12;
        while (bytes[cursor] != 0) {
            cursor++;
        }
        final int timestampOffset = cursor + 1 + 8;
        final byte[] masked = bytes.clone();
        java.util.Arrays.fill(masked, timestampOffset, timestampOffset + 8, (byte) 0);
        System.out.printf("idx\t%s\t%d\t%s%n", label, timestampOffset,
                Base64.getEncoder().encodeToString(masked));

        // The header, read back from the bytes rather than from the object, because the on-disk
        // layout is what a port has to reproduce.
        try (final LittleEndianInputStream dis =
                new LittleEndianInputStream(new BufferedInputStream(new ByteArrayInputStream(bytes)))) {
            final int magic = dis.readInt();
            final int typeId = dis.readInt();
            final int version = dis.readInt();
            final String path = dis.readString();
            final long size = dis.readLong();
            final long timestamp = dis.readLong();
            final String md5 = dis.readString();
            final int flags = dis.readInt();
            final int properties = dis.readInt();
            // The path, size and timestamp are of the run that produced the file, so they are
            // reported as shapes rather than values: a golden carrying them would be unstable.
            System.out.printf("header\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s%n", label, magic, typeId,
                    version, flags, properties,
                    path.endsWith(source.getFileName().toString()) ? "path-ends-with-name" : path,
                    size > 0 ? "size-positive" : "size-" + size,
                    md5.isEmpty() ? "md5-empty" : "md5-present");
        }

        // The per-contig numbers a linear query depends on, read back out of the bytes. Only the
        // linear layout is parsed here: the interval-tree chromosome record is a different shape,
        // and this suite's claim is about the linear index the query walks.
        if (index instanceof LinearIndex) {
            try (final LittleEndianInputStream dis = new LittleEndianInputStream(
                    new BufferedInputStream(new ByteArrayInputStream(bytes)))) {
                dis.readInt();
                dis.readInt();
                dis.readInt();
                dis.readString();
                dis.readLong();
                dis.readLong();
                dis.readString();
                dis.readInt();
                int properties = dis.readInt();
                while (properties-- > 0) {
                    dis.readString();
                    dis.readString();
                }
                int contigs = dis.readInt();
                while (contigs-- > 0) {
                    final String name = dis.readString();
                    final int binWidth = dis.readInt();
                    final int nBins = dis.readInt();
                    final int longestFeature = dis.readInt();
                    final int unused = dis.readInt();
                    final int nFeatures = dis.readInt();
                    final StringBuilder positions = new StringBuilder();
                    for (int i = 0; i <= nBins; i++) {
                        if (i > 0) {
                            positions.append(',');
                        }
                        positions.append(dis.readLong());
                    }
                    // N blocks are written as N+1 longs; the sizes are the differences.
                    System.out.printf("chr\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s%n", label, name,
                            binWidth, nBins, longestFeature, unused, nFeatures, positions);
                }
            }
        } else {
            // The interval-tree chromosome record: a name, a count, then that many
            // (start, end, position, size). Sizes are STORED here, unlike the linear layout where
            // they are the differences between consecutive positions.
            try (final LittleEndianInputStream dis = new LittleEndianInputStream(
                    new BufferedInputStream(new ByteArrayInputStream(bytes)))) {
                dis.readInt();
                dis.readInt();
                dis.readInt();
                dis.readString();
                dis.readLong();
                dis.readLong();
                dis.readString();
                dis.readInt();
                int properties = dis.readInt();
                while (properties-- > 0) {
                    dis.readString();
                    dis.readString();
                }
                int contigs = dis.readInt();
                while (contigs-- > 0) {
                    final String name = dis.readString();
                    final int count = dis.readInt();
                    final StringBuilder intervals = new StringBuilder();
                    for (int i = 0; i < count; i++) {
                        if (i > 0) {
                            intervals.append(';');
                        }
                        intervals.append(dis.readInt()).append(',').append(dis.readInt())
                                .append(',').append(dis.readLong()).append(',')
                                .append(dis.readInt());
                    }
                    System.out.printf("chr\t%s\t%s\t%d\t%s%n", label, name, count, intervals);
                }
            }
        }

        // The queries. The first three walk the longest-feature back-off, the fourth is off the
        // end of the bin list, and the last is a contig the index does not hold.
        query(label, index, "chr1", 100, 120);
        query(label, index, "chr1", 350, 360);
        // Starts after the long feature but within its reach, so the back-off decides.
        query(label, index, "chr1", 880, 890);
        query(label, index, "chr1", 1000000, 1000010);
        query(label, index, "chr2", 50, 60);
        query(label, index, "chrX", 1, 10);
    }

    static void query(final String label, final Index index, final String contig, final int start,
                      final int end) {
        final String interval = contig + ":" + start + "-" + end;
        List<Block> blocks;
        try {
            blocks = index.getBlocks(contig, start, end);
        } catch (final Exception e) {
            System.out.printf("query\t%s\t%s\tE:%s%n", label, interval, e.getClass().getName());
            return;
        }
        if (blocks == null || blocks.isEmpty()) {
            System.out.printf("query\t%s\t%s\tnone%n", label, interval);
            return;
        }
        final List<String> rendered = new ArrayList<>();
        for (final Block block : blocks) {
            rendered.add(block.getStartPosition() + "," + block.getSize());
        }
        System.out.printf("query\t%s\t%s\t%s%n", label, interval, String.join("|", rendered));
    }

    static void emptyDirectory(final Path dir) throws Exception {
        if (!Files.isDirectory(dir)) {
            return;
        }
        try (final var entries = Files.list(dir)) {
            for (final Path entry : entries.toList()) {
                Files.deleteIfExists(entry);
            }
        }
    }
}
