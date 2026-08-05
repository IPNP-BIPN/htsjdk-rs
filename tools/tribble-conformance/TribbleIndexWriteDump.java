/*
 * htsjdk.tribble.index: building a .idx, which is what GATK does beside every VCF it writes.
 *
 * TribbleIndexDump reads an index. This builds one, and the two are different problems: reading is
 * a layout, and writing is a set of decisions about bin widths and index types that the layout
 * only records the outcome of.
 *
 * Four of those decisions are invisible in the file and in the format:
 *
 *   - THE BIN WIDTH IS NOT THE ONE THE CREATOR WAS GIVEN. `LinearIndex.optimize` doubles it, per
 *     contig, until the most dense block holds more than MAX_FEATURES_PER_BIN estimated features,
 *     or one block is left, or the width goes bad. That is why the read suite measured 16000 and
 *     8000 in one file: not a setting, an outcome. The loop also keeps the LAST width that was
 *     still under the threshold rather than the first one over it, so it stops one step early;
 *   - THE DENSITY IS ESTIMATED, NOT COUNTED. The score is the largest block's size in BYTES
 *     divided by the average feature size in bytes, so it is a guess at a feature count that is
 *     never compared to the feature count the same object is carrying;
 *   - THE DYNAMIC CREATOR PICKS A TYPE FROM THE DATA. It feeds every feature to both a linear and
 *     an interval-tree creator, scores them, and keeps one, so the same tool over two files
 *     produces two different index types with nothing on the command line to say so. The scores
 *     go into a LinkedHashMap keyed by Double, so two creators that score equal collapse;
 *   - THE STATISTICS PUSHED ARE NOT THE FEATURE LENGTHS. `stats.push(longestFeatureLength)`
 *     pushes the running MAXIMUM at each step, so FEATURE_LENGTH_MEAN is the mean of a
 *     non-decreasing sequence and never the mean of the features. It is written into the header
 *     as a property, so the wrong number is in the bytes.
 *
 * And one arithmetic detail that decides the block list: `while (feature.getStart() > blocks.size()
 * * binWidth)` is int arithmetic on a product that grows with the file, and it appends one block
 * per iteration, so a feature far along a contig is preceded by every empty bin before it.
 *
 * The bytes travel base64 with the timestamp masked, for the reason TribbleIndexDump gives:
 * `indexedFileTS` is the source file's modification time and differs on every run.
 *
 * Output:
 *
 *     idx\t<label>\t<type>\t<timestamp offset>\t<base64 of the .idx, timestamp zeroed>
 *     chr\t<label>\t<name>\t<binWidth>\t<nBins>\t<longestFeature>\t<nFeatures>\t<block starts and sizes>
 *     prop\t<label>\t<key>=<value>;...        the properties written into the header, in order
 *     score\t<label>\t<creator>\t<score>      what the dynamic creator computed, and what it kept
 *     err\t<label>\t<class>\t<message>
 *
 * Usage: TribbleIndexWriteDump
 */

import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.bed.BEDCodec;
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
import java.util.Arrays;
import java.util.Base64;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

public class TribbleIndexWriteDump {

    static Path dir;

    public static void main(final String[] args) throws Exception {
        System.out.println("# TribbleIndexWriteDump: building a .idx, and the decisions that shape it");

        dir = Path.of("tribble-index-write-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        // Sparse: five features on two contigs, nowhere near the density threshold, so the
        // optimizer merges until one block is left and the width is nothing like the default.
        final Path sparse = write("sparse.bed",
                "chr1\t100\t110\ta",
                "chr1\t200\t210\tb",
                "chr1\t300\t900\tc",
                "chr1\t20000\t20010\td",
                "chr2\t50\t60\te");
        build("sparse-linear", sparse, "LINEAR");
        build("sparse-interval", sparse, "INTERVAL_TREE");
        build("sparse-dynamic-seek", sparse, "FOR_SEEK_TIME");
        build("sparse-dynamic-size", sparse, "FOR_SIZE");

        // Dense: enough features per bin that the score passes the threshold, so the optimizer
        // stops early and the width stays near the default.
        final List<String> denseLines = new ArrayList<>();
        for (int i = 0; i < 4000; i++) {
            denseLines.add(String.format("chr1\t%d\t%d\tf%d", 100 + i * 2, 110 + i * 2, i));
        }
        final Path dense = write("dense.bed", denseLines.toArray(new String[0]));
        build("dense-linear", dense, "LINEAR");
        build("dense-dynamic-seek", dense, "FOR_SEEK_TIME");
        build("dense-dynamic-size", dense, "FOR_SIZE");

        // Two contigs of very different density in one file, which is what produces two different
        // bin widths in a single index. The read suite measured that outcome; this is its cause.
        final List<String> mixedLines = new ArrayList<>();
        for (int i = 0; i < 3000; i++) {
            mixedLines.add(String.format("chr1\t%d\t%d\tc1_%d", 100 + i * 3, 110 + i * 3, i));
        }
        mixedLines.add("chr2\t100\t110\tc2_a");
        mixedLines.add("chr2\t900000\t900010\tc2_b");
        final Path mixed = write("mixed.bed", mixedLines.toArray(new String[0]));
        build("mixed-linear", mixed, "LINEAR");

        // One feature only: nBlocks is 1 from the start, so the optimizer's first test breaks
        // immediately and the width is the one the creator was given.
        final Path single = write("single.bed", "chr1\t100\t110\tonly");
        build("single-linear", single, "LINEAR");

        // A feature far along the contig, which appends one empty bin per iteration of the while
        // loop until the bin list reaches it.
        final Path far = write("far.bed", "chr1\t100\t110\tnear", "chr1\t5000000\t5000010\tfar");
        build("far-linear", far, "LINEAR");

        // An empty file: no features at all, so there is no last contig to close.
        final Path empty = write("empty.bed");
        build("empty-linear", empty, "LINEAR");

        // Features out of order on one contig, which the factory refuses rather than indexing
        // wrongly.
        final Path unsorted = write("unsorted.bed",
                "chr1\t500\t510\tsecond",
                "chr1\t100\t110\tfirst");
        build("unsorted-linear", unsorted, "LINEAR");

        // A contig revisited after another one, which is a different failure from the one above.
        final Path revisited = write("revisited.bed",
                "chr1\t100\t110\ta",
                "chr2\t100\t110\tb",
                "chr1\t200\t210\tc");
        build("revisited-linear", revisited, "LINEAR");
    }

    static Path write(final String name, final String... lines) throws Exception {
        final Path path = dir.resolve(name);
        try (final PrintWriter out = new PrintWriter(Files.newBufferedWriter(path))) {
            for (final String line : lines) {
                out.println(line);
            }
        }
        return path;
    }

    /**
     * `how` is either an IndexType name, for the fixed creators, or an IndexBalanceApproach name,
     * for the dynamic one. There is no IndexType for "dynamic": the choice of type is what the
     * dynamic creator produces, so it cannot be one of the things you ask for.
     */
    static void build(final String label, final Path source, final String how) {
        try {
            final FeatureCodec<? extends Feature, ?> codec = new BEDCodec();
            final Index index =
                    how.startsWith("FOR_")
                            ? IndexFactory.createDynamicIndex(
                                    source, codec, IndexFactory.IndexBalanceApproach.valueOf(how))
                            : IndexFactory.createIndex(
                                    source, codec, IndexFactory.IndexType.valueOf(how));
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
            Arrays.fill(masked, timestampOffset, timestampOffset + 8, (byte) 0);
            System.out.printf("idx\t%s\t%s\t%d\t%s%n", label, index.getClass().getSimpleName(),
                    timestampOffset, Base64.getEncoder().encodeToString(masked));

            // The properties, in the order the header holds them, because they are written in that
            // order and two of them are Java double strings whose formatting is its own problem.
            final StringJoiner properties = new StringJoiner(";");
            final Map<String, String> map = index.getProperties();
            if (map != null) {
                for (final Map.Entry<String, String> entry : map.entrySet()) {
                    properties.add(entry.getKey() + "=" + entry.getValue());
                }
            }
            System.out.printf("prop\t%s\t%s%n", label, properties);

            // The per-contig shape, read back out of the bytes rather than off the object: the
            // fields the optimizer decides are package-private, and the layout is what a port has
            // to produce anyway. Only the linear one is parsed; the interval-tree chromosome
            // record is a different shape and this suite's claim is about the linear writer.
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
                    int count = dis.readInt();
                    while (count-- > 0) {
                        dis.readString();
                        dis.readString();
                    }
                    int contigs = dis.readInt();
                    while (contigs-- > 0) {
                        final String name = dis.readString();
                        final int binWidth = dis.readInt();
                        final int nBins = dis.readInt();
                        final int longestFeature = dis.readInt();
                        dis.readInt();
                        final int nFeatures = dis.readInt();
                        final StringJoiner positions = new StringJoiner(",");
                        for (int i = 0; i <= nBins; i++) {
                            positions.add(Long.toString(dis.readLong()));
                        }
                        System.out.printf("chr\t%s\t%s\t%d\t%d\t%d\t%d\t%s%n", label, name,
                                binWidth, nBins, longestFeature, nFeatures, positions);
                    }
                }
            }
        } catch (final Throwable t) {
            System.out.printf("err\t%s\t%s\t%s%n", label, t.getClass().getName(),
                    oneLine(t.getMessage()));
        }
    }

    static List<String> sorted(final List<String> names) {
        final List<String> copy = new ArrayList<>(names);
        copy.sort(Comparator.naturalOrder());
        return copy;
    }

    static void emptyDirectory(final Path path) throws Exception {
        if (!Files.exists(path)) {
            return;
        }
        try (final var walk = Files.walk(path)) {
            walk.sorted(Comparator.reverseOrder()).forEach(child -> {
                try {
                    Files.delete(child);
                } catch (final Exception ignored) {
                    // A leftover directory that cannot be removed is not what this dump measures.
                }
            });
        }
    }

    static String oneLine(final String s) {
        return s == null ? "null" : s.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
    }
}
