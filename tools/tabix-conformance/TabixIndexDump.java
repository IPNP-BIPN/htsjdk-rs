import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.tribble.Feature;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndex;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
import htsjdk.tribble.util.LittleEndianOutputStream;

import java.io.ByteArrayOutputStream;
import java.util.ArrayList;
import java.util.List;

/**
 * Dumps the `.tbi` `TabixIndexCreator` builds for a stream of features, as bytes.
 *
 * The creator is fed features and their file positions directly rather than through a codec: what
 * is measured is the index, and a VCF reader in front of it would measure the reader as well. The
 * positions are the virtual offsets a BGZF stream would have produced, chosen so that the cases
 * cover the places the layout is decided rather than described.
 *
 * Six behaviours this is built to catch.
 *
 *   - A FEATURE IS INDEXED ONE FEATURE LATE: a chunk needs both ends and the end of one feature is
 *     the START of the next, so nothing is indexed until the next arrives and the last one waits
 *     for finalizeIndex's own position;
 *   - THE BIN IS COMPUTED FROM A ZERO-BASED HALF-OPEN REGION, regionToBin(start - 1, end), because
 *     TabixFeature.getIndexingBin returns null;
 *   - A FEATURE WITH NO END IS ONE BASE for the bin AND shifts a window boundary for the linear
 *     index, which are two different rules reached by the same feature;
 *   - THE LINEAR INDEX FILLS ITS GAPS with the last non-empty offset, so an empty 16K window is
 *     not zero and not -1;
 *   - A SEQUENCE WITH NO FEATURES WRITES A ZERO BIN COUNT and nothing else, which is how a
 *     truncated tail differs from an absent one;
 *   - AND THE NAME BLOCK COUNTS ITS NULL TERMINATORS, so its declared size is the names plus one
 *     byte each.
 *
 * Output:
 *
 *     body\t&lt;case&gt;\t&lt;the little-endian index, hex&gt;
 *     file\t&lt;case&gt;\t&lt;the same index inside a BGZF stream, hex&gt;
 *     error\t&lt;case&gt;\t&lt;exception class&gt;: &lt;message&gt;
 *
 * Usage: TabixIndexDump
 */
public class TabixIndexDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t').append(payload).append('\n');
    }

    static String hex(final byte[] bytes) {
        final StringBuilder out = new StringBuilder(bytes.length * 2);
        for (final byte b : bytes) out.append(String.format("%02x", b));
        return out.toString();
    }

    /** One feature and the position it starts at. */
    record Row(String contig, int start, int end, long position) implements Feature {
        @Override public String getContig() { return contig; }
        @Override public int getStart() { return start; }
        @Override public int getEnd() { return end; }
    }

    static List<Row> rows(final Row... rows) {
        return new ArrayList<>(List.of(rows));
    }

    /**
     * Build one index and emit it twice: the body the port composes, and the file that body
     * becomes once it is block compressed, which is what lands beside the feature file.
     */
    static void run(final String name, final TabixFormat format, final long finalPosition,
                    final List<Row> features) {
        final TabixIndexCreator creator = new TabixIndexCreator(format);
        try {
            for (final Row row : features) {
                creator.addFeature(row, row.position());
            }
            final Index index = creator.finalizeIndex(finalPosition);
            final ByteArrayOutputStream body = new ByteArrayOutputStream();
            try (final LittleEndianOutputStream los = new LittleEndianOutputStream(body)) {
                ((TabixIndex) index).write(los);
            }
            emit("body", name, hex(body.toByteArray()));

            final ByteArrayOutputStream file = new ByteArrayOutputStream();
            try (final LittleEndianOutputStream los =
                         new LittleEndianOutputStream(new BlockCompressedOutputStream(file, (java.io.File) null))) {
                ((TabixIndex) index).write(los);
            }
            emit("file", name, hex(file.toByteArray()));
        } catch (final Exception e) {
            emit("error", name, e.getClass().getName() + ": " + e.getMessage());
        }
    }

    public static void main(final String[] args) throws Exception {
        // One feature on one contig: the smallest index that is not empty, and the case where the
        // only chunk is closed by finalizeIndex rather than by a following feature.
        run("one-feature", TabixFormat.VCF, 4096L, rows(new Row("chr1", 100, 100, 512L)));

        // Two features in the same 16K window: one bin, one chunk pair, and the linear index's
        // first entry is the EARLIER offset.
        run("two-in-one-window", TabixFormat.VCF, 8192L,
                rows(new Row("chr1", 100, 100, 512L), new Row("chr1", 200, 200, 1024L)));

        // Features far enough apart to leave empty windows between them, which is where the
        // linear index's gap filling shows.
        run("gap-in-the-linear-index", TabixFormat.VCF, 1L << 20,
                rows(new Row("chr1", 1, 1, 512L), new Row("chr1", 200_000, 200_000, 65536L)));

        // A feature spanning several windows, so one chunk start reaches more than one entry.
        run("spans-many-windows", TabixFormat.VCF, 1L << 20,
                rows(new Row("chr1", 1, 100_000, 512L), new Row("chr1", 150_000, 150_000, 70000L)));

        // A feature whose end is unset, which is one base for the bin and a shifted window for the
        // linear index. 16385 sits one past a window boundary, which is where the two rules differ.
        run("unset-end", TabixFormat.VCF, 8192L,
                rows(new Row("chr1", 16385, 0, 512L), new Row("chr1", 16385, 0, 1024L)));

        // Two contigs, so the reference index advances and a second block is written.
        run("two-contigs", TabixFormat.VCF, 8192L,
                rows(new Row("chr1", 100, 100, 512L), new Row("chr2", 100, 100, 2048L)));

        // Bins at several levels: a feature confined to 16K, one to 128K, one to 1M and one that
        // reaches bin 0, so the ladder in regionToBin is exercised rather than assumed.
        run("every-bin-level", TabixFormat.VCF, 1L << 28,
                rows(
                        new Row("chr1", 1, 1000, 512L),
                        new Row("chr1", 20_000, 100_000, 1024L),
                        new Row("chr1", 200_000, 900_000, 2048L),
                        new Row("chr1", 1_000_000, 9_000_000, 4096L),
                        new Row("chr1", 10_000_000, 200_000_000, 8192L)));

        // The four other format specs, whose six header integers are the only difference.
        run("format-gff", TabixFormat.GFF, 4096L, rows(new Row("chr1", 100, 200, 512L)));
        run("format-bed", TabixFormat.BED, 4096L, rows(new Row("chr1", 100, 200, 512L)));
        run("format-sam", TabixFormat.SAM, 4096L, rows(new Row("chr1", 100, 200, 512L)));
        run("format-psltbl", TabixFormat.PSLTBL, 4096L, rows(new Row("chr1", 100, 200, 512L)));

        // A contig name long enough that the name block's declared size is worth reading, and one
        // with a character above ASCII, since StringUtil.stringToBytes takes the low byte.
        run("long-and-wide-names", TabixFormat.VCF, 8192L,
                rows(
                        new Row("a_very_long_contig_name_0123456789", 100, 100, 512L),
                        new Row("chré", 100, 100, 2048L)));

        // No features at all: the header, no names, and no blocks.
        run("no-features", TabixFormat.VCF, 0L, rows());

        // A virtual offset with a real block address, which is what a BGZF stream actually hands
        // the creator: the chunk is not a small integer and the linear entry is the same number.
        run("virtual-offsets", TabixFormat.VCF, (30000L << 16),
                rows(
                        new Row("chr1", 100, 100, (12345L << 16) | 678L),
                        new Row("chr1", 20_000, 20_000, (23456L << 16) | 90L)));

        // The three refusals, each of them an IllegalArgumentException.
        run("sequence-out-of-order", TabixFormat.VCF, 8192L,
                rows(
                        new Row("chr1", 100, 100, 512L),
                        new Row("chr2", 100, 100, 1024L),
                        new Row("chr1", 100, 100, 2048L)));
        run("features-out-of-order", TabixFormat.VCF, 8192L,
                rows(new Row("chr1", 200, 200, 512L), new Row("chr1", 100, 100, 1024L)));
        run("position-did-not-advance", TabixFormat.VCF, 8192L,
                rows(new Row("chr1", 100, 100, 512L), new Row("chr1", 200, 200, 512L)));
        // The same refusal reached by finalizeIndex rather than by the next feature.
        run("final-position-did-not-advance", TabixFormat.VCF, 512L,
                rows(new Row("chr1", 100, 100, 512L)));
        // Equal starts are NOT out of order: compareTo does not look at the end.
        run("equal-starts", TabixFormat.VCF, 8192L,
                rows(new Row("chr1", 100, 500, 512L), new Row("chr1", 100, 200, 1024L)));

        System.out.print(buf);
    }
}
