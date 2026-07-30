/*
 * BEDCodec, taken from the reference.
 *
 * This is what GATK's `-L regions.bed` runs through, so its coordinate convention is the one a
 * whole interval argument inherits, and five of its decisions are not what "a BED parser"
 * suggests:
 *
 *   - the start is SHIFTED by a constructor argument. StartOffset.ONE is the default, so a BED
 *     file's 0-based start becomes 1-based on the way in, and the feature carries the shifted
 *     number without carrying the shift;
 *   - a two-token line is a POINT: with no end column, end = start, which is one base at a
 *     coordinate the file never mentions;
 *   - the separator is "\t|( +)": one tab, or a run of spaces however long, and the limit is -1
 *     so a trailing separator produces an empty field rather than being dropped;
 *   - a blank line and a header line decode to NULL rather than to an error, and the header
 *     prefixes (#, track, browser) are matched on the raw line, so " track" is data;
 *   - a bad score returns the feature EARLY, so a line whose score is "." keeps its name and
 *     silently loses its strand, its colour and its exons even though those columns are present.
 *
 * Output:
 *
 *     bed\t<escaped line>\t<offset>\t<null|contig:start-end|name|score|strand|color|exons>
 *     bed\t<escaped line>\t<offset>\tE:<class>:<message>
 *     split\t<escaped line>\t<field count>\t<fields, pipe-separated and escaped>
 *     candecode\t<path>\t<true|false>
 *
 * Usage: BedCodecDump
 */

import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.bed.FullBEDFeature;

import java.util.List;
import java.util.StringJoiner;

public class BedCodecDump {

    /** Every line the codec is asked to decode, at both start offsets. */
    static final String[] LINES = {
        // The ordinary shapes.
        "chr1\t100\t200",
        "chr1\t0\t1",
        "chr1\t100\t200\tname",
        "chr1\t100\t200\tname\t500",
        "chr1\t100\t200\tname\t500\t+",
        "chr1\t100\t200\tname\t500\t-",
        "chr1\t100\t200\tname\t500\t.",
        // Two tokens: the end is the shifted start.
        "chr1\t100",
        "chr1\t0",
        // One token, which is not a feature at all.
        "chr1",
        "",
        // Blank and header lines, which are null rather than errors.
        "   ",
        "\t",
        "# a comment",
        "track name=whatever",
        "browser position chr1",
        // The prefixes are matched raw, so a leading space makes them data.
        " track name=whatever",
        // Separators: a run of spaces is one, a tab beside a space is two.
        "chr1 100 200",
        "chr1   100   200",
        "chr1 \t100\t200",
        "chr1\t100\t200\t",
        "chr1\t100\t200\t\t",
        // A quoted name, whose quotes are stripped wherever they are.
        "chr1\t100\t200\t\"na\"me\"",
        // A score that does not parse: the feature comes back early.
        "chr1\t100\t200\tname\t.\t+\t100\t200\t255,0,0\t2\t10,10\t0,50",
        "chr1\t100\t200\tname\tabc\t-",
        // A score Java parses and Rust does not, and the reverse.
        "chr1\t100\t200\tname\t1.5f\t+",
        "chr1\t100\t200\tname\tInfinity\t+",
        "chr1\t100\t200\tname\tinf\t+",
        // Colours: RGB, hex, a name, and something that is none of them.
        "chr1\t100\t200\tname\t500\t+\t100\t200\t255,0,0",
        "chr1\t100\t200\tname\t500\t+\t100\t200\t#00ff00",
        "chr1\t100\t200\tname\t500\t+\t100\t200\tred",
        "chr1\t100\t200\tname\t500\t+\t100\t200\tnonsense",
        "chr1\t100\t200\tname\t500\t+\t100\t200\t300,0,0",
        // Exons, on both strands, and with mismatched lists.
        "chr1\t100\t200\tname\t500\t+\t100\t200\t255,0,0\t2\t10,20\t0,50",
        "chr1\t100\t200\tname\t500\t-\t100\t200\t255,0,0\t2\t10,20\t0,50",
        "chr1\t100\t200\tname\t500\t+\t100\t200\t255,0,0\t3\t10,20\t0,50",
        "chr1\t100\t200\tname\t500\t+\t100\t200\t255,0,0\t0\t\t",
        // Numbers that do not parse, which throw rather than returning null.
        "chr1\tabc\t200",
        "chr1\t100\tabc",
        "chr1\t \t200",
        "chr1\t+100\t200",
        "chr1\t-100\t-50",
    };

    static final String[] PATHS = {
        "regions.bed", "regions.BED", "regions.bed.gz", "regions.BED.GZ", "regions.bed.bgz",
        "regions.bed.gz.gz", "regions.bedgraph", "regions.txt", "bed", ".bed", "regions.bed.tbi",
    };

    public static void main(final String[] args) {
        System.out.println("# BedCodecDump: BEDCodec, its split, and canDecode");

        for (final String line : LINES) {
            for (final BEDCodec.StartOffset offset : BEDCodec.StartOffset.values()) {
                emit(line, offset);
            }
            split(line);
        }

        final BEDCodec codec = new BEDCodec();
        for (final String path : PATHS) {
            System.out.printf("candecode\t%s\t%b%n", path, codec.canDecode(path));
        }
    }

    static void emit(final String line, final BEDCodec.StartOffset offset) {
        try {
            final BEDFeature feature = new BEDCodec(offset).decode(line);
            System.out.printf("bed\t%s\t%s\t%s%n", escape(line), offset, show(feature));
        } catch (final Exception | AssertionError e) {
            System.out.printf("bed\t%s\t%s\tE:%s:%s%n", escape(line), offset,
                    e.getClass().getName(),
                    e.getMessage() == null ? "" : e.getMessage().replace('\n', ' '));
        }
    }

    /** The split on its own, because the field count is what decides which branch runs. */
    static void split(final String line) {
        final String[] tokens = line.split("\\t|( +)", -1);
        final StringJoiner joiner = new StringJoiner("|");
        for (final String token : tokens) {
            joiner.add(escape(token));
        }
        System.out.printf("split\t%s\t%d\t%s%n", escape(line), tokens.length, joiner);
    }

    static String show(final BEDFeature feature) {
        if (feature == null) {
            return "null";
        }
        final StringJoiner joiner = new StringJoiner("|");
        joiner.add(String.format("%s:%d-%d", feature.getContig(), feature.getStart(),
                feature.getEnd()));
        joiner.add(feature.getName() == null ? "null" : escape(feature.getName()));
        joiner.add(Float.toString(feature.getScore()));
        joiner.add(feature.getStrand() == null ? "null" : feature.getStrand().toString());
        joiner.add(feature.getColor() == null ? "null" : String.format("%d,%d,%d",
                feature.getColor().getRed(), feature.getColor().getGreen(),
                feature.getColor().getBlue()));

        final StringJoiner exons = new StringJoiner(";");
        if (feature instanceof FullBEDFeature) {
            final List<FullBEDFeature.Exon> list = ((FullBEDFeature) feature).getExons();
            if (list != null) {
                for (final FullBEDFeature.Exon exon : list) {
                    // Exon.start and Exon.end are package-private fields with no accessor, so
                    // reflection is the only way to observe what addExon recorded. Everything
                    // else about the exon has a getter.
                    exons.add(String.format("%d-%d#%d#%d-%d", intField(exon, "start"),
                            intField(exon, "end"), exon.getNumber(), exon.getCdStart(),
                            exon.getCdEnd()));
                }
            }
        }
        joiner.add(exons.toString());
        return joiner.toString();
    }

    static int intField(final Object target, final String name) {
        try {
            final java.lang.reflect.Field field = target.getClass().getDeclaredField(name);
            field.setAccessible(true);
            return field.getInt(target);
        } catch (final ReflectiveOperationException e) {
            throw new IllegalStateException("no field " + name + " on " + target.getClass(), e);
        }
    }

    /** Tabs and spaces are the subject here, so they travel escaped. */
    static String escape(final String text) {
        final StringBuilder out = new StringBuilder();
        for (final char c : text.toCharArray()) {
            if (c == '\t') {
                out.append("\\t");
            } else if (c == ' ') {
                out.append("\\s");
            } else if (c < 0x20 || c > 0x7e) {
                out.append(String.format("\\u%04x", (int) c));
            } else if (c == '\\') {
                out.append("\\\\");
            } else {
                out.append(c);
            }
        }
        return out.toString();
    }
}
