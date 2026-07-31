/*
 * IntervalListCodec, taken from the reference.
 *
 * This is the codec `-L regions.interval_list` runs through in GATK, and it is a stricter parser
 * than the reader `IntervalList.fromReader` uses on the same file. Five of its decisions are not
 * what "the interval list parser" suggests:
 *
 *   - the field count is checked against 5 EXACTLY, and the count comes from `split("\t")`, which
 *     drops trailing empty fields, so a line ending in a tab has fewer fields than separators;
 *   - the strand column is decoded by `Strand.decode`, which answers NONE for anything that is not
 *     a single "+" or "-", and NONE is then REJECTED, so "." is an error here while the writer
 *     never emits one;
 *   - an interval on a contig the dictionary does not hold is NULL, not an error: the line is
 *     dropped with a warning and the file still loads;
 *   - an interval past the end of a contig IS an error, unless the contig's length is 0;
 *   - `start == end + 1` is legal and describes an empty interval, while `start == end + 2` is
 *     refused with a message that quotes WarGames.
 *
 * Output:
 *
 *     interval\t<escaped line>\t<dict>\t<null|contig:start-end|strand|name>
 *     interval\t<escaped line>\t<dict>\tE:<class>:<message>
 *     candecode\t<path>\t<true|false>
 *
 * Usage: IntervalListCodecDump
 */

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Interval;
import htsjdk.tribble.IntervalList.IntervalListCodec;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class IntervalListCodecDump {

    /** Every line the codec is asked to decode, against every dictionary. */
    static final String[] LINES = {
        // The ordinary shapes.
        "chr1\t1\t10\t+\tname",
        "chr1\t1\t10\t-\tname",
        "chr1\t1\t10\t+\t.",
        "chr1\t1\t10\t+\tmy name",
        "chr1\t200\t200\t+\tname",
        // The empty interval, and the first illegal one past it.
        "chr1\t11\t10\t+\tname",
        "chr1\t12\t10\t+\tname",
        "chr1\t1\t0\t+\tname",
        // Coordinates below 1.
        "chr1\t0\t10\t+\tname",
        "chr1\t-5\t10\t+\tname",
        // Past the end of the contig, and the same on a contig of declared length 0.
        "chr1\t1\t201\t+\tname",
        "chr3\t1\t999999\t+\tname",
        // A contig the dictionary does not hold, which is dropped rather than refused.
        "chrX\t1\t10\t+\tname",
        "CHR1\t1\t10\t+\tname",
        "\t1\t10\t+\tname",
        // Field counts. split("\t") drops trailing empties, so the last two are short.
        "chr1\t1\t10\t+\tname\textra",
        "chr1\t1\t10\t+",
        "chr1\t1\t10\t+\t",
        "chr1\t1\t10\t+\t\t",
        "chr1\t1\t10",
        "chr1",
        // An empty strand field is present, so the count is 5 and the strand is what fails.
        "chr1\t1\t10\t\tname",
        "chr1\t1\t10\t.\tname",
        "chr1\t1\t10\t++\tname",
        "chr1\t1\t10\tplus\tname",
        "chr1\t1\t10\tF\tname",
        // Numbers Java takes and numbers it refuses.
        "chr1\t+1\t10\t+\tname",
        "chr1\t01\t10\t+\tname",
        "chr1\t 1\t10\t+\tname",
        "chr1\t1.0\t10\t+\tname",
        "chr1\t1\t1e1\t+\tname",
        "chr1\t2147483648\t10\t+\tname",
        "chr1\t\t10\t+\tname",
        // Header and blank lines, which are null rather than errors.
        "@HD\tVN:1.6\tSO:coordinate",
        "@SQ\tSN:chr1\tLN:200",
        "",
        "   ",
        "\t",
        "\t\t",
        // A leading space is not a header prefix here, because there is no prefix: the line is
        // trimmed only for the emptiness test, and split sees the untrimmed line.
        " chr1\t1\t10\t+\tname",
    };

    static final String[] PATHS = {
        "a.interval_list", "a.INTERVAL_LIST", "a.interval_list.gz", "a.interval_list.GZ",
        "a.interval_list.bgz", "a.interval_list.gz.gz", "a.intervals", "a.list", "a.bed",
        "interval_list", ".interval_list", "a.interval_list.tbi",
    };

    public static void main(final String[] args) {
        System.out.println("# IntervalListCodecDump: IntervalListCodec and canDecode");

        for (final String line : LINES) {
            for (final Dict dict : dictionaries()) {
                emit(line, dict);
            }
        }

        final IntervalListCodec codec = new IntervalListCodec();
        for (final String path : PATHS) {
            System.out.printf("candecode\t%s\t%b%n", path, codec.canDecode(path));
        }
    }

    /** A named dictionary, so the golden says which one a row was decoded against. */
    static final class Dict {
        final String label;
        final SAMSequenceDictionary dictionary;

        Dict(final String label, final SAMSequenceDictionary dictionary) {
            this.label = label;
            this.dictionary = dictionary;
        }
    }

    static List<Dict> dictionaries() {
        final List<Dict> dicts = new ArrayList<>();
        dicts.add(new Dict("two-contigs", new SAMSequenceDictionary(Arrays.asList(
                new SAMSequenceRecord("chr1", 200),
                new SAMSequenceRecord("chr2", 200),
                // A length of 0 turns off the end-of-contig check, so this contig accepts
                // anything.
                new SAMSequenceRecord("chr3", 0)))));
        dicts.add(new Dict("empty", new SAMSequenceDictionary()));
        // No dictionary at all: decode refuses before it parses anything, which is why a header
        // is mandatory for this codec and optional for the BED one.
        dicts.add(new Dict("null", null));
        return dicts;
    }

    static void emit(final String line, final Dict dict) {
        try {
            final Interval interval = new IntervalListCodec(dict.dictionary).decode(line);
            System.out.printf("interval\t%s\t%s\t%s%n", escape(line), dict.label, show(interval));
        } catch (final Exception | AssertionError e) {
            System.out.printf("interval\t%s\t%s\tE:%s:%s%n", escape(line), dict.label,
                    e.getClass().getName(),
                    e.getMessage() == null ? "" : e.getMessage().replace('\n', ' '));
        }
    }

    static String show(final Interval interval) {
        if (interval == null) {
            return "null";
        }
        return String.format("%s:%d-%d|%s|%s", escape(interval.getContig()), interval.getStart(),
                interval.getEnd(), interval.isNegativeStrand() ? "-" : "+",
                interval.getName() == null ? "null" : escape(interval.getName()));
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
