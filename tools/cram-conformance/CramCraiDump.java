/*
 * The CRAM index: six numbers a line, and the arithmetic that queries them.
 *
 * A `.crai` is not a structure, it is a sorted text file: one line per slice, six tab-separated
 * integers, gzipped. What is worth pinning is therefore not a layout but four decisions.
 *
 *   - UNMAPPED-UNPLACED SORTS LAST, whatever its alignment start says, and its alignment start is
 *     not consulted at all. Everything else sorts by reference, then start, then container offset,
 *     then slice offset;
 *   - AN UNMAPPED ENTRY NEVER INTERSECTS, not even with itself, which is stated in the code as a
 *     special case rather than falling out of the arithmetic;
 *   - THE OVERLAP TEST IS A MIDPOINT COMPARISON: |a0 + b0 - a1 - b1| < span0 + span1, which is not
 *     the same expression as a0 < b1 && a1 < b0 and does not agree with it on a zero span;
 *   - A QUERY WITH A START OR A SPAN BELOW ONE MATCHES THE WHOLE SEQUENCE, so 0 and -1 are not
 *     out-of-range values here, they are a wildcard.
 *
 * Output:
 *
 *     entry\t<sequence>\t<start>\t<span>\t<container>\t<slice offset>\t<slice size>\t<serialized>
 *     parse\t<line>\t<serialized again>
 *     sort\t<label>\t<entries in, semicolon separated>\t<entries out>
 *     intersect\t<first>\t<second>\t<result>
 *     find\t<label>\t<sequence>\t<start>\t<span>\t<matches>
 *     leftmost\t<label>\t<result>
 *     err\t<what>\t<detail>\t<class>\t<message>
 *
 * An entry is written `sequence:start:span:container:offset:size` where a list is needed. A
 * serialized line is shown with its tabs as spaces and its newlines as `|`, because a row here
 * is itself tab-separated.
 *
 * Usage: CramCraiDump
 */

import htsjdk.samtools.cram.CRAIEntry;
import htsjdk.samtools.cram.CRAIIndex;

import java.io.ByteArrayOutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.StringJoiner;

public class CramCraiDump {

    public static void main(final String[] args) {
        System.out.println("# CramCraiDump: six numbers a line, and the arithmetic that queries them");

        // The line one entry makes.
        entry(0, 100, 50, 1000L, 20, 300);
        entry(0, 1, 1, 0L, 0, 1);
        entry(3, 2147483647, 2147483647, 9223372036854775807L, 2147483647, 2147483647);
        entry(-1, 0, 0, 500L, 10, 20);
        entry(-1, -1, -1, 500L, 10, 20);

        // And the entry a line makes.
        parse("0\t100\t50\t1000\t20\t300");
        parse("-1\t0\t0\t500\t10\t20");
        parse("2\t-5\t-5\t-1\t-2\t-3");

        // Sorting, which is what writing an index does to it.
        sort("unmapped-last", entries("1:100:10:0:0:1", "-1:1:1:0:0:1", "0:100:10:0:0:1"));
        sort("by-start", entries("0:300:10:0:0:1", "0:100:10:0:0:1", "0:200:10:0:0:1"));
        sort("by-container-then-slice", entries("0:100:10:20:5:1", "0:100:10:10:9:1",
                "0:100:10:10:2:1"));
        sort("unmapped-ignores-its-start", entries("-1:900:10:20:0:1", "-1:100:10:10:0:1"));
        sort("across-references", entries("2:1:1:0:0:1", "-1:1:1:0:0:1", "0:1:1:0:0:1",
                "1:1:1:0:0:1"));

        // The overlap test, at every boundary it has.
        intersect("0:100:10:0:0:1", "0:105:10:0:0:1");
        intersect("0:100:10:0:0:1", "0:110:10:0:0:1");
        intersect("0:100:10:0:0:1", "0:109:1:0:0:1");
        intersect("0:100:10:0:0:1", "0:90:10:0:0:1");
        intersect("0:100:10:0:0:1", "0:100:10:0:0:1");
        intersect("0:100:0:0:0:1", "0:100:0:0:0:1");
        intersect("0:100:10:0:0:1", "0:100:0:0:0:1");
        intersect("0:100:10:0:0:1", "1:100:10:0:0:1");
        intersect("-1:100:10:0:0:1", "-1:100:10:0:0:1");

        // Querying a list, including the two wildcards.
        final List<CRAIEntry> list = entries("0:100:10:0:0:1", "0:200:10:0:1:1", "0:300:10:0:2:1",
                "1:100:10:0:3:1", "-1:1:1:0:4:1");
        find("in-range", list, 0, 195, 10);
        find("out-of-range", list, 0, 400, 10);
        find("whole-sequence-by-start", list, 0, 0, 10);
        find("whole-sequence-by-span", list, 0, 100, 0);
        find("another-reference", list, 1, 0, 0);
        find("unmapped", list, -1, 0, 0);
        leftmost("the-list", list);
        leftmost("empty", new ArrayList<>());

        // What it refuses.
        errMultiRef();
        errLine("0\t100\t50", "too few columns");
        errLine("0\t100\t50\t1000\t20\t300\t400", "too many columns");
        errLine("x\t100\t50\t1000\t20\t300", "not a number");
        errLine("", "an empty line");
    }

    static void entry(final int sequenceId, final int start, final int span, final long container,
            final int sliceOffset, final int sliceSize) {
        final CRAIEntry entry = new CRAIEntry(sequenceId, start, span, container, sliceOffset,
                sliceSize);
        // The serialized form is itself tab-separated, so it is shown with spaces: a row of this
        // dump is split on tabs and could not carry it otherwise.
        System.out.printf("entry\t%d\t%d\t%d\t%d\t%d\t%d\t%s%n", sequenceId, start, span, container,
                sliceOffset, sliceSize, entry.toString().replace('\t', ' '));
    }

    static void parse(final String line) {
        final CRAIEntry entry = new CRAIEntry(line);
        System.out.printf("parse\t%s\t%s%n", line.replace('\t', ' '),
                entry.toString().replace('\t', ' '));
    }

    /** `sequence:start:span:container:offset:size`, which is how a list is written on one row. */
    static List<CRAIEntry> entries(final String... specifications) {
        final List<CRAIEntry> entries = new ArrayList<>();
        for (final String specification : specifications) {
            final String[] parts = specification.split(":");
            entries.add(new CRAIEntry(Integer.parseInt(parts[0]), Integer.parseInt(parts[1]),
                    Integer.parseInt(parts[2]), Long.parseLong(parts[3]),
                    Integer.parseInt(parts[4]), Integer.parseInt(parts[5])));
        }
        return entries;
    }

    static void sort(final String label, final List<CRAIEntry> entries) {
        final String before = show(entries);
        final List<CRAIEntry> sorted = new ArrayList<>(entries);
        Collections.sort(sorted);
        System.out.printf("sort\t%s\t%s\t%s%n", label, before, show(sorted));

        // And the bytes an index of them makes, which is the sort plus a newline each.
        final CRAIIndex index = new CRAIIndex();
        index.addEntries(entries);
        final ByteArrayOutputStream out = new ByteArrayOutputStream();
        index.writeIndex(out);
        // Lines separated by `|` and fields by spaces, for the same reason.
        System.out.printf("index\t%s\t%s%n", label,
                new String(out.toByteArray()).replace("\n", "|").replace('\t', ' ').trim());
    }

    static void intersect(final String first, final String second) {
        final List<CRAIEntry> pair = entries(first, second);
        System.out.printf("intersect\t%s\t%s\t%b%n", first, second,
                CRAIEntry.intersect(pair.get(0), pair.get(1)));
    }

    static void find(final String label, final List<CRAIEntry> list, final int sequenceId,
            final int start, final int span) {
        final List<CRAIEntry> found = CRAIIndex.find(list, sequenceId, start, span);
        System.out.printf("find\t%s\t%d\t%d\t%d\t%s%n", label, sequenceId, start, span,
                show(found));
    }

    static void leftmost(final String label, final List<CRAIEntry> list) {
        final CRAIEntry entry = CRAIIndex.getLeftmost(list);
        System.out.printf("leftmost\t%s\t%s%n", label, entry == null ? "-" : show(entry));
    }

    static void errMultiRef() {
        try {
            new CRAIEntry(-2, 1, 1, 0L, 0, 1);
            System.out.printf("err\tmulti-ref\tsequence=-2\t-\t-%n");
        } catch (final Throwable t) {
            System.out.printf("err\tmulti-ref\tsequence=-2\t%s\t%s%n", t.getClass().getSimpleName(),
                    String.valueOf(t.getMessage()));
        }
    }

    static void errLine(final String line, final String what) {
        try {
            new CRAIEntry(line);
            System.out.printf("err\t%s\t%s\t-\t-%n", what, line.replace('\t', ' '));
        } catch (final Throwable t) {
            System.out.printf("err\t%s\t%s\t%s\t%s%n", what, line.replace('\t', ' '),
                    t.getClass().getSimpleName(), String.valueOf(t.getMessage()));
        }
    }

    static String show(final List<CRAIEntry> entries) {
        if (entries.isEmpty()) {
            return "-";
        }
        final StringJoiner joiner = new StringJoiner(";");
        for (final CRAIEntry entry : entries) {
            joiner.add(show(entry));
        }
        return joiner.toString();
    }

    static String show(final CRAIEntry entry) {
        return String.format("%d:%d:%d:%d:%d:%d", entry.getSequenceId(), entry.getAlignmentStart(),
                entry.getAlignmentSpan(), entry.getContainerStartByteOffset(),
                entry.getSliceByteOffsetFromCompressionHeaderStart(), entry.getSliceByteSize());
    }
}
