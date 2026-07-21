/*
 * Writes .interval_list files with htsjdk's IntervalList and emits them.
 *
 * Output: <case name> <TAB> <the interval lines only, escaped>
 *
 * Only the body is emitted. The SAM header that precedes it is written by SAMTextHeaderCodec,
 * which has its own conformance suite already; what is under test here is the interval lines and
 * above all their *order*, which comes from IntervalCoordinateComparator and therefore from the
 * sequence dictionary rather than from the alphabet.
 */

import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import java.io.File;
import java.nio.file.Files;
import java.util.*;

public class IntervalListDump {

    static String esc(String s) {
        return s.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
    }

    /** A dictionary whose index order is not its lexicographic order. */
    static SAMFileHeader header(String... contigs) {
        SAMFileHeader h = new SAMFileHeader();
        SAMSequenceDictionary d = new SAMSequenceDictionary();
        for (String c : contigs) d.addSequence(new SAMSequenceRecord(c, 100000));
        h.setSequenceDictionary(d);
        return h;
    }

    static void emit(String name, IntervalList list) throws Exception {
        File f = File.createTempFile("ilist", ".interval_list");
        list.write(f);
        StringBuilder body = new StringBuilder();
        for (String line : new String(Files.readAllBytes(f.toPath())).split("\n")) {
            if (!line.startsWith("@")) body.append(line).append("\n");
        }
        System.out.println(name + "\t" + esc(body.toString()));
        f.delete();
    }

    static Interval iv(String c, int s, int e) { return new Interval(c, s, e); }
    static Interval iv(String c, int s, int e, boolean neg, String n) {
        return new Interval(c, s, e, neg, n);
    }

    public static void main(String[] args) throws Exception {
        SAMFileHeader h = header("chr1", "chr2", "chr10", "chrX");

        IntervalList plain = new IntervalList(h.clone());
        plain.add(iv("chr1", 100, 200));
        emit("one_interval", plain);

        // The order test. Inserted so that neither the input order nor the alphabet is the
        // answer: only the dictionary index is.
        IntervalList unsorted = new IntervalList(h.clone());
        unsorted.add(iv("chr10", 1, 10));
        unsorted.add(iv("chrX", 1, 10));
        unsorted.add(iv("chr2", 1, 10));
        unsorted.add(iv("chr1", 1, 10));
        emit("unsorted_as_given", unsorted);
        emit("unsorted_sorted", unsorted.sorted());

        // Strand and name rendering, including the null name.
        IntervalList named = new IntervalList(h.clone());
        named.add(iv("chr1", 1, 10, false, "plus_named"));
        named.add(iv("chr1", 20, 30, true, "minus_named"));
        named.add(iv("chr1", 40, 50, false, null));
        named.add(iv("chr1", 60, 70, true, null));
        emit("strands_and_names", named.sorted());

        // The comparator's tail: same coordinates, differing only by strand then by name.
        IntervalList tail = new IntervalList(h.clone());
        tail.add(iv("chr1", 1, 10, true, "zzz"));
        tail.add(iv("chr1", 1, 10, false, "zzz"));
        tail.add(iv("chr1", 1, 10, false, "aaa"));
        tail.add(iv("chr1", 1, 10, false, null));
        emit("comparator_tail", tail.sorted());

        // uniqued: overlapping, abutting, and separated by one base.
        IntervalList overlapping = new IntervalList(h.clone());
        overlapping.add(iv("chr1", 1, 10, false, "first"));
        overlapping.add(iv("chr1", 5, 20, false, "second"));
        overlapping.add(iv("chr1", 21, 30, false, "abutting"));
        overlapping.add(iv("chr1", 32, 40, false, "separated"));
        emit("uniqued", overlapping.uniqued());
        emit("uniqued_concatenated", overlapping.uniqued(true));
        emit("sorted_not_uniqued", overlapping.sorted());

        // padded, including the clamp at 1.
        IntervalList toPad = new IntervalList(h.clone());
        toPad.add(iv("chr1", 5, 10, false, "near_start"));
        toPad.add(iv("chr2", 500, 600, false, "middle"));
        emit("padded", toPad.padded(100, 50));

        // A single-base interval, where start == end.
        IntervalList single = new IntervalList(h.clone());
        single.add(iv("chr1", 42, 42, false, "point"));
        emit("single_base", single);

        // An empty list writes a header and no interval lines at all.
        emit("empty", new IntervalList(h.clone()));
    }
}
