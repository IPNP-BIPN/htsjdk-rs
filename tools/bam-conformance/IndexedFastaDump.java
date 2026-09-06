/*
 * What htsjdk's indexed FASTA reader answers, so the port's arithmetic is measured rather than
 * argued.
 *
 * `IndexedFastaSequenceFile.getSubsequenceAt` reads a query by seeking: the first base's byte
 * offset comes from the `.fai`'s bases-per-line and bytes-per-line, and every line boundary the
 * query crosses is a jump over the terminator. Nothing scans for a newline, so a file whose `.fai`
 * disagrees with its own line lengths is read wrongly and silently -- which is why the port has to
 * do the same arithmetic rather than a reasonable one.
 *
 * The cases are the boundaries: inside one line, exactly one line, across a boundary, the whole
 * contig, the last base, a second contig (whose offset the index carries), a two-byte terminator,
 * and the two refusals.
 *
 * Output:
 *     index\t<file>\t<name>\t<size>\t<location>\t<basesPerLine>\t<bytesPerLine>
 *     query\t<file>\t<contig>:<start>-<stop>\t<bases>
 *     error\t<file>\t<contig>:<start>-<stop>\t<exception class>\t<message>
 */

import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.FastaSequenceIndexCreator;
import htsjdk.samtools.reference.FastaSequenceIndexEntry;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

public class IndexedFastaDump {

    static String escape(final String s) {
        return s.replace("\\", "\\\\").replace("\t", "\\t").replace("\r", "\\r").replace("\n", "\\n");
    }

    public static void main(final String[] args) throws Exception {
        final Path dir = Files.createTempDirectory("indexed-fasta");

        // The `.fai` is written by htsjdk's own creator rather than by hand. Writing one by hand
        // is how the first version of this dump measured a shifted file: an offset that is one
        // byte wrong makes the reader return a TERMINATOR byte as a base, and htsjdk reports it as
        // an answer rather than as an error. The index is part of the golden below, so the
        // creator's arithmetic is pinned too.
        //
        // LF terminators, six bases per line, two contigs of different shapes.
        write(dir, "lf.fasta", ">chr1\nACGTAC\nGTACGT\nAC\n>chr2\nTTTTTT\nGG\n");
        // CRLF: the terminator is two bytes, and only the index says so.
        write(dir, "crlf.fasta", ">chr1\r\nACGTAC\r\nGTACGT\r\n");
        // One line that is the whole contig, which is what a small reference looks like.
        write(dir, "single.fasta", ">chrS\nACGTACGTAC\n");

        for (final String name : new String[] {"lf.fasta", "crlf.fasta", "single.fasta"}) {
            final Path fasta = dir.resolve(name);
            final FastaSequenceIndex index = new FastaSequenceIndex(dir.resolve(name + ".fai"));
            for (final FastaSequenceIndexEntry entry : index) {
                System.out.printf("index\t%s\t%s\t%d\t%d\t%d\t%d%n", name, entry.getContig(),
                        entry.getSize(), entry.getLocation(), entry.getBasesPerLine(),
                        entry.getBytesPerLine());
            }
            try (final IndexedFastaSequenceFile reader =
                         new IndexedFastaSequenceFile(fasta, index)) {
                for (final long[] q : queries(name)) {
                    query(reader, name, contig(name, q[2]), q[0], q[1]);
                }
            }
        }
    }

    static String contig(final String file, final long which) {
        if (file.startsWith("single")) {
            return "chrS";
        }
        return which == 0 ? "chr1" : "chr2";
    }

    /** start, stop, which contig. */
    static long[][] queries(final String file) {
        if (file.startsWith("single")) {
            return new long[][] {{1, 10, 0}, {3, 7, 0}, {10, 10, 0}, {1, 11, 0}};
        }
        if (file.startsWith("crlf")) {
            return new long[][] {{1, 12, 0}, {6, 7, 0}, {7, 12, 0}, {2, 2, 0}};
        }
        return new long[][] {
            {1, 6, 0}, {2, 4, 0}, {1, 14, 0}, {5, 9, 0}, {7, 12, 0}, {14, 14, 0},
            {5, 4, 0}, {6, 4, 0}, {1, 15, 0}, {1, 8, 1}, {7, 8, 1},
        };
    }

    static void query(final IndexedFastaSequenceFile reader, final String file,
                      final String contig, final long start, final long stop) {
        final String label = contig + ":" + start + "-" + stop;
        try {
            final ReferenceSequence sequence = reader.getSubsequenceAt(contig, start, stop);
            // Escaped: a wrong index makes the reader answer with a terminator byte, and a raw
            // newline in the middle of a value would break the row rather than reporting it.
            System.out.printf("query\t%s\t%s\t%s%n", file, label,
                    escape(new String(sequence.getBases(), StandardCharsets.UTF_8)));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s\t%s\t%s%n", file, label,
                    e.getClass().getName(), e.getMessage());
        }
    }

    static void write(final Path dir, final String name, final String fasta) throws Exception {
        final Path path = dir.resolve(name);
        Files.writeString(path, fasta, StandardCharsets.UTF_8);
        FastaSequenceIndexCreator.create(path, true);
    }
}
