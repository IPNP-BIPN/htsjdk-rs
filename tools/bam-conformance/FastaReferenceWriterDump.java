/*
 * FastaReferenceWriter: the FASTA, the .fai and the .dict it writes, taken from the reference.
 *
 * Three GATK tools write a reference through this class -- FastaReferenceMaker,
 * FastaAlternateReferenceMaker and ShiftFasta -- and none of them can be ported until it is. What
 * they need is not "a FASTA writer" but this one: the index it writes beside the sequence has to
 * agree with the reader that will open it, and the offsets in it are byte counts of the file this
 * writer produced.
 *
 * Seven behaviours this is built to catch.
 *
 *   - THE INDEX OFFSET IS TAKEN AFTER THE HEADER LINE, so it counts the `>name description\n` bytes
 *     of THIS sequence and every byte of every sequence before it, newlines included;
 *   - BYTES-PER-LINE IS BASES-PER-LINE PLUS ONE, always, even for a sequence shorter than one line
 *     and even for the last line of a sequence that does not fill it;
 *   - A SEQUENCE IS TERMINATED BY EXACTLY ONE NEWLINE, written when the sequence is closed rather
 *     than after each line, so a sequence whose length is a multiple of the line width does not get
 *     a blank line and one that is not does not get a missing one;
 *   - THE LINE WIDTH IS PER SEQUENCE, not per file: startSequence takes its own, and the index
 *     records the one that sequence used;
 *   - APPENDING IN CHUNKS IS THE SAME FILE as appending at once, because the line breaks are
 *     decided by a running count rather than by the chunk boundaries;
 *   - THE MD5 IS OF THE UPPER-CASED BASES, computed chunk by chunk as they are written, so a
 *     lower-case sequence keeps its case in the FASTA and hashes as though it did not;
 *   - AND THE REFUSALS ARE PART OF THE CONTRACT: a blank in a name, a control character, an empty
 *     name, a sequence with no bases, a repeated name, a base that is not IUPAC, and a line width
 *     of zero. Each is a different exception class or message and a port that refused them all the
 *     same way would be wrong in a way no valid input reveals.
 *
 * Output:
 *
 *     fasta\t<label>\t<the FASTA text, escaped>
 *     fai\t<label>\t<the .fai text, escaped>
 *     dict\t<label>\t<the .dict text, escaped>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: FastaReferenceWriterDump
 */

import htsjdk.samtools.reference.FastaReferenceWriter;
import htsjdk.samtools.reference.FastaReferenceWriterBuilder;

import java.io.ByteArrayOutputStream;
import java.nio.charset.StandardCharsets;

public class FastaReferenceWriterDump {

    /** One writer's three outputs, kept as byte arrays so nothing touches a filesystem. */
    static class Outputs {
        final ByteArrayOutputStream fasta = new ByteArrayOutputStream();
        final ByteArrayOutputStream fai = new ByteArrayOutputStream();
        final ByteArrayOutputStream dict = new ByteArrayOutputStream();

        FastaReferenceWriter writer(final int basesPerLine, final boolean md5) throws Exception {
            return new FastaReferenceWriterBuilder()
                    .setFastaOutput(fasta)
                    .setIndexOutput(fai)
                    .setDictOutput(dict)
                    .setBasesPerLine(basesPerLine)
                    .setEmitMd5(md5)
                    .build();
        }
    }

    public static void main(final String[] args) throws Exception {
        System.out.println("# FastaReferenceWriterDump: the FASTA, the .fai and the .dict");

        // One sequence at the default width, whose length is not a multiple of it.
        emit("default-width", 60, false, writer -> {
            writer.startSequence("chr1");
            writer.appendBases(repeat("ACGT", 25).getBytes());   // 100 bases
        });

        // A width that divides the length exactly, so the last line is full.
        emit("exact-multiple", 10, false, writer -> {
            writer.startSequence("chr1");
            writer.appendBases(repeat("ACGT", 10).getBytes());   // 40 bases, four lines
        });

        // Shorter than one line.
        emit("short", 60, false, writer -> {
            writer.startSequence("chr1");
            writer.appendBases("ACGT".getBytes());
        });

        // Two sequences, the second with a description and its own width.
        emit("two-sequences", 12, false, writer -> {
            writer.startSequence("chr1");
            writer.appendBases(repeat("ACGT", 7).getBytes());    // 28 bases at 12
            writer.startSequence("chr2", "the second one", 5);
            writer.appendBases(repeat("TTGC", 3).getBytes());    // 12 bases at 5
        });

        // The same file, appended in chunks that do not line up with the line width.
        emit("chunked", 12, false, writer -> {
            writer.startSequence("chr1");
            writer.appendBases("ACG".getBytes());
            writer.appendBases("TACGTACG".getBytes());
            writer.appendBases("TACGTACGTACGTACGT".getBytes());  // 28 bases in three pieces
        });

        // Lower case and IUPAC codes, which the writer keeps as they are.
        emit("mixed-case", 10, false, writer -> {
            writer.startSequence("chr1");
            writer.appendBases("acgtRYKMSWacgtNNNN".getBytes());
        });

        // The md5, which is of the upper-cased bases.
        emit("md5", 10, true, writer -> {
            writer.startSequence("chr1");
            writer.appendBases("acgtacgtac".getBytes());
        });
        emit("md5-uppercase", 10, true, writer -> {
            writer.startSequence("chr1");
            writer.appendBases("ACGTACGTAC".getBytes());
        });

        // And the refusals.
        error("empty-name", 60, writer -> writer.startSequence(""));
        error("blank-in-name", 60, writer -> writer.startSequence("chr 1"));
        error("control-in-name", 60, writer -> writer.startSequence("chr\u0001"));
        error("control-in-description", 60, writer -> writer.startSequence("chr1", "a\u0001b", 60));
        error("tab-in-description", 60, writer -> {
            // A tab IS allowed in a description, so this one is not an error and the row records
            // the header line it produces.
            writer.startSequence("chr1", "a\tb", 60);
            writer.appendBases("ACGT".getBytes());
        });
        error("no-bases", 60, writer -> {
            writer.startSequence("chr1");
            writer.startSequence("chr2");
        });
        error("duplicate-name", 60, writer -> {
            writer.startSequence("chr1");
            writer.appendBases("ACGT".getBytes());
            writer.startSequence("chr1");
        });
        error("bad-base", 60, writer -> {
            writer.startSequence("chr1");
            writer.appendBases("ACGZ".getBytes());
        });
        error("bases-before-sequence", 60, writer -> writer.appendBases("ACGT".getBytes()));
        error("zero-width", 60, writer -> writer.startSequence("chr1", "", 0));
        error("close-with-no-bases", 60, writer -> {
            writer.startSequence("chr1");
            writer.close();
        });
    }

    interface Body {
        void run(FastaReferenceWriter writer) throws Exception;
    }

    /** Runs a body to completion and prints the three outputs. */
    static void emit(final String label, final int basesPerLine, final boolean md5, final Body body)
            throws Exception {
        final Outputs outputs = new Outputs();
        try (final FastaReferenceWriter writer = outputs.writer(basesPerLine, md5)) {
            body.run(writer);
        }
        System.out.printf("fasta\t%s\t%s%n", label, escape(outputs.fasta.toString("UTF-8")));
        System.out.printf("fai\t%s\t%s%n", label, escape(outputs.fai.toString("UTF-8")));
        System.out.printf("dict\t%s\t%s%n", label, escape(outputs.dict.toString("UTF-8")));
    }

    /**
     * Runs a body that is expected to fail, and prints whatever came out.
     *
     * A body that does not fail prints its outputs instead, which is how the tab-in-description
     * case records that a tab is allowed where other control characters are not.
     */
    static void error(final String label, final int basesPerLine, final Body body) {
        final Outputs outputs = new Outputs();
        try {
            final FastaReferenceWriter writer = outputs.writer(basesPerLine, false);
            body.run(writer);
            writer.close();
            System.out.printf("fasta\t%s\t%s%n", label,
                    escape(new String(outputs.fasta.toByteArray(), StandardCharsets.UTF_8)));
            System.out.printf("fai\t%s\t%s%n", label,
                    escape(new String(outputs.fai.toByteArray(), StandardCharsets.UTF_8)));
            System.out.printf("dict\t%s\t%s%n", label,
                    escape(new String(outputs.dict.toByteArray(), StandardCharsets.UTF_8)));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(),
                    escape(String.valueOf(e.getMessage())));
        }
    }

    static String repeat(final String unit, final int times) {
        final StringBuilder builder = new StringBuilder(unit.length() * times);
        for (int i = 0; i < times; i++) {
            builder.append(unit);
        }
        return builder.toString();
    }

    static String escape(final String text) {
        return text.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n")
                .replace("\u0001", "\\u0001");
    }
}
