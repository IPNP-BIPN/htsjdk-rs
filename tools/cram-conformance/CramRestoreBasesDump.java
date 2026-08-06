/*
 * The bases restored: the last piece of the record model's reverse direction.
 *
 * The cigar comes back from the feature positions. The bases come back from three sources at once:
 * the features say what differs, the REFERENCE supplies everything else, and the SUBSTITUTION
 * MATRIX turns a substitution's code back into a base. Nothing in the record carries the matching
 * bases themselves.
 *
 * Seven things here are decisions rather than layout.
 *
 *   - IT IS TWO PASSES, AND THE SECOND ONE OVERWRITES. ReadBase and Bases are skipped in the main
 *     loop with a comment saying "defer until after the reference bases are retrieved", then
 *     applied at the end straight into the array. So they win over whatever the reference fill put
 *     there, and a Bases feature can overwrite bases the features before it produced;
 *   - THE COMMENT ABOVE THE LOOP IS WRONG. It says "ReadFeatures use a 0-based feature position",
 *     and every position in it is used one-based. The forward direction measured the same thing
 *     from the other side;
 *   - THE TRAILING FILL STOPS AT THE END OF THE REFERENCE, and the array it leaves behind is a
 *     fresh byte[], so the bases past that point are NUL rather than N. What they become is decided
 *     later, by the lookup table;
 *   - toBamReadBasesInPlace INDEXES A 127-BYTE TABLE WITH A SIGNED BYTE. Every base goes through
 *     it: NUL becomes N, a lower-case base is folded up, and a byte of 127 or above is an index
 *     out of that table;
 *   - getByteOrDefault GUARDS THE UPPER BOUND AND NOT THE LOWER, the same asymmetry the forward
 *     direction has;
 *   - A SUBSTITUTION IS RESOLVED THROUGH THE MATRIX BY CODE, against a NORMALIZED reference base,
 *     so the read base restored depends on a table stored once per container;
 *   - UNKNOWN BASES OR A READ LENGTH OF ZERO RETURN THE NULL SEQUENCE, which is a shared empty
 *     array rather than a run of Ns.
 *
 * Output:
 *
 *     reference\t<the reference bases>\t<the region's zero-based start>
 *     matrix\t<the 20 encoded matrix bytes, hex>
 *     case\t<label>\t<alignment start>\t<read length>\t<unknown bases>\t<features>\t<restored bases>
 *     lookup\t<input byte>\t<what toBamReadBasesInPlace makes of it>
 *     err\t<label>\t<alignment start>\t<read length>\t<unknown bases>\t<features>\t<class>\t<message>
 *
 * Usage: CramRestoreBasesDump
 */

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.cram.build.CRAMReferenceRegion;
import htsjdk.samtools.cram.encoding.readfeatures.BaseQualityScore;
import htsjdk.samtools.cram.encoding.readfeatures.Bases;
import htsjdk.samtools.cram.encoding.readfeatures.Deletion;
import htsjdk.samtools.cram.encoding.readfeatures.HardClip;
import htsjdk.samtools.cram.encoding.readfeatures.InsertBase;
import htsjdk.samtools.cram.encoding.readfeatures.Insertion;
import htsjdk.samtools.cram.encoding.readfeatures.Padding;
import htsjdk.samtools.cram.encoding.readfeatures.ReadBase;
import htsjdk.samtools.cram.encoding.readfeatures.ReadFeature;
import htsjdk.samtools.cram.encoding.readfeatures.RefSkip;
import htsjdk.samtools.cram.encoding.readfeatures.Scores;
import htsjdk.samtools.cram.encoding.readfeatures.SoftClip;
import htsjdk.samtools.cram.encoding.readfeatures.Substitution;
import htsjdk.samtools.cram.ref.CRAMReferenceSource;
import htsjdk.samtools.cram.structure.CRAMRecordReadFeatures;
import htsjdk.samtools.cram.structure.SubstitutionMatrix;
import htsjdk.samtools.util.SequenceUtil;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.StringJoiner;

public class CramRestoreBasesDump {

    // Aperiodic on purpose: against ACGTACGT... a reference cursor off by four reads the same
    // bases as a correct one, and every case below would pass anyway.
    private static final String REFERENCE = "ACGTTGCAAGCTMRWSTTACGGCA";

    /**
     * The matrix as the reference encodes it: five bytes, one per reference base, each packing the
     * four substitutions for that base two bits at a time. This one is the identity ordering, so
     * code 0 is the first substitute of ACGTN order for that reference base.
     */
    private static final byte[] MATRIX = new byte[] {0x1B, 0x1B, 0x1B, 0x1B, 0x1B};

    public static void main(final String[] args) {
        System.out.println("# CramRestoreBasesDump: the bases restored from the features and the reference");

        final SubstitutionMatrix matrix = new SubstitutionMatrix(MATRIX);
        System.out.printf("reference\t%s\t%d%n", REFERENCE, region().getRegionStart());
        System.out.printf("matrix\t%s%n", hex(matrix.getEncodedMatrix()));

        // What the matrix answers, so the substitution cases below are readable.
        for (final byte refBase : "ACGTN".getBytes()) {
            final StringJoiner joiner = new StringJoiner(",");
            for (byte code = 0; code < 4; code++) {
                joiner.add(String.format("%d=%c", code, matrix.base(refBase, code)));
            }
            System.out.printf("matrixrow\t%c\t%s%n", refBase, joiner.toString());
        }

        // Nothing but the reference.
        emit("no-features", 1, 8, false, Collections.emptyList(), matrix);
        emit("no-features-mid-reference", 5, 8, false, Collections.emptyList(), matrix);

        // One substitution, resolved through the matrix.
        emit("substitution-code-0", 1, 8, false, Arrays.asList(new Substitution(4, (byte) 0)),
                matrix);
        emit("substitution-code-2", 1, 8, false, Arrays.asList(new Substitution(4, (byte) 2)),
                matrix);
        emit("substitution-no-code", 1, 8, false,
                Arrays.asList(new Substitution(4, (byte) -1)), matrix);
        emit("substitution-code-out-of-range", 1, 8, false,
                Arrays.asList(new Substitution(4, (byte) 9)), matrix);
        emit("substitution-against-iupac", 13, 4, false,
                Arrays.asList(new Substitution(1, (byte) 0)), matrix);

        // The features that put their own bases in.
        emit("insert-base", 1, 8, false, Arrays.asList(new InsertBase(3, (byte) 'T')), matrix);
        emit("insertion", 1, 8, false, Arrays.asList(new Insertion(3, "TTT".getBytes())), matrix);
        emit("soft-clip", 1, 8, false, Arrays.asList(new SoftClip(1, "TTT".getBytes())), matrix);

        // The features that move the reference cursor and put nothing in.
        emit("deletion", 1, 8, false, Arrays.asList(new Deletion(5, 2)), matrix);
        emit("ref-skip", 1, 8, false, Arrays.asList(new RefSkip(5, 2)), matrix);
        emit("padding", 1, 8, false, Arrays.asList(new Padding(5, 2)), matrix);
        emit("hard-clip", 1, 8, false, Arrays.asList(new HardClip(1, 2)), matrix);

        // The deferred pass, which overwrites what the reference fill wrote.
        emit("read-base", 1, 8, false, Arrays.asList(new ReadBase(4, (byte) 'M', (byte) 40)),
                matrix);
        emit("bases", 1, 8, false, Arrays.asList(new Bases(1, "TTTT".getBytes())), matrix);
        emit("bases-over-an-insertion", 1, 8, false,
                Arrays.asList(new Insertion(1, "GGG".getBytes()), new Bases(1, "TTT".getBytes())),
                matrix);
        emit("read-base-over-a-substitution", 1, 8, false,
                Arrays.asList(new Substitution(4, (byte) 0), new ReadBase(4, (byte) 'M',
                        (byte) 40)), matrix);

        // The features restoreReadBases does not act on at all.
        emit("scores", 1, 8, false, Arrays.asList(new Scores(1, "IIII".getBytes())), matrix);
        emit("base-quality-score", 1, 8, false,
                Arrays.asList(new BaseQualityScore(4, (byte) 40)), matrix);

        // Past the end of the reference, where the trailing fill stops.
        emit("past-the-end", 20, 8, false, Collections.emptyList(), matrix);
        emit("entirely-past-the-end", 30, 4, false, Collections.emptyList(), matrix);
        emit("past-the-end-with-a-feature", 20, 8, false,
                Arrays.asList(new InsertBase(1, (byte) 'T')), matrix);

        // The two shortcuts at the top.
        emit("unknown-bases", 1, 8, true, Collections.emptyList(), matrix);
        emit("zero-read-length", 1, 0, false, Collections.emptyList(), matrix);

        // A feature that writes past the end of the read.
        emit("insert-base-past-the-read", 1, 8, false,
                Arrays.asList(new InsertBase(9, (byte) 'T')), matrix);
        emit("insertion-past-the-read", 1, 8, false,
                Arrays.asList(new Insertion(7, "TTTT".getBytes())), matrix);

        // A base outside the lookup table, which only a feature can put there.
        emit("lower-case-inserted-base", 1, 8, false,
                Arrays.asList(new InsertBase(3, (byte) 'a')), matrix);
        emit("high-inserted-base", 1, 8, false,
                Arrays.asList(new InsertBase(3, (byte) 0xE9)), matrix);
        emit("del-inserted-base", 1, 8, false,
                Arrays.asList(new InsertBase(3, (byte) 127)), matrix);

        // The table itself, at every boundary that matters.
        for (final int value : new int[] {0, 'A', 'a', 'C', 'c', '=', ']', 'M', 'm', 'Y', 'y', '.',
                '*', 126, 127, 128, 233, 255}) {
            lookup(value);
        }
    }

    static void emit(final String label, final int alignmentStart, final int readLength,
            final boolean unknownBases, final List<ReadFeature> features,
            final SubstitutionMatrix matrix) {
        try {
            final CRAMReferenceRegion region = region();
            final byte[] bases = CRAMRecordReadFeatures.restoreReadBases(features, unknownBases,
                    alignmentStart, readLength, region, matrix);
            System.out.printf("case\t%s\t%d\t%d\t%b\t%s\t%s%n", label, alignmentStart, readLength,
                    unknownBases, describe(features), escape(new String(bases)));
        } catch (final Throwable t) {
            System.out.printf("err\t%s\t%d\t%d\t%b\t%s\t%s\t%s%n", label, alignmentStart,
                    readLength, unknownBases, describe(features), t.getClass().getSimpleName(),
                    escape(String.valueOf(t.getMessage())));
        }
    }

    /** One byte through the lookup every restored base goes through. */
    static void lookup(final int value) {
        try {
            final byte[] one = new byte[] {(byte) value};
            SequenceUtil.toBamReadBasesInPlace(one);
            System.out.printf("lookup\t%d\t%s%n", value, escape(new String(one)));
        } catch (final Throwable t) {
            System.out.printf("lookup\t%d\t%s%n", value, t.getClass().getSimpleName());
        }
    }

    /** A region over the whole of REFERENCE, already fetched. */
    static CRAMReferenceRegion region() {
        final SAMSequenceRecord sequence = new SAMSequenceRecord("chr1", REFERENCE.length());
        final SAMSequenceDictionary dictionary =
                new SAMSequenceDictionary(Collections.singletonList(sequence));
        final CRAMReferenceSource source = new CRAMReferenceSource() {
            @Override
            public byte[] getReferenceBases(final SAMSequenceRecord record,
                    final boolean tryNameVariants) {
                return REFERENCE.getBytes();
            }

            @Override
            public byte[] getReferenceBasesByRegion(final SAMSequenceRecord record,
                    final int zeroBasedStart, final int requestedRegionLength) {
                return Arrays.copyOfRange(REFERENCE.getBytes(), zeroBasedStart,
                        zeroBasedStart + requestedRegionLength);
            }
        };
        final CRAMReferenceRegion region = new CRAMReferenceRegion(source, dictionary);
        region.fetchReferenceBases(0);
        return region;
    }

    /** The feature list as text, so a golden row carries the input it came from. */
    static String describe(final List<ReadFeature> features) {
        if (features.isEmpty()) {
            return "-";
        }
        final StringJoiner joiner = new StringJoiner(",");
        for (final ReadFeature feature : features) {
            final List<String> parts = new ArrayList<>();
            parts.add(String.format("%c@%d", feature.getOperator(), feature.getPosition()));
            if (feature instanceof Deletion) {
                parts.add("len=" + ((Deletion) feature).getLength());
            } else if (feature instanceof RefSkip) {
                parts.add("len=" + ((RefSkip) feature).getLength());
            } else if (feature instanceof Padding) {
                parts.add("len=" + ((Padding) feature).getLength());
            } else if (feature instanceof HardClip) {
                parts.add("len=" + ((HardClip) feature).getLength());
            } else if (feature instanceof SoftClip) {
                parts.add("bases=" + new String(((SoftClip) feature).getSequence()));
            } else if (feature instanceof Insertion) {
                parts.add("bases=" + new String(((Insertion) feature).getSequence()));
            } else if (feature instanceof Bases) {
                parts.add("bases=" + new String(((Bases) feature).getBases()));
            } else if (feature instanceof Scores) {
                parts.add("scores=" + new String(((Scores) feature).getScores()));
            } else if (feature instanceof InsertBase) {
                parts.add("base=" + (((InsertBase) feature).getBase() & 0xFF));
            } else if (feature instanceof BaseQualityScore) {
                parts.add("quality=" + ((BaseQualityScore) feature).getQualityScore());
            } else if (feature instanceof ReadBase) {
                parts.add("base=" + (char) ((ReadBase) feature).getBase());
                parts.add("quality=" + ((ReadBase) feature).getQualityScore());
            } else if (feature instanceof Substitution) {
                parts.add("code=" + ((Substitution) feature).getCode());
            }
            joiner.add(String.join(" ", parts));
        }
        return joiner.toString();
    }

    /** Non-printable characters as \\uXXXX, so a golden stays a text file. */
    static String escape(final String text) {
        final StringBuilder b = new StringBuilder();
        for (int i = 0; i < text.length(); i++) {
            final char c = text.charAt(i);
            if (c < 0x20 || c > 0x7E) {
                b.append(String.format("\\u%04X", (int) c));
            } else {
                b.append(c);
            }
        }
        return b.length() == 0 ? "-" : b.toString();
    }

    static String hex(final byte[] bytes) {
        final StringBuilder b = new StringBuilder(bytes.length * 2);
        for (final byte value : bytes) {
            b.append(String.format("%02x", value));
        }
        return b.length() == 0 ? "-" : b.toString();
    }
}
