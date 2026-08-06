/*
 * The cigar rebuilt from read features: the first half of the way back.
 *
 * The forward direction is pinned: a record becomes read features, and everything matching the
 * reference becomes nothing at all. Going back, the cigar is not stored anywhere. It is rebuilt
 * from the feature positions and the read length by `getCigarForReadFeatures`, and the matches come
 * back as the GAPS BETWEEN FEATURES.
 *
 * Seven things here are decisions rather than layout.
 *
 *   - THE MATCHES ARE THE GAPS. `gap = feature.getPosition() - (lastOpPos + lastOpLen)` is the only
 *     source of M in the output. Nothing in the features says a base matched;
 *   - A FEATURE THAT CONSUMES NO READ BASES WINDS THE READ CURSOR BACK. `lastOpPos -=
 *     readFeatureLength` after a D, N or P, so the position bookkeeping is in read space and a
 *     reference-only operator must be undone from it;
 *   - THE SWITCH SILENTLY IGNORES WHAT IT DOES NOT NAME. `default: continue` drops
 *     BaseQualityScore, Scores and Bases, and Bases carries read bases: a feature list holding one
 *     produces a cigar that does not account for it;
 *   - A SUBSTITUTION AND A READBASE ARE BOTH M, so the rebuilt cigar cannot distinguish a match
 *     from a mismatch and never emits X or EQ, whatever the record started with;
 *   - THE TAIL HAS THREE BRANCHES, and only one of them appends the trailing M. Whether it appends
 *     anything depends on `readLength >= lastOpPos + lastOpLen`, in read coordinates that the
 *     previous rule has already been adjusting;
 *   - AN EMPTY FEATURE LIST IS ONE M OF THE WHOLE READ, and so is a list whose features all fell
 *     through the default branch: the empty-list check at the end catches both;
 *   - A READ LENGTH OF 0 IS A SPECIAL CASE inside the tail, taking the accumulated length rather
 *     than the read length.
 *
 * Output:
 *
 *     reference\t<the reference the round trips are measured against>
 *     case\t<label>\t<read length>\t<features>\t<rebuilt cigar>
 *     roundtrip\t<label>\t<original cigar>\t<bases>\t<start>\t<rebuilt cigar>\t<same|changed>
 *     err\t<label>\t<class>\t<message>
 *
 * Usage: CramRecordCigarDump
 */

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.TextCigarCodec;
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
import htsjdk.samtools.cram.structure.CRAMRecordReadFeatures;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.StringJoiner;

public class CramRecordCigarDump {

    private static final String REFERENCE = "ACGTACGTACGTMRWSACGTACGT";

    public static void main(final String[] args) {
        System.out.println("# CramRecordCigarDump: the cigar rebuilt from read features");

        System.out.printf("reference\t%s%n", REFERENCE);

        // Hand-built feature lists, so the rebuild is exercised on shapes the writer never emits
        // as well as on shapes it does.
        emit("empty", 8, Collections.emptyList());
        emit("one-substitution", 8, Arrays.asList(sub(4)));
        emit("substitution-at-one", 8, Arrays.asList(sub(1)));
        emit("substitution-at-end", 8, Arrays.asList(sub(8)));
        emit("two-substitutions", 8, Arrays.asList(sub(2), sub(7)));
        emit("adjacent-substitutions", 8, Arrays.asList(sub(4), sub(5)));
        emit("read-base-is-also-m", 8, Arrays.asList(new ReadBase(4, (byte) 'M', (byte) 40)));

        emit("insert-base", 8, Arrays.asList(new InsertBase(3, (byte) 'T')));
        emit("three-insert-bases", 8,
                Arrays.asList(new InsertBase(3, (byte) 'T'), new InsertBase(4, (byte) 'T'),
                        new InsertBase(5, (byte) 'T')));
        emit("insertion", 8, Arrays.asList(new Insertion(3, "TTT".getBytes())));
        emit("soft-clip-left", 8, Arrays.asList(new SoftClip(1, "TTT".getBytes())));
        emit("soft-clip-right", 8, Arrays.asList(new SoftClip(6, "TTT".getBytes())));
        emit("soft-clip-both", 8,
                Arrays.asList(new SoftClip(1, "TT".getBytes()), new SoftClip(7, "TT".getBytes())));
        emit("hard-clip", 8, Arrays.asList(new HardClip(1, 2)));

        // The features that consume reference and not read, which wind the read cursor back.
        emit("deletion", 8, Arrays.asList(new Deletion(5, 2)));
        emit("deletion-then-substitution", 8, Arrays.asList(new Deletion(5, 2), sub(5)));
        emit("ref-skip", 8, Arrays.asList(new RefSkip(5, 2)));
        emit("padding", 8, Arrays.asList(new Padding(5, 2)));
        emit("deletion-at-one", 8, Arrays.asList(new Deletion(1, 2)));

        // Everything the switch does not name.
        emit("base-quality-score-only", 8,
                Arrays.asList(new BaseQualityScore(4, (byte) 40)));
        emit("scores-only", 8, Arrays.asList(new Scores(1, "IIII".getBytes())));
        emit("bases-only", 8, Arrays.asList(new Bases(1, "ACGT".getBytes())));
        emit("bases-then-substitution", 8, Arrays.asList(new Bases(1, "ACGT".getBytes()), sub(6)));

        // The read length, which is the only thing that says where the read ends.
        emit("zero-read-length", 0, Arrays.asList(sub(4)));
        emit("read-length-shorter-than-features", 3, Arrays.asList(sub(6)));
        emit("read-length-one", 1, Arrays.asList(sub(1)));

        // And the round trip: a record's own cigar, through the features, and back.
        roundTrip("match", "8M", "ACGTACGT", 1);
        roundTrip("mismatch", "8M", "ACGTACGA", 1);
        roundTrip("all-mismatch", "8M", "TTTTTTTT", 1);
        roundTrip("insertion", "2M3I3M", "ACTTTTAC", 1);
        roundTrip("soft-clip", "3S5M", "TTTACGTA", 4);
        roundTrip("soft-clip-both-ends", "2S4M2S", "TTACGTTT", 3);
        roundTrip("deletion", "4M2D4M", "ACGTCGTA", 1);
        roundTrip("ref-skip", "4M2N4M", "ACGTCGTA", 1);
        roundTrip("padding", "4M2P4M", "ACGTACGT", 1);
        roundTrip("hard-clip", "2H8M", "ACGTACGT", 1);
        roundTrip("x-operator", "8X", "ACGTACGT", 1);
        roundTrip("eq-operator", "8=", "TTTTTTTT", 1);
        roundTrip("insertion-and-deletion", "2M2I2M2D2M", "ACTTACGT", 1);
    }

    static Substitution sub(final int position) {
        return new Substitution(position, (byte) 'T', (byte) 'A');
    }

    static void emit(final String label, final int readLength, final List<ReadFeature> features) {
        final CRAMRecordReadFeatures readFeatures = new CRAMRecordReadFeatures(features);
        try {
            final Cigar cigar = readFeatures.getCigarForReadFeatures(readLength);
            System.out.printf("case\t%s\t%d\t%s\t%s%n", label, readLength, describe(features),
                    cigar.toString());
        } catch (final Throwable t) {
            System.out.printf("err\t%s\t%s\t%s%n", label, t.getClass().getSimpleName(),
                    String.valueOf(t.getMessage()));
        }
    }

    /** A record's own cigar, through the features it produces, and back. */
    static void roundTrip(final String label, final String cigarText, final String bases,
            final int alignmentStart) {
        final SAMFileHeader header = new SAMFileHeader();
        header.addSequence(new SAMSequenceRecord("chr1", REFERENCE.length()));
        final SAMRecord record = new SAMRecord(header);
        record.setReadName("read");
        record.setReferenceIndex(0);
        record.setAlignmentStart(alignmentStart);
        record.setCigar(TextCigarCodec.decode(cigarText));
        record.setReadBases(bases.getBytes());
        record.setBaseQualityString("IIIIIIII".substring(0, bases.length()));

        final CRAMRecordReadFeatures features =
                new CRAMRecordReadFeatures(record, record.getReadBases(), REFERENCE.getBytes());
        final int readLength = Cigar.getReadLength(record.getCigar().getCigarElements());
        final String rebuilt = features.getCigarForReadFeatures(readLength).toString();
        System.out.printf("roundtrip\t%s\t%s\t%s\t%d\t%s\t%s%n", label, cigarText, bases,
                alignmentStart, rebuilt, cigarText.equals(rebuilt) ? "same" : "changed");
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
                parts.add("base=" + (char) ((InsertBase) feature).getBase());
            } else if (feature instanceof BaseQualityScore) {
                parts.add("quality=" + ((BaseQualityScore) feature).getQualityScore());
            } else if (feature instanceof ReadBase) {
                parts.add("base=" + (char) ((ReadBase) feature).getBase());
                parts.add("quality=" + ((ReadBase) feature).getQualityScore());
            } else if (feature instanceof Substitution) {
                parts.add("base=" + (char) ((Substitution) feature).getBase());
                parts.add("ref=" + (char) ((Substitution) feature).getReferenceBase());
            }
            joiner.add(String.join(" ", parts));
        }
        return joiner.toString();
    }
}
