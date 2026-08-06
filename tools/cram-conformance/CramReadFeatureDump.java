/*
 * How a record becomes read features: the first half of the CRAM record model.
 *
 * The frames are pinned through the slice header. A slice's records are not stored as bases and a
 * cigar; they are stored as an alignment start, a read length, and a list of READ FEATURES, each a
 * one-letter operator and a payload. Everything that matches the reference is stored as nothing at
 * all, which is where CRAM's compression actually comes from.
 *
 * Eight things here are decisions rather than layout.
 *
 *   - THE POSITIONS ARE ONE-BASED, AND THE INTERFACE SAYS THEY ARE NOT. Every construction site
 *     passes `zeroBasedPositionInRead + 1`, while `ReadFeature.getPosition`'s javadoc says
 *     "zero-based position in the read". The doc is wrong about all twelve implementations;
 *   - AN INSERTION OF n BASES BECOMES n InsertBase FEATURES, not one Insertion. htsjdk's own
 *     comment says it should use a Bases feature and does not, because that would need a
 *     ByteArrayLenEncoding and therefore a frequency distribution over lengths. So the Insertion
 *     feature exists and the writer never emits it;
 *   - A SOFT CLIP OF n BASES BECOMES ONE SoftClip FEATURE carrying all n. The opposite decision to
 *     the insertion, in the same loop, five lines apart;
 *   - A MISMATCH SPLITS TWO WAYS, and the split is on the ALPHABET, not on the cigar. An
 *     ACGTN-to-ACGTN mismatch is a Substitution; anything else is a ReadBase, which carries the
 *     quality score a second time even when the scores are already preserved as an array;
 *   - THE CIGAR'S OWN CLAIM IS IGNORED. M, X and EQ all go through the same comparison, so an X
 *     over bases that match emits NOTHING and an EQ over bases that differ emits a substitution.
 *     The bases decide, and the cigar operator only says how far to walk;
 *   - A READ RUNNING PAST THE END OF THE REFERENCE READS 'N' rather than failing, so every base
 *     out there mismatches into a substitution against a reference base of N;
 *   - SEQ="*" MANUFACTURES 'N's, one per read base the cigar consumes, and those Ns then mismatch
 *     the reference like any other base: a record with no sequence at all produces a substitution
 *     per position;
 *   - THE MISSING-QUALITY TEST IS AN IDENTITY TEST. `baseQualities.equals(SAMRecord.NULL_QUALS)`
 *     is Object.equals on a byte[], so it is true only for that one array instance. A record whose
 *     qualities are an equal but distinct empty array takes the other branch and indexes it.
 *
 * Output:
 *
 *     op\t<class>\t<operator letter>
 *     consts\t<missing quality score>
 *     case\t<label>\t<cigar>\t<bases>\t<quals>\t<alignment start>\t<feature count>
 *     feature\t<label>\t<index>\t<operator>\t<position>\t<payload>
 *     alignend\t<label>\t<alignment start>\t<read length>\t<alignment end>
 *     err\t<label>\t<cigar>\t<bases>\t<quals>\t<alignment start>\t<class>\t<message>
 *
 * Usage: CramReadFeatureDump
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
import htsjdk.samtools.cram.structure.CRAMCompressionRecord;
import htsjdk.samtools.cram.structure.CRAMRecordReadFeatures;

import java.util.List;

public class CramReadFeatureDump {

    /** A reference long enough for every case, and deliberately not all ACGT. */
    private static final String REFERENCE = "ACGTACGTACGTMRWSACGTACGT";

    public static void main(final String[] args) {
        System.out.println("# CramReadFeatureDump: how a record becomes read features");

        // The operator letters, taken from the classes rather than from the specification.
        System.out.printf("op\tBaseQualityScore\t%c%n", BaseQualityScore.operator);
        System.out.printf("op\tBases\t%c%n", Bases.operator);
        System.out.printf("op\tDeletion\t%c%n", Deletion.operator);
        System.out.printf("op\tHardClip\t%c%n", HardClip.operator);
        System.out.printf("op\tInsertBase\t%c%n", InsertBase.operator);
        System.out.printf("op\tInsertion\t%c%n", Insertion.operator);
        System.out.printf("op\tPadding\t%c%n", Padding.operator);
        System.out.printf("op\tReadBase\t%c%n", ReadBase.operator);
        System.out.printf("op\tRefSkip\t%c%n", RefSkip.operator);
        System.out.printf("op\tScores\t%c%n", Scores.operator);
        System.out.printf("op\tSoftClip\t%c%n", SoftClip.operator);
        System.out.printf("op\tSubstitution\t%c%n", Substitution.operator);

        System.out.printf("consts\t%d%n", CRAMCompressionRecord.MISSING_QUALITY_SCORE);
        System.out.printf("reference\t%s%n", REFERENCE);

        // Everything that matches is stored as nothing.
        emit("perfect-match", "8M", "ACGTACGT", "IIIIIIII", 1);

        // An ACGTN mismatch on both sides is a substitution; the position is one-based.
        emit("one-substitution", "8M", "ACGTACGA", "IIIIIIII", 1);
        emit("first-base-substitution", "8M", "CCGTACGT", "IIIIIIII", 1);
        emit("every-base-substitution", "4M", "TTTT", "IIII", 1);

        // A mismatch where either side is outside ACGTN is a ReadBase, which carries the quality
        // score a second time.
        emit("read-base-not-acgtn", "4M", "ACGM", "IIII", 1);
        emit("reference-not-acgtn", "4M", "ACGT", "IIII", 13);

        // An insertion becomes one feature per base; a soft clip becomes one feature for all.
        emit("insertion-of-three", "2M3I3M", "ACTTTTAC", "IIIIIIII", 1);
        emit("insertion-of-one", "2M1I5M", "ACTTACGT", "IIIIIIII", 1);
        emit("soft-clip-of-three", "3S5M", "TTTACGTA", "IIIIIIII", 4);
        emit("soft-clip-both-ends", "2S4M2S", "TTACGTTT", "IIIIIIII", 3);

        // The features that carry a length rather than bases.
        emit("deletion", "4M2D4M", "ACGTCGTA", "IIIIIIII", 1);
        emit("ref-skip", "4M2N4M", "ACGTCGTA", "IIIIIIII", 1);
        emit("padding", "4M2P4M", "ACGTACGT", "IIIIIIII", 1);
        emit("hard-clip", "2H8M", "ACGTACGT", "IIIIIIII", 1);

        // The cigar's own claim about match and mismatch is not consulted.
        emit("x-over-matching-bases", "8X", "ACGTACGT", "IIIIIIII", 1);
        emit("eq-over-mismatching-bases", "8=", "TTTTTTTT", "IIIIIIII", 1);

        // Past the end of the reference every base is compared against 'N'.
        emit("past-the-end", "8M", "ACGTACGT", "IIIIIIII", 20);
        emit("n-past-the-end", "4M", "NNNN", "IIII", 22);

        // SEQ="*": the bases are manufactured, and then they mismatch.
        emitNoBases("no-sequence", "4M", 1);

        // An alignment start of 0 puts the reference cursor before the start of the reference.
        emit("alignment-start-zero", "4M", "ACGT", "IIII", 0);

        // The missing-quality test compares array identity, not contents.
        emitNullQuals("null-quals-non-acgtn", "4M", "ACGM", 1);
        emitEmptyQuals("empty-quals-non-acgtn", "4M", "ACGM", 1);
    }

    static SAMRecord record(final String cigar, final String bases, final int alignmentStart) {
        final SAMFileHeader header = new SAMFileHeader();
        header.addSequence(new SAMSequenceRecord("chr1", REFERENCE.length()));
        final SAMRecord record = new SAMRecord(header);
        record.setReadName("read");
        record.setReferenceIndex(0);
        record.setAlignmentStart(alignmentStart);
        record.setCigar(TextCigarCodec.decode(cigar));
        if (bases != null) {
            record.setReadBases(bases.getBytes());
        }
        return record;
    }

    static void emit(final String label, final String cigar, final String bases,
            final String quals, final int alignmentStart) {
        final SAMRecord record = record(cigar, bases, alignmentStart);
        record.setBaseQualityString(quals);
        run(label, record, cigar, bases, quals, alignmentStart);
    }

    /** A record with SEQ="*", whose bases htsjdk manufactures as 'N's. */
    static void emitNoBases(final String label, final String cigar, final int alignmentStart) {
        final SAMRecord record = record(cigar, null, alignmentStart);
        record.setBaseQualityString("IIII");
        run(label, record, cigar, "*", "IIII", alignmentStart);
    }

    /** QUAL="*", which leaves the qualities as the NULL_QUALS singleton. */
    static void emitNullQuals(final String label, final String cigar, final String bases,
            final int alignmentStart) {
        final SAMRecord record = record(cigar, bases, alignmentStart);
        record.setBaseQualityString("*");
        run(label, record, cigar, bases, "*", alignmentStart);
    }

    /** An empty quality array that is not the singleton: equal contents, different identity. */
    static void emitEmptyQuals(final String label, final String cigar, final String bases,
            final int alignmentStart) {
        final SAMRecord record = record(cigar, bases, alignmentStart);
        record.setBaseQualities(new byte[0]);
        run(label, record, cigar, bases, "<empty array>", alignmentStart);
    }

    static void run(final String label, final SAMRecord record, final String cigar,
            final String bases, final String quals, final int alignmentStart) {
        final byte[] readBases = record.getReadBases();
        final CRAMRecordReadFeatures features;
        try {
            features = new CRAMRecordReadFeatures(record, readBases, REFERENCE.getBytes());
        } catch (final Throwable t) {
            System.out.printf("err\t%s\t%s\t%s\t%s\t%d\t%s\t%s%n", label, cigar, bases, quals,
                    alignmentStart, t.getClass().getSimpleName(), String.valueOf(t.getMessage()));
            return;
        }

        final List<ReadFeature> list = features.getReadFeaturesList();
        System.out.printf("case\t%s\t%s\t%s\t%s\t%d\t%d%n", label, cigar, bases, quals,
                alignmentStart, list.size());
        for (int i = 0; i < list.size(); i++) {
            final ReadFeature feature = list.get(i);
            System.out.printf("feature\t%s\t%d\t%c\t%d\t%s%n", label, i, feature.getOperator(),
                    feature.getPosition(), payload(feature));
        }

        final int readLength = Cigar.getReadLength(record.getCigar().getCigarElements());
        System.out.printf("alignend\t%s\t%d\t%d\t%d%n", label, alignmentStart, readLength,
                features.getAlignmentEnd(alignmentStart, readLength));
    }

    /** Everything the feature carries beyond its position. */
    static String payload(final ReadFeature feature) {
        if (feature instanceof Substitution) {
            final Substitution substitution = (Substitution) feature;
            // The code is assigned later, from the substitution matrix; here it is still unset.
            return String.format("base=%c ref=%c code=%d", substitution.getBase(),
                    substitution.getReferenceBase(), substitution.getCode());
        }
        if (feature instanceof ReadBase) {
            final ReadBase readBase = (ReadBase) feature;
            return String.format("base=%c quality=%d", readBase.getBase(),
                    readBase.getQualityScore());
        }
        if (feature instanceof InsertBase) {
            return String.format("base=%c", ((InsertBase) feature).getBase());
        }
        if (feature instanceof SoftClip) {
            return "bases=" + new String(((SoftClip) feature).getSequence());
        }
        if (feature instanceof Insertion) {
            return "bases=" + new String(((Insertion) feature).getSequence());
        }
        if (feature instanceof Bases) {
            return "bases=" + new String(((Bases) feature).getBases());
        }
        if (feature instanceof Scores) {
            return "scores=" + new String(((Scores) feature).getScores());
        }
        if (feature instanceof BaseQualityScore) {
            return "quality=" + ((BaseQualityScore) feature).getQualityScore();
        }
        if (feature instanceof Deletion) {
            return "length=" + ((Deletion) feature).getLength();
        }
        if (feature instanceof RefSkip) {
            return "length=" + ((RefSkip) feature).getLength();
        }
        if (feature instanceof Padding) {
            return "length=" + ((Padding) feature).getLength();
        }
        if (feature instanceof HardClip) {
            return "length=" + ((HardClip) feature).getLength();
        }
        return "?";
    }
}
