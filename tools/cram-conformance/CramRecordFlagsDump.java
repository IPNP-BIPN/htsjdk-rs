/*
 * The three flag words a CRAM record carries, and the mate chain that fills two of them.
 *
 * A BAM record has one flag word. A CRAM record has three: the BAM flags it will become, a CRAM
 * flags byte the format adds, and a mate flags byte that duplicates two bits of the first one at
 * different positions. Two of the three are one byte wide and masked on the way out.
 *
 * Seven things here are decisions rather than layout.
 *
 *   - MATE UNMAPPED AND MATE REVERSE ARE STORED TWICE, at 0x2 and 0x1 in the mate flags and at 0x8
 *     and 0x20 in the BAM flags. Setting either through the record sets both, and nothing keeps
 *     them in step if the two words are set independently;
 *   - THE TWO SMALL WORDS ARE MASKED TO A BYTE on the way out, so a value above 255 stored in them
 *     comes back truncated rather than refused;
 *   - RESTORING MATE INFO WALKS A RING. Each record takes its mate's position, reference and two
 *     flags from the next; the last takes them from the first, so a chain of N becomes a cycle;
 *   - A MATE ON NO REFERENCE LOSES ITS POSITION. If the mate's reference index is -1 the mate
 *     alignment start is forced to 0, whatever the mate's own start was;
 *   - THE TEMPLATE LENGTH IS COMPUTED ONCE AND NEGATED ONCE, on the first and last of the chain.
 *     Every record between them keeps whatever it was constructed with;
 *   - THE LENGTH IS ZERO IF EITHER END IS UNMAPPED OR THE TWO ARE ON DIFFERENT REFERENCES, and
 *     otherwise it is the distance between 5' ends plus a sign, so the shortest possible template
 *     is 1 and never 0;
 *   - THE 5' END DEPENDS ON THE STRAND: a negative-strand record measures from its alignment end,
 *     which is derived from its read features rather than stored.
 *
 * Output:
 *
 *     flags\t<bam>\t<cram>\t<mate>\t<predicates>
 *     mask\t<field>\t<stored>\t<returned>
 *     chain\t<label>\t<index>\t<bam>\t<mate flags>\t<mate ref>\t<mate start>\t<template size>
 *     insert\t<label>\t<first template size>\t<last template size>
 *     detach\t<cram before>\t<records to next before>\t<cram after>\t<records to next after>
 *
 * Usage: CramRecordFlagsDump
 */

import htsjdk.samtools.cram.structure.CRAMCompressionRecord;

import java.util.ArrayList;
import java.util.List;
import java.util.StringJoiner;

public class CramRecordFlagsDump {

    public static void main(final String[] args) {
        System.out.println("# CramRecordFlagsDump: three flag words and the mate chain");

        // Every bit of the CRAM flags byte, and every bit of the mate flags byte.
        for (final int cram : new int[] {0, 1, 2, 4, 8, 15, 3, 5, 6}) {
            flags(0, cram, 0);
        }
        for (final int mate : new int[] {0, 1, 2, 3}) {
            flags(0, 0, mate);
        }
        // And the BAM flags the record reads for itself.
        for (final int bam : new int[] {0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048,
                0x900, 0xFFF}) {
            flags(bam, 0, 0);
        }

        // The two narrow words are masked to a byte on the way out.
        mask(0x1FF, 0x1FF);
        mask(0x100, 0x100);
        mask(-1, -1);

        // A chain of records, each pointing at the next, restored into a ring.
        chain("pair-forward", record(0, 100, 0, 50, false), record(0, 200, 1, 50, false));
        chain("pair-reverse-second", record(0, 100, 0, 50, false), record(0, 200, 1, 50, true));
        chain("pair-same-start", record(0, 100, 0, 50, false), record(0, 100, 1, 50, false));
        chain("pair-second-before", record(0, 200, 0, 50, false), record(0, 100, 1, 50, false));
        chain("pair-other-reference", record(0, 100, 0, 50, false), record(1, 200, 1, 50, false));
        chain("pair-unmapped-mate", record(0, 100, 0, 50, false), unmapped(1));
        chain("pair-mate-no-reference", record(0, 100, 0, 50, false), noReference(1));
        chain("triple", record(0, 100, 0, 50, false), record(0, 200, 1, 50, false),
                record(0, 300, 2, 50, true));
        chain("single", record(0, 100, 0, 50, false));

        // What detaching does to the CRAM flags and to the distance to the next fragment.
        detach(0, 5);
        detach(6, 5);
        detach(2, -1);
    }

    static void flags(final int bam, final int cram, final int mate) {
        final CRAMCompressionRecord record = new CRAMCompressionRecord(0, bam, cram, "r", 10, 0,
                100, 0, 30, new byte[10], new byte[10], null, null, -1, mate, -1, 0, -1);
        final StringJoiner joiner = new StringJoiner(",");
        joiner.add("detached=" + record.isDetached());
        joiner.add("mateDownstream=" + record.isHasMateDownStream());
        joiner.add("forceQuality=" + record.isForcePreserveQualityScores());
        joiner.add("unknownBases=" + record.isUnknownBases());
        joiner.add("paired=" + record.isReadPaired());
        joiner.add("unmapped=" + record.isSegmentUnmapped());
        joiner.add("first=" + record.isFirstSegment());
        joiner.add("last=" + record.isLastSegment());
        joiner.add("secondary=" + record.isSecondaryAlignment());
        System.out.printf("flags\t%d\t%d\t%d\t%s%n", bam, cram, mate, joiner);
    }

    static void mask(final int cram, final int mate) {
        final CRAMCompressionRecord record = new CRAMCompressionRecord(0, 0, cram, "r", 10, 0, 100,
                0, 30, new byte[10], new byte[10], null, null, -1, mate, -1, 0, -1);
        System.out.printf("mask\tcram\t%d\t%d%n", cram, record.getCRAMFlags());
        System.out.printf("mask\tmate\t%d\t%d%n", mate, record.getMateFlags());
    }

    /** A mapped record on a reference, with no read features so its end follows from its length. */
    static CRAMCompressionRecord record(final int referenceIndex, final int alignmentStart,
            final int index, final int readLength, final boolean negativeStrand) {
        final int bamFlags = 0x1 | (negativeStrand ? 0x10 : 0);
        return new CRAMCompressionRecord(index, bamFlags, 0, "r" + index, readLength,
                referenceIndex, alignmentStart, 0, 30, new byte[readLength], new byte[readLength],
                null, null, -1, 0, -1, 0, -1);
    }

    /** Unmapped, which is what makes a template length zero however the two are placed. */
    static CRAMCompressionRecord unmapped(final int index) {
        return new CRAMCompressionRecord(index, 0x1 | 0x4, 0, "r" + index, 50, 0, 200, 0, 0,
                new byte[50], new byte[50], null, null, -1, 0, -1, 0, -1);
    }

    /** Mapped by its flags but on no reference, which is what empties the mate's start. */
    static CRAMCompressionRecord noReference(final int index) {
        return new CRAMCompressionRecord(index, 0x1, 0, "r" + index, 50, -1, 200, 0, 30,
                new byte[50], new byte[50], null, null, -1, 0, -1, 0, -1);
    }

    static void chain(final String label, final CRAMCompressionRecord... records) {
        final List<CRAMCompressionRecord> list = new ArrayList<>();
        for (final CRAMCompressionRecord record : records) {
            list.add(record);
        }
        for (int i = 0; i + 1 < list.size(); i++) {
            list.get(i).setNextSegment(list.get(i + 1));
            list.get(i + 1).setPreviousSegment(list.get(i));
        }
        list.get(0).restoreMateInfo();

        for (int i = 0; i < list.size(); i++) {
            final CRAMCompressionRecord record = list.get(i);
            System.out.printf("chain\t%s\t%d\t%d\t%d\t%d\t%d\t%d%n", label, i, record.getBAMFlags(),
                    record.getMateFlags(), record.getMateReferenceIndex(),
                    record.getMateAlignmentStart(), record.getTemplateSize());
        }
        System.out.printf("insert\t%s\t%d\t%d%n", label, list.get(0).getTemplateSize(),
                list.get(list.size() - 1).getTemplateSize());
    }

    static void detach(final int cramFlags, final int recordsToNext) {
        final CRAMCompressionRecord record = new CRAMCompressionRecord(0, 0x1, cramFlags, "r", 10,
                0, 100, 0, 30, new byte[10], new byte[10], null, null, -1, 0, -1, 0, recordsToNext);
        final int before = record.getCRAMFlags();
        final int beforeNext = record.getRecordsToNextFragment();
        record.setToDetachedState();
        System.out.printf("detach\t%d\t%d\t%d\t%d%n", before, beforeNext, record.getCRAMFlags(),
                record.getRecordsToNextFragment());
    }
}
