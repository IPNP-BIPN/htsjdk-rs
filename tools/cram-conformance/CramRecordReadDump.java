/*
 * Reading a CRAM record: the data series, in the order the specification prescribes.
 *
 * Every field of a record comes from its own data series, and the series share streams. So the
 * order the reads happen in is not a detail of the implementation: read two series in the wrong
 * order and both come back wrong, with nothing to say so. This dump writes records through the
 * reference's own writer, hands back the bytes, and records what its own reader made of them.
 *
 * Six things here are decisions rather than layout.
 *
 *   - THE ORDER IS PRESCRIBED. Flags, then the reference for a multi-reference slice, then read
 *     length, alignment start, read group, read name, the mate block, the tag list, and only then
 *     the read features;
 *   - THE ALIGNMENT START IS A DELTA when the compression header says so, and the delta may be
 *     negative. The first record's previous start is the one the caller passes in;
 *   - THE READ NAME MOVES. Preserved names are read before the mate block, generated ones after
 *     the mate flags inside it, and the specification says so explicitly;
 *   - THE MATE FLAGS ARE PROPAGATED INTO THE BAM FLAGS, because a writer is not required to have
 *     put them there. Two bits, at different positions in the two words;
 *   - THE TAG LIST IS AN INDEX INTO THE HEADER'S DICTIONARY, and the tags are then read in the
 *     dictionary's order, one series per tag id;
 *   - AN UNMAPPED RECORD READS NO READ FEATURES AT ALL, and its bases and scores come from
 *     elsewhere.
 *
 * Output:
 *
 *     header\t<hex of the compression header block>
 *     block\t<label>\t<core|content id>\t<compression method>\t<uncompressed hex>
 *     start\t<label>\t<the slice's alignment start>\t<its reference context>
 *     record\t<label>\t<index>\t<bam>\t<cram>\t<mate>\t<name>\t<length>\t<ref>\t<start>\t<mateRef>\t<mateStart>\t<template>\t<mq>\t<rg>\t<features>\t<tags>\t<bases>\t<scores>
 *     err\t<what>\t<class>\t<message>
 *
 * Usage: CramRecordReadDump
 */

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.cram.build.CompressionHeaderFactory;
import htsjdk.samtools.cram.common.CramVersions;
import htsjdk.samtools.cram.encoding.reader.CramRecordReader;
import htsjdk.samtools.cram.encoding.readfeatures.Deletion;
import htsjdk.samtools.cram.encoding.readfeatures.InsertBase;
import htsjdk.samtools.cram.encoding.readfeatures.Insertion;
import htsjdk.samtools.cram.encoding.readfeatures.ReadFeature;
import htsjdk.samtools.cram.encoding.readfeatures.SoftClip;
import htsjdk.samtools.cram.structure.CRAMCompressionRecord;
import htsjdk.samtools.cram.structure.CRAMEncodingStrategy;
import htsjdk.samtools.cram.structure.CompressionHeader;
import htsjdk.samtools.cram.structure.CompressorCache;
import htsjdk.samtools.cram.structure.Slice;

import java.io.ByteArrayOutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.StringJoiner;

public class CramRecordReadDump {

    public static void main(final String[] args) {
        System.out.println("# CramRecordReadDump: the data series a record is read from");

        // A single mapped record with no read features at all.
        roundTrip("plain", Arrays.asList(mapped(0, 100, "r0", 10, null)));

        // Two records on the same reference, which is where the alignment start delta shows.
        roundTrip("delta", Arrays.asList(mapped(0, 100, "r0", 10, null),
                mapped(0, 140, "r1", 10, null)));

        // A start that goes backwards, which is legal and gives a negative delta.
        roundTrip("delta-negative", Arrays.asList(mapped(0, 200, "r0", 10, null),
                mapped(0, 100, "r1", 10, null)));

        // The read features the writer can encode without a reference to derive a substitution
        // from. A Substitution carries a code rather than a base, and building the compression
        // header's substitution matrix from one refuses: "Attempt to generate a substitution code
        // for invalid reference base with value '-1'".
        roundTrip("features", Arrays.asList(mapped(0, 100, "r0", 10, Arrays.asList(
                new InsertBase(4, (byte) 'A'),
                new Insertion(5, new byte[] {'C', 'G'}),
                new Deletion(7, 3),
                new SoftClip(8, new byte[] {'T', 'T'})))));

        // An unmapped record, which reads no read features whatever it carries.
        roundTrip("unmapped", Arrays.asList(unmapped("r0", 10)));

        // Quality scores kept as an array, which is a data series the other cases never touch.
        roundTrip("scores-preserved", Arrays.asList(withQualityArray(0, 100, "r0", 10)));

        // An unmapped record reads its bases one at a time from the same series a following
        // record's InsertBase reads from, so a reader that skips them desynchronises here and
        // nowhere else.
        roundTrip("unmapped-then-mapped", Arrays.asList(unmapped("r0", 10),
                mapped(0, 150, "r1", 12, Arrays.asList(new InsertBase(3, (byte) 'G')))));

        // A mixture, which is what a real slice holds.
        roundTrip("mixed", Arrays.asList(mapped(0, 100, "r0", 10, null), unmapped("r1", 10),
                mapped(0, 150, "r2", 12, Arrays.asList(new InsertBase(3, (byte) 'G')))));
    }

    /** A mapped, unpaired record. Detached, which is how every record is written. */
    static CRAMCompressionRecord mapped(final int referenceIndex, final int alignmentStart,
            final String name, final int readLength, final List<ReadFeature> features) {
        final byte[] bases = new byte[readLength];
        Arrays.fill(bases, (byte) 'A');
        final byte[] scores = new byte[readLength];
        Arrays.fill(scores, (byte) 30);
        return new CRAMCompressionRecord(0, 0, CRAMCompressionRecord.CF_DETACHED, name, readLength,
                referenceIndex, alignmentStart, 0, 40, scores, bases, null, features, -1, 0,
                SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX, 0, -1);
    }

    /** Mapped, with CF_QS_PRESERVED_AS_ARRAY set so the scores go out as one array. */
    static CRAMCompressionRecord withQualityArray(final int referenceIndex,
            final int alignmentStart, final String name, final int readLength) {
        final byte[] bases = new byte[readLength];
        Arrays.fill(bases, (byte) 'A');
        final byte[] scores = new byte[readLength];
        for (int i = 0; i < readLength; i++) {
            scores[i] = (byte) (10 + i);
        }
        return new CRAMCompressionRecord(0, 0,
                CRAMCompressionRecord.CF_DETACHED | CRAMCompressionRecord.CF_QS_PRESERVED_AS_ARRAY,
                name, readLength, referenceIndex, alignmentStart, 0, 40, scores, bases, null, null,
                -1, 0, SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX, 0, -1);
    }

    static CRAMCompressionRecord unmapped(final String name, final int readLength) {
        final byte[] bases = new byte[readLength];
        Arrays.fill(bases, (byte) 'C');
        final byte[] scores = new byte[readLength];
        Arrays.fill(scores, (byte) 20);
        return new CRAMCompressionRecord(0, 0x4, CRAMCompressionRecord.CF_DETACHED, name,
                readLength, SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX, 0, 0, 0, scores, bases, null,
                null, -1, 0, SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX, 0, -1);
    }

    /** Write the records through the reference's writer, then read them back through its reader. */
    static void roundTrip(final String label, final List<CRAMCompressionRecord> records) {
        try {
            final CompressionHeader compressionHeader =
                    new CompressionHeaderFactory(new CRAMEncodingStrategy())
                            .createCompressionHeader(records, true);
            final Slice slice = new Slice(records, compressionHeader, 0L, 0L);

            final ByteArrayOutputStream headerOut = new ByteArrayOutputStream();
            compressionHeader.write(CramVersions.CRAM_v3, headerOut);
            System.out.printf("header\t%s\t%s%n", label, hex(headerOut.toByteArray()));

            // Each block's content after its compression is undone, which is what the codecs see.
            // The compression itself is another suite's business, and carrying the compressed
            // stream here would make this suite depend on it.
            final CompressorCache cache = new CompressorCache();
            System.out.printf("block\t%s\tcore\tRAW\t%s%n", label,
                    hex(slice.getSliceBlocks().getCoreBlock().getRawContent()));
            for (final Integer contentId : slice.getSliceBlocks().getExternalContentIDs()) {
                final htsjdk.samtools.cram.structure.block.Block block =
                        slice.getSliceBlocks().getExternalBlock(contentId);
                System.out.printf("block\t%s\t%d\t%s\t%s%n", label, contentId,
                        block.getCompressionMethod().name(),
                        hex(block.getUncompressedContent(cache)));
            }

            // The reader is handed the slice's own alignment start as the previous one, which is
            // what Slice.getRecords does. The row carries it so nothing has to be assumed.
            final CramRecordReader reader = new CramRecordReader(slice, new CompressorCache(),
                    ValidationStringency.SILENT);
            int previousStart = slice.getAlignmentContext().getAlignmentStart();
            System.out.printf("start\t%s\t%d\t%s%n", label, previousStart,
                    slice.getAlignmentContext().getReferenceContext().toString());
            for (int i = 0; i < records.size(); i++) {
                final CRAMCompressionRecord back = reader.readCRAMRecord(i, previousStart);
                previousStart = back.getAlignmentStart();
                System.out.printf(
                        "record\t%s\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s%n",
                        label, i, back.getBAMFlags(), back.getCRAMFlags(), back.getMateFlags(),
                        String.valueOf(back.getReadName()), back.getReadLength(),
                        back.getReferenceIndex(), back.getAlignmentStart(),
                        back.getMateReferenceIndex(), back.getMateAlignmentStart(),
                        back.getTemplateSize(), back.getMappingQuality(), back.getReadGroupID(),
                        features(back.getReadFeatures()), tags(back), hex(back.getReadBases()),
                        hex(back.getQualityScores()));
            }
        } catch (final Throwable t) {
            System.out.printf("err\t%s\t%s\t%s%n", label, t.getClass().getSimpleName(),
                    String.valueOf(t.getMessage()));
        }
    }

    static String features(final List<ReadFeature> features) {
        if (features == null || features.isEmpty()) {
            return "-";
        }
        final StringJoiner joiner = new StringJoiner(";");
        for (final ReadFeature feature : features) {
            joiner.add(String.format("%c@%d", feature.getOperator(), feature.getPosition()));
        }
        return joiner.toString();
    }

    static String tags(final CRAMCompressionRecord record) {
        if (record.getTags() == null || record.getTags().isEmpty()) {
            return "-";
        }
        final List<String> names = new ArrayList<>();
        record.getTags().forEach(tag -> names.add(tag.getKey()));
        return String.join(";", names);
    }

    static String hex(final byte[] bytes) {
        final StringBuilder builder = new StringBuilder(bytes.length * 2);
        for (final byte value : bytes) {
            builder.append(String.format("%02x", value));
        }
        return builder.length() == 0 ? "-" : builder.toString();
    }
}
