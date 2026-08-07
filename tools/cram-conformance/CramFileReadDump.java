/*
 * A whole CRAM file, from its first byte to its records.
 *
 * Every piece has been measured on its own: the file definition, the container header, the
 * compression header's three maps, the slice header, the blocks, the codecs, the record reader.
 * This walks a real file through all of them at once, which is the only thing that shows whether
 * the pieces fit.
 *
 * Four things here are decisions rather than layout.
 *
 *   - A CONTAINER'S BLOCKS ARE COMPRESSED AND THE COMPRESSION HEADER'S IS NOT. The compression
 *     header block is raw by definition; a slice's external blocks carry whichever compressor the
 *     writer chose, so reading records means undoing three or four different ones;
 *   - THE SLICE'S ALIGNMENT START IS THE FIRST RECORD'S PREVIOUS ONE, so a container read out of
 *     order gives every position in it a different value;
 *   - THE EOF CONTAINER IS A CONTAINER, not a marker: it parses as one whose record count is zero,
 *     and a reader that stops on a byte pattern rather than on that count is reading the wrong
 *     thing;
 *   - THE RECORDS OF A SLICE ARE READ IN ONE PASS over shared streams, so a file cannot be read
 *     record by record out of order at all.
 *
 * Output:
 *
 *     file\t<name>\t<bytes>\t<version>\t<id>
 *     samheader\t<name>\t<sequences>\t<read groups>\t<sort order>
 *     container\t<name>\t<index>\t<byte offset>\t<length>\t<reference>\t<start>\t<span>\t<records>\t<blocks>
 *     block\t<name>\t<container>\t<index>\t<type>\t<method>\t<content id>\t<uncompressed hex>
 *     slice\t<name>\t<container>\t<reference>\t<start>\t<span>\t<records>\t<blocks>
 *     record\t<name>\t<index>\t<readName>\t<flags>\t<ref>\t<start>\t<cigar>\t<bases>\t<quals>
 *     err\t<what>\t<class>\t<message>
 *
 * Usage: CramFileReadDump [<cram file> <name> <reference fasta>]
 *
 * The three default to ce5.cram, ce5 and ce.fa, which are copied into the harness beside this
 * file: ce#5.2.1.cram from htsjdk's own test resources, renamed because a `#` in a path is not
 * worth the trouble, and the reference its reads were compressed against.
 */

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.cram.common.CramVersions;
import htsjdk.samtools.cram.structure.Container;
import htsjdk.samtools.cram.structure.CompressorCache;
import htsjdk.samtools.cram.structure.CramHeader;
import htsjdk.samtools.cram.structure.Slice;
import htsjdk.samtools.cram.structure.block.Block;
import htsjdk.samtools.cram.io.CountingInputStream;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

public class CramFileReadDump {

    public static void main(final String[] args) throws Exception {
        System.out.println("# CramFileReadDump: a whole CRAM file, from its first byte to its records");
        final String path = args.length > 0 ? args[0] : "/harness/ce5.cram";
        final String name = args.length > 1 ? args[1] : "ce5";

        final byte[] bytes = Files.readAllBytes(Paths.get(path));
        final CountingInputStream stream = new CountingInputStream(new ByteArrayInputStream(bytes));
        final CramHeader header = htsjdk.samtools.cram.build.CramIO.readCramHeader(stream);
        System.out.printf("file\t%s\t%d\t%s\t%s%n", name, bytes.length,
                header.getCRAMVersion().toString(), new String(header.getId()).trim());

        // The first container is the SAM header's, whose block is a FILE_HEADER rather than a
        // compression header. A reader that treats it as an ordinary container refuses it, which
        // is the first thing composing the pieces shows.
        final htsjdk.samtools.SAMFileHeader samHeader =
                Container.readSAMFileHeaderContainer(header.getCRAMVersion(), stream, name);
        System.out.printf("samheader\t%s\t%d\t%d\t%s%n", name,
                samHeader.getSequenceDictionary().size(), samHeader.getReadGroups().size(),
                samHeader.getSortOrder().name());

        final CompressorCache cache = new CompressorCache();
        int index = 0;
        while (true) {
            final long offset = stream.getCount();
            final Container container = new Container(header.getCRAMVersion(), stream, offset);
            System.out.printf("container\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d%n", name, index,
                    offset, container.getContainerHeader().getContainerBlocksByteSize(),
                    container.getAlignmentContext().getReferenceContext().getReferenceContextID(),
                    container.getAlignmentContext().getAlignmentStart(),
                    container.getAlignmentContext().getAlignmentSpan(),
                    container.getContainerHeader().getNumberOfRecords(),
                    container.getContainerHeader().getBlockCount());

            if (container.getContainerHeader().isEOF()) {
                break;
            }

            for (final Slice slice : container.getSlices()) {
                System.out.printf("slice\t%s\t%d\t%d\t%d\t%d\t%d\t%d%n", name, index,
                        slice.getAlignmentContext().getReferenceContext().getReferenceContextID(),
                        slice.getAlignmentContext().getAlignmentStart(),
                        slice.getAlignmentContext().getAlignmentSpan(),
                        slice.getNumberOfRecords(),
                        slice.getSliceBlocks().getNumberOfExternalBlocks() + 1);

                final Block core = slice.getSliceBlocks().getCoreBlock();
                System.out.printf("block\t%s\t%d\tcore\t%s\t%s\t%d\t%s%n", name, index,
                        core.getContentType().name(), core.getCompressionMethod().name(),
                        core.getContentId(), hex(core.getUncompressedContent(cache)));
                for (final Integer contentId : slice.getSliceBlocks().getExternalContentIDs()) {
                    final Block block = slice.getSliceBlocks().getExternalBlock(contentId);
                    System.out.printf("block\t%s\t%d\t%d\t%s\t%s\t%d\t%s%n", name, index, contentId,
                            block.getContentType().name(), block.getCompressionMethod().name(),
                            block.getContentId(), hex(block.getUncompressedContent(cache)));
                }
            }
            index++;
        }

        // And the records, through the reader that composes all of it.
        final String reference = args.length > 2 ? args[2] : "/harness/ce.fa";
        SamReaderFactory factory = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT);
        if (reference != null) {
            factory = factory.referenceSequence(Paths.get(reference));
        }
        try (final SamReader reader = factory.open(new File(path))) {
            int recordIndex = 0;
            for (final SAMRecord record : reader) {
                System.out.printf("record\t%s\t%d\t%s\t%d\t%d\t%d\t%s\t%s\t%s%n", name,
                        recordIndex++, record.getReadName(), record.getFlags(),
                        record.getReferenceIndex(), record.getAlignmentStart(),
                        record.getCigarString(), record.getReadString(),
                        record.getBaseQualityString().replace('\t', ' '));
            }
        }
    }

    static String hex(final byte[] bytes) {
        final StringBuilder builder = new StringBuilder(bytes.length * 2);
        for (final byte value : bytes) {
            builder.append(String.format("%02x", value));
        }
        return builder.length() == 0 ? "-" : builder.toString();
    }
}
