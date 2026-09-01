/*
 * Dumps each .bai twice: as bytes, and as the text TextualBAMIndexWriter prints for it.
 *
 * This is the read-side index, which differs from the write-side index (setCreateIndex): a chunk
 * that ends on a BGZF block boundary is recorded as (nextBlockAddress, 0) by the reading
 * getFilePointer, but as (blockAddress, blockLength) by the writing one. build_bam_index reproduces
 * the read-side form, so this is its oracle.
 *
 * Output, one line per case:
 *     bam <TAB> <case name> <TAB> <hex of the BAM>
 *     bai <TAB> <case name> <TAB> <hex of the read-side BAI>
 *
 * The JDK deflater is pinned per the oracle contract.
 */

import htsjdk.samtools.*;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;
import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

public class TextualIndexDump {

    static String hex(final byte[] bytes) {
        final StringBuilder sb = new StringBuilder();
        for (final byte b : bytes) sb.append(String.format("%02x", b));
        return sb.toString();
    }

    static SAMFileHeader header(final int... lengths) {
        final SAMFileHeader h = new SAMFileHeader();
        final SAMSequenceDictionary d = new SAMSequenceDictionary();
        for (int i = 0; i < lengths.length; i++) {
            d.addSequence(new SAMSequenceRecord("chr" + (i + 1), lengths[i]));
        }
        h.setSequenceDictionary(d);
        h.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        return h;
    }

    static SAMRecord read(final SAMFileHeader h, final String name, final int ref,
                          final int start, final String cigar) {
        final SAMRecord r = new SAMRecord(h);
        r.setReadName(name);
        r.setReferenceIndex(ref);
        if (ref >= 0) {
            r.setAlignmentStart(start);
            r.setCigarString(cigar);
        } else {
            r.setReadUnmappedFlag(true);
        }
        r.setMappingQuality(ref >= 0 ? 60 : 0);
        final int len = ref >= 0 ? TextCigarCodec.decode(cigar).getReadLength() : 4;
        final byte[] bases = new byte[len];
        java.util.Arrays.fill(bases, (byte) 'A');
        r.setReadBases(bases);
        final byte[] quals = new byte[len];
        java.util.Arrays.fill(quals, (byte) 30);
        r.setBaseQualities(quals);
        return r;
    }

    static void emit(final String name, final SAMFileHeader header, final List<SAMRecord> records)
            throws Exception {
        final File dir = Files.createTempDirectory("buildindex").toFile();
        final File bam = new File(dir, "in.bam");
        final SAMFileWriter w = new SAMFileWriterFactory()
                .setCreateIndex(false).setCreateMd5File(false).setUseAsyncIo(false)
                .makeBAMWriter(header, true, bam);
        for (final SAMRecord r : records) w.addAlignment(r);
        w.close();

        // The read side, exactly as picard BuildBamIndex does it: open with source-in-records so the
        // file pointers are available, then BAMIndexer.createIndex.
        final File bai = new File(dir, "in.bai");
        final SamReader reader = SamReaderFactory.makeDefault()
                .disable(SamReaderFactory.Option.EAGERLY_DECODE)
                .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS)
                .open(bam);
        BAMIndexer.createIndex(reader, bai);
        reader.close();

        // The text form, through the public door: createAndWriteIndex(input, output, true) is what
        // GATK's PrintFileDiagnostics reaches, and TextualBAMIndexWriter itself is package-private.
        final File text = new File(dir, "in.bai.txt");
        BAMIndexer.createAndWriteIndex(bai, text, true);

        System.out.println("bai\t" + name + "\t" + hex(Files.readAllBytes(bai.toPath())));
        for (final String line : Files.readAllLines(text.toPath())) {
            // A blank line is part of the format (an empty bin prints one), so it is emitted as an
            // explicit marker rather than lost to the line splitting on the other side.
            System.out.println("text\t" + name + "\t" + (line.isEmpty() ? "<blank>" : line));
        }
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        {
            final SAMFileHeader h = header(250000000);
            emit("one_read", h, List.of(read(h, "r1", 0, 100, "4M")));
        }
        {
            final SAMFileHeader h = header(250000000);
            emit("same_window", h, List.of(
                    read(h, "r1", 0, 100, "50M"),
                    read(h, "r2", 0, 200, "50M")));
        }
        {
            final SAMFileHeader h = header(250000000);
            emit("sparse_windows", h, List.of(
                    read(h, "r1", 0, 1, "50M"),
                    read(h, "r2", 0, 5 * 16384 + 1, "50M"),
                    read(h, "r3", 0, 40 * 16384 + 1, "50M")));
        }
        {
            final SAMFileHeader h = header(250000000, 200000000, 100000000);
            emit("multi_reference", h, List.of(
                    read(h, "a", 0, 100, "50M"),
                    read(h, "b", 2, 500, "50M")));
        }
        {
            final SAMFileHeader h = header(250000000);
            final List<SAMRecord> rs = new ArrayList<>();
            for (int i = 0; i < 4000; i++) rs.add(read(h, "r" + i, 0, 1 + i * 100, "50M"));
            emit("many_blocks", h, rs);
        }
        {
            // placed reads then unplaced (no-coordinate) reads, which bump the no-coord counter.
            final SAMFileHeader h = header(250000000);
            emit("with_unmapped", h, List.of(
                    read(h, "m", 0, 100, "50M"),
                    read(h, "u1", -1, 0, null),
                    read(h, "u2", -1, 0, null)));
        }
    }
}
