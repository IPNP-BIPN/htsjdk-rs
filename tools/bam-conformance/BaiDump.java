/*
 * Dumps BAI indices produced by htsjdk alongside the BAM files they index.
 *
 * Output, one line per case:
 *
 *     bam <TAB> <case name> <TAB> <hex of the BAM>
 *     bai <TAB> <case name> <TAB> <hex of the BAI>
 *
 * The BAM is emitted too so the Rust side can confirm it is indexing the same bytes: an index
 * that matches for a file that does not is a coincidence, not a result.
 */

import htsjdk.samtools.*;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

public class BaiDump {

    static String hex(final byte[] bytes) {
        final StringBuilder sb = new StringBuilder();
        for (final byte b : bytes) sb.append(String.format("%02x", b));
        return sb.toString();
    }

    static void emit(final String name, final SAMFileHeader header, final List<SAMRecord> records)
            throws Exception {
        // setCreateIndex needs a real path, so write to a temp file and read both back.
        final File bam = File.createTempFile("bai-dump-", ".bam");
        bam.deleteOnExit();
        final SAMFileWriter w = new SAMFileWriterFactory()
                .setCreateIndex(true)
                .setCreateMd5File(false)
                .setUseAsyncIo(false)
                .makeBAMWriter(header, true, bam);
        for (final SAMRecord r : records) w.addAlignment(r);
        w.close();

        final File bai = new File(bam.getPath().replaceAll("\\.bam$", ".bai"));
        bai.deleteOnExit();
        System.out.println("bam\t" + name + "\t" + hex(Files.readAllBytes(bam.toPath())));
        System.out.println("bai\t" + name + "\t" + hex(Files.readAllBytes(bai.toPath())));
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
        r.setAlignmentStart(start);
        r.setMappingQuality(60);
        r.setCigarString(cigar);
        final int len = TextCigarCodec.decode(cigar).getReadLength();
        final byte[] bases = new byte[len];
        java.util.Arrays.fill(bases, (byte) 'A');
        r.setReadBases(bases);
        final byte[] quals = new byte[len];
        java.util.Arrays.fill(quals, (byte) 30);
        r.setBaseQualities(quals);
        return r;
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // A reference with no reads at all: exercises writeNullContent.
        emit("empty", header(250000000), new ArrayList<>());

        // One read: the smallest index that has a real bin and the pseudo-bin.
        {
            final SAMFileHeader h = header(250000000);
            emit("one_read", h, List.of(read(h, "r1", 0, 100, "4M")));
        }

        // Two reads in the same 16 kb window and the same block: chunks must coalesce.
        {
            final SAMFileHeader h = header(250000000);
            emit("same_window", h, List.of(
                    read(h, "r1", 0, 100, "50M"),
                    read(h, "r2", 0, 200, "50M")));
        }

        // Reads far apart, leaving empty linear windows that must be back-filled.
        {
            final SAMFileHeader h = header(250000000);
            emit("sparse_windows", h, List.of(
                    read(h, "r1", 0, 1, "50M"),
                    read(h, "r2", 0, 5 * 16384 + 1, "50M"),
                    read(h, "r3", 0, 40 * 16384 + 1, "50M")));
        }

        // Reads at every bin level, so several bins appear with different chunk lists.
        {
            final SAMFileHeader h = header(250000000);
            final List<SAMRecord> rs = new ArrayList<>();
            rs.add(read(h, "lvl5", 0, 1, "100M"));
            rs.add(read(h, "lvl4", 0, 100, "200000M"));
            rs.add(read(h, "lvl3", 0, 300000, "2000000M"));
            rs.add(read(h, "lvl2", 0, 3000000, "20000000M"));
            rs.add(read(h, "lvl1", 0, 30000000, "100000000M"));
            emit("all_levels", h, rs);
        }

        // Several references, one of them with no reads in the middle.
        {
            final SAMFileHeader h = header(250000000, 200000000, 100000000);
            emit("multi_reference_with_gap", h, List.of(
                    read(h, "a", 0, 100, "50M"),
                    read(h, "b", 2, 500, "50M")));
        }

        // Unmapped reads: placed ones are indexed and counted as unaligned; unplaced ones only
        // bump the no-coordinate counter.
        {
            final SAMFileHeader h = header(250000000);
            final SAMRecord placed = read(h, "placed", 0, 100, "50M");
            placed.setReadUnmappedFlag(true);
            final SAMRecord unplaced = read(h, "unplaced", 0, 100, "50M");
            unplaced.setReadUnmappedFlag(true);
            unplaced.setReferenceIndex(-1);
            unplaced.setAlignmentStart(0);
            unplaced.setCigarString("*");
            emit("unmapped", h, List.of(
                    read(h, "mapped", 0, 50, "50M"), placed, unplaced, unplaced));
        }

        // Enough reads to cross BGZF block boundaries, so chunks stop coalescing.
        {
            final SAMFileHeader h = header(250000000);
            final List<SAMRecord> rs = new ArrayList<>();
            for (int i = 0; i < 20000; i++) {
                rs.add(read(h, "r" + i, 0, 1 + i * 11, "50M"));
            }
            emit("many_blocks", h, rs);
        }

        // Reads on a window boundary, where the -1 in the no-end path matters.
        {
            final SAMFileHeader h = header(250000000);
            emit("window_boundary", h, List.of(
                    read(h, "before", 0, 16384, "1M"),
                    read(h, "on", 0, 16385, "1M"),
                    read(h, "after", 0, 16386, "1M")));
        }
    }
}
