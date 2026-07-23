/*
 * Dumps htsjdk's block-copy gather (BamFileIoUtils.gatherWithBlockCopying) alongside its inputs.
 *
 * Output, one line per case:
 *
 *     input     <TAB> <case name> <TAB> <hex of one input BAM>   (repeated, in order)
 *     gathered  <TAB> <case name> <TAB> <hex of the gathered BAM>
 *
 * The inputs are emitted so the Rust side gathers exactly the same bytes htsjdk did. The JDK
 * deflater is pinned per the oracle contract, so the BGZF bytes match htsjdk-rs's writer. All inputs
 * in a case share one header, as GatherBamFiles' block-copy fast path requires.
 */

import htsjdk.samtools.*;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;
import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

public class GatherDump {

    static String hex(final byte[] bytes) {
        final StringBuilder sb = new StringBuilder();
        for (final byte b : bytes) sb.append(String.format("%02x", b));
        return sb.toString();
    }

    static SAMFileHeader header() {
        final SAMFileHeader h = new SAMFileHeader();
        final SAMSequenceDictionary d = new SAMSequenceDictionary();
        d.addSequence(new SAMSequenceRecord("chr1", 1000000));
        h.setSequenceDictionary(d);
        h.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        return h;
    }

    static SAMRecord read(final SAMFileHeader h, final String name, final int start) {
        final SAMRecord r = new SAMRecord(h);
        r.setReadName(name);
        r.setReferenceIndex(0);
        r.setAlignmentStart(start);
        r.setMappingQuality(60);
        r.setCigarString("50M");
        final byte[] bases = new byte[50];
        java.util.Arrays.fill(bases, (byte) 'A');
        r.setReadBases(bases);
        final byte[] quals = new byte[50];
        java.util.Arrays.fill(quals, (byte) 30);
        r.setBaseQualities(quals);
        return r;
    }

    /** Writes one BAM file from a list of read specs (name,start) and returns it. */
    static File writeBam(final List<int[]> starts, final int base) throws Exception {
        final SAMFileHeader h = header();
        final File f = File.createTempFile("gather-in-", ".bam");
        f.deleteOnExit();
        final SAMFileWriter w = new SAMFileWriterFactory()
                .setCreateIndex(false).setCreateMd5File(false).setUseAsyncIo(false)
                .makeBAMWriter(h, true, f);
        for (final int[] s : starts) w.addAlignment(read(h, "r" + base + "_" + s[0], s[1]));
        w.close();
        return f;
    }

    static void emit(final String name, final List<File> inputs) throws Exception {
        final File out = File.createTempFile("gather-out-", ".bam");
        out.deleteOnExit();
        BamFileIoUtils.gatherWithBlockCopying(inputs, out, false, false);
        for (final File in : inputs) {
            System.out.println("input\t" + name + "\t" + hex(Files.readAllBytes(in.toPath())));
        }
        System.out.println("gathered\t" + name + "\t" + hex(Files.readAllBytes(out.toPath())));
    }

    static List<int[]> reads(final int... starts) {
        final List<int[]> out = new ArrayList<>();
        for (int i = 0; i < starts.length; i++) out.add(new int[]{i, starts[i]});
        return out;
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // Two files, a couple of reads each.
        emit("two_files", List.of(
                writeBam(reads(10, 20), 0),
                writeBam(reads(30, 40), 1)));

        // Three files.
        emit("three_files", List.of(
                writeBam(reads(10), 0),
                writeBam(reads(20, 25), 1),
                writeBam(reads(30), 2)));

        // First file empty (header only): its header leads, the rest follow.
        emit("first_empty", List.of(
                writeBam(reads(), 0),
                writeBam(reads(10, 20), 1)));

        // Middle file empty.
        emit("middle_empty", List.of(
                writeBam(reads(10), 0),
                writeBam(reads(), 1),
                writeBam(reads(20), 2)));

        // Many reads per file, so records span several BGZF blocks.
        {
            final List<int[]> a = new ArrayList<>();
            for (int i = 0; i < 3000; i++) a.add(new int[]{i, 1 + i * 100});
            final List<int[]> b = new ArrayList<>();
            for (int i = 0; i < 3000; i++) b.add(new int[]{i, 2 + i * 100});
            emit("many_blocks", List.of(writeBam(a, 0), writeBam(b, 1)));
        }
    }
}
