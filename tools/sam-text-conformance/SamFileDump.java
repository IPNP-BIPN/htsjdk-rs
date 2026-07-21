/*
 * Writes complete SAM text files with htsjdk and emits them, alongside the BAM built from the
 * same records.
 *
 * Output:
 *   sam <TAB> <case> <TAB> <the whole SAM file, \n and \t escaped>
 *   bam <TAB> <case> <TAB> <hex of the equivalent BAM>
 *
 * Both are emitted so the Rust side can check that reading the BAM and writing the SAM gives
 * htsjdk's SAM, which is the conversion a user actually performs.
 */

import htsjdk.samtools.*;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;
import java.io.ByteArrayOutputStream;
import java.nio.file.Files;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

public class SamFileDump {

    static String hex(byte[] b) {
        StringBuilder sb = new StringBuilder();
        for (byte x : b) sb.append(String.format("%02x", x));
        return sb.toString();
    }

    static String esc(String s) {
        return s.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
    }

    static void emit(String name, SAMFileHeader h, List<SAMRecord> records) throws Exception {
        File sam = File.createTempFile("sf-", ".sam");
        sam.deleteOnExit();
        SAMFileWriter w = new SAMFileWriterFactory()
                .setCreateIndex(false).setCreateMd5File(false).setUseAsyncIo(false)
                .makeSAMWriter(h, false, sam);
        for (SAMRecord r : records) w.addAlignment(r);
        w.close();
        System.out.println("sam\t" + name + "\t"
                + esc(new String(Files.readAllBytes(sam.toPath()))));

        ByteArrayOutputStream bam = new ByteArrayOutputStream();
        SAMFileWriter bw = new SAMFileWriterFactory()
                .setCreateIndex(false).setCreateMd5File(false).setUseAsyncIo(false)
                .makeBAMWriter(h, false, bam);
        for (SAMRecord r : records) bw.addAlignment(r);
        bw.close();
        System.out.println("bam\t" + name + "\t" + hex(bam.toByteArray()));
    }

    static SAMFileHeader header(String... names) {
        SAMFileHeader h = new SAMFileHeader();
        SAMSequenceDictionary d = new SAMSequenceDictionary();
        int len = 250000000;
        for (String n : names) d.addSequence(new SAMSequenceRecord(n, len -= 10000000));
        h.setSequenceDictionary(d);
        return h;
    }

    static SAMRecord read(SAMFileHeader h, String name, int ref, int start, String cigar,
                          String bases, int flags) {
        SAMRecord r = new SAMRecord(h);
        r.setReadName(name);
        r.setReferenceIndex(ref);
        r.setAlignmentStart(start);
        r.setMappingQuality(60);
        r.setCigarString(cigar);
        r.setReadBases(bases.getBytes());
        byte[] q = new byte[bases.length()];
        for (int i = 0; i < q.length; i++) q[i] = (byte) (20 + (i % 25));
        r.setBaseQualities(q);
        r.setFlags(flags);
        r.setMateReferenceIndex(-1);
        r.setMateAlignmentStart(0);
        return r;
    }

    public static void main(String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        SAMFileHeader h1 = header("chr1");
        h1.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        emit("header_only", h1, new ArrayList<>());

        SAMFileHeader h2 = header("chr1", "chr2");
        h2.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        List<SAMRecord> few = new ArrayList<>();
        few.add(read(h2, "r1", 0, 100, "4M", "ACGT", 0));
        few.add(read(h2, "r2", 1, 200, "10M", "ACGTACGTAC", 0));
        emit("two_references", h2, few);

        // Read groups, programs and comments, so the header text has every section.
        SAMFileHeader h3 = header("chr1");
        h3.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        SAMReadGroupRecord rg = new SAMReadGroupRecord("rg1");
        rg.setSample("s1");
        rg.setLibrary("lib1");
        h3.addReadGroup(rg);
        SAMProgramRecord pg = new SAMProgramRecord("p1");
        pg.setProgramName("test");
        h3.addProgramRecord(pg);
        h3.addComment("a comment");
        List<SAMRecord> tagged = new ArrayList<>();
        for (int i = 0; i < 30; i++) {
            SAMRecord r = read(h3, "t" + i, 0, 100 + i * 11, "6M", "ACGTNA", 0);
            r.setAttribute("NM", i % 5);
            r.setAttribute("RG", "rg1");
            r.setAttribute("XF", 0.5f * i);
            r.setAttribute("XB", new int[]{i, i + 1});
            tagged.add(r);
        }
        emit("full_header_and_tags", h3, tagged);

        // Unmapped and unplaced records, where the mandatory-field rules bite.
        SAMFileHeader h4 = header("chr1");
        List<SAMRecord> mixed = new ArrayList<>();
        mixed.add(read(h4, "mapped", 0, 100, "4M", "ACGT", 0));
        SAMRecord placedUnmapped = read(h4, "placed", 0, 100, "4M", "ACGT", 0x4);
        mixed.add(placedUnmapped);
        SAMRecord unplaced = read(h4, "unplaced", -1, 0, "*", "ACGT", 0x4);
        mixed.add(unplaced);
        emit("unmapped", h4, mixed);

        // Enough records that the file is not trivially short.
        SAMFileHeader h5 = header("chr1");
        h5.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        List<SAMRecord> many = new ArrayList<>();
        for (int i = 0; i < 500; i++) {
            many.add(read(h5, "m" + i, 0, 1 + i * 13, "20M", "ACGTNACGTNACGTNACGTN", 0));
        }
        emit("many_records", h5, many);
    }
}
