/*
 * Dumps whole BAM files produced by htsjdk, and the SAM text header alone, for a range of
 * header shapes.
 *
 * Output, one line per case:
 *
 *     header <TAB> <case name> <TAB> <SAM text header, newlines as \n>
 *     file   <TAB> <case name> <TAB> <hex of the complete BAM file>
 *
 * The BAM writer pins the JDK deflater, matching the oracle contract, so the compressed bytes
 * are the ones htsjdk-rs must reproduce.
 */

import htsjdk.samtools.*;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.BinaryCodec;
import htsjdk.samtools.util.zip.DeflaterFactory;
import java.io.ByteArrayOutputStream;
import java.io.StringWriter;
import java.util.Arrays;

public class BamFileDump {

    static void emitHeader(final String name, final SAMFileHeader header) {
        final StringWriter sw = new StringWriter();
        new SAMTextHeaderCodec().encode(sw, header, true);
        System.out.println("header\t" + name + "\t" + sw.toString().replace("\n", "\\n"));
    }

    static void emitFile(final String name, final SAMFileHeader header, final SAMRecord... records) {
        final ByteArrayOutputStream out = new ByteArrayOutputStream();
        final SAMFileWriter w = new SAMFileWriterFactory()
                .setCreateIndex(false)
                .setCreateMd5File(false)
                .setUseAsyncIo(false)
                .makeBAMWriter(header, true, out);
        for (final SAMRecord r : records) w.addAlignment(r);
        w.close();
        final byte[] bytes = out.toByteArray();
        final StringBuilder sb = new StringBuilder();
        for (final byte b : bytes) sb.append(String.format("%02x", b));
        System.out.println("file\t" + name + "\t" + sb);
    }

    static SAMFileHeader minimal() {
        final SAMFileHeader h = new SAMFileHeader();
        final SAMSequenceDictionary d = new SAMSequenceDictionary();
        d.addSequence(new SAMSequenceRecord("chr1", 250000000));
        h.setSequenceDictionary(d);
        return h;
    }

    public static void main(final String[] args) {
        // The oracle contract pins the JDK deflater, so BGZF blocks come from java.util.zip.
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // ---- header shapes ---------------------------------------------------------------
        emitHeader("empty_dict", new SAMFileHeader());
        emitHeader("minimal", minimal());

        final SAMFileHeader sorted = minimal();
        sorted.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        emitHeader("sorted", sorted);

        final SAMFileHeader queryname = minimal();
        queryname.setSortOrder(SAMFileHeader.SortOrder.queryname);
        queryname.setGroupOrder(SAMFileHeader.GroupOrder.query);
        emitHeader("queryname_grouped", queryname);

        // Several sequences, with the optional SQ attributes.
        final SAMFileHeader manySq = new SAMFileHeader();
        final SAMSequenceDictionary d = new SAMSequenceDictionary();
        final SAMSequenceRecord s1 = new SAMSequenceRecord("chr1", 250000000);
        s1.setAssembly("GRCh38");
        s1.setMd5("d41d8cd98f00b204e9800998ecf8427e");
        s1.setSpecies("Homo sapiens");
        s1.setAttribute("UR", "file:/ref/chr1.fa");
        d.addSequence(s1);
        d.addSequence(new SAMSequenceRecord("chr2", 200000000));
        final SAMSequenceRecord s3 = new SAMSequenceRecord("chrM", 16571);
        s3.setAttribute("M5", "c68f52674c9fb33aef52dcf399755519");
        d.addSequence(s3);
        manySq.setSequenceDictionary(d);
        emitHeader("many_sq", manySq);

        // Read groups.
        final SAMFileHeader withRg = minimal();
        final SAMReadGroupRecord rg1 = new SAMReadGroupRecord("rg1");
        rg1.setSample("sample1");
        rg1.setLibrary("lib1");
        rg1.setPlatform("ILLUMINA");
        rg1.setPlatformUnit("unit1");
        withRg.addReadGroup(rg1);
        final SAMReadGroupRecord rg2 = new SAMReadGroupRecord("rg2");
        rg2.setSample("sample2");
        withRg.addReadGroup(rg2);
        emitHeader("read_groups", withRg);

        // Attribute insertion order: set, overwrite, then set another. LinkedHashMap keeps the
        // overwritten key in its ORIGINAL position, which a naive remove-then-append loses.
        final SAMFileHeader reorder = minimal();
        final SAMReadGroupRecord rg = new SAMReadGroupRecord("rgx");
        rg.setAttribute("SM", "first");
        rg.setAttribute("LB", "lib");
        rg.setAttribute("SM", "second");   // overwrite: must stay before LB
        rg.setAttribute("PL", "ILLUMINA");
        reorder.addReadGroup(rg);
        emitHeader("attribute_overwrite_keeps_position", reorder);

        // Program records and comments.
        final SAMFileHeader withPg = minimal();
        final SAMProgramRecord pg1 = new SAMProgramRecord("prog1");
        pg1.setProgramName("MarkDuplicates");
        pg1.setProgramVersion("3.4.0");
        pg1.setCommandLine("MarkDuplicates I=in.bam O=out.bam");
        withPg.addProgramRecord(pg1);
        final SAMProgramRecord pg2 = new SAMProgramRecord("prog2");
        pg2.setPreviousProgramGroupId("prog1");
        withPg.addProgramRecord(pg2);
        withPg.addComment("a comment");
        withPg.addComment("another comment");
        emitHeader("programs_and_comments", withPg);

        // Everything at once, in the order htsjdk emits: HD, SQ, RG, PG, CO.
        final SAMFileHeader full = manySq;
        full.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        full.addReadGroup(rg1);
        full.addProgramRecord(pg1);
        full.addComment("full header");
        emitHeader("full", full);

        // ---- whole files ------------------------------------------------------------------
        emitFile("file_empty", sorted);

        final SAMFileHeader fh = minimal();
        fh.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMRecord r1 = new SAMRecord(fh);
        r1.setReadName("read1");
        r1.setReferenceIndex(0);
        r1.setAlignmentStart(100);
        r1.setMappingQuality(60);
        r1.setCigarString("4M");
        r1.setReadBases("ACGT".getBytes());
        r1.setBaseQualities(new byte[]{30, 30, 30, 30});
        r1.setAttribute("NM", 0);
        emitFile("file_one_record", fh, r1);

        final SAMRecord[] many = new SAMRecord[500];
        for (int i = 0; i < many.length; i++) {
            final SAMRecord r = new SAMRecord(fh);
            r.setReadName("read" + i);
            r.setReferenceIndex(0);
            r.setAlignmentStart(100 + i * 37);
            r.setMappingQuality(60);
            r.setCigarString("10M");
            final byte[] bases = new byte[10];
            for (int j = 0; j < 10; j++) bases[j] = (byte) "ACGT".charAt((i + j) % 4);
            r.setReadBases(bases);
            final byte[] quals = new byte[10];
            Arrays.fill(quals, (byte) (20 + i % 20));
            r.setBaseQualities(quals);
            r.setAttribute("NM", i % 5);
            r.setAttribute("RG", "rg1");
            many[i] = r;
        }
        emitFile("file_500_records", fh, many);

        // Enough records to spill past a single BGZF block, so block boundaries are exercised.
        final SAMRecord[] lots = new SAMRecord[20000];
        for (int i = 0; i < lots.length; i++) {
            final SAMRecord r = new SAMRecord(fh);
            r.setReadName("r" + i);
            r.setReferenceIndex(0);
            r.setAlignmentStart(1 + i * 11);
            r.setMappingQuality(60);
            r.setCigarString("50M");
            final byte[] bases = new byte[50];
            for (int j = 0; j < 50; j++) bases[j] = (byte) "ACGTN".charAt((i * 7 + j) % 5);
            r.setReadBases(bases);
            final byte[] quals = new byte[50];
            for (int j = 0; j < 50; j++) quals[j] = (byte) ((i + j) % 60);
            r.setBaseQualities(quals);
            lots[i] = r;
        }
        emitFile("file_20000_records", fh, lots);
    }
}
