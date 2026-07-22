/*
 * Oracle dump harness for SamPairUtil.setMateInfo conformance in htsjdk-rs.
 *
 * For each case it builds a pair with deliberately-absent mate info, writes it as a 2-record SAM
 * (the `input` row), runs SamPairUtil.setMateInfo(rec1, rec2, true), and writes the result (the
 * `output` row). Every mate field lands in the SAM columns (flags, RNEXT/PNEXT/TLEN) and the MC/MQ
 * tags, so comparing the output SAM validates all of setMateInfo at once. The Rust port reads the
 * input SAM, applies set_mate_info to the two records, and must reproduce the output SAM.
 *
 *   java -cp htsjdk.jar:. SetMateInfoDump
 */
import htsjdk.samtools.*;

import java.io.File;
import java.nio.file.Files;

public class SetMateInfoDump {
    static SAMFileHeader header() {
        final SAMFileHeader h = new SAMFileHeader();
        h.addSequence(new SAMSequenceRecord("chr1", 1000));
        h.addSequence(new SAMSequenceRecord("chr2", 1000));
        h.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        return h;
    }

    static SAMRecord read(final SAMFileHeader h, final String name, final boolean first,
                          final boolean unmapped, final String rname, final int pos,
                          final boolean negative, final int mapq, final String cigar) {
        final SAMRecord r = new SAMRecord(h);
        r.setReadName(name);
        r.setReadPairedFlag(true);
        r.setFirstOfPairFlag(first);
        r.setSecondOfPairFlag(!first);
        r.setReadBases("ACGTACGTAC".getBytes());
        r.setBaseQualities("IIIIIIIIII".getBytes());
        if (unmapped) {
            r.setReadUnmappedFlag(true);
        } else {
            r.setReferenceName(rname);
            r.setAlignmentStart(pos);
            r.setReadNegativeStrandFlag(negative);
            r.setMappingQuality(mapq);
            r.setCigarString(cigar);
        }
        return r;
    }

    public static void main(final String[] args) throws Exception {
        emit("bothMapped",   read(header(), "p", true,  false, "chr1", 100, false, 60, "10M"),
                             read(header(), "p", false, false, "chr1", 200, true,  50, "10M"));
        emit("samePos",      read(header(), "p", true,  false, "chr1", 100, false, 60, "10M"),
                             read(header(), "p", false, false, "chr1", 100, true,  60, "10M"));
        emit("crossContig",  read(header(), "p", true,  false, "chr1", 100, false, 60, "10M"),
                             read(header(), "p", false, false, "chr2", 300, false, 40, "10M"));
        emit("oneUnmapped",  read(header(), "p", true,  false, "chr1", 500, true,  60, "10M"),
                             read(header(), "p", false, true,  null,   0,   false, 0,  null));
        emit("bothUnmapped", read(header(), "p", true,  true,  null, 0, false, 0, null),
                             read(header(), "p", false, true,  null, 0, true,  0, null));
    }

    static void emit(final String kase, final SAMRecord r1, final SAMRecord r2) throws Exception {
        final SAMFileHeader h = r1.getHeader();
        final File in = File.createTempFile("mi-in-", ".sam"); in.deleteOnExit();
        writeSam(h, in, r1, r2);

        SamPairUtil.setMateInfo(r1, r2, true);

        final File out = File.createTempFile("mi-out-", ".sam"); out.deleteOnExit();
        writeSam(h, out, r1, r2);

        emitRow("input", kase, new String(Files.readAllBytes(in.toPath())));
        emitRow("output", kase, new String(Files.readAllBytes(out.toPath())));
    }

    static void writeSam(final SAMFileHeader h, final File f, final SAMRecord r1, final SAMRecord r2) {
        final SAMFileWriter w = new SAMFileWriterFactory().makeSAMWriter(h, false, f);
        w.addAlignment(r1);
        w.addAlignment(r2);
        w.close();
    }

    static void emitRow(final String kind, final String kase, final String payload) {
        final String esc = payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
        System.out.println(kind + "\t" + kase + "\t" + esc);
    }
}
