/*
 * Oracle dump harness for SAMUtils.calculateReadGroupRecordChecksum conformance in htsjdk-rs.
 *
 * Emits an escaped TSV to stdout: a `sam` row (the header-only SAM the checksum is computed from)
 * and a `checksum` row (the 32-char MD5 hex htsjdk produced). The committed corpus is
 * `java ... RgChecksumDump | gzip > tests/data/rg_checksum.txt.gz`, regenerated and compared in CI.
 *
 * The read groups are chosen to exercise the sort: a record missing the leading PU tag (which sorts
 * first), records that tie on PU and break on a later tag, attributes given out of key order (the
 * digest must visit them in key order), and a value containing a space.
 *
 *   java -cp htsjdk.jar:. RgChecksumDump
 */
import htsjdk.samtools.*;

import java.io.File;
import java.nio.file.Files;

public class RgChecksumDump {
    public static void main(final String[] args) throws Exception {
        final SAMFileHeader header = new SAMFileHeader();
        header.addSequence(new SAMSequenceRecord("chr1", 10000));

        // Attributes are set in a deliberately non-sorted order to prove the digest sorts by key.
        final SAMReadGroupRecord a = new SAMReadGroupRecord("A");
        a.setPlatformUnit("unit1");
        a.setSample("sampleX");
        a.setLibrary("lib2");
        a.setPlatform("ILLUMINA");

        final SAMReadGroupRecord b = new SAMReadGroupRecord("B");
        b.setPlatformUnit("unit2");
        b.setLibrary("lib1");
        b.setSample("sampleX");

        // No PU: a null leading tag sorts this record before the others.
        final SAMReadGroupRecord c = new SAMReadGroupRecord("C");
        c.setSample("sampleY");
        c.setDescription("a description with spaces");

        // Ties A on PU=unit1; breaks later. Also carries a date and a center.
        final SAMReadGroupRecord d = new SAMReadGroupRecord("D");
        d.setPlatformUnit("unit1");
        d.setLibrary("lib3");
        d.setSample("sampleX");
        d.setSequencingCenter("centerZ");

        // Added in an order unrelated to the sort order.
        header.addReadGroup(b);
        header.addReadGroup(d);
        header.addReadGroup(a);
        header.addReadGroup(c);

        final File sam = File.createTempFile("rgck-", ".sam");
        sam.deleteOnExit();
        final SAMFileWriter w = new SAMFileWriterFactory().makeSAMWriter(header, false, sam);
        w.close(); // header only

        final String checksum = SAMUtils.calculateReadGroupRecordChecksum(sam, null);

        emit("sam", "multi", new String(Files.readAllBytes(sam.toPath())));
        emit("checksum", "multi", checksum);
    }

    private static void emit(final String kind, final String kase, final String payload) {
        final String esc = payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
        System.out.println(kind + "\t" + kase + "\t" + esc);
    }
}
