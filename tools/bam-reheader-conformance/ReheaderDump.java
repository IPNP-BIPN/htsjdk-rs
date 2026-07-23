/*
 * Dumps htsjdk's block-copy reheader (BamFileIoUtils.reheaderBamFile) alongside the input BAM.
 *
 * Output, one line per case:
 *
 *     input      <TAB> <case name> <TAB> <hex of the input BAM>
 *     comment    <TAB> <case name> <TAB> <one comment>   (repeated, in order)
 *     reheadered <TAB> <case name> <TAB> <hex of the reheadered BAM>
 *
 * The input BAM is emitted so the Rust side reheaders exactly the same bytes htsjdk did. The JDK
 * deflater is pinned per the oracle contract, so the BGZF bytes match htsjdk-rs's writer.
 */

import htsjdk.samtools.*;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;
import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

public class ReheaderDump {

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

    static void emit(final String name, final SAMFileHeader header, final List<SAMRecord> records,
                     final List<String> comments) throws Exception {
        final File in = File.createTempFile("reheader-in-", ".bam");
        in.deleteOnExit();
        final SAMFileWriter w = new SAMFileWriterFactory()
                .setCreateIndex(false)
                .setCreateMd5File(false)
                .setUseAsyncIo(false)
                .makeBAMWriter(header, true, in);
        for (final SAMRecord r : records) w.addAlignment(r);
        w.close();

        // Read the header back the way AddCommentsToBam does, add the comments, reheader.
        final SAMFileHeader newHeader = SamReaderFactory.makeDefault().getFileHeader(in);
        for (final String c : comments) newHeader.addComment(c);

        final File out = File.createTempFile("reheader-out-", ".bam");
        out.deleteOnExit();
        BamFileIoUtils.reheaderBamFile(newHeader, in.toPath(), out.toPath(), false, false);

        System.out.println("input\t" + name + "\t" + hex(Files.readAllBytes(in.toPath())));
        for (final String c : comments) System.out.println("comment\t" + name + "\t" + c);
        System.out.println("reheadered\t" + name + "\t" + hex(Files.readAllBytes(out.toPath())));
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // One comment on a small coordinate-sorted file.
        {
            final SAMFileHeader h = header(1000);
            emit("one_comment", h, List.of(
                    read(h, "a", 0, 10, "4M"),
                    read(h, "b", 0, 20, "4M")), List.of("a comment"));
        }

        // Several comments at once.
        {
            final SAMFileHeader h = header(1000);
            emit("three_comments", h, List.of(read(h, "a", 0, 10, "4M")),
                    List.of("first", "second with spaces", "third"));
        }

        // Multiple references.
        {
            final SAMFileHeader h = header(1000, 2000, 500);
            emit("multi_ref", h, List.of(
                    read(h, "a", 0, 10, "10M"),
                    read(h, "b", 2, 100, "10M")), List.of("multiref comment"));
        }

        // No records at all: the first-record offset lands at the end of the header, so the
        // re-compressed tail is empty and the copy is just the terminator.
        {
            final SAMFileHeader h = header(1000);
            emit("no_records", h, new ArrayList<>(), List.of("comment on empty"));
        }

        // The input already carries a comment; a second is appended.
        {
            final SAMFileHeader h = header(1000);
            h.addComment("pre-existing");
            emit("existing_comment", h, List.of(read(h, "a", 0, 10, "4M")),
                    List.of("appended"));
        }

        // Many reads, so records span several BGZF blocks and the raw copy covers real blocks.
        {
            final SAMFileHeader h = header(1000000);
            final List<SAMRecord> rs = new ArrayList<>();
            for (int i = 0; i < 4000; i++) rs.add(read(h, "r" + i, 0, 1 + i * 100, "50M"));
            emit("many_blocks", h, rs, List.of("spanning comment"));
        }
    }
}
