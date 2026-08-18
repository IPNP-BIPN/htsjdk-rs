/*
 * Which @HD version a written BAM carries, taken from the reference.
 *
 * `BAMFileWriter` has two header paths and they disagree. The ordinary writer, the one
 * `SAMFileWriterFactory` builds, goes through `SAMFileWriterImpl.writeHeader`, which encodes with
 * keepExistingVersionNumber = FALSE: the @HD line is rebuilt from a fresh `SAMFileHeader`, so the
 * version becomes the current one whatever the input said. The static
 * `BAMFileWriter.writeHeader(BinaryCodec, SAMFileHeader)` passes TRUE and keeps it, and that method
 * is only reachable from the block-copy reheader path.
 *
 * Three things this is built to catch.
 *
 *   - A WRITTEN BAM CARRIES THE CURRENT VERSION, not the input's. A header handed to the factory
 *     with VN:1.5 produces a file whose header text says VN:1.6;
 *   - THE KEPT PATH KEEPS IT, so the same header encoded with keepExistingVersionNumber = true is
 *     still VN:1.5. Both are printed here, because a port with one function for both would pass
 *     whichever half it implemented;
 *   - AND VN CANNOT BE ANYWHERE BUT FIRST. `SAMFileHeader`'s constructor sets it before anything is
 *     parsed, and `setAttribute` overwrites in place, so even a header parsed from text whose @HD
 *     line reads `SO` before `VN` comes back with VN first. That is worth a row rather than an
 *     assumption: it says the attribute order needs no fixing.
 *
 * Output:
 *
 *     header\t<name>\t<the encoded text, newlines escaped>
 *     file\t<name>\t<the whole BAM, hex>
 *
 * Usage: BamHeaderVersionDump
 */

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.BufferedLineReader;
import htsjdk.samtools.util.zip.DeflaterFactory;

import java.io.ByteArrayOutputStream;
import java.io.StringWriter;

public class BamHeaderVersionDump {

    /** The keepExistingVersionNumber = true path, which is the block-copy reheader's. */
    static void emitHeader(final String name, final SAMFileHeader header) {
        final StringWriter sw = new StringWriter();
        new SAMTextHeaderCodec().encode(sw, header, true);
        System.out.println("header\t" + name + "\t" + sw.toString().replace("\n", "\\n"));
    }

    /** The ordinary writer, which is the path that replaces the version. */
    static void emitFile(final String name, final SAMFileHeader header) {
        final ByteArrayOutputStream out = new ByteArrayOutputStream();
        final SAMFileWriter w = new SAMFileWriterFactory()
                .setCreateIndex(false)
                .setCreateMd5File(false)
                .setUseAsyncIo(false)
                .makeBAMWriter(header, true, out);
        w.close();
        final StringBuilder sb = new StringBuilder();
        for (final byte b : out.toByteArray()) sb.append(String.format("%02x", b));
        System.out.println("file\t" + name + "\t" + sb);
    }

    static SAMFileHeader minimal() {
        final SAMFileHeader h = new SAMFileHeader();
        final SAMSequenceDictionary d = new SAMSequenceDictionary();
        d.addSequence(new SAMSequenceRecord("chr1", 250000000));
        h.setSequenceDictionary(d);
        return h;
    }

    /** A header parsed from text, which is the only way to try putting VN anywhere but first. */
    static SAMFileHeader decode(final String text) {
        return new SAMTextHeaderCodec().decode(BufferedLineReader.fromString(text),
                "BamHeaderVersionDump");
    }

    public static void main(final String[] args) {
        // The oracle contract pins the JDK deflater, so BGZF blocks come from java.util.zip.
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // The header as constructed: already the current version, so the two paths agree.
        final SAMFileHeader current = minimal();
        emitHeader("current_kept", current);
        emitFile("current_written", current);

        // A header from an older file, which is where the two paths part.
        final SAMFileHeader old = minimal();
        old.setAttribute("VN", "1.5");
        emitHeader("old_kept", old);
        emitFile("old_written", old);

        // The same version, parsed from text whose @HD line puts SO before VN.
        final SAMFileHeader parsed = decode(
                "@HD\tSO:coordinate\tVN:1.5\n@SQ\tSN:chr1\tLN:250000000\n");
        emitHeader("parsed_kept", parsed);
        emitFile("parsed_written", parsed);

        // And a sorted header, so the replacement is seen beside another attribute rather than
        // alone on the line.
        final SAMFileHeader sorted = minimal();
        sorted.setAttribute("VN", "1.4");
        sorted.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        emitHeader("sorted_kept", sorted);
        emitFile("sorted_written", sorted);
    }
}
