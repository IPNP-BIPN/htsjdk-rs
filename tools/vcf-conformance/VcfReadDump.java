/*
 * Reading a whole VCF file: the header frame, then every line after it, through one codec.
 *
 * The header frame (VcfHeaderParseDump), one data line (VcfRecordParseDump) and one genotype block
 * (VcfGenotypeParseDump) each have a suite already. What none of them measures is the loop that
 * joins them, and the loop is where the surprises are, because the codec is stateful and the state
 * it carries between lines is not visible in any single line.
 *
 *   - the line counter is shared and it is incremented in two places, so the SAME malformed line
 *     reports two different line numbers depending on which check refuses it. decodeLine's column
 *     count check runs before parseVCFLine's lineNo++, and generateException runs after it;
 *   - a line beginning with '#' anywhere in the body decodes to null rather than to a record or to
 *     a refusal, so a reader that pushes every decode result onto a list gets a NullPointerException
 *     at some later and unrelated point, and a reader that keeps them silently loses records;
 *   - the version is read once, from the header, and it selects a text transformer that is then
 *     applied to every INFO and every genotype value for the rest of the file. So the same data
 *     line means different things in a v4.2 file and a v4.3 one, and nothing on the line says so;
 *   - that transformer's decoder is Integer.parseInt(s, 16), which accepts a sign. '%+1' is a
 *     control character and '%-1' is left alone, and neither is percent-encoding;
 *   - the header handed back is not the header the file contains. doOnTheFlyModifications defaults
 *     to true, so every INFO and FORMAT line whose ID is one htsjdk holds a standard for is
 *     REWRITTEN when the two disagree on count or type, and the rebuilt header keeps its version
 *     only when that version is v4.3 or later. So a v4.2 file read back has forgotten which
 *     version it is, and the codec is the only thing that still knows.
 *
 * The round trip is measured rather than assumed. Reading a file and writing it back is not the
 * identity and is not close to it: the writer substitutes its own fileformat line and sorts the
 * metadata, so the bytes differ even for a file htsjdk itself produced.
 *
 * Output:
 *
 *     file\t<label>\t<codec version>\t<header version>\t<samples>\t<record count>
 *     rec\t<label>\t<index>\t<chr>\t<start>\t<stop>\t<id>\t<alleles>\t<qual>\t<filters>\t<attrs>\t<genotypes>
 *     null\t<label>\t<index>          a line that decoded to null, with the index it would have had
 *     hdr\t<label>\t<the ID'd header lines as the codec hands them back, in input order>
 *     err\t<label>\t<class>\t<message>
 *     trip\t<label>\t<same|differs>\t<first differing offset or ->\t<rewritten file, escaped>
 *     pct\t<label>\t<raw>\t<decoded under 4.3>\t<decoded under 4.2>
 *
 * Usage: VcfReadDump
 */

import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.SynchronousLineReader;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;

import java.io.ByteArrayOutputStream;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.TreeMap;
import java.util.TreeSet;

public class VcfReadDump {

    /** The metadata every case shares, so a difference in a row is a difference in the body. */
    static final String META =
              "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
            + "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Frequency\">\n"
            + "##INFO=<ID=NOTE,Number=1,Type=String,Description=\"A note\">\n"
            + "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"In dbSNP\">\n"
            + "##FILTER=<ID=LowQual,Description=\"Low quality\">\n"
            + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            + "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n"
            + "##FORMAT=<ID=SB,Number=1,Type=String,Description=\"A string\">\n"
            + "##contig=<ID=chr1,length=100000>\n"
            + "##contig=<ID=chr2,length=100000>\n";

    static final String COLUMNS_SITES = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    static final String COLUMNS_SAMPLES =
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA1\tNA2\n";

    static String sites(final String version, final String body) {
        return "##fileformat=" + version + "\n" + META + COLUMNS_SITES + body;
    }

    static String samples(final String version, final String body) {
        return "##fileformat=" + version + "\n" + META + COLUMNS_SAMPLES + body;
    }

    public static void main(final String[] args) {
        System.out.println("# VcfReadDump: reading a whole VCF file through one stateful codec");

        // The ordinary shapes, so the loop itself is under test before its edges are.
        read("sites-only", sites("VCFv4.2",
                "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\n"
              + "chr1\t200\trs1\tC\tG\t.\t.\tDP=20;AF=0.5\n"
              + "chr2\t1\t.\tGG\tG\t30\tLowQual\tDB\n"));

        read("genotyped", samples("VCFv4.2",
                "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\tGT:GQ\t0/1:30\t1|1:40\n"
              + "chr1\t200\t.\tC\tG\t50\tPASS\tDP=20\tGT\t./.\t0/0\n"));

        read("header-only", sites("VCFv4.2", ""));

        // A body line that decodes to null rather than to a record or to a refusal. Nothing on the
        // way in says a record was dropped.
        read("comment-in-body", sites("VCFv4.2",
                "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\n"
              + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
              + "chr1\t200\t.\tC\tG\t50\tPASS\tDP=20\n"));
        read("hash-comment-in-body", sites("VCFv4.2",
                "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\n"
              + "# a comment\n"
              + "chr1\t200\t.\tC\tG\t50\tPASS\tDP=20\n"));

        // A blank line, which is neither a comment nor a record.
        read("blank-line-in-body", sites("VCFv4.2",
                "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\n"
              + "\n"
              + "chr1\t200\t.\tC\tG\t50\tPASS\tDP=20\n"));

        // The two line numbers. Both files fail on their fourth data line, and the two failures
        // report different numbers for it.
        read("short-line-reports-its-number", sites("VCFv4.2",
                "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\n"
              + "chr1\t200\t.\tC\tG\t50\tPASS\tDP=20\n"
              + "chr1\t300\t.\tC\tG\t50\tPASS\tDP=30\n"
              + "chr1\t400\t.\tC\tG\t50\tPASS\n"));
        read("bad-qual-reports-its-number", sites("VCFv4.2",
                "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\n"
              + "chr1\t200\t.\tC\tG\t50\tPASS\tDP=20\n"
              + "chr1\t300\t.\tC\tG\t50\tPASS\tDP=30\n"
              + "chr1\t400\t.\tC\tG\tx\tPASS\tDP=40\n"));
        // The same short line, first rather than fourth, so the offset is visible and not inferred.
        read("short-line-first", sites("VCFv4.2", "chr1\t100\t.\tA\tT\t50\tPASS\n"));

        // Line endings. The reader is asked whether it strips a carriage return or leaves it on the
        // last column, where it would become part of an INFO value and never fail.
        read("crlf", sites("VCFv4.2",
                "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\r\n"
              + "chr1\t200\t.\tC\tG\t50\tPASS\tNOTE=x\r\n"));
        read("no-trailing-newline", sites("VCFv4.2", "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10"));

        // Header refusals reached through the whole-file path rather than through the frame reader.
        read("no-fileformat", META + COLUMNS_SITES + "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\n");
        read("no-chrom-line", "##fileformat=VCFv4.2\n" + META
                + "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\n");
        read("empty-file", "");
        // split("=") with no limit: three fields is not two, so this is not a version line at all
        // and the file is refused for having no version rather than for having a bad one.
        read("fileformat-with-two-equals",
                "##fileformat=VCFv4.2=x\n" + META + COLUMNS_SITES);
        read("unsupported-version", "##fileformat=VCFv3.3\n" + META + COLUMNS_SITES);

        // The version selects a transformer that applies to every later line. Same body, two files.
        final String encodedBody =
                "chr1\t100\t.\tA\tT\t50\tPASS\tNOTE=a%3Ab%3Bc\n"
              + "chr1\t200\t.\tC\tG\t50\tPASS\tNOTE=100%25\n";
        read("percent-under-4-2", sites("VCFv4.2", encodedBody));
        read("percent-under-4-3", sites("VCFv4.3", encodedBody));
        // And to genotype values, which arrive through a different call site.
        final String encodedGenotypes =
                "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\tGT:SB\t0/1:a%3Ab\t1/1:c%3Ad\n";
        read("percent-genotype-4-2", samples("VCFv4.2", encodedGenotypes));
        read("percent-genotype-4-3", samples("VCFv4.3", encodedGenotypes));
        // The decode runs BEFORE the flag test, so under 4.3 a declared Flag written as DB=%30 is
        // dropped exactly as DB=0 would be, and under 4.2 it is kept as the string "%30".
        final String encodedFlag = "chr1\t100\t.\tA\tT\t50\tPASS\tDB=%30\n";
        read("percent-flag-zero-4-2", sites("VCFv4.2", encodedFlag));
        read("percent-flag-zero-4-3", sites("VCFv4.3", encodedFlag));

        // The header lines the reader hands back, which are not always the ones the file declared.
        // DP, GQ and AF are IDs htsjdk holds a standard for; XX is not.
        headerLines("repair-none", "##fileformat=VCFv4.2\n"
                + "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth; some reads may have been filtered\">\n"
                + COLUMNS_SITES);
        headerLines("repair-wrong-type", "##fileformat=VCFv4.2\n"
                + "##INFO=<ID=DP,Number=1,Type=Float,Description=\"Depth\">\n"
                + COLUMNS_SITES);
        headerLines("repair-wrong-count", "##fileformat=VCFv4.2\n"
                + "##INFO=<ID=DP,Number=2,Type=Integer,Description=\"Depth\">\n"
                + COLUMNS_SITES);
        headerLines("repair-wrong-count-type", "##fileformat=VCFv4.2\n"
                + "##INFO=<ID=DP,Number=A,Type=Integer,Description=\"Depth\">\n"
                + COLUMNS_SITES);
        // Only the description differs, and REPAIR_BAD_DESCRIPTIONS is false.
        headerLines("repair-wrong-description-only", "##fileformat=VCFv4.2\n"
                + "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"my own depth\">\n"
                + COLUMNS_SITES);
        headerLines("repair-format-gq", "##fileformat=VCFv4.2\n"
                + "##FORMAT=<ID=GQ,Number=1,Type=String,Description=\"Genotype Quality\">\n"
                + COLUMNS_SITES);
        headerLines("repair-not-a-standard-id", "##fileformat=VCFv4.2\n"
                + "##INFO=<ID=XX,Number=2,Type=Float,Description=\"Mine\">\n"
                + COLUMNS_SITES);
        // The repair rebuilds the header, and the rebuilt one keeps its version only from 4.3.
        headerLines("repair-under-4-3", "##fileformat=VCFv4.3\n"
                + "##INFO=<ID=DP,Number=1,Type=Float,Description=\"Depth\">\n"
                + COLUMNS_SITES);

        // Which fileformat line comes back, for a file that declared each of the four versions and
        // for one that declared a source line as well, so its position among the rest is visible.
        headerLines("fileformat-4-0", "##fileformat=VCFv4.0\n##source=mine\n" + COLUMNS_SITES);
        headerLines("fileformat-4-1", "##fileformat=VCFv4.1\n##source=mine\n" + COLUMNS_SITES);
        headerLines("fileformat-4-2", "##fileformat=VCFv4.2\n##source=mine\n" + COLUMNS_SITES);
        headerLines("fileformat-4-3", "##fileformat=VCFv4.3\n##source=mine\n" + COLUMNS_SITES);
        // A second fileformat line, further down, which is not the version-setting one.
        headerLines("fileformat-twice",
                "##fileformat=VCFv4.2\n##source=mine\n##fileformat=VCFv4.1\n" + COLUMNS_SITES);

        // Duplicate sample names on the column line, which the header deduplicates and the record
        // then has to be split against.
        read("duplicate-samples", "##fileformat=VCFv4.2\n" + META
                + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA1\tNA1\n"
                + "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\tGT\t0/1\t1/1\n");

        // Records that are not in coordinate order, and one contig the header never declared. A
        // reader that validated either would refuse a file htsjdk reads.
        read("unsorted-and-undeclared-contig", sites("VCFv4.2",
                "chr2\t500\t.\tA\tT\t50\tPASS\tDP=10\n"
              + "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\n"
              + "chrX\t9\t.\tA\tT\t50\tPASS\tDP=10\n"));

        // The percent decoder on its own, because its edges are arithmetic rather than format. The
        // sign cases are the ones no reading of the specification produces.
        for (final String raw : new String[] {
                "plain", "%41", "%3D%41", "a%3Ab", "%", "x%", "%4", "x%4", "%4G", "%G4",
                "%+1", "%-1", "%%41", "%09", "%00", "%7e", "%7E", "100%25", "%zz", "% 1" }) {
            percent(raw);
        }

        // The round trip, on a file htsjdk itself would accept. Reading and writing back is not the
        // identity: the writer declares its own version and sorts the metadata.
        trip("trip-sites", sites("VCFv4.2",
                "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\n"
              + "chr1\t200\trs1\tC\tG\t.\t.\tDP=20;AF=0.5\n"));
        trip("trip-genotyped", samples("VCFv4.2",
                "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\tGT:GQ\t0/1:30\t1|1:40\n"));
        trip("trip-4-3", sites("VCFv4.3",
                "chr1\t100\t.\tA\tT\t50\tPASS\tNOTE=a%3Ab\n"));
    }

    /** Read a whole file the way a feature reader does: the frame, then every remaining line. */
    static void read(final String label, final String text) {
        final VCFCodec codec = new VCFCodec();
        final List<VariantContext> records = new ArrayList<>();
        try {
            final LineIterator it =
                    new LineIteratorImpl(new SynchronousLineReader(new StringReader(text)));
            final VCFHeader header = (VCFHeader) codec.readActualHeader(it);
            int index = 0;
            while (it.hasNext()) {
                final VariantContext vc = codec.decode(it.next());
                if (vc == null) {
                    System.out.printf("null\t%s\t%d%n", label, index);
                    continue;
                }
                records.add(vc);
                index++;
            }
            System.out.printf("file\t%s\t%s\t%s\t%s\t%d%n", label,
                    version(codec.getVersion()), version(header.getVCFHeaderVersion()),
                    header.getGenotypeSamples().isEmpty()
                            ? "-" : String.join(",", header.getGenotypeSamples()),
                    records.size());
            for (int i = 0; i < records.size(); i++) {
                emitRecord(label, i, records.get(i));
            }
        } catch (final Throwable t) {
            // The records read before the refusal are reported too: a file that fails halfway is
            // not the same as one that fails at the first line, and a reader that returns nothing
            // in both cases loses that.
            System.out.printf("file\t%s\t%s\t%s\t%s\t%d%n", label, "aborted", "aborted", "-",
                    records.size());
            for (int i = 0; i < records.size(); i++) {
                emitRecord(label, i, records.get(i));
            }
            System.out.printf("err\t%s\t%s\t%s%n", label, t.getClass().getName(),
                    oneLine(t.getMessage()));
        }
    }

    static String version(final htsjdk.variant.vcf.VCFHeaderVersion v) {
        return v == null ? "none" : v.getVersionString();
    }

    /**
     * Every metadata line as the codec hands it back, in input order, which is not what the file
     * said. Two rewrites are visible in this row and neither is in the format:
     *
     *   - the file's own fileformat line is REMOVED from the stored metadata by the VCFHeader
     *     constructor, and getMetaDataInInputOrder puts a SYNTHESIZED one back at the front. The
     *     synthesized one says VCFv4.2 for everything below v4.3, so a v4.0 file comes back
     *     claiming to be v4.2;
     *   - a standard INFO or FORMAT ID whose count or type disagrees with htsjdk's own is replaced
     *     wholesale, description included.
     */
    static void headerLines(final String label, final String text) {
        try {
            final VCFCodec codec = new VCFCodec();
            final LineIterator it =
                    new LineIteratorImpl(new SynchronousLineReader(new StringReader(text)));
            final VCFHeader header = (VCFHeader) codec.readActualHeader(it);
            final StringJoiner lines = new StringJoiner(" | ");
            header.getMetaDataInInputOrder().forEach(line -> lines.add(line.toString()));
            System.out.printf("hdr\t%s\t%s\t%s\t%s%n", label,
                    version(codec.getVersion()), version(header.getVCFHeaderVersion()), lines);
        } catch (final Throwable t) {
            System.out.printf("err\t%s\t%s\t%s%n", label, t.getClass().getName(),
                    oneLine(t.getMessage()));
        }
    }

    static void emitRecord(final String label, final int index, final VariantContext vc) {
        final StringJoiner alleles = new StringJoiner(",");
        for (final Allele allele : vc.getAlleles()) {
            alleles.add(allele.getDisplayString() + (allele.isReference() ? "*" : ""));
        }

        // Sorted: the attributes live in a HashMap upstream and its iteration order is not a
        // property of the file.
        final StringJoiner attributes = new StringJoiner(";");
        for (final Map.Entry<String, Object> entry : new TreeMap<>(vc.getAttributes()).entrySet()) {
            attributes.add(entry.getKey() + "=" + render(entry.getValue()));
        }

        final String filters;
        if (!vc.filtersWereApplied()) {
            filters = "unfiltered";
        } else if (vc.getFilters().isEmpty()) {
            filters = "PASS";
        } else {
            filters = String.join(",", new TreeSet<>(vc.getFilters()));
        }

        final StringJoiner genotypes = new StringJoiner(" ");
        for (final Genotype g : vc.getGenotypes()) {
            final StringJoiner ga = new StringJoiner(",");
            for (final Map.Entry<String, Object> entry
                    : new TreeMap<>(g.getExtendedAttributes()).entrySet()) {
                ga.add(entry.getKey() + "=" + render(entry.getValue()));
            }
            final StringJoiner called = new StringJoiner(g.isPhased() ? "|" : "/");
            for (final Allele allele : g.getAlleles()) {
                called.add(allele.getDisplayString());
            }
            genotypes.add(g.getSampleName() + ":" + called
                    + ":GQ=" + g.getGQ() + ":" + (ga.length() == 0 ? "-" : ga.toString()));
        }

        System.out.printf("rec\t%s\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s%n",
                label, index, vc.getContig(), vc.getStart(), vc.getEnd(), vc.getID(), alleles,
                vc.hasLog10PError() ? Double.toString(vc.getLog10PError()) : "none",
                filters, attributes.length() == 0 ? "-" : attributes.toString(),
                genotypes.length() == 0 ? "-" : genotypes.toString());
    }

    /**
     * One INFO value carrying the raw text, read under both versions.
     *
     * Going through a whole file rather than calling the transformer directly, because the
     * transformer is chosen by the codec and choosing it is half of the behaviour.
     */
    static void percent(final String raw) {
        System.out.printf("pct\t%s\t%s\t%s%n", escape(raw),
                escape(noteValue("VCFv4.3", raw)), escape(noteValue("VCFv4.2", raw)));
    }

    static String noteValue(final String version, final String raw) {
        try {
            final VCFCodec codec = new VCFCodec();
            final String text = sites(version, "chr1\t100\t.\tA\tT\t50\tPASS\tNOTE=" + raw + "\n");
            final LineIterator it =
                    new LineIteratorImpl(new SynchronousLineReader(new StringReader(text)));
            codec.readActualHeader(it);
            final VariantContext vc = codec.decode(it.next());
            return String.valueOf(vc.getAttribute("NOTE"));
        } catch (final Throwable t) {
            return t.getClass().getSimpleName() + ": " + oneLine(t.getMessage());
        }
    }

    /** Read a file, write it back, and say where the two first differ. */
    static void trip(final String label, final String text) {
        try {
            final VCFCodec codec = new VCFCodec();
            final LineIterator it =
                    new LineIteratorImpl(new SynchronousLineReader(new StringReader(text)));
            final VCFHeader header = (VCFHeader) codec.readActualHeader(it);
            final List<VariantContext> records = new ArrayList<>();
            while (it.hasNext()) {
                final VariantContext vc = codec.decode(it.next());
                if (vc != null) {
                    records.add(vc);
                }
            }
            final ByteArrayOutputStream out = new ByteArrayOutputStream();
            final VariantContextWriter writer = new VariantContextWriterBuilder()
                    .setOutputVCFStream(out)
                    .unsetOption(Options.INDEX_ON_THE_FLY)
                    .build();
            writer.writeHeader(header);
            for (final VariantContext vc : records) {
                writer.add(vc);
            }
            writer.close();

            final String rewritten = out.toString();
            int at = -1;
            for (int i = 0; i < Math.min(text.length(), rewritten.length()); i++) {
                if (text.charAt(i) != rewritten.charAt(i)) {
                    at = i;
                    break;
                }
            }
            if (at == -1 && text.length() != rewritten.length()) {
                at = Math.min(text.length(), rewritten.length());
            }
            System.out.printf("trip\t%s\t%s\t%s\t%s%n", label,
                    at == -1 ? "same" : "differs", at == -1 ? "-" : Integer.toString(at),
                    escape(rewritten));
        } catch (final Throwable t) {
            System.out.printf("err\t%s\t%s\t%s%n", label, t.getClass().getName(),
                    oneLine(t.getMessage()));
        }
    }

    static String render(final Object value) {
        if (value instanceof List) {
            final StringJoiner joined = new StringJoiner(",");
            for (final Object item : (List<?>) value) {
                joined.add(String.valueOf(item));
            }
            return "[" + joined + "]";
        }
        return String.valueOf(value);
    }

    static String escape(final String s) {
        if (s == null) {
            return "null";
        }
        final StringBuilder b = new StringBuilder(s.length());
        for (int i = 0; i < s.length(); i++) {
            final char c = s.charAt(i);
            switch (c) {
                case '\\': b.append("\\\\"); break;
                case '\t': b.append("\\t"); break;
                case '\n': b.append("\\n"); break;
                case '\r': b.append("\\r"); break;
                default:
                    if (c < 0x20 || c > 0x7e) {
                        b.append(String.format("\\u%04x", (int) c));
                    } else {
                        b.append(c);
                    }
            }
        }
        return b.toString();
    }

    static String oneLine(final String s) {
        return s == null ? "null" : escape(s);
    }
}
