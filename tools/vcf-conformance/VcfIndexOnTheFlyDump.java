/*
 * Writing a VCF and its .idx in one pass, which is what every GATK tool that emits a VCF does.
 *
 * The pieces exist separately already: VcfFileDump measures the file, TribbleIndexWriteDump
 * measures an index built from a feature list. What neither can see is the join, and the join is
 * where the file's bytes and the index's numbers have to agree.
 *
 * Four things decide that agreement and none of them is in either format:
 *
 *   - INDEX_ON_THE_FLY IS ON BY DEFAULT. VariantContextWriterBuilder.DEFAULT_OPTIONS is
 *     EnumSet.of(Options.INDEX_ON_THE_FLY), so a caller who asks for nothing gets an index, and a
 *     caller who asks for nothing AND supplies no dictionary gets an exception rather than a file;
 *   - THE POSITION RECORDED IS THE ONE BEFORE THE RECORD. IndexingVariantContextWriter.add feeds
 *     the indexer locationSource.getPosition() and only then does VCFWriter.add write the line, so
 *     a record's block starts where the record starts. The final position, handed to
 *     finalizeIndex, is the whole file's length INCLUDING the header;
 *   - THE HEADER IS COUNTED. Positions are absolute in the output stream, so the first record's
 *     position is the header's length, and an index built from a feature list that forgot the
 *     header is uniformly off by it;
 *   - THE SEQUENCE DICTIONARY BECOMES PROPERTIES, NOT A FLAG. setIndexSequenceDictionary writes
 *     one DICT:<contig> = <length> property per sequence, in dictionary order, and the flags field
 *     stays zero: the SEQUENCE_DICTIONARY_FLAG is only read for version < 3. Those properties go in
 *     BEFORE the dynamic creator's four statistics, because finalizeIndex copies the creator's own
 *     map first and appends the statistics after.
 *
 * The index type is not stated by the caller either: the writer always uses a DynamicIndexCreator
 * with FOR_SEEK_TIME, so which layout lands beside a VCF depends on the variants in it.
 *
 * Output:
 *
 *     vcf\t<label>\t<file length>\t<header length>\t<record count>\t<record positions>
 *     idx\t<label>\t<index class>\t<timestamp offset>\t<base64 of the .idx, timestamp zeroed>
 *     prop\t<label>\t<key>=<value>;...
 *     chr\t<label>\t<name>\t<binWidth>\t<nBins>\t<longestFeature>\t<nFeatures>\t<positions>
 *     err\t<label>\t<class>\t<message>
 *
 * Usage: VcfIndexOnTheFlyDump
 */

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.tribble.index.linear.LinearIndex;
import htsjdk.tribble.util.LittleEndianInputStream;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import java.io.BufferedInputStream;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Base64;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;

public class VcfIndexOnTheFlyDump {

    static Path dir;

    public static void main(final String[] args) throws Exception {
        System.out.println("# VcfIndexOnTheFlyDump: a VCF and its .idx written in one pass");

        dir = Path.of("vcf-index-dump");
        emptyDirectory(dir);
        Files.createDirectories(dir);

        // Two contigs, records on both, so the index has two chromosome records and the second
        // one's positions depend on how long the first one's lines were.
        write("two-contigs", header(), List.of(
                vc("chr1", 100, "A", "T"),
                vc("chr1", 20000, "C", "G"),
                vc("chr2", 50, "GG", "G")));

        // One record, so the index has a single block and the optimizer breaks immediately.
        write("one-record", header(), List.of(vc("chr1", 100, "A", "T")));

        // Enough records that the linear score passes the threshold and the dynamic creator's
        // choice may go the other way. Which one wins is the measurement.
        final List<VariantContext> many = new ArrayList<>();
        for (int i = 0; i < 2000; i++) {
            many.add(vc("chr1", 100 + i * 5, "A", "T"));
        }
        write("many-records", header(), many);

        // No records at all: the header is written, the index is finalized at the header's length,
        // and there is no contig to close.
        write("header-only", header(), List.of());

        // A record on a contig the dictionary does not declare. The dictionary decides the DICT:
        // properties; the index's chromosome records come from the features.
        write("undeclared-contig", header(), List.of(vc("chrX", 9, "A", "T")));

        // A long INFO value, so one record's line is much longer than its neighbours and the block
        // positions are visibly uneven.
        write("uneven-lines", header(), List.of(
                vc("chr1", 100, "A", "T"),
                longInfo("chr1", 200),
                vc("chr1", 300, "A", "T")));

        // The refusals. Indexing on the fly needs a dictionary, and it needs a file.
        noDictionary("no-dictionary");
        toStream("to-a-stream");

        // And the same file written with indexing off, which is the byte comparison that says the
        // index costs the VCF nothing.
        write("indexing-off", header(), List.of(vc("chr1", 100, "A", "T")), false);
    }

    static VCFHeader header() {
        final Set<VCFHeaderLine> lines = new LinkedHashSet<>();
        lines.add(new VCFInfoHeaderLine("DP", 1, VCFHeaderLineType.Integer, "Depth"));
        lines.add(new VCFInfoHeaderLine("NOTE", 1, VCFHeaderLineType.String, "A note"));
        for (int i = 0; i < 2; i++) {
            final String id = "chr" + (i + 1);
            final Map<String, String> fields = new LinkedHashMap<>();
            fields.put("ID", id);
            fields.put("length", String.valueOf(100000 * (i + 1)));
            lines.add(new VCFContigHeaderLine(fields, i));
        }
        return new VCFHeader(lines);
    }

    static SAMSequenceDictionary dictionary() {
        final SAMSequenceDictionary dict = new SAMSequenceDictionary();
        dict.addSequence(new SAMSequenceRecord("chr1", 100000));
        dict.addSequence(new SAMSequenceRecord("chr2", 200000));
        return dict;
    }

    static VariantContext vc(final String contig, final int start, final String ref, final String alt) {
        return new VariantContextBuilder("src", contig, start, start + ref.length() - 1,
                Arrays.asList(Allele.create(ref, true), Allele.create(alt, false)))
                .attribute("DP", 10)
                .make();
    }

    static VariantContext longInfo(final String contig, final int start) {
        return new VariantContextBuilder("src", contig, start, start,
                Arrays.asList(Allele.create("A", true), Allele.create("T", false)))
                .attribute("NOTE", "x".repeat(400))
                .make();
    }

    static void write(final String label, final VCFHeader header, final List<VariantContext> records) {
        write(label, header, records, true);
    }

    static void write(final String label, final VCFHeader header, final List<VariantContext> records,
                      final boolean index) {
        try {
            final Path vcf = dir.resolve(label + ".vcf");
            final VariantContextWriterBuilder builder = new VariantContextWriterBuilder()
                    .setOutputPath(vcf)
                    .setReferenceDictionary(dictionary());
            if (!index) {
                builder.unsetOption(Options.INDEX_ON_THE_FLY);
            }
            // Deliberately NOT setting INDEX_ON_THE_FLY in the indexing case: it is the default,
            // and that being the default is half of what this row measures.
            try (final VariantContextWriter writer = builder.build()) {
                writer.writeHeader(header);
                for (final VariantContext vc : records) {
                    writer.add(vc);
                }
            }

            final byte[] text = Files.readAllBytes(vcf);
            // The header's length is where the first record starts, computed from the text rather
            // than from the header object so that it is the same number the writer counted.
            int headerLength = 0;
            int seen = 0;
            for (int i = 0; i < text.length; i++) {
                if (text[i] == '\n') {
                    seen++;
                    final String line = new String(text, headerLength, i - headerLength + 1);
                    headerLength = i + 1;
                    if (!line.startsWith("#")) {
                        headerLength = headerLength - line.length();
                        break;
                    }
                }
            }
            final StringJoiner positions = new StringJoiner(",");
            int at = headerLength;
            for (int i = 0; i < records.size(); i++) {
                positions.add(Integer.toString(at));
                while (at < text.length && text[at] != '\n') {
                    at++;
                }
                at++;
            }
            System.out.printf("vcf\t%s\t%d\t%d\t%d\t%s%n", label, text.length, headerLength,
                    records.size(), positions.length() == 0 ? "-" : positions.toString());

            final Path idx = dir.resolve(label + ".vcf.idx");
            if (!Files.exists(idx)) {
                System.out.printf("idx\t%s\tnone\t-\t-%n", label);
                return;
            }
            emitIndex(label, Files.readAllBytes(idx));
        } catch (final Throwable t) {
            System.out.printf("err\t%s\t%s\t%s%n", label, t.getClass().getName(),
                    oneLine(t.getMessage()));
        }
    }

    static void emitIndex(final String label, final byte[] bytes) throws Exception {
        int cursor = 12;
        while (bytes[cursor] != 0) {
            cursor++;
        }
        final int timestampOffset = cursor + 1 + 8;
        final byte[] masked = bytes.clone();
        Arrays.fill(masked, timestampOffset, timestampOffset + 8, (byte) 0);

        try (final LittleEndianInputStream dis = new LittleEndianInputStream(
                new BufferedInputStream(new ByteArrayInputStream(bytes)))) {
            dis.readInt();
            final int type = dis.readInt();
            dis.readInt();
            dis.readString();
            dis.readLong();
            dis.readLong();
            dis.readString();
            final int flags = dis.readInt();
            System.out.printf("idx\t%s\t%s\t%d\t%s%n", label,
                    type == LinearIndex.INDEX_TYPE ? "LinearIndex" : "IntervalTreeIndex",
                    timestampOffset, Base64.getEncoder().encodeToString(masked));

            int count = dis.readInt();
            final StringJoiner properties = new StringJoiner(";");
            while (count-- > 0) {
                properties.add(dis.readString() + "=" + dis.readString());
            }
            // The flags travel with the properties because the sequence dictionary used to live in
            // the flags and now lives in the properties, and the row should show both.
            System.out.printf("prop\t%s\tflags=%d;%s%n", label, flags, properties);

            if (type != LinearIndex.INDEX_TYPE) {
                return;
            }
            int contigs = dis.readInt();
            while (contigs-- > 0) {
                final String name = dis.readString();
                final int binWidth = dis.readInt();
                final int nBins = dis.readInt();
                final int longestFeature = dis.readInt();
                dis.readInt();
                final int nFeatures = dis.readInt();
                final StringJoiner blocks = new StringJoiner(",");
                for (int i = 0; i <= nBins; i++) {
                    blocks.add(Long.toString(dis.readLong()));
                }
                System.out.printf("chr\t%s\t%s\t%d\t%d\t%d\t%d\t%s%n", label, name, binWidth,
                        nBins, longestFeature, nFeatures, blocks);
            }
        }
    }

    /** Indexing on the fly with no dictionary, which is refused at build time. */
    static void noDictionary(final String label) {
        try {
            final Path vcf = dir.resolve(label + ".vcf");
            new VariantContextWriterBuilder().setOutputPath(vcf).build().close();
            System.out.printf("err\t%s\tnone\tno exception%n", label);
        } catch (final Throwable t) {
            System.out.printf("err\t%s\t%s\t%s%n", label, t.getClass().getName(),
                    oneLine(t.getMessage()));
        }
    }

    /** Indexing on the fly to a stream, which has no path to write the index beside. */
    static void toStream(final String label) {
        try {
            final ByteArrayOutputStream out = new ByteArrayOutputStream();
            new VariantContextWriterBuilder()
                    .setOutputVCFStream(out)
                    .setReferenceDictionary(dictionary())
                    .build()
                    .close();
            System.out.printf("err\t%s\tnone\tno exception%n", label);
        } catch (final Throwable t) {
            System.out.printf("err\t%s\t%s\t%s%n", label, t.getClass().getName(),
                    oneLine(t.getMessage()));
        }
    }

    static void emptyDirectory(final Path path) throws Exception {
        if (!Files.exists(path)) {
            return;
        }
        try (final var walk = Files.walk(path)) {
            walk.sorted(Comparator.reverseOrder()).forEach(child -> {
                try {
                    Files.delete(child);
                } catch (final Exception ignored) {
                    // Not what this dump measures.
                }
            });
        }
    }

    static String oneLine(final String s) {
        return s == null ? "null" : s.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
    }
}
