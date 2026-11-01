/*
 * When, and with whose line, a record's genotype columns are refused.
 *
 * `VCFCodec` hands a line's genotype columns to a `LazyGenotypesContext` and parses them on first
 * access, through `createGenotypeMap`, with the codec's state AT THAT TIME: `lineNo` is the codec's
 * counter then, and `parts`, which one message quotes, is the split of the line the codec read
 * last. `VCFFileReader` iterates one record ahead, so a tool touching record N meets a codec that
 * has already read N+1. This dump measures what that does to each refusal, per way of reading:
 *
 *   - stream: `VCFFileReader`, each record touched as it is handed out;
 *   - whole:  `VCFFileReader`, every record collected first, then each touched;
 *   - codec:  `VCFCodec.decode` line by line, each record touched before the next line;
 *   - copy:   every record written untouched, which parses nothing;
 *   - unsorted: a file whose sample names are not sorted, which `parseVCFLine` decodes at the read.
 *
 * Output (the temporary file's path is written as <file>):
 *     size\t<scenario>\t<record>\t<size()>\t<isEmpty()>
 *     access\t<scenario>\t<record>\tok\t<genotype count>
 *     access\t<scenario>\t<record>\tthrow\t<class>\t<message>
 *     write\t<scenario>\tok\t<the file written, escaped>
 *     write\t<scenario>\tthrow\t<class>\t<message>
 *     read\t<scenario>\tthrow\t<class>\t<message>
 */

import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.SynchronousLineReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

import java.io.File;
import java.io.StringReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;

public class VcfDeferredGenotypesDump {

    static Path dir;

    static String escape(final String s) {
        return s.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
    }

    static String header(final String samples) {
        return "##fileformat=VCFv4.2\n"
                + "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
                + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                + "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Quality\">\n"
                + "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
                + "##contig=<ID=chr1,length=100000>\n"
                + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + samples + "\n";
    }

    // A good record after each bad one, with its own POS, so a message quoting "the line the codec
    // read last" says which line that was.
    static final String BODY =
            "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\tGT:GQ\t0/1:60\t1/1:70\n"
            // s0 has three values for two keys: "There are too many keys ...", which quotes parts.
            + "chr1\t200\t.\tC\tG\t50\tPASS\tDP=10\tGT:GQ\t0/1:60:99\t1/1:70\n"
            + "chr1\t300\t.\tG\tA\t50\tPASS\tDP=10\tGT:GQ:DP\t0/0:30:5\t0/1:40:6\n"
            // One sample column for two samples: the count check, which names chr1:400 itself.
            + "chr1\t400\t.\tT\tC\t50\tPASS\tDP=10\tGT\t0/1\n"
            // An allele index past REF/ALT: InternalCodecException, which carries no line.
            + "chr1\t500\t.\tA\tG\t50\tPASS\tDP=10\tGT\t0/3\t0/1\n"
            + "chr1\t600\t.\tC\tT\t50\tPASS\tDP=10\tGT:GQ\t1/1:80\t0/1:90\n";

    static String clean(final Throwable t) {
        final String message = String.valueOf(t.getMessage()).replace(dir.toString(), "<dir>");
        return t.getClass().getName() + "\t" + escape(message);
    }

    static File write(final String name, final String text) throws Exception {
        final File f = dir.resolve(name).toFile();
        Files.write(f.toPath(), text.getBytes(StandardCharsets.UTF_8));
        return f;
    }

    static void touch(final String scenario, final int record, final VariantContext vc) {
        System.out.printf("size\t%s\t%d\t%d\t%s%n", scenario, record,
                vc.getGenotypes().size(), vc.getGenotypes().isEmpty());
        try {
            System.out.printf("access\t%s\t%d\tok\t%d%n", scenario, record,
                    vc.getGenotypes().get(0) == null ? -1 : vc.getGenotypes().getSampleNames().size());
        } catch (final RuntimeException e) {
            System.out.printf("access\t%s\t%d\tthrow\t%s%n", scenario, record, clean(e));
        }
    }

    public static void main(final String[] args) throws Exception {
        dir = Files.createTempDirectory("deferred");
        final File sorted = write("sorted.vcf", header("s0\ts1") + BODY);

        // stream: touched as handed out, with the reader one record ahead.
        try (VCFFileReader reader = new VCFFileReader(sorted, false)) {
            int i = 0;
            for (final VariantContext vc : reader) {
                touch("stream", i++, vc);
            }
        }

        // whole: every record collected, then touched.
        final List<VariantContext> all = new ArrayList<>();
        try (VCFFileReader reader = new VCFFileReader(sorted, false)) {
            for (final VariantContext vc : reader) {
                all.add(vc);
            }
        }
        for (int i = 0; i < all.size(); i++) {
            touch("whole", i, all.get(i));
        }

        // codec: no reader, no read-ahead.
        {
            final VCFCodec codec = new VCFCodec();
            final LineIterator it = new LineIteratorImpl(
                    new SynchronousLineReader(new StringReader(header("s0\ts1") + BODY)));
            codec.readActualHeader(it);
            int i = 0;
            while (it.hasNext()) {
                final VariantContext vc = codec.decode(it.next());
                if (vc != null) {
                    touch("codec", i++, vc);
                }
            }
        }

        // copy: nothing parses.
        copy("copy", sorted);

        // unsorted: parseVCFLine decodes at the read, inside the reader.
        final File unsorted = write("unsorted.vcf", header("s1\ts0") + BODY);
        try (VCFFileReader reader = new VCFFileReader(unsorted, false)) {
            int i = 0;
            for (final VariantContext vc : reader) {
                touch("unsorted", i++, vc);
            }
        } catch (final RuntimeException e) {
            System.out.printf("read\tunsorted\tthrow\t%s%n", clean(e));
        }
    }

    static void copy(final String scenario, final File in) throws Exception {
        final File out = dir.resolve(scenario + ".out.vcf").toFile();
        try (VCFFileReader reader = new VCFFileReader(in, false)) {
            final VCFHeader header = reader.getFileHeader();
            try (VariantContextWriter writer = new VariantContextWriterBuilder()
                    .setOutputFile(out)
                    .setOptions(EnumSet.noneOf(Options.class))
                    .build()) {
                writer.writeHeader(header);
                for (final VariantContext vc : reader) {
                    writer.add(vc);
                }
            }
            final String text = new String(Files.readAllBytes(out.toPath()), StandardCharsets.UTF_8);
            System.out.printf("write\t%s\tok\t%s%n", scenario, escape(text));
        } catch (final RuntimeException e) {
            System.out.printf("write\t%s\tthrow\t%s%n", scenario, clean(e));
        }
    }
}
