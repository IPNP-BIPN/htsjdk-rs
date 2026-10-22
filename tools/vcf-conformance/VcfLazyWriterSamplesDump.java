/*
 * What a copied record's genotype columns become under a writer header whose samples are not the
 * file's.
 *
 * `VCFEncoder` writes a `LazyGenotypesContext`'s text verbatim without asking the header anything:
 *
 *     if (gc.isLazyWithData() && ... instanceof String) { append(FIELD_SEPARATOR); append(text); }
 *
 * So a record nothing has looked at keeps its own columns, in its own order and number, whatever
 * samples the writer's header names: reordered, a subset, a superset, renamed, or none. A record
 * that was looked at is decoded and written by sample name in the HEADER's order instead. And a
 * file whose samples are not sorted is decoded at the read (`parseVCFLine`), so its records take
 * the second path even untouched, which is the `MergeVcfs` case of inputs ordered differently.
 *
 * Each input is read with `VCFFileReader` and written with `VariantContextWriter` under
 * `new VCFHeader(<the file's lines>, <samples>)`, either untouched (copy) or after
 * `getGenotypes().get(0)` on every record (touch).
 *
 * Output:
 *     write\t<input>\t<writer samples>\t<mode>\tok\t<the file written, escaped>
 *     write\t<input>\t<writer samples>\t<mode>\tthrow\t<class>\t<message>
 */

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.LinkedHashSet;
import java.util.List;

public class VcfLazyWriterSamplesDump {

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
                + "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Depths\">\n"
                + "##contig=<ID=chr1,length=100000>\n"
                + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + samples + "\n";
    }

    // A FORMAT out of the computed order (GT:GQ:DP, which re-encodes as GT:DP:GQ) and a trailing
    // missing field (which re-encoding trims), so a copied line and a re-encoded one never agree.
    static final String BODY =
            "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\tGT:GQ:DP\t0/1:60:10\t1/1:70:20\n"
            + "chr1\t200\t.\tG\tA\t50\tPASS\tDP=10\tGT:GQ:AD\t0/1:60:5,6\t1/1:.:.\n";

    public static void main(final String[] args) throws Exception {
        dir = Files.createTempDirectory("lazy-writer-samples");
        final File sorted = write("sorted.vcf", header("s0\ts1") + BODY);
        for (final String samples : new String[] {"s0,s1", "s1,s0", "s1", "s0,s1,s2", "t0,t1", ""}) {
            emit("sorted", sorted, samples, false);
            emit("sorted", sorted, samples, true);
        }
        final File unsorted = write("unsorted.vcf", header("s1\ts0") + BODY);
        for (final String samples : new String[] {"s1,s0", "s0,s1"}) {
            emit("unsorted", unsorted, samples, false);
            emit("unsorted", unsorted, samples, true);
        }
    }

    static File write(final String name, final String text) throws Exception {
        final File f = dir.resolve(name).toFile();
        Files.write(f.toPath(), text.getBytes(StandardCharsets.UTF_8));
        return f;
    }

    static void emit(final String input, final File in, final String samples, final boolean touch)
            throws Exception {
        final String mode = touch ? "touch" : "copy";
        final File out = dir.resolve("out.vcf").toFile();
        final List<String> names =
                samples.isEmpty() ? new ArrayList<>() : Arrays.asList(samples.split(","));
        try (VCFFileReader reader = new VCFFileReader(in, false)) {
            final VCFHeader header = new VCFHeader(
                    new LinkedHashSet<VCFHeaderLine>(reader.getFileHeader().getMetaDataInInputOrder()),
                    names);
            try (VariantContextWriter writer = new VariantContextWriterBuilder()
                    .setOutputFile(out)
                    .setOptions(EnumSet.noneOf(Options.class))
                    .build()) {
                writer.writeHeader(header);
                for (final VariantContext vc : reader) {
                    if (touch) {
                        // A READ: `ensureSampleOrdering` decodes, and `decode()` drops the text.
                        vc.getGenotypes().get(0);
                    }
                    writer.add(vc);
                }
            }
            final String text = new String(Files.readAllBytes(out.toPath()), StandardCharsets.UTF_8);
            System.out.printf("write\t%s\t%s\t%s\tok\t%s%n", input, samples, mode, escape(text));
        } catch (final RuntimeException e) {
            final String message = String.valueOf(e.getMessage()).replace(dir.toString(), "<dir>");
            System.out.printf("write\t%s\t%s\t%s\tthrow\t%s\t%s%n",
                    input, samples, mode, e.getClass().getName(), escape(message));
        }
    }
}
