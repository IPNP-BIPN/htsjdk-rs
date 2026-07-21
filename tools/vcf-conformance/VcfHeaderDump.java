/*
 * Dumps VCF headers written by htsjdk's VCFWriter.
 *
 * Output: <case name> <TAB> <the header text, newlines escaped>
 *
 * The property under test is the line ordering. VCFHeader.getMetaDataInSortedOrder puts every
 * line into a TreeSet, and VCFHeaderLine.compareTo compares the *rendered strings*. So the
 * order is lexicographic over whole lines, payload included, which is the opposite of the SAM
 * header's insertion order (decision 0009).
 */

import htsjdk.variant.vcf.*;
import htsjdk.variant.variantcontext.writer.*;
import java.io.ByteArrayOutputStream;
import java.util.*;

public class VcfHeaderDump {

    static void emit(String name, Set<VCFHeaderLine> lines, List<String> samples) {
        VCFHeader header = samples == null ? new VCFHeader(lines) : new VCFHeader(lines, samples);
        ByteArrayOutputStream out = new ByteArrayOutputStream();
        VariantContextWriter w = new VariantContextWriterBuilder()
                .setOutputVCFStream(out)
                .unsetOption(Options.INDEX_ON_THE_FLY)
                .build();
        w.writeHeader(header);
        w.close();
        System.out.println(name + "\t"
                + out.toString().replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"));
    }

    public static void main(String[] args) {
        emit("minimal", new LinkedHashSet<>(), null);

        // One of each typed line.
        Set<VCFHeaderLine> typed = new LinkedHashSet<>();
        typed.add(new VCFInfoHeaderLine("DP", 1, VCFHeaderLineType.Integer, "Total Depth"));
        typed.add(new VCFFormatHeaderLine("GT", 1, VCFHeaderLineType.String, "Genotype"));
        typed.add(new VCFFilterHeaderLine("q10", "Quality below 10"));
        typed.add(new VCFContigHeaderLine(
                new LinkedHashMap<String, String>() {{ put("ID", "chr1"); put("length", "249250621"); }},
                0));
        emit("one_of_each", typed, null);

        // Deliberately inserted in an order that is neither sorted nor reverse-sorted, so the
        // TreeSet ordering is observable.
        Set<VCFHeaderLine> shuffled = new LinkedHashSet<>();
        shuffled.add(new VCFInfoHeaderLine("ZZ", 1, VCFHeaderLineType.Integer, "last by id"));
        shuffled.add(new VCFFilterHeaderLine("zFilter", "z"));
        shuffled.add(new VCFInfoHeaderLine("AA", 1, VCFHeaderLineType.String, "first by id"));
        shuffled.add(new VCFFormatHeaderLine("ZQ", 1, VCFHeaderLineType.Float, "z format"));
        shuffled.add(new VCFFilterHeaderLine("aFilter", "a"));
        shuffled.add(new VCFFormatHeaderLine("AQ", 1, VCFHeaderLineType.Float, "a format"));
        emit("shuffled", shuffled, null);

        // Values that trigger the quoting rule from every direction.
        Set<VCFHeaderLine> quoting = new LinkedHashSet<>();
        quoting.add(new VCFInfoHeaderLine("NOSPACE", 1, VCFHeaderLineType.String, "nospace"));
        quoting.add(new VCFInfoHeaderLine("WITHSPACE", 1, VCFHeaderLineType.String, "with space"));
        quoting.add(new VCFInfoHeaderLine("WITHCOMMA", 1, VCFHeaderLineType.String, "with,comma"));
        quoting.add(new VCFInfoHeaderLine("WITHQUOTE", 1, VCFHeaderLineType.String, "with\"quote"));
        quoting.add(new VCFHeaderLine("unstructured", "plain value"));
        quoting.add(new VCFHeaderLine("alsoUnstructured", "value,with,commas"));
        emit("quoting", quoting, null);

        // Number cardinalities, which render as A, G, R and '.'.
        Set<VCFHeaderLine> counts = new LinkedHashSet<>();
        counts.add(new VCFInfoHeaderLine("FIXED", 3, VCFHeaderLineType.Integer, "three"));
        counts.add(new VCFInfoHeaderLine("PERALT", VCFHeaderLineCount.A, VCFHeaderLineType.Float, "per alt"));
        counts.add(new VCFInfoHeaderLine("PERGT", VCFHeaderLineCount.G, VCFHeaderLineType.Float, "per genotype"));
        counts.add(new VCFInfoHeaderLine("PERALLELE", VCFHeaderLineCount.R, VCFHeaderLineType.Float, "per allele"));
        counts.add(new VCFInfoHeaderLine("UNBOUNDED", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "any"));
        counts.add(new VCFInfoHeaderLine("FLAG", 0, VCFHeaderLineType.Flag, "a flag"));
        emit("cardinalities", counts, null);

        // Samples, which add the FORMAT column and the sample names.
        Set<VCFHeaderLine> withSamples = new LinkedHashSet<>();
        withSamples.add(new VCFFormatHeaderLine("GT", 1, VCFHeaderLineType.String, "Genotype"));
        emit("samples", withSamples, Arrays.asList("NA12878", "NA12891", "NA12892"));

        // A contig list, which is what a real header is mostly made of.
        Set<VCFHeaderLine> contigs = new LinkedHashSet<>();
        for (int i = 1; i <= 12; i++) {
            final int idx = i;
            contigs.add(new VCFContigHeaderLine(
                    new LinkedHashMap<String, String>() {{
                        put("ID", "chr" + idx);
                        put("length", String.valueOf(250000000 - idx * 1000000));
                    }}, idx - 1));
        }
        emit("contigs", contigs, null);
    }
}
