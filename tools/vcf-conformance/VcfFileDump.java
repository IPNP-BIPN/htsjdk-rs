/*
 * Writes complete VCF files with htsjdk's VariantContextWriter and emits them.
 *
 * Output: <case name> <TAB> <the whole file, escaped>
 *
 * The header and the data lines each have their own conformance suite already. What this adds
 * is the join between them: the version line the writer substitutes for the header's own, the
 * absence of a separator between header and first record, and the trailing newline.
 */

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.*;
import htsjdk.variant.variantcontext.writer.*;
import java.io.ByteArrayOutputStream;
import java.util.*;

public class VcfFileDump {

    static void emit(String name, VCFHeader header, List<VariantContext> records) {
        ByteArrayOutputStream out = new ByteArrayOutputStream();
        VariantContextWriter w = new VariantContextWriterBuilder()
                .setOutputVCFStream(out)
                .unsetOption(Options.INDEX_ON_THE_FLY)
                .build();
        w.writeHeader(header);
        for (VariantContext vc : records) w.add(vc);
        w.close();
        System.out.println(name + "\t"
                + out.toString().replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"));
    }

    static Allele ref(String s) { return Allele.create(s, true); }
    static Allele alt(String s) { return Allele.create(s, false); }

    static VariantContext vc(String contig, int start, String r, String a) {
        return new VariantContextBuilder("src", contig, start, start + r.length() - 1,
                Arrays.asList(ref(r), alt(a))).make();
    }

    static VCFHeader header(boolean withSamples) {
        Set<VCFHeaderLine> lines = new LinkedHashSet<>();
        lines.add(new VCFInfoHeaderLine("DP", 1, VCFHeaderLineType.Integer, "Total Depth"));
        lines.add(new VCFFormatHeaderLine("GT", 1, VCFHeaderLineType.String, "Genotype"));
        lines.add(new VCFFormatHeaderLine("GQ", 1, VCFHeaderLineType.Integer, "Genotype Quality"));
        lines.add(new VCFFilterHeaderLine("q10", "Quality below 10"));
        for (int i = 0; i < 3; i++) {
            final String id = "chr" + (i + 1);
            lines.add(new VCFContigHeaderLine(
                    new LinkedHashMap<String, String>() {{ put("ID", id); put("length", "1000"); }}, i));
        }
        return withSamples
                ? new VCFHeader(lines, Arrays.asList("s1", "s2"))
                : new VCFHeader(lines);
    }

    public static void main(String[] args) {
        emit("header_only", header(false), List.of());
        emit("header_only_with_samples", header(true), List.of());

        emit("one_record", header(false), List.of(vc("chr1", 100, "A", "T")));

        emit("many_records", header(false), List.of(
                vc("chr1", 100, "A", "T"),
                vc("chr1", 200, "C", "G"),
                vc("chr2", 50, "GG", "G"),
                vc("chr3", 1, "T", "TTTT")));

        // With genotypes, so the FORMAT and sample columns join the header's column line.
        List<Allele> at = Arrays.asList(ref("A"), alt("T"));
        VariantContext genotyped = new VariantContextBuilder("src", "chr1", 100, 100, at)
                .genotypes(
                        new GenotypeBuilder("s1", Arrays.asList(at.get(0), at.get(1))).GQ(30).make(),
                        new GenotypeBuilder("s2", Arrays.asList(at.get(0), at.get(0))).GQ(40).make())
                .make();
        emit("genotyped", header(true), List.of(genotyped));

        // A header that declares a fileformat of its own. The writer substitutes VERSION_LINE
        // and skips the header's line rather than writing both.
        Set<VCFHeaderLine> withFormat = new LinkedHashSet<>();
        withFormat.add(new VCFInfoHeaderLine("DP", 1, VCFHeaderLineType.Integer, "Total Depth"));
        withFormat.add(new VCFHeaderLine("fileformat", "VCFv4.3"));
        emit("header_declares_its_own_fileformat", new VCFHeader(withFormat),
             List.of(vc("chr1", 100, "A", "T")));

        // Records out of coordinate order: the writer does not sort, it writes what it is given.
        emit("unsorted_records", header(false), List.of(
                vc("chr3", 500, "A", "T"),
                vc("chr1", 100, "A", "T"),
                vc("chr2", 300, "A", "T")));
    }
}
