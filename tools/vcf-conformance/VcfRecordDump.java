/*
 * Dumps VCF data lines written by htsjdk's VCFEncoder.
 *
 * Output: <case name> <TAB> <the encoded line, escaped>
 *
 * Three properties are under test here, all of them invisible in the format specification:
 *
 *  1. INFO field order. VCFEncoder.write puts the attributes into a `new TreeMap<>()` before
 *     writing, so INFO keys come out in ASCII order of the key regardless of the order they
 *     were set on the VariantContext. This is the third distinct ordering rule in the same
 *     library: SAM header attributes keep insertion order (decision 0009), VCF header lines
 *     sort by rendered whole line (decision 0016), and INFO sorts by key alone.
 *
 *  2. FORMAT key order. VariantContext.calcVCFGenotypeKeys sorts the keys and then puts GT in
 *     front, so GT is first and everything else is in ASCII order.
 *
 *  3. Number formatting. formatQualValue is "%.2f" with a trailing ".00" trimmed, and
 *     formatVCFDouble switches format string on the value's magnitude - with a signed
 *     comparison, so the sign matters as much as the magnitude.
 */

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.*;
import java.util.*;

public class VcfRecordDump {

    static VCFHeader header;
    static VCFEncoder encoder;

    static String esc(String s) {
        return s.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
    }

    static void emit(String name, String value) {
        System.out.println(name + "\t" + esc(value));
    }

    static void emitVc(String name, VariantContext vc) {
        emit(name, encoder.encode(vc));
    }

    static Allele ref(String s) { return Allele.create(s, true); }
    static Allele alt(String s) { return Allele.create(s, false); }

    static VariantContextBuilder base(String... alleles) {
        List<Allele> as = new ArrayList<>();
        as.add(ref(alleles[0]));
        for (int i = 1; i < alleles.length; i++) as.add(alt(alleles[i]));
        return new VariantContextBuilder("src", "chr1", 100, 99 + alleles[0].length(), as);
    }

    public static void main(String[] args) {
        Set<VCFHeaderLine> lines = new LinkedHashSet<>();
        for (String k : new String[] {"AA", "DP", "ZZ", "MM", "AF", "STR", "LIST", "NEG"}) {
            lines.add(new VCFInfoHeaderLine(k, 1, VCFHeaderLineType.String, k));
        }
        lines.add(new VCFInfoHeaderLine("FLAG", 0, VCFHeaderLineType.Flag, "a flag"));
        lines.add(new VCFFormatHeaderLine("GT", 1, VCFHeaderLineType.String, "Genotype"));
        lines.add(new VCFFormatHeaderLine("GQ", 1, VCFHeaderLineType.Integer, "Genotype Quality"));
        lines.add(new VCFFormatHeaderLine("DP", 1, VCFHeaderLineType.Integer, "Depth"));
        lines.add(new VCFFormatHeaderLine("AD", VCFHeaderLineCount.R, VCFHeaderLineType.Integer, "Allele Depths"));
        lines.add(new VCFFormatHeaderLine("PL", VCFHeaderLineCount.G, VCFHeaderLineType.Integer, "Likelihoods"));
        lines.add(new VCFFormatHeaderLine("FT", 1, VCFHeaderLineType.String, "Genotype Filter"));
        lines.add(new VCFFormatHeaderLine("XX", 1, VCFHeaderLineType.String, "extended"));
        lines.add(new VCFFilterHeaderLine("q10", "Quality below 10"));
        lines.add(new VCFFilterHeaderLine("aFilter", "a"));
        lines.add(new VCFFilterHeaderLine("zFilter", "z"));

        header = new VCFHeader(lines, Arrays.asList("SAMPLE1", "SAMPLE2"));
        encoder = new VCFEncoder(header, false, false);

        // ---- the skeleton ------------------------------------------------------------------
        emitVc("minimal", base("A", "T").make());
        emitVc("with_id", base("A", "T").id("rs123").make());
        emitVc("no_alt", base("A").make());
        emitVc("multi_alt", base("A", "T", "C", "GG").make());
        emitVc("span_del", base("A", "T", "*").make());
        emitVc("symbolic", base("A", "<DEL>").make());
        emitVc("breakend", base("A", "A]chr2:456]").make());

        // Bases are uppercased on construction, so a lowercase input is not round-tripped.
        emitVc("lowercase_bases", base("acgt", "a").make());

        // ---- QUAL --------------------------------------------------------------------------
        // VariantContextBuilder takes log10 p-error; the encoder prints -10 * that.
        double[] quals = {0.0, 1.0, 10.0, 29.5, 30.0, 30.004, 30.005, 30.015, 0.125, 2.675,
                          1234.5678, 1e-3, 99999.999};
        for (double q : quals) {
            emitVc("qual_" + q, base("A", "T").log10PError(q / -10.0).make());
        }

        // ---- FILTER ------------------------------------------------------------------------
        emitVc("filter_none", base("A", "T").make());
        emitVc("filter_pass", base("A", "T").passFilters().make());
        emitVc("filter_one", base("A", "T").filter("q10").make());
        emitVc("filter_sorted", base("A", "T").filters("zFilter", "q10", "aFilter").make());

        // ---- INFO --------------------------------------------------------------------------
        Map<String, Object> attrs = new LinkedHashMap<>();
        attrs.put("ZZ", "last");
        attrs.put("AA", "first");
        attrs.put("MM", "middle");
        emitVc("info_sorted", base("A", "T").attributes(attrs).make());

        emitVc("info_flag_true", base("A", "T").attribute("FLAG", true).make());
        emitVc("info_flag_false", base("A", "T").attribute("FLAG", false).make());
        // Number=0 is rejected for every type but Flag, so the encoder's `getCount() != 0`
        // test is reachable only through a real Flag line. A non-flag key whose value is the
        // empty string therefore prints a bare '=' rather than being treated as a flag.
        emitVc("info_empty_string", base("A", "T").attribute("STR", "").make());
        emitVc("info_null", base("A", "T").attribute("STR", null).make());
        emitVc("info_int", base("A", "T").attribute("DP", 42).make());
        emitVc("info_list", base("A", "T").attribute("LIST", Arrays.asList(1, 2, 3)).make());
        emitVc("info_empty_list", base("A", "T").attribute("LIST", new ArrayList<Integer>()).make());
        emitVc("info_int_array", base("A", "T").attribute("LIST", new int[] {4, 5}).make());
        emitVc("info_double", base("A", "T").attribute("AF", 0.5).make());
        emitVc("info_double_small", base("A", "T").attribute("AF", 0.001).make());
        emitVc("info_double_negative", base("A", "T").attribute("NEG", -0.5).make());
        emitVc("info_double_list",
               base("A", "T").attribute("LIST", Arrays.asList(0.5, 0.001, 1e-30)).make());

        // ---- genotypes ---------------------------------------------------------------------
        List<Allele> at = base("A", "T").make().getAlleles();
        Allele a = at.get(0), t = at.get(1);

        emitVc("gt_simple", base("A", "T").genotypes(
                new GenotypeBuilder("SAMPLE1", Arrays.asList(a, t)).make(),
                new GenotypeBuilder("SAMPLE2", Arrays.asList(a, a)).make()).make());

        emitVc("gt_phased", base("A", "T").genotypes(
                new GenotypeBuilder("SAMPLE1", Arrays.asList(a, t)).phased(true).make(),
                new GenotypeBuilder("SAMPLE2", Arrays.asList(a, a)).make()).make());

        emitVc("gt_nocall", base("A", "T").genotypes(
                new GenotypeBuilder("SAMPLE1", Arrays.asList(Allele.NO_CALL, Allele.NO_CALL)).make(),
                new GenotypeBuilder("SAMPLE2", Arrays.asList(a, t)).make()).make());

        emitVc("gt_haploid", base("A", "T").genotypes(
                new GenotypeBuilder("SAMPLE1", Collections.singletonList(t)).make(),
                new GenotypeBuilder("SAMPLE2", Arrays.asList(a, t)).make()).make());

        emitVc("gt_triploid", base("A", "T").genotypes(
                new GenotypeBuilder("SAMPLE1", Arrays.asList(a, t, t)).make(),
                new GenotypeBuilder("SAMPLE2", Arrays.asList(a, t)).make()).make());

        // Only one sample given; the other is synthesised by GenotypeBuilder.createMissing at
        // the record's max ploidy.
        emitVc("gt_one_sample_only", base("A", "T").genotypes(
                new GenotypeBuilder("SAMPLE1", Arrays.asList(a, t, t)).make()).make());

        // Keys are sorted and GT is forced in front, so this comes out GT:AD:DP:GQ:PL.
        emitVc("gt_all_int_fields", base("A", "T").genotypes(
                new GenotypeBuilder("SAMPLE1", Arrays.asList(a, t))
                        .PL(new int[] {0, 10, 100}).AD(new int[] {5, 6}).DP(11).GQ(30).make(),
                new GenotypeBuilder("SAMPLE2", Arrays.asList(a, a))
                        .PL(new int[] {0, 20, 200}).AD(new int[] {9, 0}).DP(9).GQ(40).make())
                .make());

        // One sample has a field the other lacks: the missing one prints '.'.
        emitVc("gt_ragged", base("A", "T").genotypes(
                new GenotypeBuilder("SAMPLE1", Arrays.asList(a, t)).GQ(30).DP(11).make(),
                new GenotypeBuilder("SAMPLE2", Arrays.asList(a, a)).make()).make());

        // Trailing all-missing fields are stripped per sample, so the two samples can end up
        // with a different number of colons.
        emitVc("gt_trailing_stripped", base("A", "T").genotypes(
                new GenotypeBuilder("SAMPLE1", Arrays.asList(a, t)).GQ(30).make(),
                new GenotypeBuilder("SAMPLE2", Arrays.asList(a, a)).DP(9).make()).make());

        emitVc("gt_filtered", base("A", "T").genotypes(
                new GenotypeBuilder("SAMPLE1", Arrays.asList(a, t)).filter("q10").make(),
                new GenotypeBuilder("SAMPLE2", Arrays.asList(a, a)).make()).make());

        emitVc("gt_extended", base("A", "T").genotypes(
                new GenotypeBuilder("SAMPLE1", Arrays.asList(a, t))
                        .attribute("XX", "hello").make(),
                new GenotypeBuilder("SAMPLE2", Arrays.asList(a, a)).make()).make());

        emitVc("gt_extended_double", base("A", "T").genotypes(
                new GenotypeBuilder("SAMPLE1", Arrays.asList(a, t))
                        .attribute("XX", 0.001).make(),
                new GenotypeBuilder("SAMPLE2", Arrays.asList(a, a)).make()).make());

        // ---- formatVCFDouble, directly -----------------------------------------------------
        // The branch structure is `d < 1` then `d < 0.01`, both signed, so every negative value
        // takes the exponent branch no matter how large its magnitude.
        double[] ds = {
            0.0, -0.0, 1.0, -1.0, 0.5, -0.5, 0.01, 0.009999, 0.001, -0.001,
            1e-19, 1e-20, 1e-21, -1e-19, -1e-20, -1e-21, 0.125, 2.675, 0.0005,
            100.0, 123.456, 1e10, -1e10, 1.0 / 3.0, 2.0 / 3.0, 1e-100, 1e300,
            Double.MIN_VALUE, Double.MAX_VALUE,
        };
        for (double d : ds) {
            emit("formatVCFDouble_" + Double.toHexString(d), VCFEncoder.formatVCFDouble(d));
        }
    }
}
