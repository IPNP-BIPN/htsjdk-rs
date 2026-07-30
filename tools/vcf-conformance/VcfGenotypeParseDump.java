/*
 * Reading the genotype columns of a VCF data line, taken from the reference.
 *
 * VcfRecordParseDump stops at the eight site columns. This is the rest: FORMAT, the per-sample
 * values, and the GT field, which htsjdk defers through LazyGenotypesContext until something asks
 * for it.
 *
 * Four behaviours decide what a genotype contains and none is in the VCF specification.
 *
 *   - the GT separators are three characters, not two. VCFConstants.PHASING_TOKENS is "/|\", and
 *     the split is a StringTokenizer, which drops empty tokens. So a backslash separates alleles
 *     like a slash does, and 0//1, /0/1 and 0/1/ all yield the same two alleles as 0/1;
 *   - a malformed AD or PL is silently dropped. decodeInts catches NumberFormatException and
 *     returns null, which the builder stores as "no AD". Two lines further down the same method,
 *     DP goes through a bare Integer.parseInt whose failure is not caught at all. The same
 *     malformed integer therefore either disappears or aborts the record, decided by its key;
 *   - GT must be at position 0 when present, and before VCF 4.1 it must be present at all, so the
 *     same record is valid or invalid depending on a header line;
 *   - a key with no value is skipped, and a value of "." is also skipped, but only the second one
 *     had to be looked at: the difference shows in FT, where "." means unfiltered rather than
 *     absent.
 *
 * The genotype is dumped field by field rather than through toString, so a difference names which
 * field diverged.
 *
 * Output:
 *
 *     gt\t<label>\t<sample>\t<alleles>\t<phased>\t<gq>\t<dp>\t<ad>\t<pl>\t<filters>\t<extended>
 *     gterror\t<label>\t<class>\t<message>
 *
 * Usage: VcfGenotypeParseDump
 */

import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.SynchronousLineReader;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

import java.io.StringReader;
import java.util.Map;
import java.util.StringJoiner;
import java.util.TreeMap;

public class VcfGenotypeParseDump {

    static final String COMMON_META = "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
            + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            + "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Quality\">\n"
            + "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
            + "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Depths\">\n"
            + "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Likelihoods\">\n"
            + "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Genotype filter\">\n"
            + "##FORMAT=<ID=XX,Number=1,Type=String,Description=\"Anything\">\n"
            + "##FILTER=<ID=LowQual,Description=\"Low quality\">\n";

    /** Two samples, VCF 4.2. */
    static final String HEADER = "##fileformat=VCFv4.2\n" + COMMON_META
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA1\tNA2\n";

    /** One sample, VCF 4.0, where a record with no GT is a refusal. */
    static final String HEADER_V40 = "##fileformat=VCFv4.0\n" + COMMON_META
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA1\n";

    /** One sample, VCF 4.1, where the same record is fine. */
    static final String HEADER_V41 = "##fileformat=VCFv4.1\n" + COMMON_META
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA1\n";

    public static void main(final String[] args) {
        System.out.println("# VcfGenotypeParseDump: reading the genotype columns");

        // The ordinary shapes.
        two("plain", "GT\t0/1\t1/1");
        two("phased", "GT\t0|1\t1|1");
        two("haploid", "GT\t0\t1");
        two("no-call", "GT\t./.\t0/1");
        two("full-format", "GT:GQ:DP:AD:PL\t0/1:99:30:10,20:100,0,200\t1/1:50:10:0,10:255,30,0");
        // The three separators, and the tokenizer that drops empty tokens.
        two("backslash-separator", "GT\t0\\1\t1\\1");
        two("doubled-separator", "GT\t0//1\t1/1");
        two("leading-separator", "GT\t/0/1\t1/1");
        two("trailing-separator", "GT\t0/1/\t1/1");
        // Missing values and missing keys.
        two("missing-gq", "GT:GQ\t0/1:.\t1/1:50");
        two("short-value-list", "GT:GQ:DP\t0/1:99\t1/1:50:10");
        two("gq-minus-one", "GT:GQ\t0/1:-1\t1/1:50");
        // Math.round is floor(x + 0.5), so a .5 rounds up and a negative .5 rounds towards zero.
        two("gq-half", "GT:GQ\t0/1:2.5\t1/1:3.5");
        two("gq-negative-half", "GT:GQ\t0/1:-2.5\t1/1:-1.5");
        // The asymmetry between decodeInts and Integer.parseInt.
        two("ad-not-a-number", "GT:AD\t0/1:1,x\t1/1:0,10");
        two("pl-not-a-number", "GT:PL\t0/1:1,x,3\t1/1:0,10,20");
        two("dp-not-a-number", "GT:DP\t0/1:x\t1/1:10");
        // FT, which goes through the record's filter rules.
        two("ft-pass", "GT:FT\t0/1:PASS\t1/1:LowQual");
        two("ft-missing", "GT:FT\t0/1:.\t1/1:PASS");
        two("ft-two", "GT:FT\t0/1:LowQual;q10\t1/1:PASS");
        // An unreserved key, which becomes a plain string attribute.
        two("extended-key", "GT:XX\t0/1:hello\t1/1:world");
        // No GT at all, at 4.2.
        two("no-gt", "GQ\t99\t50");
        // The refusals.
        two("gt-not-first", "GQ:GT\t99:0/1\t50:1/1");
        two("too-many-values", "GT\t0/1:99\t1/1");
        two("too-few-columns", "GT\t0/1");
        two("allele-index-out-of-range", "GT\t0/2\t1/1");
        two("allele-index-not-a-number", "GT\tx/1\t1/1");

        // The version gate: the same record, no GT, at 4.0 and at 4.1.
        one("no-gt-v40", HEADER_V40, "GQ\t99");
        one("no-gt-v41", HEADER_V41, "GQ\t99");
    }

    static void two(final String label, final String genotypeBlock) {
        one(label, HEADER, genotypeBlock);
    }

    static void one(final String label, final String headerText, final String genotypeBlock) {
        final String line = "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\t" + genotypeBlock;
        try {
            final VCFCodec codec = new VCFCodec();
            codec.readActualHeader(
                    new LineIteratorImpl(new SynchronousLineReader(new StringReader(headerText))));
            final VariantContext vc = codec.decode(line);
            for (final Genotype genotype : vc.getGenotypes()) {
                emit(label, genotype);
            }
        } catch (final Throwable t) {
            System.out.printf("gterror\t%s\t%s\t%s%n", label, t.getClass().getName(),
                    oneLine(t.getMessage()));
        }
    }

    static void emit(final String label, final Genotype genotype) {
        final StringJoiner alleles = new StringJoiner(",");
        for (final Allele allele : genotype.getAlleles()) {
            alleles.add(allele.getDisplayString());
        }

        // Sorted: the extended attributes live in a map whose order is not a property of the file.
        final StringJoiner extended = new StringJoiner(";");
        for (final Map.Entry<String, Object> entry
                : new TreeMap<>(genotype.getExtendedAttributes()).entrySet()) {
            extended.add(entry.getKey() + "=" + entry.getValue());
        }

        System.out.printf("gt\t%s\t%s\t%s\t%b\t%s\t%s\t%s\t%s\t%s\t%s%n",
                label, genotype.getSampleName(), alleles, genotype.isPhased(),
                genotype.hasGQ() ? Integer.toString(genotype.getGQ()) : "none",
                genotype.hasDP() ? Integer.toString(genotype.getDP()) : "none",
                genotype.hasAD() ? ints(genotype.getAD()) : "none",
                genotype.hasPL() ? ints(genotype.getPL()) : "none",
                genotype.isFiltered() ? genotype.getFilters()
                        : (genotype.getFilters() == null ? "unfiltered" : "PASS"),
                extended);
    }

    static String ints(final int[] values) {
        final StringJoiner joined = new StringJoiner(",");
        for (final int value : values) {
            joined.add(Integer.toString(value));
        }
        return joined.toString();
    }

    static String oneLine(final String message) {
        return message == null ? "" : message.replace('\n', ' ').replace('\t', ' ');
    }
}
