/*
 * htsjdk.variant.variantcontext.VariantContextComparator, taken from the reference.
 *
 * One of the two things gatk-rs's MultiVariantDataSource needs, and therefore one of the two
 * gatk-rs G1.6 handed over: the multi-input walkers cannot be ported until this and
 * VCFUtils.smartMergeHeaders exist.
 *
 * Three behaviours this is built to catch, and the first is the one a port would get wrong by
 * being tidy:
 *
 *   - THE TWO CONSTRUCTORS REFUSE DIFFERENT THINGS. From a contig list the index is the position,
 *     so only a repeated name can go wrong. From header lines the index is carried by the line, so
 *     a repeated INDEX is a second, separately-worded error that a contig list cannot even express.
 *     A port that modelled one constructor and derived the other would lose that;
 *   - AN UNKNOWN CONTIG THROWS. htsjdk's own comment says "will throw NullPointerException --
 *     happily --", so it is a decision rather than a lapse: a caller sorting a variant whose contig
 *     is not in the dictionary gets an exception rather than an arbitrary order;
 *   - COMPARE RETURNS A SUBTRACTION, not a normalised -1/0/1. The magnitude is observable to any
 *     caller that inspects the value, so the dump records the number rather than its sign.
 *
 * Output:
 *
 *     ctor\t<label>\tok
 *     ctor\t<label>\tE:<class>:<message>
 *     cmp\t<label>\t<first>\t<second>\t<value>
 *     cmp\t<label>\t<first>\t<second>\tE:<class>
 *     compat\t<label>\t<true|false>
 *
 * Usage: VariantComparatorDump
 */

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFHeaderVersion;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public class VariantComparatorDump {

    public static void main(final String[] args) {
        System.out.println("# VariantComparatorDump: contig order, then start");

        // The constructors, including the four ways they refuse.
        ctorFromList("one", List.of("chr1"));
        ctorFromList("three", List.of("chr1", "chr2", "chr3"));
        ctorFromList("empty", List.of());
        ctorFromList("duplicate-name", List.of("chr1", "chr1"));
        ctorFromLines("lines-two", lines(new String[] {"chr1", "chr2"}, new int[] {0, 1}));
        // Indexes that are not 0..n-1, which a contig list can never produce.
        ctorFromLines("lines-sparse", lines(new String[] {"chr1", "chr2"}, new int[] {5, 9}));
        // Out of order: the index decides, not the position in the collection.
        ctorFromLines("lines-reversed", lines(new String[] {"chr1", "chr2"}, new int[] {1, 0}));
        ctorFromLines("lines-empty", new ArrayList<>());
        ctorFromLines("lines-duplicate-name", lines(new String[] {"chr1", "chr1"}, new int[] {0, 1}));
        ctorFromLines("lines-shared-index", lines(new String[] {"chr1", "chr2"}, new int[] {0, 0}));

        // compare, over a comparator built from a list.
        final VariantContextComparator listComparator =
                new VariantContextComparator(List.of("chr1", "chr2", "chr3", "chr4"));
        final Map<String, VariantContext> variants = new LinkedHashMap<>();
        variants.put("chr1:1", variant("chr1", 1));
        variants.put("chr1:100", variant("chr1", 100));
        variants.put("chr1:500", variant("chr1", 500));
        variants.put("chr2:1", variant("chr2", 1));
        variants.put("chr4:1", variant("chr4", 1));
        // A start far enough out to show the subtraction is a subtraction.
        variants.put("chr1:2000000000", variant("chr1", 2000000000));
        variants.put("chrX:1", variant("chrX", 1));

        for (final Map.Entry<String, VariantContext> first : variants.entrySet()) {
            for (final Map.Entry<String, VariantContext> second : variants.entrySet()) {
                compare("list", listComparator, first, second);
            }
        }

        // And over one built from lines whose indexes are sparse, so the difference between
        // "position in the collection" and "the line's own index" is visible in the value.
        final VariantContextComparator lineComparator = new VariantContextComparator(
                lines(new String[] {"chr1", "chr2"}, new int[] {5, 9}));
        for (final String a : new String[] {"chr1:1", "chr2:1"}) {
            for (final String b : new String[] {"chr1:1", "chr2:1"}) {
                compare("lines-sparse", lineComparator,
                        Map.entry(a, variants.get(a)), Map.entry(b, variants.get(b)));
            }
        }

        // isCompatible: the same name is not enough, the index must match too.
        compat("same", lineComparator, lines(new String[] {"chr1"}, new int[] {5}));
        compat("other-index", lineComparator, lines(new String[] {"chr1"}, new int[] {0}));
        compat("unknown-name", lineComparator, lines(new String[] {"chrX"}, new int[] {5}));
        compat("empty", lineComparator, new ArrayList<>());
    }

    static List<VCFContigHeaderLine> lines(final String[] names, final int[] indexes) {
        final List<VCFContigHeaderLine> result = new ArrayList<>();
        for (int i = 0; i < names.length; i++) {
            result.add(new VCFContigHeaderLine(
                    "<ID=" + names[i] + ",length=1000>", VCFHeaderVersion.VCF4_2, "contig",
                    indexes[i]));
        }
        return result;
    }

    static VariantContext variant(final String contig, final int start) {
        return new VariantContextBuilder("src", contig, start, start,
                List.of(Allele.create("A", true), Allele.create("C", false))).make();
    }

    static void ctorFromList(final String label, final List<String> contigs) {
        try {
            new VariantContextComparator(contigs);
            System.out.printf("ctor\t%s\tok%n", label);
        } catch (final Exception e) {
            System.out.printf("ctor\t%s\tE:%s:%s%n", label, e.getClass().getName(), e.getMessage());
        }
    }

    static void ctorFromLines(final String label, final List<VCFContigHeaderLine> lines) {
        try {
            new VariantContextComparator(lines);
            System.out.printf("ctor\t%s\tok%n", label);
        } catch (final Exception e) {
            System.out.printf("ctor\t%s\tE:%s:%s%n", label, e.getClass().getName(), e.getMessage());
        }
    }

    static void compare(final String label, final VariantContextComparator comparator,
                        final Map.Entry<String, VariantContext> first,
                        final Map.Entry<String, VariantContext> second) {
        String value;
        try {
            value = Integer.toString(comparator.compare(first.getValue(), second.getValue()));
        } catch (final Exception e) {
            value = "E:" + e.getClass().getName();
        }
        System.out.printf("cmp\t%s\t%s\t%s\t%s%n", label, first.getKey(), second.getKey(), value);
    }

    static void compat(final String label, final VariantContextComparator comparator,
                       final List<VCFContigHeaderLine> lines) {
        System.out.printf("compat\t%s\t%s%n", label, comparator.isCompatible(lines));
    }
}
