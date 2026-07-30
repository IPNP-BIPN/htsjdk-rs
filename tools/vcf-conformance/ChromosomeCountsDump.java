/*
 * AC, AF and AN, taken from the reference.
 *
 * GATK's ChromosomeCounts annotation is two lines that delegate straight to
 * VariantContextUtils.calculateChromosomeCounts, so this is where the three commonest INFO fields in
 * any VCF are actually decided. Four behaviours are measured and none follows from the field names.
 *
 *   - AC and AF change TYPE with the number of alternate alleles: one alt gives a scalar, two or
 *     more give a list. A consumer that fetched AC gets an Integer in one case and an ArrayList in
 *     the other, even though both render the same way;
 *   - a genotype whose FT is set contributes NOTHING, and so does a no-call allele, for different
 *     reasons. AN is therefore not the ploidy times the sample count;
 *   - AF is not AC / AN. The numerator is the founders' count for that allele and the denominator
 *     is the founders' total called count, while the AC reported alongside it is the WHOLE
 *     cohort's. With no founders the two coincide, which is why the difference only shows on a
 *     pedigree;
 *   - the function REMOVES keys as well as adding them. A site with no alternate allele has its AC
 *     and AF deleted; a site where nobody is called has all three deleted, but only when
 *     removeStaleValues is on.
 *
 * AF is dumped as raw bits, because it is a division and a decimal rendering would hide a
 * divergence in the last place.
 *
 * Output:
 *
 *     an\t<label>\t<AN, or absent>
 *     ac\t<label>\t<scalar|list, or absent>
 *     af\t<label>\t<raw bits, comma-separated, or absent>
 *     type\t<label>\t<the Java class of the AC value, or absent>
 *     called\t<label>\t<getCalledChrCount()>\t<getCalledChrCount(founders)>
 *
 * Usage: ChromosomeCountsDump
 */

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.VariantContextUtils;
import htsjdk.variant.vcf.VCFConstants;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;

public class ChromosomeCountsDump {

    static final Allele REF = Allele.create("A", true);
    static final Allele ALT1 = Allele.create("C", false);
    static final Allele ALT2 = Allele.create("G", false);
    static final Allele NO_CALL = Allele.NO_CALL;

    public static void main(final String[] args) {
        System.out.println("# ChromosomeCountsDump: AC, AF and AN, from the reference");

        // No genotypes at all: the annotation returns an empty map before reaching here, but the
        // utility is still defined and this is what it does.
        emit("no-genotypes", build(List.of(REF, ALT1)), true, Set.of());

        // One alternate allele: AC and AF are scalars.
        emit("one-alt-het", build(List.of(REF, ALT1),
                gt("s1", REF, ALT1)), true, Set.of());
        emit("one-alt-hom-var", build(List.of(REF, ALT1),
                gt("s1", ALT1, ALT1)), true, Set.of());
        emit("one-alt-hom-ref", build(List.of(REF, ALT1),
                gt("s1", REF, REF)), true, Set.of());

        // Two alternate alleles: AC and AF become lists.
        emit("two-alts", build(List.of(REF, ALT1, ALT2),
                gt("s1", REF, ALT1), gt("s2", ALT1, ALT2)), true, Set.of());

        // No alternate allele: AC and AF are removed, AN stays.
        emit("ref-only", build(List.of(REF), gt("s1", REF, REF)), true, Set.of());

        // Everybody a no-call, with and without stale removal: the same site gives three keys or
        // none depending on a boolean the annotation hard-codes to true.
        emit("all-no-call-remove", build(List.of(REF, ALT1),
                gt("s1", NO_CALL, NO_CALL)), true, Set.of());
        emit("all-no-call-keep", build(List.of(REF, ALT1),
                gt("s1", NO_CALL, NO_CALL)), false, Set.of());

        // A partial no-call: one chromosome counts and the other does not.
        emit("half-no-call", build(List.of(REF, ALT1),
                gt("s1", ALT1, NO_CALL)), true, Set.of());

        // A filtered genotype, which contributes nothing at all.
        emit("filtered-genotype", build(List.of(REF, ALT1),
                filtered(gt("s1", ALT1, ALT1), "LowGQ"), gt("s2", REF, ALT1)), true, Set.of());
        // A genotype whose FT is the empty string, which is NOT filtered.
        emit("empty-filter", build(List.of(REF, ALT1),
                filtered(gt("s1", ALT1, ALT1), ""), gt("s2", REF, ALT1)), true, Set.of());
        // Every genotype filtered, so AN is zero even though nobody is a no-call.
        emit("all-filtered", build(List.of(REF, ALT1),
                filtered(gt("s1", ALT1, ALT1), "LowGQ")), true, Set.of());

        // Founders: AC is the whole cohort's and AF is the founders'.
        final VariantContext pedigree = build(List.of(REF, ALT1),
                gt("founder", REF, ALT1), gt("child1", ALT1, ALT1), gt("child2", ALT1, ALT1));
        emit("pedigree-no-founders", pedigree, true, Set.of());
        emit("pedigree-one-founder", pedigree, true, Set.of("founder"));
        emit("pedigree-two-founders", pedigree, true, Set.of("founder", "child1"));
        // A founder set naming a sample that is not there: the denominator is zero and the
        // division is by zero rather than an error.
        emit("pedigree-absent-founder", pedigree, true, Set.of("nobody"));

        // Mixed ploidy, and a haploid call.
        emit("haploid", build(List.of(REF, ALT1), gt("s1", ALT1)), true, Set.of());
        emit("mixed-ploidy", build(List.of(REF, ALT1),
                gt("s1", ALT1), gt("s2", REF, ALT1), gt("s3", REF, ALT1, ALT1)), true, Set.of());

        // Many samples, so AF is a division that does not land on a representable decimal.
        final List<Genotype> many = new ArrayList<>();
        for (int i = 0; i < 3; i++) {
            many.add(gt("m" + i, REF, ALT1));
        }
        emit("three-hets", build(List.of(REF, ALT1), many.toArray(new Genotype[0])), true, Set.of());
    }

    static Genotype gt(final String sample, final Allele... alleles) {
        return new GenotypeBuilder(sample, Arrays.asList(alleles)).make();
    }

    static Genotype filtered(final Genotype genotype, final String filter) {
        return new GenotypeBuilder(genotype).filters(filter).make();
    }

    static VariantContext build(final List<Allele> alleles, final Genotype... genotypes) {
        final VariantContextBuilder builder = new VariantContextBuilder()
                .chr("chr1").start(100).stop(100).alleles(alleles);
        if (genotypes.length > 0) {
            builder.genotypes(Arrays.asList(genotypes));
        }
        return builder.make();
    }

    static void emit(final String label, final VariantContext vc, final boolean removeStale,
                     final Set<String> founders) {
        try {
            final Map<String, Object> attributes = VariantContextUtils.calculateChromosomeCounts(
                    vc, new LinkedHashMap<>(), removeStale, new HashSet<>(founders));

            System.out.printf("an\t%s\t%s%n", label, show(attributes.get(VCFConstants.ALLELE_NUMBER_KEY)));
            System.out.printf("ac\t%s\t%s%n", label, show(attributes.get(VCFConstants.ALLELE_COUNT_KEY)));
            System.out.printf("af\t%s\t%s%n", label, bits(attributes.get(VCFConstants.ALLELE_FREQUENCY_KEY)));
            final Object ac = attributes.get(VCFConstants.ALLELE_COUNT_KEY);
            System.out.printf("type\t%s\t%s%n", label,
                    ac == null ? "absent" : ac.getClass().getName());
            System.out.printf("called\t%s\t%d\t%d%n", label,
                    vc.getCalledChrCount(), vc.getCalledChrCount(new HashSet<>(founders)));
        } catch (final Exception | AssertionError e) {
            System.out.printf("error\t%s\t%s\t%s%n", label, e.getClass().getName(),
                    e.getMessage() == null ? "" : e.getMessage().replace('\n', ' '));
        }
    }

    static String show(final Object value) {
        if (value == null) {
            return "absent";
        }
        if (value instanceof List) {
            final StringJoiner joiner = new StringJoiner(",");
            for (final Object element : (List<?>) value) {
                joiner.add(String.valueOf(element));
            }
            return joiner.toString();
        }
        return String.valueOf(value);
    }

    /** AF is a division, so it travels as raw bits rather than as a rendering. */
    static String bits(final Object value) {
        if (value == null) {
            return "absent";
        }
        final StringJoiner joiner = new StringJoiner(",");
        if (value instanceof List) {
            for (final Object element : (List<?>) value) {
                joiner.add(Long.toString(Double.doubleToRawLongBits((Double) element)));
            }
        } else {
            joiner.add(Long.toString(Double.doubleToRawLongBits((Double) value)));
        }
        return joiner.toString();
    }
}
