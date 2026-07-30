/*
 * The genotype type ladder, taken from the reference.
 *
 * Genotype.determineType decides what every downstream consumer means by "het" and "called", and
 * three of its answers contradict the javadoc sitting above them.
 *
 *   - HET is "two called alleles that are not equal". The javadoc says "at least one ref and at
 *     least one alt", but C/G, with no reference allele in it at all, comes back HET. That is why
 *     isHetNonRef exists as a separate question;
 *   - equality is Allele.equals, which is bases AND the reference flag. A genotype holding the ref
 *     A and a non-ref A is therefore HET, and prints as A/A;
 *   - ploidy is not two. A haploid call is HOM_REF or HOM_VAR, a triploid A/A/C is HET, and MIXED
 *     is anything with at least one call and at least one no-call whatever the ploidy.
 *
 * VariantContext.isMonomorphicInSamples is dumped alongside, because it is the guard that decides
 * whether the SampleList annotation emits anything at all, and it is not "every genotype is
 * hom-ref":
 *
 *     monomorphic = !isVariant() || (hasGenotypes() && getCalledChrCount(getReference()) == getCalledChrCount())
 *
 * A site with alternate alleles and NO genotypes is therefore NOT monomorphic, because the second
 * disjunct needs genotypes and the first needs there to be no alternate allele.
 *
 * And the sample ordering, because getGenotypesOrderedByName sorts with Collections.sort on the
 * names, which is String.compareTo, which is UTF-16 code units: uppercase before lowercase, digits
 * before letters, and "10" before "2".
 *
 * Output:
 *
 *     type\t<label>\t<GenotypeType>\t<called>\t<hom>\t<homRef>\t<homVar>\t<het>\t<hetNonRef>\t<mixed>\t<noCall>\t<available>
 *     mono\t<label>\t<monomorphic>\t<polymorphic>\t<calledChrCount(ref)>\t<calledChrCount()>
 *     order\t<label>\t<escaped sample names, comma-separated, in order>
 *
 * Usage: GenotypeTypeDump
 */

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.StringJoiner;

public class GenotypeTypeDump {

    static final Allele REF = Allele.create("A", true);
    /** The same bases as the reference, but not flagged reference: a different allele to equals. */
    static final Allele REF_BASES_AS_ALT = Allele.create("A", false);
    static final Allele ALT1 = Allele.create("C", false);
    static final Allele ALT2 = Allele.create("G", false);
    static final Allele NO_CALL = Allele.NO_CALL;

    public static void main(final String[] args) {
        System.out.println("# GenotypeTypeDump: determineType, isMonomorphicInSamples, and the name order");

        // Ploidy zero, which is the only way to reach UNAVAILABLE.
        type("no-alleles");
        // Diploid, every combination that changes the answer.
        type("hom-ref", REF, REF);
        type("hom-var", ALT1, ALT1);
        type("het", REF, ALT1);
        type("het-reversed", ALT1, REF);
        // Two different alternates and no reference: HET, which the javadoc denies.
        type("het-non-ref", ALT1, ALT2);
        // Same bases as the reference but not flagged reference: HET, and it prints as A/A.
        type("het-by-ref-flag", REF, REF_BASES_AS_ALT);
        type("hom-var-unflagged-ref", REF_BASES_AS_ALT, REF_BASES_AS_ALT);
        // No-calls, whole and partial.
        type("no-call", NO_CALL, NO_CALL);
        type("mixed-ref", REF, NO_CALL);
        type("mixed-alt", ALT1, NO_CALL);
        // Haploid: a single allele is HOM, not HET and not MIXED.
        type("haploid-ref", REF);
        type("haploid-alt", ALT1);
        type("haploid-no-call", NO_CALL);
        // Triploid and beyond.
        type("triploid-hom-ref", REF, REF, REF);
        type("triploid-het", REF, REF, ALT1);
        type("triploid-hom-var", ALT1, ALT1, ALT1);
        type("triploid-two-alts", ALT1, ALT1, ALT2);
        type("triploid-mixed", REF, ALT1, NO_CALL);
        type("tetraploid-hom-var", ALT1, ALT1, ALT1, ALT1);

        // isMonomorphicInSamples, which is the SampleList guard.
        mono("no-genotypes", build(List.of(REF, ALT1)));
        mono("ref-only-site", build(List.of(REF), gt("s1", REF, REF)));
        mono("ref-only-site-no-genotypes", build(List.of(REF)));
        mono("all-hom-ref", build(List.of(REF, ALT1), gt("s1", REF, REF), gt("s2", REF, REF)));
        mono("one-het", build(List.of(REF, ALT1), gt("s1", REF, REF), gt("s2", REF, ALT1)));
        mono("all-no-call", build(List.of(REF, ALT1), gt("s1", NO_CALL, NO_CALL)));
        mono("half-no-call-ref", build(List.of(REF, ALT1), gt("s1", REF, NO_CALL)));
        mono("half-no-call-alt", build(List.of(REF, ALT1), gt("s1", ALT1, NO_CALL)));
        // A filtered genotype contributes to neither count, so the alt it carries is invisible.
        mono("filtered-het", build(List.of(REF, ALT1),
                filtered(gt("s1", REF, ALT1), "LowGQ"), gt("s2", REF, REF)));
        mono("hom-var-only", build(List.of(REF, ALT1), gt("s1", ALT1, ALT1)));

        // The name order, which is String.compareTo and not any locale's collation.
        order("ascii", "b", "A", "a", "B");
        order("digits", "10", "2", "1", "20");
        order("underscore-vs-letters", "S_1", "SA", "Sa", "S1");
        order("shared-prefix", "sample", "sample1", "sampl", "sample10", "sample2");
        // Non-ASCII, where the order is by UTF-16 code unit and not by any alphabet.
        order("non-ascii", "é", "z", "Z", "É");
        // A supplementary character, whose surrogate pair sorts above every BMP character.
        order("supplementary", "😀", "￿", "a");
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

    static void type(final String label, final Allele... alleles) {
        try {
            final Genotype genotype = gt("s1", alleles);
            System.out.printf("type\t%s\t%s\t%b\t%b\t%b\t%b\t%b\t%b\t%b\t%b\t%b%n", label,
                    genotype.getType(), genotype.isCalled(), genotype.isHom(), genotype.isHomRef(),
                    genotype.isHomVar(), genotype.isHet(), hetNonRef(genotype), genotype.isMixed(),
                    genotype.isNoCall(), genotype.isAvailable());
        } catch (final Exception | AssertionError e) {
            System.out.printf("type\t%s\tE:%s:%s%n", label, e.getClass().getName(),
                    e.getMessage() == null ? "" : e.getMessage().replace('\n', ' '));
        }
    }

    /** isHetNonRef indexes allele 1 unconditionally, so it is only asked where that exists. */
    static boolean hetNonRef(final Genotype genotype) {
        return genotype.getPloidy() >= 2 && genotype.isHetNonRef();
    }

    static void mono(final String label, final VariantContext vc) {
        System.out.printf("mono\t%s\t%b\t%b\t%d\t%d%n", label, vc.isMonomorphicInSamples(),
                vc.isPolymorphicInSamples(), vc.getCalledChrCount(vc.getReference()),
                vc.getCalledChrCount());
    }

    static void order(final String label, final String... samples) {
        final List<Genotype> genotypes = new ArrayList<>();
        for (final String sample : samples) {
            genotypes.add(gt(sample, REF, ALT1));
        }
        final VariantContext vc = build(List.of(REF, ALT1), genotypes.toArray(new Genotype[0]));
        final StringJoiner joiner = new StringJoiner(",");
        for (final Genotype genotype : vc.getGenotypesOrderedByName()) {
            joiner.add(escape(genotype.getSampleName()));
        }
        System.out.printf("order\t%s\t%s%n", label, joiner);
    }

    /** Sample names deliberately contain non-ASCII, so they travel escaped. */
    static String escape(final String text) {
        final StringBuilder out = new StringBuilder();
        for (final char c : text.toCharArray()) {
            if (c < 0x20 || c > 0x7e) {
                out.append(String.format("\\u%04x", (int) c));
            } else if (c == '\\') {
                out.append("\\\\");
            } else {
                out.append(c);
            }
        }
        return out.toString();
    }
}
