/*
 * htsjdk.variant.vcf.VCFUtils.smartMergeHeaders, taken from the reference.
 *
 * The second of the two things gatk-rs's MultiVariantDataSource needs, alongside
 * VariantContextComparator. With both, the multi-input walkers gatk-rs G1.6 handed over become
 * portable.
 *
 * Eighty lines of Java, and most of what matters in them is not the merge:
 *
 *   - THE OUTPUT ORDER IS FIRST-SEEN ACROSS SOURCES AND SORTED WITHIN ONE. The loop reads
 *     getMetaDataInSortedOrder(), which is a TreeSet, so a source's lines do not arrive in file
 *     order. The method's own comment says this is what keeps contig lines from being scrambled,
 *     so the ordering is the point rather than a detail;
 *   - FIRST ONE WINS, EXCEPT THREE WAYS. A Number difference promotes the STORED line to
 *     unbounded, in place. Integer against Float keeps the Float in both directions — and one of
 *     the two arms does it with a put whose value is already what the map holds, so it is a no-op
 *     that emits a message naming the other line. A type collision that is neither throws;
 *   - THE VERSION POLICY THROWS A DIFFERENT CLASS. IllegalArgumentException, where the rest of the
 *     method throws IllegalStateException, and only when a 4.3 header meets any other version.
 *
 * Output:
 *
 *     merged\t<label>\t<count>\t<rendered lines, joined by |>
 *     warn\t<label>\t<message>
 *     err\t<label>\t<class>:<message>
 *
 * Usage: SmartMergeDump
 */

import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFHeaderVersion;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFUtils;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

public class SmartMergeDump {

    public static void main(final String[] args) {
        System.out.println("# SmartMergeDump: what survives a header merge, and in what order");

        final VCFInfoHeaderLine dpInt =
                new VCFInfoHeaderLine("DP", 1, VCFHeaderLineType.Integer, "depth");
        final VCFInfoHeaderLine dpIntOtherDesc =
                new VCFInfoHeaderLine("DP", 1, VCFHeaderLineType.Integer, "read depth");
        final VCFInfoHeaderLine afOne =
                new VCFInfoHeaderLine("AF", 1, VCFHeaderLineType.Float, "af");
        final VCFInfoHeaderLine afTwo =
                new VCFInfoHeaderLine("AF", 2, VCFHeaderLineType.Float, "af");
        final VCFInfoHeaderLine xInt =
                new VCFInfoHeaderLine("X", 1, VCFHeaderLineType.Integer, "x");
        final VCFInfoHeaderLine xFloat =
                new VCFInfoHeaderLine("X", 1, VCFHeaderLineType.Float, "x");
        final VCFInfoHeaderLine xString =
                new VCFInfoHeaderLine("X", 1, VCFHeaderLineType.String, "x");
        final VCFFormatHeaderLine gq =
                new VCFFormatHeaderLine("GQ", 1, VCFHeaderLineType.Integer, "gq");
        final VCFFilterHeaderLine lowQual = new VCFFilterHeaderLine("LowQual", "low quality");
        final VCFFilterHeaderLine lowQualOther = new VCFFilterHeaderLine("LowQual", "poor");

        // The base cases.
        merge("identical", headers(lines(dpInt), lines(dpInt)));
        merge("disjoint", headers(lines(dpInt), lines(gq)));

        // Number differs, same type: the stored line is promoted to `.` in place.
        merge("number-differs", headers(lines(afOne), lines(afTwo)));
        merge("number-differs-reversed", headers(lines(afTwo), lines(afOne)));

        // Integer against Float, both orders. Both keep the Float; one of the Java's arms is a
        // no-op, and the message names the stored line either way.
        merge("int-then-float", headers(lines(xInt), lines(xFloat)));
        merge("float-then-int", headers(lines(xFloat), lines(xInt)));

        // A collision that promotes to nothing.
        merge("int-then-string", headers(lines(xInt), lines(xString)));

        // Descriptions differ: a warning, and the stored one wins.
        merge("description-differs", headers(lines(dpInt), lines(dpIntOtherDesc)));
        // The same on a FILTER, whose own arm in the Java is unreachable because the ID is part of
        // the key. This is the branch that catches it instead.
        merge("filter-description-differs", headers(lines(lowQual), lines(lowQualOther)));

        // An unstructured line under one key, twice, with different values: the generic branch.
        merge("unstructured-conflict", headers(
                lines(new VCFHeaderLine("source", "one")),
                lines(new VCFHeaderLine("source", "two"))));

        // ORDER. The first source's lines are read sorted, so B and A come out A then B.
        final VCFInfoHeaderLine a = new VCFInfoHeaderLine("A", 1, VCFHeaderLineType.Integer, "a");
        final VCFInfoHeaderLine b = new VCFInfoHeaderLine("B", 1, VCFHeaderLineType.Integer, "b");
        final VCFInfoHeaderLine c = new VCFInfoHeaderLine("C", 1, VCFHeaderLineType.Integer, "c");
        merge("order-sorted-within-source", headers(lines(b, a), lines(c)));
        // Contigs, which is what the method's comment says the ordering exists to protect.
        merge("contigs", headers(
                lines(contig("chr2", 1), contig("chr1", 0)),
                lines(contig("chr3", 2))));

        // The version policy, which throws a different class from everything else.
        merge("version-43-and-42", List.of(withVersion(VCFHeaderVersion.VCF4_3, dpInt),
                withVersion(VCFHeaderVersion.VCF4_2, dpInt)));
        merge("version-43-and-43", List.of(withVersion(VCFHeaderVersion.VCF4_3, dpInt),
                withVersion(VCFHeaderVersion.VCF4_3, dpInt)));
        merge("version-42-and-41", List.of(withVersion(VCFHeaderVersion.VCF4_2, dpInt),
                withVersion(VCFHeaderVersion.VCF4_1, dpInt)));
        // An empty header never reaches the policy, because the check sits inside the per-line
        // loop rather than the per-header one.
        merge("version-43-and-empty", List.of(withVersion(VCFHeaderVersion.VCF4_3, dpInt),
                withVersion(VCFHeaderVersion.VCF4_2)));
    }

    static VCFContigHeaderLine contig(final String id, final int index) {
        return new VCFContigHeaderLine(
                "<ID=" + id + ",length=1000>", VCFHeaderVersion.VCF4_2, "contig", index);
    }

    static Set<VCFHeaderLine> lines(final VCFHeaderLine... lines) {
        return new LinkedHashSet<>(List.of(lines));
    }

    static List<VCFHeader> headers(final Set<VCFHeaderLine> first, final Set<VCFHeaderLine> second) {
        return List.of(new VCFHeader(first), new VCFHeader(second));
    }

    /** A header whose `fileformat` version is set, which is what the policy reads. */
    static VCFHeader withVersion(final VCFHeaderVersion version, final VCFHeaderLine... lines) {
        final Set<VCFHeaderLine> set = new LinkedHashSet<>(List.of(lines));
        set.add(new VCFHeaderLine(version.getFormatString(), version.getVersionString()));
        return new VCFHeader(set);
    }

    static void merge(final String label, final List<VCFHeader> headers) {
        // The warnings go to a logger, so they are captured off stderr rather than returned.
        final PrintStream realErr = System.err;
        final ByteArrayOutputStream captured = new ByteArrayOutputStream();
        System.setErr(new PrintStream(captured));
        Set<VCFHeaderLine> merged = null;
        String error = null;
        try {
            merged = VCFUtils.smartMergeHeaders(headers, true);
        } catch (final Exception e) {
            error = e.getClass().getName() + ":" + e.getMessage();
        } finally {
            System.setErr(realErr);
        }

        if (error != null) {
            System.out.printf("err\t%s\t%s%n", label, error);
            return;
        }
        final List<String> rendered = new ArrayList<>();
        for (final VCFHeaderLine line : merged) {
            rendered.add(line.toString());
        }
        System.out.printf("merged\t%s\t%d\t%s%n", label, rendered.size(),
                rendered.stream().collect(Collectors.joining("|")));
    }
}
