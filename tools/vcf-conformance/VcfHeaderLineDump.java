/*
 * One ## line into a typed header line, taken from the reference.
 *
 * VcfHeaderParseDump measured the scanner: which pairs a line carries. This measures what they
 * mean, and it is where the refusals are. Four different exception classes come out of this layer,
 * and which one a file gets is a property of the field that is wrong:
 *
 *   Number=x     java.lang.NumberFormatException, uncaught: nothing wraps Integer.parseInt
 *   Number=-1    TribbleException$InvalidHeader, "Count < 0 for fixed size VCF header field"
 *   Number=0     java.lang.IllegalArgumentException, from validate(), for a non-Flag type
 *   Type=integer plain TribbleException, "not a valid type ... types are case-sensitive"
 *
 * Two of those are unchecked Java exceptions no catch in the codec touches, so a malformed Number
 * does not report "malformed header" at all. A port that funnelled every failure into one error
 * type would answer a different question from the reference, and the class is dumped for that
 * reason.
 *
 * Three asymmetries are probed because they make the same line valid or invalid depending on
 * context rather than on the line:
 *
 *   - INFO allows Type=Flag and FORMAT does not, so one line is a header line under one key and an
 *     IllegalArgumentException under the other;
 *   - a Flag with a non-zero count is silently rewritten to count 0 rather than refused, which is a
 *     value change and not an error;
 *   - Source and Version are read only from VCF 4.2, and under 4.1 they are not recommended tags
 *     either, so the tag-order check rejects the same line the 4.2 codec accepts.
 *
 * The rendered line is dumped rather than the fields, because the rendering is what a writer emits
 * and it carries the quoting rule as well as the values.
 *
 * Output:
 *
 *     hline\t<label>\t<class>\t<rendered line>
 *     hlineerror\t<label>\t<class>\t<message>
 *
 * Usage: VcfHeaderLineDump
 */

import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderVersion;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class VcfHeaderLineDump {

    public static void main(final String[] args) {
        System.out.println("# VcfHeaderLineDump: one ## line into a typed header line");

        // INFO, the ordinary shapes first.
        info("info-int", "<ID=DP,Number=1,Type=Integer,Description=\"Depth\">", VCFHeaderVersion.VCF4_2);
        info("info-a", "<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">", VCFHeaderVersion.VCF4_2);
        info("info-r", "<ID=AD,Number=R,Type=Integer,Description=\"Depths\">", VCFHeaderVersion.VCF4_2);
        info("info-g", "<ID=PL,Number=G,Type=Integer,Description=\"Likelihoods\">", VCFHeaderVersion.VCF4_2);
        info("info-unbounded", "<ID=X,Number=.,Type=String,Description=\"Any\">", VCFHeaderVersion.VCF4_2);
        info("info-flag", "<ID=DB,Number=0,Type=Flag,Description=\"In dbSNP\">", VCFHeaderVersion.VCF4_2);
        // A Flag whose count is not zero, which is rewritten rather than refused.
        info("info-flag-count-2", "<ID=DB,Number=2,Type=Flag,Description=\"In dbSNP\">", VCFHeaderVersion.VCF4_2);
        // Source and Version, admitted from 4.2 and rejected before it.
        info("info-source-42", "<ID=DP,Number=1,Type=Integer,Description=\"d\",Source=\"s\",Version=\"3\">", VCFHeaderVersion.VCF4_2);
        info("info-source-41", "<ID=DP,Number=1,Type=Integer,Description=\"d\",Source=\"s\">", VCFHeaderVersion.VCF4_1);
        // No Description at all, which is a default rather than a refusal.
        info("info-no-description", "<ID=DP,Number=1,Type=Integer>", VCFHeaderVersion.VCF4_2);
        // The four failures.
        info("info-number-not-a-number", "<ID=DP,Number=x,Type=Integer,Description=\"d\">", VCFHeaderVersion.VCF4_2);
        info("info-number-negative", "<ID=DP,Number=-1,Type=Integer,Description=\"d\">", VCFHeaderVersion.VCF4_2);
        info("info-number-zero", "<ID=DP,Number=0,Type=Integer,Description=\"d\">", VCFHeaderVersion.VCF4_2);
        info("info-type-lowercase", "<ID=DP,Number=1,Type=integer,Description=\"d\">", VCFHeaderVersion.VCF4_2);
        // No Number tag at all: the codec dereferences it without a null check.
        info("info-no-number", "<ID=DP>", VCFHeaderVersion.VCF4_2);
        // An ID carrying characters the validation refuses. They survive the scanner only because
        // they are quoted; unquoted, the scanner drops them first.
        info("info-id-angle", "<ID=\"a<b\",Number=1,Type=Integer,Description=\"d\">", VCFHeaderVersion.VCF4_2);
        info("info-id-equals", "<ID=\"a=b\",Number=1,Type=Integer,Description=\"d\">", VCFHeaderVersion.VCF4_2);

        // FORMAT: the same lines, with the Flag asymmetry.
        format("format-int", "<ID=GQ,Number=1,Type=Integer,Description=\"Quality\">", VCFHeaderVersion.VCF4_2);
        format("format-flag", "<ID=DB,Number=0,Type=Flag,Description=\"In dbSNP\">", VCFHeaderVersion.VCF4_2);

        // FILTER, whose expected tags are ID and Description and which admits no recommended tags.
        filter("filter-plain", "<ID=LowQual,Description=\"Low quality\">", VCFHeaderVersion.VCF4_2);
        filter("filter-no-description", "<ID=LowQual>", VCFHeaderVersion.VCF4_2);
        filter("filter-wrong-order", "<Description=\"d\",ID=LowQual>", VCFHeaderVersion.VCF4_2);
        filter("filter-no-id", "<Description=\"d\">", VCFHeaderVersion.VCF4_2);
        filter("filter-extra-tag", "<ID=LowQual,Description=\"d\",Extra=1>", VCFHeaderVersion.VCF4_2);

        // contig, whose expected tag order is null, so its fields may come in any order.
        contig("contig-plain", "<ID=chr1,length=1000>", VCFHeaderVersion.VCF4_2, 0);
        contig("contig-reordered", "<length=1000,ID=chr1>", VCFHeaderVersion.VCF4_2, 3);
        contig("contig-extra", "<ID=chr1,length=1000,assembly=b37,md5=abc>", VCFHeaderVersion.VCF4_2, 7);
        contig("contig-no-id", "<length=1000>", VCFHeaderVersion.VCF4_2, 0);
        contig("contig-negative-index", "<ID=chr1,length=1000>", VCFHeaderVersion.VCF4_2, -1);
    }

    static void info(final String label, final String value, final VCFHeaderVersion version) {
        emit(label, () -> new VCFInfoHeaderLine(value, version));
    }

    static void format(final String label, final String value, final VCFHeaderVersion version) {
        emit(label, () -> new VCFFormatHeaderLine(value, version));
    }

    static void filter(final String label, final String value, final VCFHeaderVersion version) {
        emit(label, () -> new VCFFilterHeaderLine(value, version));
    }

    static void contig(final String label, final String value, final VCFHeaderVersion version,
                       final int index) {
        emit(label, () -> new VCFContigHeaderLine(value, version, "contig", index));
    }

    interface Build {
        VCFHeaderLine get();
    }

    static void emit(final String label, final Build build) {
        try {
            final VCFHeaderLine line = build.get();
            // toString is the sort key and the written form at once, so it is what a divergence
            // would reach the output through.
            System.out.printf("hline\t%s\t%s\t%s%n", label, line.getClass().getName(),
                    oneLine(line.toString()));
        } catch (final Throwable t) {
            System.out.printf("hlineerror\t%s\t%s\t%s%n", label, t.getClass().getName(),
                    oneLine(t.getMessage()));
        }
    }

    static String oneLine(final String message) {
        return message == null ? "" : message.replace('\n', ' ').replace('\t', ' ');
    }
}
