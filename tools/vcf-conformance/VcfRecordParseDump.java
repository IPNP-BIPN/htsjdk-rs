/*
 * Reading one VCF data line, taken from the reference: the eight site columns.
 *
 * VcfRecordDump measures writing a record. This measures reading one, which is what GATK's
 * VariantWalker needs and what nothing here has done before. Genotypes are left out: htsjdk defers
 * them through LazyGenotypesContext and the deferral is itself observable, so it gets its own
 * slice; this stops at the site columns and reports the genotype block as the splitter left it.
 *
 * Five behaviours decide what a record means and none of them is in the VCF specification.
 *
 *   - the splitter is not String.split. ParsingUtils.split fills a fixed-size array, skips a
 *     delimiter at position 0 rather than producing an empty first token, and condenses every
 *     trailing column into the last slot. So a line starting with a tab loses its first character,
 *     and every genotype column after the first arrives joined back together in parts[8];
 *   - END overrides the stop and nothing checks it, so a record can end before it starts;
 *   - a declared Flag written as KEY=0 is dropped from the attributes entirely rather than stored
 *     as false, and a bare key whose header type is not Flag becomes the string "." rather than a
 *     flag. Both mean the header changes what a record contains;
 *   - REF is upper-cased silently, so a lower-case REF is a rewrite and not an error;
 *   - an ALT of "." is checked and then not added, so the record has one allele rather than two.
 *
 * The refusals are dumped with their classes and messages, because parseQual throws a
 * NumberFormatException that nothing catches while the neighbouring failures are TribbleExceptions
 * carrying a line number.
 *
 * Every record is read against the same header unless the label says otherwise, so a difference in
 * a row is a difference in the line and not in the setup.
 *
 * Output:
 *
 *     rec\t<label>\t<chr>\t<start>\t<stop>\t<id>\t<alleles>\t<qual>\t<filters>\t<attributes>\t<genotype block>
 *     recerror\t<label>\t<class>\t<message>
 *     recnull\t<label>
 *
 * Usage: VcfRecordParseDump
 */

import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.SynchronousLineReader;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

import java.io.StringReader;
import java.util.Map;
import java.util.StringJoiner;
import java.util.TreeMap;

public class VcfRecordParseDump {

    /** A header declaring the types the INFO parsing consults, with two samples. */
    static final String HEADER = "##fileformat=VCFv4.2\n"
            + "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"In dbSNP\">\n"
            + "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
            + "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Frequency\">\n"
            + "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n"
            + "##FILTER=<ID=LowQual,Description=\"Low quality\">\n"
            + "##contig=<ID=chr1,length=100000>\n"
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA1\tNA2\n";

    /** The same header with no samples at all, which changes the column count check. */
    static final String SITES_ONLY_HEADER = "##fileformat=VCFv4.2\n"
            + "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"In dbSNP\">\n"
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    public static void main(final String[] args) {
        System.out.println("# VcfRecordParseDump: reading one VCF data line, site columns only");

        // The ordinary shapes.
        record("snp", "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\tGT\t0/1\t1/1");
        record("multiallelic", "chr1\t100\trs1\tA\tT,C\t50\tPASS\tDP=10;AF=0.5,0.25\tGT\t0/1\t1/2");
        record("deletion", "chr1\t100\t.\tACGT\tA\t50\t.\tDP=10\tGT\t0/1\t0/0");
        record("no-qual", "chr1\t100\t.\tA\tT\t.\tPASS\tDP=10\tGT\t0/1\t1/1");
        // The VCF 3 encoding of a missing quality, which is silently treated as missing.
        record("qual-minus-one", "chr1\t100\t.\tA\tT\t-1\tPASS\tDP=10\tGT\t0/1\t1/1");
        record("qual-minus-one-point-zero", "chr1\t100\t.\tA\tT\t-1.0\tPASS\tDP=10\tGT\t0/1\t1/1");
        // Filters: unfiltered, passed, and two of them.
        record("unfiltered", "chr1\t100\t.\tA\tT\t50\t.\tDP=10\tGT\t0/1\t1/1");
        record("two-filters", "chr1\t100\t.\tA\tT\t50\tLowQual;q10\tDP=10\tGT\t0/1\t1/1");
        // END, which overrides the stop, including backwards.
        record("end-key", "chr1\t100\t.\tA\t<DEL>\t50\tPASS\tEND=200\tGT\t0/1\t1/1");
        record("end-before-start", "chr1\t100\t.\tA\t<DEL>\t50\tPASS\tEND=50\tGT\t0/1\t1/1");
        // The INFO cases the header decides.
        record("flag-bare", "chr1\t100\t.\tA\tT\t50\tPASS\tDB\tGT\t0/1\t1/1");
        record("flag-equals-zero", "chr1\t100\t.\tA\tT\t50\tPASS\tDB=0\tGT\t0/1\t1/1");
        record("flag-equals-one", "chr1\t100\t.\tA\tT\t50\tPASS\tDB=1\tGT\t0/1\t1/1");
        // A bare key the header types as Integer, which becomes the string "." rather than a flag.
        record("bare-non-flag", "chr1\t100\t.\tA\tT\t50\tPASS\tDP\tGT\t0/1\t1/1");
        // A key the header does not declare at all, bare and with =0.
        record("undeclared-bare", "chr1\t100\t.\tA\tT\t50\tPASS\tXX\tGT\t0/1\t1/1");
        record("undeclared-zero", "chr1\t100\t.\tA\tT\t50\tPASS\tXX=0\tGT\t0/1\t1/1");
        // key= with nothing after it, which is the missing value rather than an empty string.
        record("empty-value", "chr1\t100\t.\tA\tT\t50\tPASS\tDP=\tGT\t0/1\t1/1");
        record("empty-info", "chr1\t100\t.\tA\tT\t50\tPASS\t.\tGT\t0/1\t1/1");
        // A repeated key, which a HashMap resolves to the last value.
        record("repeated-key", "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10;DP=20\tGT\t0/1\t1/1");
        // Silent rewrites.
        record("lowercase-ref", "chr1\t100\t.\ta\tT\t50\tPASS\tDP=10\tGT\t0/1\t1/1");
        record("alt-missing", "chr1\t100\t.\tA\t.\t50\tPASS\tDP=10\tGT\t0/1\t1/1");
        record("alt-star", "chr1\t100\t.\tA\tT,*\t50\tPASS\tDP=10\tGT\t0/1\t1/1");
        record("symbolic-alt", "chr1\t100\t.\tA\t<NON_REF>\t50\tPASS\tDP=10\tGT\t0/1\t1/1");
        // A line starting with a tab, where the splitter steps over the first character.
        record("leading-tab", "\tchr1\t100\t.\tA\tT\t50\tPASS\tDP=10\tGT\t0/1\t1/1");
        // A header line handed to the record decoder, which is skipped rather than refused.
        record("header-line", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");

        // The refusals.
        record("pos-not-a-number", "chr1\tx\t.\tA\tT\t50\tPASS\tDP=10\tGT\t0/1\t1/1");
        record("empty-id", "chr1\t100\t\tA\tT\t50\tPASS\tDP=10\tGT\t0/1\t1/1");
        record("qual-not-a-number", "chr1\t100\t.\tA\tT\tx\tPASS\tDP=10\tGT\t0/1\t1/1");
        record("filter-zero", "chr1\t100\t.\tA\tT\t50\t0\tDP=10\tGT\t0/1\t1/1");
        record("info-with-space", "chr1\t100\t.\tA\tT\t50\tPASS\tDP=1 0\tGT\t0/1\t1/1");
        record("info-empty-string", "chr1\t100\t.\tA\tT\t50\tPASS\t\tGT\t0/1\t1/1");
        record("end-not-a-number", "chr1\t100\t.\tA\t<DEL>\t50\tPASS\tEND=x\tGT\t0/1\t1/1");
        record("ref-missing", "chr1\t100\t.\t.\tT\t50\tPASS\tDP=10\tGT\t0/1\t1/1");
        record("ref-symbolic", "chr1\t100\t.\t<DEL>\tT\t50\tPASS\tDP=10\tGT\t0/1\t1/1");
        record("ref-bad-base", "chr1\t100\t.\tQ\tT\t50\tPASS\tDP=10\tGT\t0/1\t1/1");
        record("alt-breakend", "chr1\t100\t.\tA\tA[chr2:200[\t50\tPASS\tDP=10\tGT\t0/1\t1/1");
        record("vcf3-deletion", "chr1\t100\t.\tA\tD2\t50\tPASS\tDP=10\tGT\t0/1\t1/1");
        record("too-few-columns", "chr1\t100\t.\tA\tT\t50\tPASS");

        // The same lines against a header with no samples, where nine columns is one too many.
        sitesOnly("sites-only-eight", "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10");
        sitesOnly("sites-only-nine", "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\tGT");
        // Against that header DB is still a Flag, but DP is not declared at all.
        sitesOnly("sites-only-undeclared-dp", "chr1\t100\t.\tA\tT\t50\tPASS\tDP");
    }

    static void record(final String label, final String line) {
        emit(label, HEADER, line);
    }

    static void sitesOnly(final String label, final String line) {
        emit(label, SITES_ONLY_HEADER, line);
    }

    static void emit(final String label, final String headerText, final String line) {
        try {
            final VCFCodec codec = new VCFCodec();
            codec.readActualHeader(
                    new LineIteratorImpl(new SynchronousLineReader(new StringReader(headerText))));
            final VariantContext vc = codec.decode(line);
            if (vc == null) {
                System.out.printf("recnull\t%s%n", label);
                return;
            }

            final StringJoiner alleles = new StringJoiner(",");
            for (final Allele allele : vc.getAlleles()) {
                alleles.add(allele.getDisplayString() + (allele.isReference() ? "*" : ""));
            }

            // Sorted, because the attributes live in a HashMap upstream and its iteration order is
            // not a property of the file. Sorting is the honest comparison; the order is not one.
            final StringJoiner attributes = new StringJoiner(";");
            for (final Map.Entry<String, Object> entry
                    : new TreeMap<>(vc.getAttributes()).entrySet()) {
                attributes.add(entry.getKey() + "=" + render(entry.getValue()));
            }

            // Three states, not two: never applied, applied and passed, applied and failed. The
            // first two are different files and most ports collapse them.
            final String filters;
            if (!vc.filtersWereApplied()) {
                filters = "unfiltered";
            } else if (vc.getFilters().isEmpty()) {
                filters = "PASS";
            } else {
                filters = String.join(",", new java.util.TreeSet<>(vc.getFilters()));
            }

            System.out.printf("rec\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s%n",
                    label, vc.getContig(), vc.getStart(), vc.getEnd(), vc.getID(), alleles,
                    vc.hasLog10PError() ? Double.toString(vc.getLog10PError()) : "none",
                    filters, attributes, sampleNames(codec));
        } catch (final Throwable t) {
            System.out.printf("recerror\t%s\t%s\t%s%n", label, t.getClass().getName(),
                    oneLine(t.getMessage()));
        }
    }

    /**
     * The sample names from the header, not from the record.
     *
     * Asking the record forces LazyGenotypesContext to decode, which drags the genotype layer into
     * a suite that is meant to stop at the site columns: an ALT of "." made GT 0/1 refer to an
     * allele that does not exist and the row became an InternalCodecException instead of a record.
     * That failure is real and belongs to the genotype slice, so it is measured there.
     */
    static String sampleNames(final VCFCodec codec) {
        return String.join(",", codec.getHeader().getGenotypeSamples());
    }

    static String render(final Object value) {
        if (value instanceof java.util.List) {
            final StringJoiner joined = new StringJoiner(",");
            for (final Object item : (java.util.List<?>) value) {
                joined.add(String.valueOf(item));
            }
            return "[" + joined + "]";
        }
        return String.valueOf(value);
    }

    static String oneLine(final String message) {
        return message == null ? "" : message.replace('\n', ' ').replace('\t', ' ');
    }
}
