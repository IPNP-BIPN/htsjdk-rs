/*
 * Reading a VCF header, taken from the reference: the structured-line scanner and the frame the
 * codec establishes before any line is given a type.
 *
 * Everything so far in this repository writes VCF; nothing reads one. GATK's VariantWalker cannot
 * exist without a reader, and the reader's first job is this. Two layers are probed.
 *
 * VCFHeaderLineTranslator's VCF4Parser is a hand-written scanner whose switch falls through:
 *
 *     case ('<') : if (index == 0) break;   // no break when index != 0: falls into '>'
 *     case ('>') : if (index == valueLine.length()-1) ret.put(key,builder.toString().trim()); break;
 *     case ('=') : key = builder.toString().trim(); builder = new StringBuilder(); break;
 *     case (',') : ret.put(key,builder.toString().trim()); builder = new StringBuilder(); break;
 *     default: builder.append(c);
 *
 * so an unquoted '<' anywhere but position 0 is dropped, an unquoted '>' anywhere is dropped, and a
 * line that does not end in '>' loses its last field entirely and silently. None of that is in the
 * VCF specification, and a port written from the specification would keep those characters and
 * store that field.
 *
 * Quotes are a toggle rather than a delimiter, because the c == '"' test comes before the inQuote
 * test: a quote in the middle of an unquoted value opens a quoted region. Inside one, \" gives ",
 * \\ gives \, and a backslash before anything else is kept together with the character, so \n stays
 * two characters. An unclosed quote is a refusal.
 *
 * VCFCodec.readActualHeader is the second layer. It splits the version line on every '=' and
 * records a version only when that yields exactly two fields, so ##fileformat=VCFv4.2=x records
 * nothing and fails later with a different message. Its refusals are dumped with their messages,
 * because the messages are what separate two failures of the same file: "we never saw a header line
 * specifying VCF version" and "we never saw the required CHROM header line".
 *
 * Only unstructured meta lines are fed to readActualHeader here. Turning ##INFO=<...> into a typed
 * compound header line is the next layer and is measured on its own; mixing the two would leave a
 * failure ambiguous between them.
 *
 * Output:
 *
 *     line\t<label>\t<key=value pairs, separated by |>
 *     lineerror\t<label>\t<class>\t<message>
 *     frame\t<label>\t<version>\t<meta line keys, comma-separated>\t<samples, comma-separated>
 *     frameerror\t<label>\t<class>\t<message>
 *
 * Usage: VcfHeaderParseDump
 */

import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.SynchronousLineReader;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineTranslator;
import htsjdk.variant.vcf.VCFHeaderVersion;

import java.io.StringReader;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

public class VcfHeaderParseDump {

    /** The tags an INFO line is required to carry, in order, which is what drives the validation. */
    static final List<String> INFO_TAGS = List.of("ID", "Number", "Type", "Description");
    static final List<String> INFO_RECOMMENDED = List.of("Source", "Version");

    public static void main(final String[] args) {
        System.out.println("# VcfHeaderParseDump: reading a VCF header, scanner and frame");

        // The scanner, with no tag validation at all (Java's null), so the machine is measured
        // before the rules on top of it.
        line("plain", "<ID=DP,Number=1,Type=Integer,Description=\"Approximate depth\">", null);
        line("comma-in-quotes", "<ID=A,Description=\"one, two\">", null);
        line("escaped-quote", "<ID=A,Description=\"say \\\"hi\\\"\">", null);
        // A backslash before anything but a quote or a backslash is kept with the character.
        line("backslash-n", "<ID=A,Description=\"path\\nnext\">", null);
        line("double-backslash", "<ID=A,Description=\"a\\\\b\">", null);
        line("unclosed-quote", "<ID=A,Description=\"unterminated>", null);
        // The fall-through: '<' away from position 0, and '>' in the middle.
        line("angle-open-inside", "<ID=A,Foo=x<y>", null);
        line("angle-close-inside", "<ID=A,Foo=x>y>", null);
        // No trailing '>', so the last field is never stored.
        line("no-trailing-angle", "<ID=A,Number=1", null);
        // Whitespace either side of keys and values.
        line("spaces", "< ID = A , Number = 1 >", null);
        line("empty-value", "<ID=,Number=1>", null);
        // A quote opening in the middle of an unquoted value, which swallows the comma after it.
        line("quote-mid-value", "<ID=A,Desc=a\"b,c\"d>", null);
        // A repeated key: LinkedHashMap keeps the first position and the last value.
        line("repeated-key", "<ID=A,Number=1,ID=B>", null);
        // Nothing at all between the brackets.
        line("empty-brackets", "<>", null);
        // No brackets at all, which the machine does not require.
        line("no-brackets", "ID=A,Number=1", null);

        // The same scanner with the INFO tag order enforced.
        line("info-ok", "<ID=DP,Number=1,Type=Integer,Description=\"d\">", INFO_TAGS);
        line("info-wrong-order", "<Number=1,ID=DP,Type=Integer,Description=\"d\">", INFO_TAGS);
        line("info-unexpected-tag", "<ID=DP,Foo=1,Type=Integer,Description=\"d\">", INFO_TAGS);
        line("info-recommended-early", "<ID=DP,Source=\"s\",Number=1,Type=Integer,Description=\"d\">",
                INFO_TAGS);
        line("info-recommended-late",
                "<ID=DP,Number=1,Type=Integer,Description=\"d\",Source=\"s\">", INFO_TAGS);
        line("info-no-tags", "<>", INFO_TAGS);
        // More tags than expected: the validation only walks as far as the expected list.
        line("info-extra-trailing",
                "<ID=DP,Number=1,Type=Integer,Description=\"d\",Whatever=1>", INFO_TAGS);

        // The frame: version detection, the column line, and the refusals.
        frame("minimal", "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        frame("with-samples", "##fileformat=VCFv4.2\n##source=probe\n"
                + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA1\tNA2\n");
        frame("v40", "##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        frame("v43", "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        frame("v33", "##fileformat=VCFv3.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        frame("v44", "##fileformat=VCFv4.4\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        frame("unknown-version",
                "##fileformat=VCFv9.9\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        // Three equals signs, so the version line parses as three fields and records nothing.
        frame("version-line-extra-equals",
                "##fileformat=VCFv4.2=x\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        frame("no-version", "##source=probe\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        frame("no-column-line", "##fileformat=VCFv4.2\n");
        frame("data-before-column-line", "##fileformat=VCFv4.2\nchr1\t1\t.\tA\tT\t.\t.\t.\n");
        frame("too-few-columns", "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\n");
        frame("misspelled-column",
                "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTR\tINFO\n");
        frame("swapped-columns",
                "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tINFO\tFILTER\n");
        frame("format-without-samples",
                "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n");
        frame("ninth-column-not-format",
                "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tNA1\n");
        // A repeated sample name, which a LinkedHashSet collapses onto one column's worth.
        frame("repeated-sample", "##fileformat=VCFv4.2\n"
                + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA1\tNA1\tNA2\n");
    }

    static void line(final String label, final String value, final List<String> expected) {
        try {
            final Map<String, String> parsed = expected == null
                    ? VCFHeaderLineTranslator.parseLine(VCFHeaderVersion.VCF4_2, value, null)
                    : VCFHeaderLineTranslator.parseLine(VCFHeaderVersion.VCF4_2, value, expected,
                            INFO_RECOMMENDED);
            final StringJoiner pairs = new StringJoiner("|");
            for (final Map.Entry<String, String> entry : parsed.entrySet()) {
                pairs.add(entry.getKey() + "=" + entry.getValue());
            }
            System.out.printf("line\t%s\t%s%n", label, pairs);
        } catch (final Exception e) {
            System.out.printf("lineerror\t%s\t%s\t%s%n", label, e.getClass().getName(),
                    oneLine(e.getMessage()));
        }
    }

    static void frame(final String label, final String text) {
        try {
            final VCFCodec codec = new VCFCodec();
            final Object header = codec.readActualHeader(
                    new LineIteratorImpl(new SynchronousLineReader(new StringReader(text))));
            final VCFHeader vcfHeader = (VCFHeader) header;
            // The keys rather than the count, so a disagreement names which line went missing. The
            // fileformat line is not among them: VCFHeader.removeVCFVersionLines strips it from the
            // metadata it keeps, and the writer puts a constant one back.
            final StringJoiner keys = new StringJoiner(",");
            vcfHeader.getMetaDataInInputOrder().forEach(metaLine -> keys.add(metaLine.getKey()));
            System.out.printf("frame\t%s\t%s\t%s\t%s%n", label,
                    codec.getVersion().getVersionString(), keys,
                    String.join(",", vcfHeader.getGenotypeSamples()));
        } catch (final Exception e) {
            System.out.printf("frameerror\t%s\t%s\t%s%n", label, e.getClass().getName(),
                    oneLine(e.getMessage()));
        }
    }

    /** Messages are compared as bytes, so a newline in one would break the row shape. */
    static String oneLine(final String message) {
        return message == null ? "" : message.replace('\n', ' ').replace('\t', ' ');
    }
}
