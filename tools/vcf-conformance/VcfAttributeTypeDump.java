/*
 * The typed attribute accessors: the layer between the strings a VCF is parsed into and the
 * numbers GATK's annotations do arithmetic on.
 *
 * The codec stores every INFO value as a String or a List of Strings, whatever the header declared
 * its Type to be. Nothing converts them on the way in. So the declared type does not decide what a
 * record holds, it decides which accessor a caller reaches for, and the conversion happens there,
 * once per call, with a default the caller supplies.
 *
 * Four behaviours of that conversion are not what the signatures suggest.
 *
 *   - THE MISSING-VALUE TEST IS A REFERENCE COMPARISON. getAttributeAsInt reads
 *     `x == VCFConstants.MISSING_VALUE_v4`, with ==, on a String. It is true only when the stored
 *     object is literally that constant. The codec assigns the constant for a bare key and for
 *     `KEY=`, so those return the default; a value written as `KEY=.` arrives as a substring of the
 *     line, a different reference, and reaches Integer.parseInt(".") instead. Three spellings of
 *     "missing", two outcomes, and the two are a number and an exception;
 *   - getAttributeAsDouble DOES NOT HAVE THAT TEST AT ALL. Only null is checked, so every spelling
 *     of missing reaches parseVcfDouble and throws. The Int and Double accessors disagree about
 *     what missing means;
 *   - parseVcfDouble IS NOT Double.parseDouble. On failure it retries against a pattern accepting
 *     "inf", "Infinity" and "nan" in several spellings, so a VCF may carry infinities that
 *     Double.parseDouble refuses and this accepts;
 *   - getAttributeAsList OF A SCALAR IS A LIST OF ONE, and of an absent key is empty, so a caller
 *     cannot tell "one value" from "one value because there was only one" without asking twice.
 *
 * getAttributeAsBoolean is Boolean.valueOf, which is `"true".equalsIgnoreCase(s)`: every other
 * string, "1" and "TRUE " included, is false, and nothing reports that it was not a boolean.
 *
 * Output:
 *
 *     attr\t<label>\t<key>\t<stored class>\t<stored value>
 *     as\t<label>\t<key>\t<accessor>\t<result or E:class:message>
 *     dbl\t<raw>\t<parseVcfDouble result or E:class:message>
 *
 * Usage: VcfAttributeTypeDump
 */

import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.SynchronousLineReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFUtils;

import java.io.StringReader;
import java.util.List;
import java.util.StringJoiner;

public class VcfAttributeTypeDump {

    /** One declaration per type the format has, so the header's Type is a variable and not a fixture. */
    static final String HEADER = "##fileformat=VCFv4.2\n"
            + "##INFO=<ID=I1,Number=1,Type=Integer,Description=\"one integer\">\n"
            + "##INFO=<ID=IA,Number=A,Type=Integer,Description=\"per alt integer\">\n"
            + "##INFO=<ID=F1,Number=1,Type=Float,Description=\"one float\">\n"
            + "##INFO=<ID=FR,Number=R,Type=Float,Description=\"per allele float\">\n"
            + "##INFO=<ID=S1,Number=1,Type=String,Description=\"one string\">\n"
            + "##INFO=<ID=C1,Number=1,Type=Character,Description=\"one character\">\n"
            + "##INFO=<ID=CU,Number=.,Type=Character,Description=\"characters\">\n"
            + "##INFO=<ID=B1,Number=0,Type=Flag,Description=\"a flag\">\n"
            + "##contig=<ID=chr1,length=100000>\n"
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    public static void main(final String[] args) {
        System.out.println("# VcfAttributeTypeDump: what a declared Type actually does");

        // One value of each declared type, read back through every accessor. The stored class is
        // reported first, because "the header says Integer" and "the record holds an Integer" are
        // different claims and only one of them is true.
        probe("integer", "I1=42");
        probe("integer-list", "IA=1,2,3");
        probe("float", "F1=0.5");
        probe("float-list", "FR=0.5,0.25");
        probe("string", "S1=hello");
        probe("character", "C1=x");
        probe("character-list", "CU=a,b,c");
        probe("flag", "B1");

        // The three spellings of missing, which the header cannot distinguish and the accessors do.
        probe("missing-dot", "I1=.");
        probe("missing-empty", "I1=");
        probe("missing-bare", "I1");
        probe("missing-dot-float", "F1=.");
        probe("missing-in-a-list", "IA=1,.,3");

        // Values that do not match their declared type, which nothing checks on the way in.
        probe("integer-holding-a-float", "I1=1.5");
        probe("integer-holding-text", "I1=abc");
        probe("float-holding-an-integer", "F1=7");
        probe("float-holding-text", "F1=abc");
        probe("character-holding-a-word", "C1=word");
        probe("string-holding-true", "S1=true");
        probe("string-holding-TRUE", "S1=TRUE");
        probe("string-holding-one", "S1=1");

        // The infinities and NaN, which Double.parseDouble refuses in most of these spellings.
        probe("float-inf", "F1=inf");
        probe("float-minus-inf", "F1=-inf");
        probe("float-nan", "F1=nan");

        // An absent key, which every accessor answers with its default and getAttributeAsList
        // answers with an empty list rather than a list of the default.
        probe("absent", "S1=present");

        // parseVcfDouble on its own, because its pattern is the thing and a record only reaches a
        // few of its branches.
        for (final String raw : new String[] {
                "1", "1.5", "-1.5", "1e3", "1E3", "+1", " 1", "1 ", "1f", "1d", "0x1p3",
                "Infinity", "-Infinity", "+Infinity", "inf", "-inf", "+inf", "INF", "Inf",
                "infinity", "nan", "NaN", "NAN", "-nan", "+nan", ".", "", "1,5" }) {
            System.out.printf("dbl\t%s\t%s%n", escape(raw), attempt(() -> {
                final double value = VCFUtils.parseVcfDouble(raw);
                return Double.toString(value);
            }));
        }
    }

    static void probe(final String label, final String info) {
        final VCFCodec codec = new VCFCodec();
        codec.readActualHeader(
                new LineIteratorImpl(new SynchronousLineReader(new StringReader(HEADER))));
        final VariantContext vc =
                codec.decode("chr1\t100\t.\tA\tT\t50\tPASS\t" + info);

        // The key under test is the one the case wrote; "absent" asks about a key it did not.
        final String key = label.equals("absent") ? "MISSING" : info.split("[=;]")[0];

        final Object stored = vc.getAttribute(key);
        System.out.printf("attr\t%s\t%s\t%s\t%s%n", label, key,
                stored == null ? "null" : stored.getClass().getSimpleName(),
                escape(String.valueOf(stored)));

        emit(label, key, "asString", () -> vc.getAttributeAsString(key, "DEFAULT"));
        emit(label, key, "asInt", () -> Integer.toString(vc.getAttributeAsInt(key, -1)));
        emit(label, key, "asDouble", () -> Double.toString(vc.getAttributeAsDouble(key, -1.0)));
        emit(label, key, "asBoolean", () -> Boolean.toString(vc.getAttributeAsBoolean(key, false)));
        emit(label, key, "asList", () -> render(vc.getAttributeAsList(key)));
        emit(label, key, "asStringList", () -> render(vc.getAttributeAsStringList(key, "D")));
        emit(label, key, "asIntList", () -> render(vc.getAttributeAsIntList(key, -1)));
        emit(label, key, "asDoubleList", () -> render(vc.getAttributeAsDoubleList(key, -1.0)));
    }

    static void emit(final String label, final String key, final String accessor,
                     final Attempt attempt) {
        System.out.printf("as\t%s\t%s\t%s\t%s%n", label, key, accessor, attempt(attempt));
    }

    interface Attempt {
        String get() throws Exception;
    }

    static String attempt(final Attempt attempt) {
        try {
            return escape(attempt.get());
        } catch (final Throwable t) {
            return "E:" + t.getClass().getName() + ":" + escape(String.valueOf(t.getMessage()));
        }
    }

    static String render(final List<?> values) {
        final StringJoiner joined = new StringJoiner(",");
        for (final Object value : values) {
            joined.add(String.valueOf(value));
        }
        return "[" + joined + "]";
    }

    static String escape(final String s) {
        if (s == null) {
            return "null";
        }
        final StringBuilder b = new StringBuilder(s.length());
        for (int i = 0; i < s.length(); i++) {
            final char c = s.charAt(i);
            switch (c) {
                case '\\': b.append("\\\\"); break;
                case '\t': b.append("\\t"); break;
                case '\n': b.append("\\n"); break;
                default:
                    if (c < 0x20 || c > 0x7e) {
                        b.append(String.format("\\u%04x", (int) c));
                    } else {
                        b.append(c);
                    }
            }
        }
        return b.toString();
    }
}
