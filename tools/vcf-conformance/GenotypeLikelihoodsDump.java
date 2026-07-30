/*
 * GL to PL, and the double parser underneath it, taken from the reference.
 *
 * A genotype carrying GL and no PL gets its PL computed here, and every downstream tool reads the
 * PL. So this is the last thing between a VCF with likelihoods in it and a port that can decode it.
 *
 * Three layers are dumped, smallest first.
 *
 *   - Double.parseDouble, which is NOT Rust's f64 parser. Java accepts a trailing type suffix
 *     (1.5f, 1.5d), hexadecimal floating point (0x1p3), leading and trailing whitespace, and the
 *     spelled-out "Infinity"; Rust's parser accepts none of those and accepts "inf" and "nan",
 *     which Java does not;
 *   - VCFUtils.parseVcfDouble, which catches the failure and retries against
 *     ^(?<sign>[-+]?)((?<inf>(INF|INFINITY))|(?<nan>NAN))$, case-insensitive, so "inf" and
 *     "-INFINITY" and "nan" parse after all. The order matters: the pattern is tried only after
 *     the plain parse has failed;
 *   - GenotypeLikelihoods.fromGLField().getAsPLs(), which is one line with three decisions in it:
 *
 *         pls[i] = (int) Math.round(Math.min(-10 * (GLs[i] - adjust), MAX_PL));
 *
 *     adjust is the MAXIMUM of the likelihoods, so the conversion re-normalises and the best
 *     genotype always comes out at 0 whatever the input scale. MAX_PL is Integer.MAX_VALUE and the
 *     clamp is applied to the double BEFORE the rounding. Math.round is floor(x + 0.5), so it
 *     rounds half up rather than half away from zero, it answers 0 for a NaN, and it returns a long
 *     that is then truncated to int. Reordering any of the three agrees on every ordinary genotype.
 *
 * The doubles are dumped as raw bits, because a decimal rendering of a parse result hides exactly
 * the disagreements this suite exists to find.
 *
 * Output:
 *
 *     double\t<escaped input>\t<raw bits|E:class:message>
 *     vcfdouble\t<escaped input>\t<raw bits|E:class:message>
 *     gl\t<escaped input>\t<PLs, comma-separated|null|E:class:message>
 *     glbits\t<escaped input>\t<parsed likelihoods as raw bits, comma-separated|null>
 *     round\t<raw bits of input>\t<Math.round result>
 *
 * Usage: GenotypeLikelihoodsDump
 */

import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import htsjdk.variant.vcf.VCFUtils;

import java.util.StringJoiner;

public class GenotypeLikelihoodsDump {

    /**
     * Everything the two parsers are asked. The plain decimals are there so the two agree
     * somewhere; the rest is where they might not.
     */
    static final String[] DOUBLES = {
        // Ordinary.
        "0", "1", "-1", "0.0", "-0.0", "1.5", "-1.5", "3.14159265358979",
        "1e3", "1E3", "1e-3", "-1.5e-10", "1.7976931348623157E308", "4.9E-324",
        // Overflow and underflow, which are infinity and zero rather than errors.
        "1e400", "-1e400", "1e-400",
        // The type suffixes.
        "1.5f", "1.5F", "1.5d", "1.5D", "1f", "1d", "-1.5f",
        // Hexadecimal floating point, where the p exponent is mandatory.
        "0x1p3", "0X1P3", "0x1.8p1", "-0x1p-1", "0x1", "0x1p", "0x.8p1",
        // Whitespace, which Java trims and only up to and including the space character.
        " 1.5", "1.5 ", "  1.5  ", "\t1.5\n", "1.5 ",
        // The spelled-out specials, which Double.parseDouble takes and is case-sensitive about.
        "NaN", "Infinity", "-Infinity", "+Infinity", "-NaN",
        // The spellings only parseVcfDouble takes.
        "inf", "INF", "-inf", "+INF", "Inf", "infinity", "INFINITY", "-INFINITY",
        "nan", "NAN", "-nan", "NaNq",
        // Refusals.
        "", " ", ".", "abc", "1.5.5", "1,5", "--1", "1_000", "0x1p3.5",
        // Things that look like they should work and do not.
        "1.5e", "e3", "+", "-", "1.5ff",
    };

    /** GL fields, as they appear in a VCF's genotype column. */
    static final String[] GL_FIELDS = {
        // The ordinary shape: three likelihoods for a diploid biallelic site.
        "-0.1,-0.2,-0.3",
        "0,-1,-2",
        "-2,-1,0",
        "-1,-1,-1",
        // Already normalised, so the conversion is a pure scaling.
        "0,-3,-6",
        // Positive likelihoods, which are nonsense as log10 probabilities and are converted anyway.
        "1,2,3",
        // A single value, and many: the field's length is the genotype count and is not checked.
        "-1",
        "-0.1,-0.2,-0.3,-0.4,-0.5,-0.6",
        // Where the rounding shows: exactly half, on both signs of the shift.
        "-0.05,0,-0.15",
        "0,-0.05,-0.15",
        // Very large and very small, where the clamp to Integer.MAX_VALUE bites.
        "0,-1e9,-1e10",
        "0,-1e300,-1e-300",
        // Infinities and NaN, which Math.max propagates.
        "0,-Infinity,-1",
        "0,inf,-1",
        "0,nan,-1",
        "nan,nan,nan",
        // The missing-value rules: all missing, none missing, some missing.
        ".",
        ".,.,.",
        "-1,.,-3",
        ".,-2,-3",
        "-1,-2,.",
        // Empty and malformed elements.
        "",
        ",",
        "-1,,-3",
        "-1,-2,",
        ",-2,-3",
        "-1,abc,-3",
        // The type suffixes and whitespace, inside a GL field.
        "-1.5f,-2,-3",
        " -1, -2, -3",
    };

    /** Inputs for Math.round on its own, dumped as raw bits so the halves are exact. */
    static final double[] ROUNDS = {
        0.0, -0.0, 0.5, -0.5, 1.5, -1.5, 2.5, -2.5, 0.49999999999999994,
        -0.49999999999999994, 1e18, -1e18, 1e19, -1e19,
        Double.NaN, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY,
        Double.MAX_VALUE, -Double.MAX_VALUE, Integer.MAX_VALUE, Integer.MAX_VALUE + 0.5,
    };

    public static void main(final String[] args) {
        System.out.println("# GenotypeLikelihoodsDump: GL to PL and the double parser under it");

        for (final String text : DOUBLES) {
            emit("double", text, () -> Double.toString(Double.parseDouble(text)),
                    () -> Double.doubleToRawLongBits(Double.parseDouble(text)));
            emit("vcfdouble", text, () -> Double.toString(VCFUtils.parseVcfDouble(text)),
                    () -> Double.doubleToRawLongBits(VCFUtils.parseVcfDouble(text)));
        }

        for (final String field : GL_FIELDS) {
            try {
                final GenotypeLikelihoods likelihoods = GenotypeLikelihoods.fromGLField(field);
                final double[] vector = likelihoods.getAsVector();
                if (vector == null) {
                    System.out.printf("glbits\t%s\tnull%n", escape(field));
                } else {
                    final StringJoiner bits = new StringJoiner(",");
                    for (final double value : vector) {
                        bits.add(Long.toString(Double.doubleToRawLongBits(value)));
                    }
                    System.out.printf("glbits\t%s\t%s%n", escape(field), bits);
                }

                final int[] pls = likelihoods.getAsPLs();
                if (pls == null) {
                    System.out.printf("gl\t%s\tnull%n", escape(field));
                } else {
                    final StringJoiner joiner = new StringJoiner(",");
                    for (final int pl : pls) {
                        joiner.add(Integer.toString(pl));
                    }
                    System.out.printf("gl\t%s\t%s%n", escape(field), joiner);
                }
            } catch (final Exception | AssertionError e) {
                System.out.printf("gl\t%s\tE:%s:%s%n", escape(field), e.getClass().getName(),
                        oneLine(e.getMessage()));
            }
        }

        for (final double value : ROUNDS) {
            System.out.printf("round\t%d\t%d%n", Double.doubleToRawLongBits(value),
                    Math.round(value));
        }
    }

    interface Bits {
        long get();
    }

    interface Text {
        String get();
    }

    static void emit(final String kind, final String text, final Text ignored, final Bits bits) {
        try {
            System.out.printf("%s\t%s\t%d%n", kind, escape(text), bits.get());
        } catch (final Exception | AssertionError e) {
            System.out.printf("%s\t%s\tE:%s:%s%n", kind, escape(text), e.getClass().getName(),
                    oneLine(e.getMessage()));
        }
    }

    /** Whitespace and non-ASCII travel escaped, because the inputs deliberately contain both. */
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

    static String oneLine(final String message) {
        return message == null ? "" : message.replace('\n', ' ').replace('\t', ' ');
    }
}
