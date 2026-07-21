/*
 * Dumps FormatUtil's output for a stratified sample of doubles and longs, under the pinned
 * en-US locale.
 *
 * Output: <input bits in hex> <TAB> <formatted string>
 *
 * FormatUtil reaches java.text.DecimalFormat, which rounds the *shortest round-trip decimal*
 * of the double (via FloatingDecimal, the same code as Double.toString) rather than the exact
 * binary value. Java 17 still uses the pre-JDK19 algorithm, which is not always shortest. So
 * whether a Rust shortest-repr agrees is a question to measure, not to assume.
 */

import htsjdk.samtools.util.FormatUtil;
import java.util.Locale;
import java.util.Random;

public class FormatDump {

    static FormatUtil f;

    static void emit(final double v) {
        System.out.println(Long.toHexString(Double.doubleToRawLongBits(v)) + "\t" + f.format(v));
    }

    public static void main(final String[] args) {
        // The oracle contract pins the locale; see decision 0011.
        Locale.setDefault(Locale.US);
        f = new FormatUtil();

        // Named values and boundaries.
        final double[] named = {
                0.0, -0.0, 1.0, -1.0, 0.5, -0.5, 1.5, 2.5, 3.5, -1.5, -2.5,
                1.0 / 3.0, 2.0 / 3.0, 1.0 / 7.0, Math.PI, Math.E,
                0.1, 0.2, 0.3, 0.7, 1e-7, 1e-6, 1e-5, 1e6, 1e7, 1e15, 1e16, 1e17,
                Double.NaN, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY,
                Double.MIN_VALUE, Double.MIN_NORMAL, Double.MAX_VALUE,
                // Ties at the sixth fraction digit, where HALF_DOWN and the JDK default
                // HALF_EVEN disagree.
                0.1234565, 0.1234575, 0.1234585, 0.0000005, 0.0000015, 0.0000025,
                -0.1234565, -0.1234575, 0.9999995, 1.0000005,
                // Values that stress the shortest-repr algorithm.
                2.82879384806159e17, 1.9400994884341945e25, 3.7238001960564653e-4,
                9.007199254740992e15, 4.9e-324, 1.7976931348623157e308};
        for (final double v : named) emit(v);

        // Percentages and ratios, which is what metrics files are mostly made of.
        for (int i = 0; i <= 1000; i++) emit(i / 1000.0);
        for (int i = 1; i <= 200; i++) emit(1.0 / i);
        for (int i = 1; i <= 200; i++) emit(i / 7.0);

        // Counts and coverage-like magnitudes.
        for (int e = -12; e <= 12; e++) {
            for (int m = 1; m <= 9; m++) emit(m * Math.pow(10, e));
        }

        // A deterministic random sample across the exponent range.
        final Random rng = new Random(20260721L);
        for (int i = 0; i < 20000; i++) {
            final double v = Double.longBitsToDouble(rng.nextLong());
            if (Double.isNaN(v)) continue;
            emit(v);
        }
        // And a sample of the ranges metrics actually occupy.
        for (int i = 0; i < 20000; i++) {
            emit(rng.nextDouble() * Math.pow(10, rng.nextInt(19) - 9));
        }

        // Integers go through a different formatter with grouping disabled.
        final long[] longs = {0, 1, -1, 9, 10, 999, 1000, 1234567, -1234567,
                Long.MAX_VALUE, Long.MIN_VALUE, Integer.MAX_VALUE, Integer.MIN_VALUE};
        for (final long v : longs) {
            System.out.println("L" + Long.toHexString(v) + "\t" + f.format(v));
        }
    }
}
