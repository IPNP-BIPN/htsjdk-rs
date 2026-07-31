/*
 * Percentile, Median and FastMath.round, taken from commons-math3 3.5.
 *
 * 3.5 and not 3.6.1: GATK pins commons-math3 with `strictly '3.5'`, so that is the version whose
 * numbers reach a golden downstream.
 *
 * Four things here are not what "the median" suggests:
 *
 *   - the estimator INTERPOLATES, and the arithmetic is `lower + dif * (upper - lower)` rather
 *     than `(lower + upper) / 2`, which is a different double for values far apart;
 *   - the default NaNStrategy is REMOVED, so a NaN shortens the array instead of ranking, and an
 *     array of nothing but NaN answers NaN by way of an empty array;
 *   - a single value is returned BEFORE the NaN strategy runs, so `{NaN}` and `{NaN, NaN}` reach
 *     the same answer down two different paths;
 *   - MathUtils.median finishes with FastMath.round, which is `(long) floor(x + 0.5)`: the
 *     definition java.lang.Math.round stopped using in Java 7, and still one apart from it on
 *     0.49999999999999994.
 *
 * Output:
 *
 *     percentile\t<type>\t<quantile>\t<input, bits, comma-separated>\t<result bits>
 *     median\t<type>\t<int input, comma-separated>\t<int result>
 *     round\t<bits>\t<FastMath.round>\t<Math.round>
 *
 * Doubles travel as raw bits, so a golden cannot lose a last-bit difference to formatting.
 *
 * Usage: PercentileDump
 */

import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.commons.math3.util.FastMath;

public class PercentileDump {

    /** The arrays every estimator is asked about. */
    static final double[][] INPUTS = {
        {},
        {1},
        {Double.NaN},
        {1, 2},
        {2, 1},
        {1, 2, 3},
        {1, 2, 3, 4},
        {1, 2, 3, 4, 5},
        {1, 2, 3, 4, 5, 6},
        // Ties, which decide nothing here and would decide a lot in a rank-based estimator.
        {5, 5, 5, 5},
        {1, 1, 2, 2},
        // NaN, removed rather than ranked.
        {1, Double.NaN, 2, 3},
        {Double.NaN, Double.NaN},
        {Double.NaN, 1},
        {1, 2, Double.NaN},
        // Infinities, which the interpolation can turn into NaN.
        {Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY},
        {Double.NEGATIVE_INFINITY, 0, Double.POSITIVE_INFINITY},
        {Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY},
        // The two zeros, which sort apart and compare equal.
        {-0.0, 0.0},
        {0.0, -0.0},
        {-0.0, -0.0},
        // Values far enough apart that `lower + dif * (upper - lower)` and `(lower + upper) / 2`
        // are different doubles.
        {1e300, 1e-300},
        {Double.MAX_VALUE, Double.MAX_VALUE},
        {Double.MAX_VALUE, -Double.MAX_VALUE},
        {Double.MIN_VALUE, 0},
        // A median that lands exactly on the value Math.round and FastMath.round disagree about.
        {0.0, 0.999_999_999_999_999_88},
    };

    static final double[] QUANTILES = {50, 1, 25, 75, 99, 100};

    /** The int arrays MathUtils.median is asked about, which is where the rounding shows. */
    static final int[][] INT_INPUTS = {
        {}, {1}, {1, 2}, {1, 2, 3}, {1, 2, 3, 4}, {0, 1}, {0, 1, 2, 3},
        {-1, 0}, {-2, -1}, {-1, 2}, {Integer.MAX_VALUE, Integer.MAX_VALUE},
        {Integer.MIN_VALUE, Integer.MAX_VALUE}, {Integer.MIN_VALUE, Integer.MIN_VALUE},
        {10, 20, 30, 40, 50, 60, 70},
    };

    /** The doubles FastMath.round and Math.round are compared on. */
    static final double[] ROUNDED = {
        0.5, -0.5, 0.49999999999999994, -0.49999999999999994, 1.5, 2.5, -1.5, -2.5,
        0.0, -0.0, 1.0, -1.0, 0.1, 0.9, Double.NaN, Double.POSITIVE_INFINITY,
        Double.NEGATIVE_INFINITY, Double.MAX_VALUE, -Double.MAX_VALUE, Double.MIN_VALUE,
        9.007199254740992E15, 4.503599627370496E15,
    };

    public static void main(final String[] args) {
        System.out.println("# PercentileDump: commons-math3 3.5 Percentile, Median, FastMath.round");
        System.out.println("# commons-math3=" + Percentile.class.getPackage().getImplementationVersion());

        for (final Percentile.EstimationType type : new Percentile.EstimationType[] {
                Percentile.EstimationType.LEGACY, Percentile.EstimationType.R_1}) {
            for (final double[] input : INPUTS) {
                for (final double quantile : QUANTILES) {
                    emit(type, quantile, input);
                }
            }
        }

        for (final int[] input : INT_INPUTS) {
            emitMedian("LEGACY", input, medianOf(input, Percentile.EstimationType.LEGACY));
            emitMedian("R_1", input, medianOf(input, Percentile.EstimationType.R_1));
        }

        for (final double x : ROUNDED) {
            System.out.printf("round\t%d\t%d\t%d%n", Double.doubleToRawLongBits(x),
                    FastMath.round(x), Math.round(x));
        }
    }

    static void emit(final Percentile.EstimationType type, final double quantile,
                     final double[] input) {
        final String shown = bits(input);
        try {
            final double result = new Percentile(quantile).withEstimationType(type)
                    .evaluate(input.clone());
            System.out.printf("percentile\t%s\t%s\t%s\t%d%n", type, quantile, shown,
                    Double.doubleToRawLongBits(result));
        } catch (final Exception | AssertionError e) {
            System.out.printf("percentile\t%s\t%s\t%s\tE:%s:%s%n", type, quantile, shown,
                    e.getClass().getName(),
                    e.getMessage() == null ? "" : e.getMessage().replace('\n', ' '));
        }
    }

    /** `MathUtils.median(int[], type)`, inlined so this harness does not need GATK on the path. */
    static int medianOf(final int[] values, final Percentile.EstimationType type) {
        final double[] doubles = new double[values.length];
        for (int i = 0; i < values.length; i++) {
            doubles[i] = values[i];
        }
        return (int) FastMath.round(new Median().withEstimationType(type).evaluate(doubles));
    }

    static void emitMedian(final String type, final int[] input, final int result) {
        final StringBuilder shown = new StringBuilder();
        for (final int value : input) {
            if (shown.length() > 0) {
                shown.append(',');
            }
            shown.append(value);
        }
        System.out.printf("median\t%s\t%s\t%d%n", type, shown, result);
    }

    static String bits(final double[] values) {
        final StringBuilder out = new StringBuilder();
        for (final double value : values) {
            if (out.length() > 0) {
                out.append(',');
            }
            out.append(Double.doubleToRawLongBits(value));
        }
        return out.toString();
    }
}
