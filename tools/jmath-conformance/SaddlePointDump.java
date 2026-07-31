/*
 * SaddlePointExpansion and HypergeometricDistribution.logProbability, taken from commons-math3 3.5.
 *
 * GATK's FisherExactTest builds a HypergeometricDistribution and asks it for the log probability of
 * every point of the support, so this is the arithmetic under FS, the Fisher-strand annotation.
 *
 * Three of its decisions are not what a probability formula suggests:
 *
 *   - getStirlingError has THREE algorithms. Below 15 and on a half-integer it reads a table of
 *     literals; below 15 and off a half-integer it goes through Gamma.logGamma; at 15 and above it
 *     is a five-term asymptotic series. The boundaries are exact comparisons on a double, and the
 *     table's last entry is dead code because the guard is `z < 15.0`;
 *   - getDeviancePart iterates until `s1 != s` stops being true, so it terminates on the rounding
 *     rather than on a tolerance;
 *   - logBinomialProbability has four branches at the ends of its range, chosen by comparing p or
 *     q against 0.1, and only the middle one uses the Stirling error at all.
 *
 * The class and its methods are package-private, so this dump reaches them by reflection.
 *
 * Output:
 *
 *     saddle\t<function>\t<args, comma-separated bits or ints>\t<result bits>
 *
 * Usage: SaddlePointDump
 */

import java.lang.reflect.Method;

public class SaddlePointDump {

    public static void main(final String[] args) throws Exception {
        System.out.println("# SaddlePointDump: SaddlePointExpansion and the hypergeometric log probability");

        final Class<?> saddle =
                Class.forName("org.apache.commons.math3.distribution.SaddlePointExpansion");
        final Method stirling = saddle.getDeclaredMethod("getStirlingError", double.class);
        final Method deviance = saddle.getDeclaredMethod("getDeviancePart", double.class, double.class);
        final Method binomial = saddle.getDeclaredMethod("logBinomialProbability", int.class,
                int.class, double.class, double.class);
        stirling.setAccessible(true);
        deviance.setAccessible(true);
        binomial.setAccessible(true);

        // Every half-integer the table covers, and the two sides of every boundary.
        for (int i = 0; i <= 40; i++) {
            emit1("getStirlingError", i / 2.0, (double) stirling.invoke(null, i / 2.0));
        }
        for (final double z : new double[] {
                0.0, 0.25, 0.75, 1.1, 7.3, 14.4, 14.9, 14.999999999999998, 15.0,
                15.000000000000002, 15.5, 20.0, 100.0, 1000.0, 1e6, 1e300,
                Double.MIN_VALUE, Double.NaN, Double.POSITIVE_INFINITY}) {
            emit1("getStirlingError", z, (double) stirling.invoke(null, z));
        }

        for (final double[] pair : new double[][] {
                {1, 1}, {1, 1.05}, {1, 1.2}, {10, 10.5}, {10, 20}, {0, 1}, {1, 0},
                {100, 101}, {100, 200}, {1e-8, 1e-8}, {5, 5.0000001}, {1000, 1001}}) {
            System.out.printf("saddle\tgetDeviancePart\t%d,%d\t%d%n",
                    Double.doubleToRawLongBits(pair[0]), Double.doubleToRawLongBits(pair[1]),
                    Double.doubleToRawLongBits(
                            (double) deviance.invoke(null, pair[0], pair[1])));
        }

        for (final int n : new int[] {1, 2, 5, 20, 200}) {
            for (int x = 0; x <= n; x += Math.max(1, n / 7)) {
                for (final double p : new double[] {0.05, 0.1, 0.5, 0.9, 0.95}) {
                    System.out.printf("saddle\tlogBinomialProbability\t%d,%d,%d,%d\t%d%n", x, n,
                            Double.doubleToRawLongBits(p), Double.doubleToRawLongBits(1 - p),
                            Double.doubleToRawLongBits(
                                    (double) binomial.invoke(null, x, n, p, 1 - p)));
                }
            }
        }

        // The distribution itself, over tables the strand-bias annotations actually produce.
        for (final int[] shape : new int[][] {
                {10, 5, 5}, {20, 10, 10}, {400, 200, 200}, {100, 3, 50}, {7, 7, 7}, {2, 1, 1}}) {
            final org.apache.commons.math3.distribution.HypergeometricDistribution dist =
                    new org.apache.commons.math3.distribution.HypergeometricDistribution(
                            null, shape[0], shape[1], shape[2]);
            for (int x = -1; x <= Math.min(shape[1], shape[2]) + 1; x++) {
                System.out.printf("saddle\tlogProbability\t%d,%d,%d,%d\t%d%n", shape[0], shape[1],
                        shape[2], x, Double.doubleToRawLongBits(dist.logProbability(x)));
            }
        }
    }

    static void emit1(final String name, final double x, final double result) {
        System.out.printf("saddle\t%s\t%d\t%d%n", name, Double.doubleToRawLongBits(x),
                Double.doubleToRawLongBits(result));
    }
}
