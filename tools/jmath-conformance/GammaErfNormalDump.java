/*
 * Gamma, Erf and NormalDistribution, taken from commons-math3 3.5.
 *
 * Erf is the whole of NormalDistribution.cumulativeProbability, which is the whole of the normal
 * approximation in GATK's MannWhitneyU, which is what every rank-sum annotation reports. So these
 * three files decide the last digit of MQRankSum.
 *
 * Four things here are not what the names suggest:
 *
 *   - erf is not a polynomial approximation. It is an incomplete gamma function evaluated by
 *     series or by continued fraction, stopped by a RELATIVE tolerance of 1e-15 or by 10,000
 *     iterations, so where the iteration stops is part of the result;
 *   - regularizedGammaP delegates to regularizedGammaQ when x >= a + 1 and Q delegates back to P
 *     when x < a + 1, so the pair is one function with a switch in the middle and the two sides
 *     of x = a + 1 are computed by different algorithms;
 *   - erf and erfInv are inverses in mathematics and unrelated in code: one is a gamma integral,
 *     the other a rational approximation, and a round trip does not return its input;
 *   - both erf and cumulativeProbability shortcut beyond 40 (sigma), returning exactly 0, 1 or
 *     -1 without computing anything, so each is discontinuous at a value the author chose.
 *
 * Output:
 *
 *     gamma\t<function>\t<input bits...>\t<result bits>
 *     gamma\t<function>\t<input bits...>\tE:<class>
 *
 * Usage: GammaErfNormalDump
 */

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.special.Gamma;

public class GammaErfNormalDump {

    static final double[] SMALL = {
        0.0, -0.0, 1.0, -1.0, 0.5, -0.5, 0.25, 2.0, 2.5, 3.0, 8.0, 8.5, 40.0, 40.0000001, 41.0,
        -40.0, -41.0, 1e-8, -1e-8, 1e-300, 100.0, 170.0, 1e10,
        Double.NaN, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY,
        Double.MIN_VALUE, 0.9999999999999999, 1.0000000000000002,
    };

    public static void main(final String[] args) {
        System.out.println("# GammaErfNormalDump: Gamma, Erf and NormalDistribution from 3.5");

        for (final double x : SMALL) {
            one("logGamma", x, () -> Gamma.logGamma(x));
            one("lanczos", x, () -> Gamma.lanczos(x));
            one("erf", x, () -> Erf.erf(x));
            one("erfc", x, () -> Erf.erfc(x));
            one("erfInv", x, () -> Erf.erfInv(x));
            one("invGamma1pm1", x, () -> Gamma.invGamma1pm1(x));
            one("logGamma1p", x, () -> Gamma.logGamma1p(x));
            one("digamma", x, () -> Gamma.digamma(x));
            one("trigamma", x, () -> Gamma.trigamma(x));
        }

        // digamma and trigamma have three branches each, at 1e-5 and at 49, and the middle one
        // recurses until it reaches the upper branch. Cross both boundaries and both signs.
        for (final double x : new double[] {
                1e-6, 1e-5, 1.0000001e-5, 0.1, 1.0, 48.0, 48.999999, 49.0, 49.000001, 50.0,
                1000.0, -0.5, -1.5, -2.5, -48.5, -49.5, -0.0000001}) {
            one("digamma", x, () -> Gamma.digamma(x));
            one("trigamma", x, () -> Gamma.trigamma(x));
        }
        for (int i = 1; i <= 120; i++) {
            final double x = i / 2.0;
            one("digamma", x, () -> Gamma.digamma(x));
            one("trigamma", x, () -> Gamma.trigamma(x));
        }

        // A sweep, so the branch boundaries of each are crossed rather than sampled at round
        // numbers only.
        for (int i = -80; i <= 80; i++) {
            final double x = i / 8.0;
            one("erf", x, () -> Erf.erf(x));
            one("erfc", x, () -> Erf.erfc(x));
            one("logGamma", x, () -> Gamma.logGamma(x));
        }
        for (int i = -40; i <= 40; i++) {
            final double x = i / 41.0;
            one("erfInv", x, () -> Erf.erfInv(x));
        }

        // The two-argument gamma functions, either side of the x = a + 1 switch.
        for (final double a : new double[] {0.5, 1.0, 2.0, 5.0, 0.001, 100.0}) {
            for (final double x : new double[] {
                    0.0, 0.25, 0.5, a, a + 0.999999, a + 1.0, a + 1.000001, 2 * a, 10 * a, 1e6}) {
                two("regularizedGammaP", a, x, () -> Gamma.regularizedGammaP(a, x));
                two("regularizedGammaQ", a, x, () -> Gamma.regularizedGammaQ(a, x));
            }
        }

        final NormalDistribution normal = new NormalDistribution();
        for (final double x : SMALL) {
            one("cumulativeProbability", x, () -> normal.cumulativeProbability(x));
        }
        for (int i = -100; i <= 100; i++) {
            final double x = i / 10.0;
            one("cumulativeProbability", x, () -> normal.cumulativeProbability(x));
        }
        for (int i = 0; i <= 100; i++) {
            final double p = i / 100.0;
            one("inverseCumulativeProbability", p, () -> normal.inverseCumulativeProbability(p));
        }
        for (final double p : new double[] {
                1e-15, 1e-10, 0.5 - 1e-15, 0.5 + 1e-15, 1 - 1e-15, -0.1, 1.1, Double.NaN}) {
            one("inverseCumulativeProbability", p, () -> normal.inverseCumulativeProbability(p));
        }
    }

    interface Call {
        double get();
    }

    // Throwable, not Exception: Gamma.digamma has no NaN guard in 3.5, so digamma(NaN) and
    // digamma(-Infinity) recurse forever and raise a StackOverflowError, which is an Error. Caught
    // here so the run completes and the golden records which inputs do not terminate.
    static void one(final String name, final double x, final Call call) {
        try {
            System.out.printf("gamma\t%s\t%d\t%d%n", name, Double.doubleToRawLongBits(x),
                    Double.doubleToRawLongBits(call.get()));
        } catch (final Throwable e) {
            System.out.printf("gamma\t%s\t%d\tE:%s%n", name, Double.doubleToRawLongBits(x),
                    e.getClass().getName());
        }
    }

    static void two(final String name, final double a, final double x, final Call call) {
        try {
            System.out.printf("gamma\t%s\t%d,%d\t%d%n", name, Double.doubleToRawLongBits(a),
                    Double.doubleToRawLongBits(x), Double.doubleToRawLongBits(call.get()));
        } catch (final Throwable e) {
            System.out.printf("gamma\t%s\t%d,%d\tE:%s%n", name, Double.doubleToRawLongBits(a),
                    Double.doubleToRawLongBits(x), e.getClass().getName());
        }
    }
}
