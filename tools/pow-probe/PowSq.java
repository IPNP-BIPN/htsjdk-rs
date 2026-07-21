/* Is Math.pow(x, 2.0) bit-identical to x * x?
 *
 * Histogram.getStandardDeviation calls pow(value - mean, 2) through
 * `import static java.lang.Math.*`. Decision 0007 deferred Math.pow because its intrinsic uses
 * rcpps, an approximate instruction. If pow(x, 2) is exactly x * x, that dependency is
 * harmless and the deferral holds. If not, the calibration triple's first member is blocked
 * on a 2,220-line assembly port.
 *
 * Measured, not assumed.
 */
import java.util.Random;

public class PowSq {
    public static void main(String[] a) {
        long checked = 0, differ = 0;
        Object firstDiff = null;

        double[] named = {0.0, -0.0, 1.0, -1.0, 0.5, 3.0, 1e-160, 1e160, 1e-320,
                          Double.MIN_VALUE, Double.MAX_VALUE, Math.PI, 1.0/3.0,
                          Double.NaN, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY};
        Random rng = new Random(20260721L);
        double[][] sets = new double[2][];
        sets[0] = named;
        double[] rand = new double[2000000];
        for (int i = 0; i < rand.length; i++) {
            // Values in the range a standard deviation actually sees, plus wild ones.
            rand[i] = (i % 2 == 0)
                ? (rng.nextDouble() - 0.5) * Math.pow(10, rng.nextInt(12) - 6)
                : Double.longBitsToDouble(rng.nextLong());
        }
        sets[1] = rand;

        for (double[] set : sets) {
            for (double x : set) {
                if (Double.isNaN(x)) continue;
                double p = Math.pow(x, 2.0);
                double s = x * x;
                checked++;
                if (Double.doubleToRawLongBits(p) != Double.doubleToRawLongBits(s)) {
                    differ++;
                    if (firstDiff == null) {
                        firstDiff = String.format("x=%s pow=%s sq=%s",
                            Double.toHexString(x), Double.toHexString(p), Double.toHexString(s));
                    }
                }
            }
        }
        System.out.println("checked=" + checked + " differ=" + differ);
        if (firstDiff != null) System.out.println("first: " + firstDiff);
    }
}
