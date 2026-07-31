/*
 * FastMath's exponential tables and its exp(), taken from commons-math3 3.5.
 *
 * FastMath.exp is NOT java.lang.Math.exp. It is table-driven, pure Java, and Apache 2.0, so unlike
 * the HotSpot intrinsic it can be ported (decision 0023). Everything in commons-math3 that needs
 * an exponential goes through it: Gamma, Erf, NormalDistribution, and therefore GATK's
 * MannWhitneyU and every rank-sum annotation.
 *
 * Two things are measured here.
 *
 * The TABLES. The reference ships 6,175 lines of literal doubles and can also regenerate them with
 * FastMathCalc; the flag choosing between the two is a compile-time constant, set to "use the
 * literals". The port takes the computing branch, because 3,550 transcribed literals are 3,550
 * chances to typo one. So every table entry travels in this golden as raw bits, and the two
 * branches agreeing becomes a measurement rather than an assumption. The fields are package-
 * private statics of private nested classes, so reflection is the only way to read them.
 *
 * And exp() itself, on the boundaries the algorithm's own branches name: -746 (underflow), -709
 * (the two shifted recursive cases), 709 and 710 (overflow), plus the subnormal range.
 *
 * Output:
 *
 *     table\t<name>\t<index>\t<raw long bits>
 *     exp\t<input bits>\t<result bits>
 *
 * Usage: FastMathTablesDump
 */

import org.apache.commons.math3.util.FastMath;

import java.lang.reflect.Field;

public class FastMathTablesDump {

    public static void main(final String[] args) throws Exception {
        System.out.println("# FastMathTablesDump: commons-math3 3.5 FastMath exp tables and exp()");

        emitTable("EXP_INT_TABLE_A", "org.apache.commons.math3.util.FastMath$ExpIntTable",
                "EXP_INT_TABLE_A");
        emitTable("EXP_INT_TABLE_B", "org.apache.commons.math3.util.FastMath$ExpIntTable",
                "EXP_INT_TABLE_B");
        emitTable("EXP_FRAC_TABLE_A", "org.apache.commons.math3.util.FastMath$ExpFracTable",
                "EXP_FRAC_TABLE_A");
        emitTable("EXP_FRAC_TABLE_B", "org.apache.commons.math3.util.FastMath$ExpFracTable",
                "EXP_FRAC_TABLE_B");

        // The branch boundaries, named by the algorithm itself.
        for (final double x : new double[] {
                0.0, -0.0, 1.0, -1.0, 0.5, -0.5,
                -745.0, -745.5, -746.0, -746.1, -747.0,
                -708.0, -709.0, -709.5, -709.9, -710.0, -720.0, -744.0,
                708.0, 709.0, 709.5, 709.782712893384, 710.0, 1000.0,
                Double.MIN_VALUE, -Double.MIN_VALUE, Double.MAX_VALUE, -Double.MAX_VALUE,
                Double.NaN, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY,
                1e-300, -1e-300, 1e300, -1e300,
                0.0009765625, 0.00048828125, 1.0009765625, 40.19140625, 1.494140625,
        }) {
            emitExp(x);
        }
        // A sweep through the fractional table's own resolution, where the index arithmetic
        // decides which entry is used.
        for (int i = 0; i <= 1024; i += 37) {
            emitExp(i / 1024.0);
            emitExp(-(i / 1024.0));
            emitExp(5.0 + i / 1024.0);
        }
    }

    static void emitExp(final double x) {
        System.out.printf("exp\t%d\t%d%n", Double.doubleToRawLongBits(x),
                Double.doubleToRawLongBits(FastMath.exp(x)));
    }

    static void emitTable(final String name, final String owner, final String field)
            throws Exception {
        final Class<?> type = Class.forName(owner);
        final Field declared = type.getDeclaredField(field);
        declared.setAccessible(true);
        final double[] values = (double[]) declared.get(null);
        for (int i = 0; i < values.length; i++) {
            System.out.printf("table\t%s\t%d\t%d%n", name, i, Double.doubleToRawLongBits(values[i]));
        }
    }
}
