import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import org.apache.commons.math3.util.FastMath;

/**
 * Dumps (input bits -> output bits) for the math functions GATK and Picard put on their output
 * paths, so a Rust port can be checked against the reference bit for bit.
 *
 * <p>Both `Math` and `StrictMath` are emitted for every point. They are not the same function:
 * `StrictMath` is fdlibm and portable by specification, while `Math` is free to use a
 * platform-specific HotSpot intrinsic. GATK calls `Math`, which is exactly why
 * broadinstitute/gatk#9384 found results differing by ~1 ULP across CPU architectures. Emitting
 * both makes the size of that gap measurable rather than assumed.
 *
 * <p>Output is CSV: {@code function,input_bits_hex,math_bits_hex,strictmath_bits_hex}. Bits, not
 * decimal, because decimal rendering loses exactly the difference being measured.
 */
public final class JMathDump {

    /** Sampling is stratified, not uniform: the interesting failures live at the edges. */
    private static double[] sampleInputs() {
        final java.util.ArrayList<Double> xs = new java.util.ArrayList<>();

        // Exact special values.
        for (final double d : new double[]{
                0.0, -0.0, 1.0, -1.0, 2.0, 0.5, 10.0, Math.E, Math.PI,
                Double.MIN_VALUE, Double.MIN_NORMAL, Double.MAX_VALUE,
                Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NaN}) {
            xs.add(d);
        }

        // Probability-like range: what log10 sees in every genotype likelihood computation.
        for (int i = 1; i <= 2000; i++) {
            xs.add(i / 2000.0);
        }

        // Phred-like and log-space magnitudes: what exp sees in PairHMM and BQSR.
        for (int i = -700; i <= 700; i += 3) {
            xs.add((double) i);
            xs.add(i + 0.5);
        }

        // Dense sweep just around 1.0, where relative error is most visible.
        for (int i = -1000; i <= 1000; i++) {
            xs.add(1.0 + i * 1e-9);
        }

        // Subnormals and the boundary into them.
        long bits = 1L;
        for (int i = 0; i < 200; i++) {
            xs.add(Double.longBitsToDouble(bits));
            bits <<= 1;
            if (bits == 0) break;
        }

        // Deterministic pseudo-random bit patterns, filtered to finite values.
        long s = 88172645463325252L;
        for (int i = 0; i < 40000; i++) {
            s ^= s << 13; s ^= s >>> 7; s ^= s << 17;
            final double d = Double.longBitsToDouble(s);
            if (!Double.isNaN(d) && !Double.isInfinite(d)) {
                xs.add(d);
            }
        }

        final double[] out = new double[xs.size()];
        for (int i = 0; i < out.length; i++) out[i] = xs.get(i);
        return out;
    }

    private static void emit(BufferedWriter w, String fn, double x, double m, double s, double fm)
            throws IOException {
        w.write(fn);
        w.write(',');
        w.write(Long.toHexString(Double.doubleToRawLongBits(x)));
        w.write(',');
        w.write(Long.toHexString(Double.doubleToRawLongBits(m)));
        w.write(',');
        w.write(Long.toHexString(Double.doubleToRawLongBits(s)));
        w.write(',');
        w.write(Long.toHexString(Double.doubleToRawLongBits(fm)));
        w.write('\n');
    }

    public static void main(String[] args) throws IOException {
        final String path = args.length > 0 ? args[0] : "/out/jmath.csv";
        final double[] xs = sampleInputs();

        try (BufferedWriter w = new BufferedWriter(new FileWriter(path))) {
            w.write("# function,input_bits,math_bits,strictmath_bits,fastmath_bits\n");
            w.write("# java.version=" + System.getProperty("java.version")
                    + " os.arch=" + System.getProperty("os.arch") + "\n");

            for (final double x : xs) {
                emit(w, "exp", x, Math.exp(x), StrictMath.exp(x), FastMath.exp(x));
                emit(w, "log", x, Math.log(x), StrictMath.log(x), FastMath.log(x));
                emit(w, "log10", x, Math.log10(x), StrictMath.log10(x), FastMath.log10(x));
                emit(w, "log1p", x, Math.log1p(x), StrictMath.log1p(x), FastMath.log1p(x));
                emit(w, "expm1", x, Math.expm1(x), StrictMath.expm1(x), FastMath.expm1(x));
                emit(w, "sqrt", x, Math.sqrt(x), StrictMath.sqrt(x), FastMath.sqrt(x));
                emit(w, "cbrt", x, Math.cbrt(x), StrictMath.cbrt(x), FastMath.cbrt(x));
                emit(w, "sin", x, Math.sin(x), StrictMath.sin(x), FastMath.sin(x));
                emit(w, "cos", x, Math.cos(x), StrictMath.cos(x), FastMath.cos(x));
            }

            // pow needs two arguments; sweep a fixed set of exponents against the same inputs.
            final double[] exponents = {0.0, 1.0, 2.0, 0.5, -1.0, 10.0, -0.5, 3.0, 1.5};
            for (final double x : xs) {
                for (final double y : exponents) {
                    w.write("pow,");
                    w.write(Long.toHexString(Double.doubleToRawLongBits(x)));
                    w.write(':');
                    w.write(Long.toHexString(Double.doubleToRawLongBits(y)));
                    w.write(',');
                    w.write(Long.toHexString(Double.doubleToRawLongBits(Math.pow(x, y))));
                    w.write(',');
                    w.write(Long.toHexString(Double.doubleToRawLongBits(StrictMath.pow(x, y))));
                    w.write(',');
                    w.write(Long.toHexString(Double.doubleToRawLongBits(FastMath.pow(x, y))));
                    w.write('\n');
                }
            }
        }
        System.out.println("wrote " + path + " for " + xs.length + " inputs");
    }
}
