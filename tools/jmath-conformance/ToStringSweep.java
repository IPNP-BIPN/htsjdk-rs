/*
 * Dumps Double.toString and Float.toString over the families decision 0046 measured, so CI can
 * check jmath::to_string against the pinned JVM on every push.
 *
 * Output: d|f <TAB> <hex bits> <TAB> <toString>
 *
 * Stratified rather than uniform, because each rule in decision 0046 lives where a uniform sample
 * rarely lands: integers above 2^52, powers of two, short decimals on a rounding boundary, the
 * overflow at binary exponent 84 (85 for float), and the smallest subnormals.
 */

import java.util.SplittableRandom;

public class ToStringSweep {

    static final StringBuilder out = new StringBuilder(1 << 20);

    static void d(final double v) {
        if (Double.isNaN(v) || Double.isInfinite(v)) {
            return;
        }
        out.append("d\t").append(Long.toHexString(Double.doubleToRawLongBits(v))).append('\t')
                .append(Double.toString(v)).append('\n');
        flush();
    }

    static void f(final float v) {
        if (Float.isNaN(v) || Float.isInfinite(v)) {
            return;
        }
        out.append("f\t").append(Integer.toHexString(Float.floatToRawIntBits(v))).append('\t')
                .append(Float.toString(v)).append('\n');
        flush();
    }

    static void flush() {
        if (out.length() > (1 << 19)) {
            System.out.print(out);
            out.setLength(0);
        }
    }

    public static void main(final String[] args) {
        // Fixed seed: the sweep must be the same file on every run.
        final SplittableRandom r = new SplittableRandom(20261013L);

        // Uniform bit patterns.
        for (int i = 0; i < 500_000; i++) {
            d(Double.longBitsToDouble(r.nextLong() & 0x7fffffffffffffffL));
        }
        // Short decimals at every decimal exponent, and their neighbours: the exclusive-bound and
        // tie rules.
        for (int i = 0; i < 400_000; i++) {
            final int k = 1 + r.nextInt(17);
            long m = 0;
            for (int j = 0; j < k; j++) {
                m = m * 10 + r.nextInt(10);
            }
            final double b = Double.parseDouble(m + "e" + r.nextInt(-330, 310));
            d(b);
            d(Double.longBitsToDouble(Double.doubleToRawLongBits(b) + r.nextInt(-3, 4)));
        }
        // Few significant bits at every binary exponent: powers of two and the narrow interval.
        for (int e = -1074; e <= 1023; e++) {
            for (long odd = 1; odd < 128; odd += 2) {
                final double v = Math.scalb((double) odd, e);
                if (v != 0 && Math.scalb(v, -e) == (double) odd) {
                    d(v);
                }
            }
        }
        // The smallest subnormals, where the six exceptions are.
        for (long bits = 1; bits <= 200_000; bits++) {
            d(Double.longBitsToDouble(bits));
        }
        // Integers above 2^52, and every exponent up to 200 where the overflow could reach.
        for (int i = 0; i < 300_000; i++) {
            final long be = r.nextInt(52, 201) + 1023;
            d(Double.longBitsToDouble((be << 52) | (r.nextLong() & 0xfffffffffffffL)));
        }
        // Binary exponent 84, where the high-side test overflows.
        for (int i = 0; i < 200_000; i++) {
            d(Double.longBitsToDouble((84L + 1023) << 52 | (r.nextLong() & 0xfffffffffffffL)));
        }

        // Floats: uniform, the exact-integer range, exponents 82 and 85, and small subnormals.
        for (int i = 0; i < 500_000; i++) {
            f(Float.intBitsToFloat(r.nextInt() & 0x7fffffff));
        }
        for (int i = 0; i < 200_000; i++) {
            final int be = r.nextInt(23, 90) + 127;
            f(Float.intBitsToFloat((be << 23) | (r.nextInt() & 0x7fffff)));
        }
        for (final int exponent : new int[] {82, 85}) {
            for (int i = 0; i < 100_000; i++) {
                f(Float.intBitsToFloat(((exponent + 127) << 23) | (r.nextInt() & 0x7fffff)));
            }
        }
        for (int bits = 1; bits <= 100_000; bits++) {
            f(Float.intBitsToFloat(bits));
        }
        System.out.print(out);
        System.out.flush();
    }
}
