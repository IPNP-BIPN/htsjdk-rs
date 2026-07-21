/*
 * Measures the agreement between the JVM's %f / %e and a model of them.
 *
 * VCFEncoder routes every double through String.format(Locale.US, "%.2f" | "%.3f" | "%.3e", d).
 * Those resolve into java.util.Formatter and FloatingDecimal, which decision 0014 puts out of
 * reach for transcription. crates/htsjdk-vcf/src/jformat.rs therefore implements a *stated
 * contract* instead: round the shortest round-tripping decimal representation half-up at the
 * requested precision. This harness produces the evidence for or against that contract.
 *
 * Output: <hex bits of the double> <TAB> %.2f <TAB> %.3f <TAB> %.3e
 *
 * The sample is stratified rather than uniform because the interesting doubles are not where a
 * uniform sample lands: exact ties, values just below a rounding boundary, and the magnitudes
 * VCF actually carries (allele frequencies, phred qualities, likelihoods).
 */

import java.util.*;

public class JavaFormatSweep {

    static void emit(double d) {
        System.out.println(Long.toHexString(Double.doubleToRawLongBits(d))
                + "\t" + String.format(Locale.US, "%.2f", d)
                + "\t" + String.format(Locale.US, "%.3f", d)
                + "\t" + String.format(Locale.US, "%.3e", d));
    }

    public static void main(String[] args) {
        // Fixed seed: the corpus must be the same file on every run.
        Random rng = new Random(20260721L);

        // Exact ties at two and three decimals, where half-up is observable.
        for (int i = 0; i < 2000; i++) {
            emit(i / 200.0);   // ...005, ...015, exactly representable ties included
            emit(i / 2000.0);
            emit(-i / 200.0);
        }

        // Phred qualities and allele frequencies: the ranges VCF really carries.
        for (int i = 0; i < 20000; i++) {
            emit(rng.nextDouble() * 100.0);
            emit(rng.nextDouble());
            emit(rng.nextDouble() * 1e-3);
        }

        // Uniformly random bit patterns, finite only. This is where a shortest-representation
        // disagreement is most likely to surface, since it is not concentrated near nice values.
        int emitted = 0;
        while (emitted < 60000) {
            double d = Double.longBitsToDouble(rng.nextLong());
            if (Double.isFinite(d)) {
                emit(d);
                emitted++;
            }
        }

        // Powers of ten and their neighbours, where digit-string carries happen.
        for (int e = -300; e <= 300; e++) {
            double p = Math.pow(10, e);
            emit(p);
            emit(Math.nextUp(p));
            emit(Math.nextDown(p));
        }
    }
}
