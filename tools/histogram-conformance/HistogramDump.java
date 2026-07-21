/*
 * Dumps every Histogram statistic htsjdk computes, over a range of shapes.
 *
 * Output: <case name> <TAB> <stat> <TAB> <raw double bits in hex>
 *
 * Raw bits, not printed decimals: the point is bit-identity, and a decimal rendering would
 * hide exactly the last-bit differences that accumulation order produces.
 */

import htsjdk.samtools.util.Histogram;
import java.util.Random;

public class HistogramDump {

    static void emit(final String name, final String stat, final double v) {
        System.out.println(name + "\t" + stat + "\t" + Long.toHexString(Double.doubleToRawLongBits(v)));
    }

    static void dump(final String name, final double[][] pairs) {
        final Histogram<Double> h = new Histogram<>("bin", "count");
        for (final double[] p : pairs) h.increment(p[0], p[1]);

        emit(name, "count", h.getCount());
        emit(name, "sum", h.getSum());
        emit(name, "sumOfValues", h.getSumOfValues());
        emit(name, "mean", h.getMean());
        emit(name, "standardDeviation", h.getStandardDeviation());
        emit(name, "median", h.getMedian());
        emit(name, "medianAbsoluteDeviation", h.getMedianAbsoluteDeviation());
        emit(name, "estimateSdViaMad", h.estimateSdViaMad());
        emit(name, "meanBinSize", h.getMeanBinSize());
        emit(name, "medianBinSize", h.getMedianBinSize());
        if (h.size() > 0) {
            emit(name, "mode", h.getMode());
            emit(name, "min", h.getMin());
            emit(name, "max", h.getMax());
        }
        for (final double p : new double[]{0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99}) {
            try {
                emit(name, "percentile_" + p, h.getPercentile(p));
            } catch (final Exception e) {
                System.out.println(name + "\tpercentile_" + p + "\tERROR");
            }
        }
        for (final double v : new double[]{-1, 0, 1, 10, 100, 1000}) {
            emit(name, "cumulativeProbability_" + v, h.getCumulativeProbability(v));
        }
    }

    public static void main(final String[] args) {
        dump("single", new double[][]{{7, 1}});
        dump("two", new double[][]{{1, 1}, {3, 1}});
        dump("uniform_odd", new double[][]{{1, 1}, {2, 1}, {3, 1}});
        dump("uniform_even", new double[][]{{1, 1}, {2, 1}, {3, 1}, {4, 1}});
        dump("one_bin_many", new double[][]{{9, 4}});
        dump("weighted", new double[][]{{1, 3}, {5, 1}});
        dump("tenths", new double[][]{{0.1, 3}, {0.2, 7}, {0.3, 11}});

        // A realistic insert-size distribution: unimodal, long right tail.
        final double[][] insert = new double[600][2];
        for (int i = 0; i < insert.length; i++) {
            final double x = 100 + i;
            insert[i][0] = x;
            insert[i][1] = Math.floor(10000 * Math.exp(-Math.pow((x - 350) / 90.0, 2)));
        }
        dump("insert_size", insert);

        // Quality scores 0..45, skewed high, the shape MeanQualityByCycle sees.
        final double[][] quals = new double[46][2];
        for (int i = 0; i <= 45; i++) {
            quals[i][0] = i;
            quals[i][1] = i * i;
        }
        dump("quality_skew", quals);

        // Non-integer counts, which arise from weighted increments.
        final double[][] fractional = new double[50][2];
        final Random rng = new Random(20260721L);
        for (int i = 0; i < fractional.length; i++) {
            fractional[i][0] = i * 0.37;
            fractional[i][1] = rng.nextDouble() * 100;
        }
        dump("fractional_counts", fractional);

        // Wide magnitude range, where summation order shows up.
        final double[][] wide = new double[200][2];
        for (int i = 0; i < wide.length; i++) {
            wide[i][0] = Math.pow(10, (i - 100) / 10.0);
            wide[i][1] = 1 + (i % 7);
        }
        dump("wide_magnitudes", wide);

        // Ties in the bin values, where the mode's strict < matters.
        dump("tied_mode", new double[][]{{5, 2}, {1, 2}, {3, 2}});

        // Negative ids, which are legal.
        dump("negative_ids", new double[][]{{-5, 2}, {-1, 3}, {0, 1}, {4, 2}});
    }
}
