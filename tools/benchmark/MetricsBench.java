import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.metrics.StringHeader;
import htsjdk.samtools.util.Histogram;
import java.io.File;
import java.util.Locale;

/**
 * The reference side of `metrics-bench`: htsjdk's `MetricsFile` writing the same beans and
 * histograms.
 *
 * Usage: MetricsBench dir reps
 */
public class MetricsBench {
  public static class BenchMetric extends MetricBase {
    public long READS;
    public double MEAN;
    public String NAME;
    public boolean PASS;
  }

  static final int ROWS = 200_000;
  static final int BINS = 1_000;

  static Histogram<Integer> histogram(String valueLabel, double scale) {
    Histogram<Integer> h = new Histogram<>("insert_size", valueLabel);
    for (int k = 0; k < BINS; k++) h.increment(k, Math.floor(k * scale));
    return h;
  }

  static void write(File out) {
    MetricsFile<BenchMetric, Integer> file = new MetricsFile<>();
    file.addHeader(new StringHeader("benchmark"));
    for (int i = 0; i < ROWS; i++) {
      BenchMetric m = new BenchMetric();
      m.READS = i * 37L;
      m.MEAN = i / 7.0;
      m.NAME = "sample_" + i;
      m.PASS = i % 3 == 0;
      file.addMetric(m);
    }
    file.addHistogram(histogram("count", 1.5));
    file.addHistogram(histogram("count_dups", 0.25));
    file.write(out);
  }

  public static void main(String[] args) {
    File dir = new File(args[0]);
    int reps = args.length > 1 ? Integer.parseInt(args[1]) : 3;
    File out = new File(dir, "java_metrics_write.txt");
    write(out);
    for (int r = 0; r < reps; r++) {
      long t0 = System.nanoTime();
      write(out);
      System.out.printf(Locale.ROOT, "java_metrics_write_run%d_secs=%.6f%n", r, (System.nanoTime() - t0) / 1e9);
    }
  }
}
