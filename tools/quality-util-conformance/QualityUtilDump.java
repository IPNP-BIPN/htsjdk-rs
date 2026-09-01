import htsjdk.samtools.util.QualityUtil;

/**
 * Dumps every answer `htsjdk.samtools.util.QualityUtil` gives, as bit patterns rather than as
 * decimals, so the comparison is the double and not its rendering.
 *
 * Three families:
 *
 *   table  <score>            <bits>   the static error-probability table, all 101 entries
 *   phred  <probability bits> <int>    getPhredScoreFromErrorProbability, over a sweep that
 *                                      includes the values a metrics file actually reaches and
 *                                      the ones that break a naive port: zero, one, above one,
 *                                      and the double just below a half after scaling
 *   obs    <obs bits> <err bits> <int> getPhredScoreFromObsAndErrors, whose argument is a ratio
 *
 * The reference's own `Math.pow` builds the table and `Math.round` closes the other two, which is
 * why this is measured here rather than recomputed in the port.
 */
public class QualityUtilDump {
  static String bits(double d) {
    return String.format("%016x", Double.doubleToRawLongBits(d));
  }

  public static void main(String[] args) {
    for (int score = 0; score <= 100; score++) {
      System.out.printf("table\t%d\t%s%n", score, bits(QualityUtil.getErrorProbabilityFromPhredScore(score)));
    }

    double[] probabilities = {
      0.0, 1.0, 0.5, 0.1, 0.01, 0.001, 1e-10, 1e-100, 1e-300,
      // Above one: -10 * log10(p) is negative, and Math.round is half UP rather than half away
      // from zero, so this is where a Rust port using f64::round diverges.
      1.5, 2.0, 3.1622776601683795, 10.0, 100.0,
      // The ratios a metrics file produces: errors over observations, at the sizes tools see.
      1.0 / 3.0, 2.0 / 3.0, 1.0 / 7.0, 1.0 / 1000.0, 999.0 / 1000.0,
      // Denormal and the smallest normal, which scale to the extremes of the log.
      Double.MIN_VALUE, Double.MIN_NORMAL, Double.MAX_VALUE,
      Double.NaN, Double.POSITIVE_INFINITY,
    };
    for (double p : probabilities) {
      System.out.printf("phred\t%s\t%d%n", bits(p), QualityUtil.getPhredScoreFromErrorProbability(p));
    }

    double[][] obsAndErrors = {
      {100, 0}, {100, 1}, {100, 50}, {100, 100}, {100, 101},
      {1, 0}, {1, 1}, {0, 0}, {0, 1},
      {1e9, 1}, {1e9, 12345}, {3, 1}, {7, 2},
    };
    for (double[] pair : obsAndErrors) {
      System.out.printf(
          "obs\t%s\t%s\t%d%n",
          bits(pair[0]), bits(pair[1]), QualityUtil.getPhredScoreFromObsAndErrors(pair[0], pair[1]));
    }
  }
}
