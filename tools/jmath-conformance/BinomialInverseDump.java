/*
 * BinomialDistribution.inverseCumulativeProbability, taken from commons-math3 3.5.
 *
 * The quantile of a binomial, which GATK reaches through
 * `PowerCalculationUtils.calculateMinCountForSignal`: the count a validation pileup must reach
 * before it is signal rather than noise. It is not a formula. It is a bracket narrowed by a
 * one-sided Chebyshev inequality and then bisected, so WHERE THE BRACKET STARTS is part of the
 * answer whenever the search lands on a boundary.
 *
 * Six behaviours this is built to catch.
 *
 *   - THE ANSWER IS THE SMALLEST x WITH cdf(x) >= p, found by bisection, and the bisection returns
 *     `upper` rather than `lower`, so the invariant it maintains is `cdf(lower) < p <= cdf(upper)`;
 *   - THE LOWER BOUND IS DECREMENTED BEFORE THE SEARCH, which is what makes that invariant true at
 *     the start and is why `p` just above zero answers zero rather than minus one;
 *   - THE CHEBYSHEV NARROWING USES ceil AND THEN SUBTRACTS ONE on both ends, so the bracket is
 *     asymmetric and can exclude the true answer's neighbour without excluding the answer;
 *   - IT IS SKIPPED ENTIRELY when the variance is zero, which for a binomial means `p` of exactly
 *     zero or exactly one -- the degenerate cases take the unnarrowed bracket;
 *   - p = 0 RETURNS THE SUPPORT'S LOWER BOUND AND p = 1 ITS UPPER, both before any search;
 *   - AND A p OUTSIDE [0, 1] IS AN OutOfRangeException -- WHILE A NaN IS NOT. `NaN < 0` and
 *     `NaN > 1` are both false, so the range check passes; the two equality tests fail; the
 *     Chebyshev arithmetic produces NaNs, no comparison against them narrows anything, and the
 *     bisection returns the support's upper bound. A NaN quantile therefore answers the number of
 *     trials rather than throwing, which is the sort of thing only a measurement finds.
 *
 * Output:
 *
 *     inverse\t<trials>,<probability bits>,<p bits>=<answer>
 *     cumulative\t<trials>,<probability bits>,<x>=<bits>
 *     error\t<label>\t<exception class>:<message>
 *
 * Usage: BinomialInverseDump
 */

import org.apache.commons.math3.distribution.BinomialDistribution;

public class BinomialInverseDump {

    /** The probabilities GATK's noise model reaches, plus the two degenerate ends. */
    static final double[] PROBABILITIES = {0.0, 1e-6, 0.001, 0.01, 0.05, 0.1, 0.5, 0.9, 0.99, 1.0};

    /** The quantiles asked for, including the 0.99 `P_VALUE_FOR_NOISE` GATK uses. */
    static final double[] QUANTILES = {0.0, 1e-9, 0.01, 0.25, 0.5, 0.75, 0.99, 0.999999, 1.0};

    /** Trial counts from empty to a few hundred reads. */
    static final int[] TRIALS = {0, 1, 2, 5, 10, 30, 100, 317};

    public static void main(final String[] args) {
        System.out.println("# BinomialInverseDump: the quantile GATK's noise floor is read off");

        for (final int trials : TRIALS) {
            for (final double probability : PROBABILITIES) {
                for (final double quantile : QUANTILES) {
                    final BinomialDistribution distribution =
                            new BinomialDistribution(null, trials, probability);
                    System.out.printf("inverse\t%d,%016x,%016x=%d%n", trials,
                            Double.doubleToRawLongBits(probability),
                            Double.doubleToRawLongBits(quantile),
                            distribution.inverseCumulativeProbability(quantile));
                }
                // The cumulative probabilities the search is bisecting over, so a port that
                // disagrees can be told apart from one whose search is wrong.
                for (int x = -1; x <= Math.min(trials, 5); x++) {
                    System.out.printf("cumulative\t%d,%016x,%d=%016x%n", trials,
                            Double.doubleToRawLongBits(probability), x,
                            Double.doubleToRawLongBits(
                                    new BinomialDistribution(null, trials, probability)
                                            .cumulativeProbability(x)));
                }
            }
        }

        // The mean and the variance, which decide whether the Chebyshev narrowing applies at all.
        for (final int trials : TRIALS) {
            for (final double probability : PROBABILITIES) {
                final BinomialDistribution distribution =
                        new BinomialDistribution(null, trials, probability);
                System.out.printf("moments\t%d,%016x=%016x,%016x%n", trials,
                        Double.doubleToRawLongBits(probability),
                        Double.doubleToRawLongBits(distribution.getNumericalMean()),
                        Double.doubleToRawLongBits(distribution.getNumericalVariance()));
            }
        }

        // And the refusals.
        error("below-zero", -0.1);
        error("above-one", 1.1);
        // Not an error at all, as it turns out: the row is printed as an answer.
        error("nan", Double.NaN);
    }

    static void error(final String label, final double quantile) {
        try {
            System.out.printf("inverse\tunexpected-%s=%d%n", label,
                    new BinomialDistribution(null, 10, 0.5).inverseCumulativeProbability(quantile));
        } catch (final Exception e) {
            System.out.printf("error\t%s\t%s:%s%n", label, e.getClass().getName(), e.getMessage());
        }
    }
}
