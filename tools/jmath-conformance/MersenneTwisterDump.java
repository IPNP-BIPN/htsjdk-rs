/*
 * commons-math3's MersenneTwister and RandomDataGenerator.nextPermutation, taken from 3.5.
 *
 * A fingerprint metric randomizes: `CalculateFingerprintMetrics` builds a MersenneTwister on the
 * fixed seed 42 and permutes each site's likelihoods a hundred times. Nothing about that is
 * approximate, so the port has to produce the SAME stream, and the stream is what this measures.
 *
 * Nine behaviours this is built to catch.
 *
 *   - THE SEEDING IS NOT THE REFERENCE ALGORITHM'S: commons-math3 seeds from an int with Knuth's
 *     multiplier 1812433253 and then, unlike the original, does NOT run the array-seeding step, so
 *     the stream differs from a textbook MT19937 seeded the same way;
 *   - `nextInt()` RETURNS A SIGNED INT, the tempered word reinterpreted, so half the values are
 *     negative;
 *   - `nextInt(n)` IS NOT `nextInt() % n`: it rejects and redraws, and for a power of two it
 *     takes the HIGH bits instead;
 *   - `nextDouble()` IS BUILT FROM TWO WORDS, 26 bits and 27, which is not `nextInt()` scaled;
 *   - `nextLong()` IS TWO WORDS, the first shifted left by 32, so it consumes two draws;
 *   - `nextBoolean()` AND `nextFloat()` EACH CONSUME A WHOLE WORD, so mixing them shifts
 *     everything after them;
 *   - THE GENERATOR IS STATEFUL AND THE ORDER OF CALLS IS PART OF THE ANSWER, which is why the
 *     dump records a mixed sequence as well as the pure ones;
 *   - `RandomDataGenerator.nextPermutation(n, k)` IS A PARTIAL FISHER-YATES over `natural(n)`,
 *     taking the FIRST k of the shuffled array, and its refusals name the two bounds;
 *   - TWO SEEDS THAT DIFFER GIVE STREAMS THAT DO NOT: seeds 0 and -1 produce the same twelve
 *     values from the second draw on, differing only in the first. A port that checks itself on
 *     one seed would not notice getting the seeding wrong;
 *   - AND `MathUtil.permute` READS THE PERMUTATION AS AN INDEX MAP, `out[i] = in[perm[i]]`, which
 *     is the direction that decides where a value lands.
 *
 * Output:
 *
 *     ints\t<seed>\t<the first values of nextInt(), comma separated>
 *     bounded\t<seed>/<bound>\t<the first values of nextInt(bound)>
 *     doubles\t<seed>\t<the first values of nextDouble(), as Java renders them>
 *     longs\t<seed>\t<the first values of nextLong()>
 *     booleans\t<seed>\t<the first values of nextBoolean()>
 *     floats\t<seed>\t<the first values of nextFloat()>
 *     mixed\t<seed>\t<one of each, in a fixed order>
 *     permutation\t<seed>/<n>/<k>\t<the permutation, comma separated>
 *     permuted\t<seed>\t<MathUtil.permute's own answer over a known array>
 *     error\t<case>\t<exception class>:<message>
 *
 * Usage: MersenneTwisterDump
 */

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomDataGenerator;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

public class MersenneTwisterDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    static final int[] SEEDS = {42, 0, 1, -1, 2147483647};
    static final int DRAWS = 12;

    public static void main(final String[] args) {
        for (final int seed : SEEDS) {
            final String name = Integer.toString(seed);

            final MersenneTwister ints = new MersenneTwister(seed);
            final List<String> values = new ArrayList<>();
            for (int i = 0; i < DRAWS; i++) {
                values.add(Integer.toString(ints.nextInt()));
            }
            emit("ints", name, String.join(",", values));

            // A power of two and a value that is not one, which take two different paths.
            for (final int bound : new int[]{2, 7, 64, 1000}) {
                final MersenneTwister bounded = new MersenneTwister(seed);
                final List<String> drawn = new ArrayList<>();
                for (int i = 0; i < DRAWS; i++) {
                    drawn.add(Integer.toString(bounded.nextInt(bound)));
                }
                emit("bounded", name + "/" + bound, String.join(",", drawn));
            }

            final MersenneTwister doubles = new MersenneTwister(seed);
            final List<String> asDoubles = new ArrayList<>();
            for (int i = 0; i < DRAWS; i++) {
                asDoubles.add(Double.toString(doubles.nextDouble()));
            }
            emit("doubles", name, String.join(",", asDoubles));

            final MersenneTwister longs = new MersenneTwister(seed);
            final List<String> asLongs = new ArrayList<>();
            for (int i = 0; i < DRAWS; i++) {
                asLongs.add(Long.toString(longs.nextLong()));
            }
            emit("longs", name, String.join(",", asLongs));

            final MersenneTwister booleans = new MersenneTwister(seed);
            final List<String> asBooleans = new ArrayList<>();
            for (int i = 0; i < DRAWS; i++) {
                asBooleans.add(Boolean.toString(booleans.nextBoolean()));
            }
            emit("booleans", name, String.join(",", asBooleans));

            final MersenneTwister floats = new MersenneTwister(seed);
            final List<String> asFloats = new ArrayList<>();
            for (int i = 0; i < DRAWS; i++) {
                asFloats.add(Float.toString(floats.nextFloat()));
            }
            emit("floats", name, String.join(",", asFloats));

            // The order of calls is part of the answer: each of these consumes its own words.
            final MersenneTwister mixed = new MersenneTwister(seed);
            final List<String> sequence = new ArrayList<>();
            sequence.add(Integer.toString(mixed.nextInt()));
            sequence.add(Boolean.toString(mixed.nextBoolean()));
            sequence.add(Double.toString(mixed.nextDouble()));
            sequence.add(Integer.toString(mixed.nextInt(7)));
            sequence.add(Long.toString(mixed.nextLong()));
            sequence.add(Float.toString(mixed.nextFloat()));
            sequence.add(Integer.toString(mixed.nextInt()));
            emit("mixed", name, String.join(",", sequence));

            // The permutation, which is what a fingerprint's likelihoods go through.
            for (final int[] pair : new int[][]{{3, 3}, {4, 4}, {10, 10}, {10, 3}, {1, 1}}) {
                final RandomDataGenerator rdg =
                        new RandomDataGenerator(new MersenneTwister(seed));
                final int[] permutation = rdg.nextPermutation(pair[0], pair[1]);
                emit("permutation", name + "/" + pair[0] + "/" + pair[1],
                        java.util.Arrays.stream(permutation).mapToObj(Integer::toString)
                                .collect(Collectors.joining(",")));
            }

            // And the direction a caller reads that permutation in. Picard's `MathUtil.permute`
            // is `out[i] = in[perm[i]]`, an INDEX MAP rather than a destination map, and the two
            // give different arrays for the same permutation. The picard jar is not on this
            // oracle's class path, so the loop is written out here rather than called; the source
            // it transcribes is `picard.util.MathUtil.permute`.
            final RandomDataGenerator rdg = new RandomDataGenerator(new MersenneTwister(seed));
            final double[] array = {10.0, 20.0, 30.0, 40.0};
            final int[] permutation = rdg.nextPermutation(array.length, array.length);
            final double[] permuted = new double[array.length];
            for (int i = 0; i < array.length; i++) {
                permuted[i] = array[permutation[i]];
            }
            emit("permuted", name, java.util.Arrays.stream(permuted)
                    .mapToObj(Double::toString).collect(Collectors.joining(",")));
        }

        // The two refusals, which name the bounds they were given.
        for (final int[] pair : new int[][]{{3, 5}, {0, 0}, {5, 0}}) {
            final RandomDataGenerator rdg = new RandomDataGenerator(new MersenneTwister(42));
            try {
                rdg.nextPermutation(pair[0], pair[1]);
                emit("error", pair[0] + "/" + pair[1], "ok");
            } catch (final Exception e) {
                emit("error", pair[0] + "/" + pair[1],
                        e.getClass().getName() + ":" + e.getMessage());
            }
        }

        System.out.print(buf);
    }
}
