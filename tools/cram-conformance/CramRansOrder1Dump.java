/*
 * rANS 4x8 order 1: the same codec with a context, and a different arithmetic underneath.
 *
 * cram-rans pinned order 0. Order 1 is not order 0 with more tables. Five things change, and three
 * of them are not in the specification.
 *
 *   - THE FREQUENCY TABLE COUNTS THREE BYTES THAT ARE NOT BIGRAMS. Each of the four lanes starts
 *     with context 0, and only the first lane's first byte actually follows nothing. So
 *     calcFrequenciesOrder1 adds F[0][x]++ for the byte at each quarter boundary, three counts the
 *     input does not contain, and adds 3 to T[0]. A port that counts the bigrams it can see
 *     produces a different table on every input;
 *   - THE NORMALISATION IS FLOATING POINT, where order 0's is fixed point. p = 4096.0 / T[i], then
 *     F[i][j] *= p, which in Java is a compound assignment on an int and so TRUNCATES toward zero.
 *     Two normalisations, two arithmetics, in the same class;
 *   - THE FREQUENCY TABLE HAS TWO LEVELS OF RUN LENGTH, one over contexts and one over symbols
 *     within a context, each with the same rule as order 0: the run byte appears only once the
 *     previous entry was also present, and the reader infers it by peeking;
 *   - A FREQUENCY BYTE OF ZERO MEANS 4096 ON THE WAY IN. readStatsOrder1 has
 *     `if (D[i].frequencies[j] == 0) D[i].frequencies[j] = Constants.TOTAL_FREQ;` and order 0 has
 *     no such line. The writer never emits it, so this is a reader-only rule, and a port that
 *     treats a zero as a zero disagrees with htsjdk on a stream htsjdk accepts;
 *   - THE FOUR LANES READ FOUR QUARTERS, not four interleaved positions. Order 0's lanes take
 *     every fourth byte; order 1's take a contiguous quarter each, which is why the contexts have
 *     to be seeded and why the remainder belongs to the last lane alone.
 *
 * Output:
 *
 *     in\t<label>\t<length>\t<sha256 of the input>
 *     enc\t<label>\t<output length>\t<sha256 of the output>
 *     bytes\t<label>\t<the whole output, hex>          (only when it fits)
 *     prefix\t<label>\t<order byte>\t<compressed size field>\t<raw size field>
 *     quarter\t<label>\t<isz4>\t<q1>\t<q2>\t<q3>\t<symbol at q1>\t<symbol at q2>\t<symbol at q3>
 *     ctxtotal\t<label>\t<context>\t<raw total including the three extra counts>
 *     norm\t<label>\t<context>\t<symbol>\t<raw count>\t<normalised frequency>
 *     normdigest\t<label>\t<contexts used>\t<sha256 of the whole normalised table>
 *     freqtab\t<label>\t<size>\t<the frequency table's bytes, hex>
 *     states\t<label>\t<rans0>\t<rans1>\t<rans2>\t<rans3>
 *     roundtrip\t<label>\t<ok|MISMATCH>
 *     zerofreq\t<context>\t<symbol>\t<what the reader makes of a zero frequency byte>
 *     err\t<label>\t<class>\t<message>
 *
 * Usage: CramRansOrder1Dump
 */

import htsjdk.samtools.cram.compression.CompressionUtils;
import htsjdk.samtools.cram.compression.rans.ArithmeticDecoder;
import htsjdk.samtools.cram.compression.rans.RANSParams;
import htsjdk.samtools.cram.compression.rans.rans4x8.RANS4x8Decode;
import htsjdk.samtools.cram.compression.rans.rans4x8.RANS4x8Encode;
import htsjdk.samtools.cram.compression.rans.rans4x8.RANS4x8Params;

import java.lang.reflect.Method;
import java.nio.ByteBuffer;
import java.security.MessageDigest;
import java.util.Arrays;
import java.util.TreeSet;

public class CramRansOrder1Dump {

    /** Above this the whole output is not worth a line; the digest still fails on one wrong byte. */
    private static final int MAX_INLINE_BYTES = 4096;
    /** Above this alphabet size the per-symbol table is a digest instead of tens of thousands of
     *  rows. The small-alphabet inputs are chosen so the table itself is readable. */
    private static final int MAX_INLINE_ALPHABET = 12;

    private static final int NUMBER_OF_SYMBOLS = 256;
    private static final int PREFIX_LENGTH = 9;

    public static void main(final String[] args) throws Exception {
        System.out.println("# CramRansOrder1Dump: rANS 4x8 order 1, the context and its arithmetic");

        // Four bytes is the smallest input order 1 is permitted on, and it is the degenerate case:
        // isz4 is 1, so the lane cursors start at -1 and the main loop never runs.
        emit("four-bytes", "ACGT".getBytes());
        emit("five-bytes", "ACGTA".getBytes());
        emit("six-bytes", "ACGTAC".getBytes());
        emit("seven-bytes", "ACGTACG".getBytes());
        emit("eight-bytes", "ACGTACGT".getBytes());

        emit("acgt-1000", repeat("ACGT", 1000));
        emit("acgt-1001", repeat("ACGT", 1001));
        emit("acgt-1002", repeat("ACGT", 1002));
        emit("acgt-1003", repeat("ACGT", 1003));

        // One symbol: every context that exists is the same context, and it takes the whole 4096.
        emit("uniform-1000", repeat("A", 1000));

        // Two symbols, lopsided, so a normalised frequency needs its two-byte form.
        final byte[] two = new byte[1000];
        Arrays.fill(two, 0, 900, (byte) 'A');
        Arrays.fill(two, 900, 1000, (byte) 'B');
        emit("two-symbols", two);

        // A contiguous block of symbols, which is what both levels of run length are for.
        emit("contiguous-run", spread(5, 13, 100));

        // Something order 1 actually pays for: a repeating motif, where the next base is nearly
        // determined by the last one.
        emit("motif-1000", repeat("ACGTACGTAA", 1000));

        // The shapes rANS meets in a CRAM, and the shape it does worst on.
        final byte[] quals = new byte[5000];
        long seed = 0x5DEECE66DL;
        for (int i = 0; i < quals.length; i++) {
            seed = seed * 6364136223846793005L + 1442695040888963407L;
            quals[i] = (byte) (30 + ((seed >>> 33) % 8));
        }
        emit("quality-band", quals);

        final byte[] noise = new byte[10000];
        seed = 0x1234_5678L;
        for (int i = 0; i < noise.length; i++) {
            seed = seed * 6364136223846793005L + 1442695040888963407L;
            noise[i] = (byte) (seed >>> 33);
        }
        emit("noise-10000", noise);

        final byte[] once = new byte[256];
        for (int i = 0; i < 256; i++) {
            once[i] = (byte) i;
        }
        emit("all-256-once", once);

        // A frequency byte of zero, which the writer never emits and the reader reads as 4096.
        // The table is built by hand because no input can produce one.
        //   context 0, symbol 65, frequency byte 0, end of row, end of table.
        emitZeroFrequency(new byte[] {0x00, 0x41, 0x00, 0x00, 0x00});
    }

    static byte[] repeat(final String source, final int length) {
        final byte[] pattern = source.getBytes();
        final byte[] out = new byte[length];
        for (int i = 0; i < length; i++) {
            out[i] = pattern[i % pattern.length];
        }
        return out;
    }

    static byte[] spread(final int from, final int to, final int count) {
        final byte[] out = new byte[(to - from) * count];
        int at = 0;
        for (int n = 0; n < count; n++) {
            for (int symbol = from; symbol < to; symbol++) {
                out[at++] = (byte) symbol;
            }
        }
        return out;
    }

    static void emit(final String label, final byte[] input) throws Exception {
        System.out.printf("in\t%s\t%d\t%s%n", label, input.length, sha256(input));

        final ByteBuffer compressed = new RANS4x8Encode()
                .compress(CompressionUtils.wrap(input), new RANS4x8Params(RANSParams.ORDER.ONE));
        final byte[] out = new byte[compressed.remaining()];
        compressed.duplicate().get(out);

        System.out.printf("enc\t%s\t%d\t%s%n", label, out.length, sha256(out));
        if (out.length > 0 && out.length <= MAX_INLINE_BYTES) {
            System.out.printf("bytes\t%s\t%s%n", label, hex(out));
        }

        final ByteBuffer prefix = CompressionUtils.wrap(out);
        System.out.printf("prefix\t%s\t%d\t%d\t%d%n", label, prefix.get() & 0xFF, prefix.getInt(),
                prefix.getInt());

        // The three positions whose bytes are counted as if they followed nothing, because three of
        // the four lanes begin at them with context 0.
        final int isz4 = input.length >> 2;
        System.out.printf("quarter\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d%n", label, isz4,
                isz4, 2 * isz4, 3 * isz4,
                0xFF & input[isz4], 0xFF & input[2 * isz4], 0xFF & input[3 * isz4]);

        // The raw counts, as the reference computes them: the bigrams plus those three.
        final int[][] raw = new int[NUMBER_OF_SYMBOLS][NUMBER_OF_SYMBOLS];
        final int[] totals = new int[NUMBER_OF_SYMBOLS];
        int last = 0;
        for (final byte b : input) {
            raw[last][0xFF & b]++;
            totals[last]++;
            last = 0xFF & b;
        }
        raw[0][0xFF & input[isz4]]++;
        raw[0][0xFF & input[2 * isz4]]++;
        raw[0][0xFF & input[3 * isz4]]++;
        totals[0] += 3;

        final int[][] normalised = calcFrequenciesOrder1(input);

        final TreeSet<Integer> alphabet = new TreeSet<>();
        for (final byte b : input) {
            alphabet.add(0xFF & b);
        }
        alphabet.add(0);

        int contexts = 0;
        for (int i = 0; i < NUMBER_OF_SYMBOLS; i++) {
            if (totals[i] == 0) {
                continue;
            }
            contexts++;
            System.out.printf("ctxtotal\t%s\t%d\t%d%n", label, i, totals[i]);
        }
        if (alphabet.size() <= MAX_INLINE_ALPHABET) {
            for (int i = 0; i < NUMBER_OF_SYMBOLS; i++) {
                for (int j = 0; j < NUMBER_OF_SYMBOLS; j++) {
                    if (normalised[i][j] != 0) {
                        System.out.printf("norm\t%s\t%d\t%d\t%d\t%d%n", label, i, j, raw[i][j],
                                normalised[i][j]);
                    }
                }
            }
        }
        System.out.printf("normdigest\t%s\t%d\t%s%n", label, contexts, digestOf(normalised));

        final ByteBuffer table = CompressionUtils.allocateByteBuffer(NUMBER_OF_SYMBOLS
                * (NUMBER_OF_SYMBOLS * 3 + 3) + 1);
        final int tableSize = writeFrequenciesOrder1(table, deepCopy(normalised));
        final byte[] tableBytes = new byte[tableSize];
        table.rewind();
        table.get(tableBytes);
        System.out.printf("freqtab\t%s\t%d\t%s%n", label, tableSize,
                tableSize <= MAX_INLINE_BYTES ? hex(tableBytes) : sha256(tableBytes));

        final ByteBuffer blob = CompressionUtils.wrap(out);
        blob.position(PREFIX_LENGTH + tableSize);
        System.out.printf("states\t%s\t%d\t%d\t%d\t%d%n", label,
                0xFFFFFFFFL & blob.getInt(), 0xFFFFFFFFL & blob.getInt(),
                0xFFFFFFFFL & blob.getInt(), 0xFFFFFFFFL & blob.getInt());

        final ByteBuffer back = new RANS4x8Decode().uncompress(CompressionUtils.wrap(out));
        final byte[] decoded = new byte[back.remaining()];
        back.get(decoded);
        System.out.printf("roundtrip\t%s\t%s%n", label,
                Arrays.equals(decoded, input) ? "ok" : "MISMATCH");
    }

    /** What the reader makes of a frequency table the writer would never produce. */
    static void emitZeroFrequency(final byte[] handBuiltTable) throws Exception {
        final RANS4x8Decode decode = new RANS4x8Decode();
        final Method initialize =
                Class.forName("htsjdk.samtools.cram.compression.rans.RANSDecode")
                        .getDeclaredMethod("initializeRANSDecoder");
        initialize.setAccessible(true);
        initialize.invoke(decode);

        final Method readStats =
                RANS4x8Decode.class.getDeclaredMethod("readStatsOrder1", ByteBuffer.class);
        readStats.setAccessible(true);
        readStats.invoke(decode, CompressionUtils.wrap(handBuiltTable));

        final Method getD = Class.forName("htsjdk.samtools.cram.compression.rans.RANSDecode")
                .getDeclaredMethod("getD");
        getD.setAccessible(true);
        final ArithmeticDecoder[] decoders = (ArithmeticDecoder[]) getD.invoke(decode);
        System.out.printf("zerofreq\t%d\t%d\t%d%n", 0, 0x41, decoders[0].frequencies[0x41]);
    }

    /** `RANS4x8Encode.calcFrequenciesOrder1`, which is private and is the whole of the arithmetic. */
    static int[][] calcFrequenciesOrder1(final byte[] input) throws Exception {
        final Method method =
                RANS4x8Encode.class.getDeclaredMethod("calcFrequenciesOrder1", ByteBuffer.class);
        method.setAccessible(true);
        return (int[][]) method.invoke(null, CompressionUtils.wrap(input));
    }

    /** `RANS4x8Encode.writeFrequenciesOrder1`, likewise private, and likewise the point. */
    static int writeFrequenciesOrder1(final ByteBuffer cp, final int[][] frequencies)
            throws Exception {
        final Method method = RANS4x8Encode.class
                .getDeclaredMethod("writeFrequenciesOrder1", ByteBuffer.class, int[][].class);
        method.setAccessible(true);
        return (int) method.invoke(null, cp, (Object) frequencies);
    }

    static int[][] deepCopy(final int[][] source) {
        final int[][] copy = new int[source.length][];
        for (int i = 0; i < source.length; i++) {
            copy[i] = source[i].clone();
        }
        return copy;
    }

    /** The whole 256 by 256 table as one digest, for the inputs whose alphabet is too wide to
     *  print. It fails on one wrong frequency anywhere in it. */
    static String digestOf(final int[][] table) throws Exception {
        final StringBuilder text = new StringBuilder();
        for (int i = 0; i < NUMBER_OF_SYMBOLS; i++) {
            for (int j = 0; j < NUMBER_OF_SYMBOLS; j++) {
                if (table[i][j] != 0) {
                    text.append(i).append(':').append(j).append('=').append(table[i][j]).append(';');
                }
            }
        }
        return sha256(text.toString().getBytes());
    }

    static String hex(final byte[] bytes) {
        final StringBuilder b = new StringBuilder(bytes.length * 2);
        for (final byte value : bytes) {
            b.append(String.format("%02x", value));
        }
        return b.length() == 0 ? "-" : b.toString();
    }

    static String sha256(final byte[] bytes) throws Exception {
        return hex(MessageDigest.getInstance("SHA-256").digest(bytes));
    }
}
