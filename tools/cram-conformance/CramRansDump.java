/*
 * rANS 4x8 order 0: the entropy codec an ordinary CRAM actually uses.
 *
 * cram-block measured the methods present in a four-read file and found RAW, GZIP and rANS. This is
 * the rANS, at the order every writer reaches for first.
 *
 * The stream is a nine-byte prefix, a frequency table, and a blob. Six things are arithmetic rather
 * than layout, and arithmetic is where a port diverges in silence.
 *
 *   - AN EMPTY INPUT PRODUCES ZERO BYTES, not a nine-byte prefix with an empty table. `compress`
 *     returns EMPTY_BUFFER before it does anything else, so a decoder that reads the prefix
 *     unconditionally fails on the one input it should find easiest;
 *   - THE REQUESTED ORDER IS NOT ALWAYS THE WRITTEN ORDER. Below MINIMUM_ORDER_1_SIZE = 4 bytes,
 *     `compress` uses order 0 whatever the parameters say, and the order byte records the order it
 *     used. So a caller asking for order 1 on three bytes gets a stream that says 0;
 *   - THE FOUR FINAL STATES ARE WRITTEN BIG ENDIAN AND THE WHOLE BLOB IS THEN REVERSED, so they
 *     arrive at the head of the blob in little-endian order and in the order rans0, rans1, rans2,
 *     rans3, which is the opposite of the order they were written in. Two reversals that cancel;
 *   - THE NORMALISATION IS FIXED POINT, NOT A DIVISION. tr = (4096 << 31) / T + (1 << 30) / T, then
 *     F[j] = (F[j] * tr) >> 31, with a floor of 1 for any symbol that rounds to zero. The whole
 *     rounding error is then dumped on ONE symbol: the one with the largest RAW count, ties going
 *     to the LOWEST index, and it is adjusted up or down so the total is exactly 4096;
 *   - THE FREQUENCY TABLE'S RUN LENGTHS START AT THE SECOND CONSECUTIVE SYMBOL. The run byte is
 *     written only when F[j-1] is also non-zero, so a run of two symbols writes a run byte of 0,
 *     and the decoder infers that a run byte follows by PEEKING at whether the next symbol byte is
 *     the current symbol plus one. The marker is inferred from the data, never signalled;
 *   - A FREQUENCY OF 128 OR MORE TAKES TWO BYTES, high bit set on the first. Frequencies never
 *     exceed 4096, so the second byte's top bit is free and the decoder masks it twice.
 *
 * The encoding symbol table is dumped field by field because it is where the fixed-point reciprocal
 * lives: rcpFreq = ((1 << (shift + 31)) + freq - 1) / freq with shift = ceil(log2(freq)), and
 * rcpShift carries a +32 that has been folded in so the hot loop needs one shift instead of two.
 *
 * Output:
 *
 *     in\t<label>\t<length>\t<sha256 of the input>
 *     enc\t<label>\t<output length>\t<sha256 of the output>
 *     bytes\t<label>\t<the whole output, hex>            (only when it fits)
 *     prefix\t<label>\t<order byte>\t<compressed size field>\t<raw size field>
 *     freqtab\t<label>\t<size>\t<the frequency table's bytes, hex>
 *     norm\t<label>\t<symbol>\t<raw count>\t<normalised frequency>
 *     sym\t<label>\t<symbol>\t<start>\t<freq>\t<xMax>\t<rcpFreq>\t<bias>\t<cmplFreq>\t<rcpShift>
 *     states\t<label>\t<rans0>\t<rans1>\t<rans2>\t<rans3>
 *     roundtrip\t<label>\t<ok|MISMATCH>
 *     orderused\t<length>\t<requested>\t<written>
 *     err\t<label>\t<class>\t<message>
 *
 * Usage: CramRansDump
 */

import htsjdk.samtools.cram.compression.CompressionUtils;
import htsjdk.samtools.cram.compression.rans.RANSEncodingSymbol;
import htsjdk.samtools.cram.compression.rans.RANSParams;
import htsjdk.samtools.cram.compression.rans.rans4x8.RANS4x8Decode;
import htsjdk.samtools.cram.compression.rans.rans4x8.RANS4x8Encode;
import htsjdk.samtools.cram.compression.rans.rans4x8.RANS4x8Params;

import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.nio.ByteBuffer;
import java.security.MessageDigest;
import java.util.Arrays;

public class CramRansDump {

    /** Above this the whole output is not worth a line; the digest still fails on one wrong byte. */
    private static final int MAX_INLINE_BYTES = 4096;

    private static final int TOTAL_FREQ_SHIFT = 12;
    private static final int PREFIX_LENGTH = 9;

    public static void main(final String[] args) throws Exception {
        System.out.println("# CramRansDump: rANS 4x8 order 0, the codec an ordinary CRAM uses");

        emit("empty", new byte[0]);
        emit("one-byte", new byte[] {'A'});
        emit("two-bytes", "AB".getBytes());
        emit("three-bytes", "ACG".getBytes());
        emit("four-bytes", "ACGT".getBytes());

        // The four tail lengths, which the encoder handles in a switch whose cases fall through.
        emit("acgt-1000", repeat("ACGT", 1000));
        emit("acgt-1001", repeat("ACGT", 1001));
        emit("acgt-1002", repeat("ACGT", 1002));
        emit("acgt-1003", repeat("ACGT", 1003));

        // One symbol: its normalised frequency is the whole 4096, so cmplFreq is 0.
        emit("uniform-1000", repeat("A", 1000));

        // Two symbols, lopsided, which forces a frequency over 127 into its two-byte form.
        final byte[] two = new byte[1000];
        Arrays.fill(two, 0, 900, (byte) 'A');
        Arrays.fill(two, 900, 1000, (byte) 'B');
        emit("two-symbols", two);

        // Every symbol exactly once: the normalisation lands on 16 apiece with nothing left over.
        final byte[] once = new byte[256];
        for (int i = 0; i < 256; i++) {
            once[i] = (byte) i;
        }
        emit("all-256-once", once);

        // Every symbol, count proportional to its value: 256 non-zero frequencies, all different,
        // and a rounding error that has to land somewhere.
        final byte[] skewed = new byte[256 * 257 / 2];
        int at = 0;
        for (int i = 0; i < 256; i++) {
            for (int n = 0; n < i; n++) {
                skewed[at++] = (byte) i;
            }
        }
        emit("skewed-256", Arrays.copyOf(skewed, at));

        // A contiguous block of symbols, which is what the run-length writer is for.
        emit("contiguous-run", spread(5, 13, 100));
        // The same at the top of the alphabet, where the run scan walks off the end of the table.
        emit("high-run", spread(250, 256, 100));
        // Symbol zero carries the maximum, which is the one case the run writer guards with j != 0.
        final byte[] zeroHeavy = new byte[1000];
        Arrays.fill(zeroHeavy, (byte) 0);
        for (int i = 0; i < 40; i++) {
            zeroHeavy[i * 25] = (byte) (1 + (i % 3));
        }
        emit("zero-heavy", zeroHeavy);

        // Quality scores: a narrow band, which is the shape rANS actually meets in a CRAM.
        final byte[] quals = new byte[5000];
        long seed = 0x5DEECE66DL;
        for (int i = 0; i < quals.length; i++) {
            seed = seed * 6364136223846793005L + 1442695040888963407L;
            quals[i] = (byte) (30 + ((seed >>> 33) % 8));
        }
        emit("quality-band", quals);

        // Something with no structure at all, where the table is nearly flat.
        final byte[] noise = new byte[10000];
        seed = 0x1234_5678L;
        for (int i = 0; i < noise.length; i++) {
            seed = seed * 6364136223846793005L + 1442695040888963407L;
            noise[i] = (byte) (seed >>> 33);
        }
        emit("noise-10000", noise);

        // The order the writer actually used, against the order it was asked for. Below four bytes
        // the answer is not the question.
        for (int length = 0; length <= 8; length++) {
            final byte[] input = repeat("ACGT", length);
            for (final RANSParams.ORDER requested :
                    new RANSParams.ORDER[] {RANSParams.ORDER.ZERO, RANSParams.ORDER.ONE}) {
                final ByteBuffer out = new RANS4x8Encode()
                        .compress(CompressionUtils.wrap(input), new RANS4x8Params(requested));
                final String written = out.remaining() == 0
                        ? "-"
                        : Integer.toString(out.get(0) & 0xFF);
                System.out.printf("orderused\t%d\t%s\t%s%n", length, requested, written);
            }
        }
    }

    /** `source` repeated until exactly `length` bytes have been produced. */
    static byte[] repeat(final String source, final int length) {
        final byte[] pattern = source.getBytes();
        final byte[] out = new byte[length];
        for (int i = 0; i < length; i++) {
            out[i] = pattern[i % pattern.length];
        }
        return out;
    }

    /** `count` occurrences of each symbol in [from, to), interleaved. */
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

        final RANS4x8Encode encode = new RANS4x8Encode();
        final ByteBuffer compressed =
                encode.compress(CompressionUtils.wrap(input), new RANS4x8Params(RANSParams.ORDER.ZERO));
        final byte[] out = new byte[compressed.remaining()];
        compressed.duplicate().get(out);

        System.out.printf("enc\t%s\t%d\t%s%n", label, out.length, sha256(out));
        if (out.length > 0 && out.length <= MAX_INLINE_BYTES) {
            System.out.printf("bytes\t%s\t%s%n", label, hex(out));
        }
        if (out.length == 0) {
            // The empty input, which produces nothing at all rather than a bare prefix.
            System.out.printf("roundtrip\t%s\t%s%n", label,
                    new RANS4x8Decode().uncompress(CompressionUtils.wrap(out)).remaining() == 0
                            ? "ok" : "MISMATCH");
            return;
        }

        final ByteBuffer prefix = CompressionUtils.wrap(out);
        final int orderByte = prefix.get() & 0xFF;
        final int compressedSize = prefix.getInt();
        final int rawSize = prefix.getInt();
        System.out.printf("prefix\t%s\t%d\t%d\t%d%n", label, orderByte, compressedSize, rawSize);

        // The normalised frequencies, from the private method that computes them, beside the raw
        // counts they came from. The pair is what a port has to reproduce exactly: the counts are
        // easy and the normalisation is not.
        final int[] counts = new int[256];
        for (final byte b : input) {
            counts[0xFF & b]++;
        }
        final int[] normalised = calcFrequenciesOrder0(input);
        for (int symbol = 0; symbol < 256; symbol++) {
            if (normalised[symbol] != 0) {
                System.out.printf("norm\t%s\t%d\t%d\t%d%n", label, symbol, counts[symbol],
                        normalised[symbol]);
            }
        }

        // The frequency table as the writer writes it, and its size, which is where the blob starts.
        final ByteBuffer table = CompressionUtils.allocateByteBuffer(256 * 3 + 1);
        final int tableSize = writeFrequenciesOrder0(table, normalised.clone());
        final byte[] tableBytes = new byte[tableSize];
        table.rewind();
        table.get(tableBytes);
        System.out.printf("freqtab\t%s\t%d\t%s%n", label, tableSize, hex(tableBytes));

        // The encoding symbols, field by field, because the fixed-point reciprocal is the part a
        // port gets subtly wrong and a byte-level golden alone would not say where.
        int cumulative = 0;
        for (int symbol = 0; symbol < 256; symbol++) {
            if (normalised[symbol] == 0) {
                continue;
            }
            final RANSEncodingSymbol encodingSymbol = new RANSEncodingSymbol();
            encodingSymbol.set(cumulative, normalised[symbol], TOTAL_FREQ_SHIFT);
            System.out.printf("sym\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d%n", label, symbol,
                    cumulative, normalised[symbol],
                    field(encodingSymbol, "xMax"), field(encodingSymbol, "rcpFreq"),
                    field(encodingSymbol, "bias"), field(encodingSymbol, "cmplFreq"),
                    field(encodingSymbol, "rcpShift"));
            cumulative += normalised[symbol];
        }

        // The four states at the head of the blob, read little-endian, which is what the two
        // cancelling reversals leave behind.
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

    /** `RANS4x8Encode.calcFrequenciesOrder0`, which is private and is the whole of the arithmetic. */
    static int[] calcFrequenciesOrder0(final byte[] input) throws Exception {
        final Method method =
                RANS4x8Encode.class.getDeclaredMethod("calcFrequenciesOrder0", ByteBuffer.class);
        method.setAccessible(true);
        return (int[]) method.invoke(null, CompressionUtils.wrap(input));
    }

    /** `RANS4x8Encode.writeFrequenciesOrder0`, likewise private, and likewise the point. */
    static int writeFrequenciesOrder0(final ByteBuffer cp, final int[] frequencies)
            throws Exception {
        final Method method = RANS4x8Encode.class
                .getDeclaredMethod("writeFrequenciesOrder0", ByteBuffer.class, int[].class);
        method.setAccessible(true);
        return (int) method.invoke(null, cp, frequencies);
    }

    /** One private field of an encoding symbol, as an unsigned 32-bit value where it is an int. */
    static long field(final RANSEncodingSymbol symbol, final String name) throws Exception {
        final Field declared = RANSEncodingSymbol.class.getDeclaredField(name);
        declared.setAccessible(true);
        final Object value = declared.get(symbol);
        // rcpFreq is an int holding a value the encoder always masks to 32 unsigned bits, so the
        // sign it carries in Java is not the number it stands for.
        return value instanceof Long ? (Long) value : (0xFFFFFFFFL & (Integer) value);
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
