/*
 * The three integer codecs written on the bit stream: Beta, Gamma and Subexponential.
 *
 * The bit stream is pinned. These are the three codecs the encoding map can name that are written
 * on it rather than on an external byte block, and each is a different bargain between a fixed
 * width and a variable one.
 *
 * Six things here are decisions rather than layout.
 *
 *   - EVERY CODEC CARRIES AN OFFSET, added before storage and subtracted after. It is how a range
 *     that starts below zero is stored in bits that cannot;
 *   - BETA REFUSES WHAT DOES NOT FIT, with two messages naming the value, the offset and the
 *     limit. It is the only one of the three with an upper bound at all;
 *   - GAMMA AND SUBEXPONENTIAL COMPUTE A BIT LENGTH WITH FLOATING-POINT LOG. `Math.log(v) /
 *     Math.log(2)` decides how many bits the value takes, so the codec's output depends on a
 *     double division being exact at every power of two;
 *   - GAMMA REFUSES ZERO AND BELOW, because its length prefix cannot encode a value with no bits;
 *   - SUBEXPONENTIAL HAS TWO REGIMES, split at 2^k: below it the value is written in k bits with
 *     no prefix, above it in b = floor(log2(v)) bits with a unary prefix of b - k + 1 ones;
 *   - THE CORE BLOCK IS RAW. Whatever these three write goes into the file uncompressed, and it is
 *     written first because the specification says so.
 *
 * Output:
 *
 *     beta\t<offset>\t<bits>\t<value>\t<core block hex>\t<value read back>
 *     gamma\t<offset>\t<value>\t<core block hex>\t<value read back>
 *     subexp\t<offset>\t<k>\t<value>\t<core block hex>\t<value read back>
 *     seq\t<codec>\t<params>\t<values>\t<core block hex>\t<values read back>
 *     err\t<codec>\t<params>\t<value>\t<class>\t<message>
 *
 * Usage: CramCoreCodecDump
 */

import htsjdk.samtools.cram.io.ITF8;
import htsjdk.samtools.cram.encoding.CRAMCodec;
import htsjdk.samtools.cram.encoding.CRAMEncoding;
import htsjdk.samtools.cram.encoding.core.BetaIntegerEncoding;
import htsjdk.samtools.cram.encoding.core.GammaIntegerEncoding;
import htsjdk.samtools.cram.encoding.core.SubexponentialIntegerEncoding;
import htsjdk.samtools.cram.structure.CompressionHeader;
import htsjdk.samtools.cram.structure.CompressorCache;
import htsjdk.samtools.cram.structure.SliceBlocks;
import htsjdk.samtools.cram.structure.SliceBlocksReadStreams;
import htsjdk.samtools.cram.structure.SliceBlocksWriteStreams;

import java.util.ArrayList;
import java.util.List;
import java.util.StringJoiner;

public class CramCoreCodecDump {

    public static void main(final String[] args) {
        System.out.println("# CramCoreCodecDump: the three integer codecs written on the bit stream");

        // Beta: a fixed width, and the only codec of the three with an upper bound.
        for (final int value : new int[] {0, 1, 2, 7, 8, 15}) {
            beta(0, 4, value);
        }
        beta(0, 1, 0);
        beta(0, 1, 1);
        beta(0, 8, 255);
        beta(0, 32, Integer.MAX_VALUE);
        beta(10, 4, -10);
        beta(10, 4, 5);
        beta(-5, 4, 5);

        // Gamma: Elias gamma, whose length comes from a floating-point log.
        for (final int value : new int[] {1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64,
                65, 127, 128, 255, 256, 1023, 1024, 65535, 65536, 1048575, 1048576, 16777215,
                16777216, 1073741823, 1073741824, Integer.MAX_VALUE}) {
            gamma(0, value);
        }
        gamma(1, 0);
        gamma(5, -4);

        // Subexponential: two regimes, split at 2^k.
        for (final int k : new int[] {0, 1, 2, 3}) {
            for (final int value : new int[] {0, 1, 2, 3, 4, 7, 8, 15, 16, 31, 32}) {
                subexp(0, k, value);
            }
        }
        subexp(0, 2, 1023);
        subexp(0, 2, 1024);
        subexp(0, 2, Integer.MAX_VALUE);
        subexp(3, 2, -3);

        // Several values in a row, which is where the bit packing shows.
        sequence("beta", new BetaIntegerEncoding(0, 3), "offset=0 bits=3", new int[] {1, 2, 3, 4});
        sequence("gamma", gammaEncoding(0), "offset=0", new int[] {1, 2, 3, 4});
        sequence("subexp", subexpEncoding(0, 1), "offset=0 k=1",
                new int[] {0, 1, 2, 3});

        // What each refuses.
        errBeta(0, 4, 16);
        errBeta(0, 4, -1);
        errBeta(0, 1, 2);
        errGamma(0, 0);
        errGamma(0, -1);
        errGamma(1, -1);
        errSubexp(0, 2, -1);
    }

    static void beta(final int offset, final int bits, final int value) {
        final BetaIntegerEncoding encoding = new BetaIntegerEncoding(offset, bits);
        final Result result = round(encoding, new int[] {value});
        if (result.error != null) {
            System.out.printf("err\tbeta\toffset=%d bits=%d\t%d\t%s\t%s%n", offset, bits, value,
                    result.error.getClass().getSimpleName(),
                    String.valueOf(result.error.getMessage()));
            return;
        }
        System.out.printf("beta\t%d\t%d\t%d\t%s\t%s%n", offset, bits, value, result.hex,
                result.values);
    }

    /** Only Beta's constructor is public; the other two are reached through their factory. */
    static GammaIntegerEncoding gammaEncoding(final int offset) {
        return GammaIntegerEncoding.fromSerializedEncodingParams(itf8(offset));
    }

    static SubexponentialIntegerEncoding subexpEncoding(final int offset, final int k) {
        return SubexponentialIntegerEncoding.fromSerializedEncodingParams(itf8(offset, k));
    }

    static byte[] itf8(final int... values) {
        final java.io.ByteArrayOutputStream out = new java.io.ByteArrayOutputStream();
        for (final int value : values) {
            ITF8.writeUnsignedITF8(value, out);
        }
        return out.toByteArray();
    }

    static void gamma(final int offset, final int value) {
        final GammaIntegerEncoding encoding = gammaEncoding(offset);
        final Result result = round(encoding, new int[] {value});
        if (result.error != null) {
            System.out.printf("err\tgamma\toffset=%d\t%d\t%s\t%s%n", offset, value,
                    result.error.getClass().getSimpleName(),
                    String.valueOf(result.error.getMessage()));
            return;
        }
        System.out.printf("gamma\t%d\t%d\t%s\t%s%n", offset, value, result.hex, result.values);
    }

    static void subexp(final int offset, final int k, final int value) {
        final SubexponentialIntegerEncoding encoding = subexpEncoding(offset, k);
        final Result result = round(encoding, new int[] {value});
        if (result.error != null) {
            System.out.printf("err\tsubexp\toffset=%d k=%d\t%d\t%s\t%s%n", offset, k, value,
                    result.error.getClass().getSimpleName(),
                    String.valueOf(result.error.getMessage()));
            return;
        }
        System.out.printf("subexp\t%d\t%d\t%d\t%s\t%s%n", offset, k, value, result.hex,
                result.values);
    }

    static void sequence(final String name, final CRAMEncoding<Integer> encoding,
            final String params, final int[] values) {
        final Result result = round(encoding, values);
        final StringJoiner joiner = new StringJoiner(",");
        for (final int value : values) {
            joiner.add(Integer.toString(value));
        }
        if (result.error != null) {
            System.out.printf("err\t%s\t%s\t%s\t%s\t%s%n", name, params, joiner.toString(),
                    result.error.getClass().getSimpleName(),
                    String.valueOf(result.error.getMessage()));
            return;
        }
        System.out.printf("seq\t%s\t%s\t%s\t%s\t%s%n", name, params, joiner.toString(), result.hex,
                result.values);
    }

    static void errBeta(final int offset, final int bits, final int value) {
        beta(offset, bits, value);
    }

    static void errGamma(final int offset, final int value) {
        gamma(offset, value);
    }

    static void errSubexp(final int offset, final int k, final int value) {
        subexp(offset, k, value);
    }

    /** Write the values through the codec, take the raw core block, and read them back. */
    static Result round(final CRAMEncoding<Integer> encoding, final int[] values) {
        final Result result = new Result();
        try {
            final CompressionHeader header = new CompressionHeader();
            final SliceBlocksWriteStreams writeStreams = new SliceBlocksWriteStreams(header);
            final CRAMCodec<Integer> writer = encoding.buildCodec(null, writeStreams);
            for (final int value : values) {
                writer.write(value);
            }
            final SliceBlocks blocks = writeStreams.flushStreamsToBlocks();
            result.hex = hex(blocks.getCoreBlock().getRawContent());

            final SliceBlocksReadStreams readStreams =
                    new SliceBlocksReadStreams(blocks, new CompressorCache());
            final CRAMCodec<Integer> reader = encoding.buildCodec(readStreams, null);
            final List<String> back = new ArrayList<>();
            for (int i = 0; i < values.length; i++) {
                back.add(Integer.toString(reader.read()));
            }
            final StringJoiner joiner = new StringJoiner(",");
            back.forEach(joiner::add);
            result.values = joiner.toString();
        } catch (final Throwable t) {
            result.error = t;
        }
        return result;
    }

    static class Result {
        String hex;
        String values;
        Throwable error;
    }

    static String hex(final byte[] bytes) {
        final StringBuilder b = new StringBuilder(bytes.length * 2);
        for (final byte value : bytes) {
            b.append(String.format("%02x", value));
        }
        return b.length() == 0 ? "-" : b.toString();
    }
}
