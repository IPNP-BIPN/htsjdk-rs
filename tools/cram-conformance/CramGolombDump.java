/*
 * Golomb, Golomb-Rice and Golomb-Long: the three codecs the CRAM specification is removing.
 *
 * They are reachable all the same. The encoding factory dispatches to them by identifier, so a file
 * that names one has to be read, and a port that skips them cannot claim to read every legal file.
 * htsjdk marks them experimental and logs a warning when one is built; nothing else guards them.
 *
 * Seven things here are decisions rather than layout.
 *
 *   - THE QUOTIENT IS UNARY, so the bits a value costs grow with the value divided by m. A value of
 *     1000 with m = 2 writes five hundred ones before anything else;
 *   - GOLOMB REFUSES m < 2 AND GOLOMB-RICE DOES NOT. The same mistake is an exception in one and a
 *     silently different encoding in the other;
 *   - GOLOMB-RICE'S PARAMETER IS NOT m. The encoding calls it m and hands it to the codec as
 *     log2m, so an encoding built with 8 divides by 256;
 *   - THE REMAINDER IS WRITTEN AT ONE OF TWO WIDTHS, ceiling - 1 or ceiling, chosen by comparing it
 *     against 2^ceiling - m. That is what lets a non-power-of-two m avoid wasting a bit;
 *   - CEILING COMES FROM A FLOATING-POINT LOG, (int)(Math.log(m) / Math.log(2) + 1), so the width
 *     of every remainder depends on a double division landing on the right side of an integer;
 *   - THE COMPARISONS ARE AGAINST Math.pow, a double, so an int remainder is promoted to double
 *     before being compared and again before being subtracted;
 *   - THE TWO REFUSALS OF read(length) ARE WORDED DIFFERENTLY. Golomb and Golomb-Long say
 *     "Multi-value read method not defined.", Golomb-Rice says "Not implemented."
 *
 * Negative values are in the corpus for Golomb only. Golomb-Rice takes an unsigned right shift of
 * a long, so a negative one would emit about 2^60 unary bits: that is not a row, it is a hang.
 * For the same reason every value here is small next to its m. A single long of 2^32 with m = 10
 * writes four hundred million ones, which is a hundred-megabyte row, measured once and removed.
 *
 * Output:
 *
 *     golomb\t<offset>\t<m>\t<value>\t<core block hex>\t<value read back>
 *     rice\t<offset>\t<log2m>\t<value>\t<core block hex>\t<value read back>
 *     long\t<offset>\t<m>\t<value>\t<core block hex>\t<value read back>
 *     seq\t<codec>\t<params>\t<values>\t<core block hex>\t<values read back>
 *     ser\t<codec>\t<offset>\t<m>\t<hex>\t<reparsed hex>
 *     err\t<codec>\t<params>\t<value>\t<class>\t<message>
 *
 * Usage: CramGolombDump
 */

import htsjdk.samtools.cram.io.ITF8;
import htsjdk.samtools.cram.encoding.CRAMCodec;
import htsjdk.samtools.cram.encoding.CRAMEncoding;
import htsjdk.samtools.cram.encoding.core.experimental.GolombIntegerEncoding;
import htsjdk.samtools.cram.encoding.core.experimental.GolombLongEncoding;
import htsjdk.samtools.cram.encoding.core.experimental.GolombRiceIntegerEncoding;
import htsjdk.samtools.cram.structure.CompressionHeader;
import htsjdk.samtools.cram.structure.CompressorCache;
import htsjdk.samtools.cram.structure.SliceBlocks;
import htsjdk.samtools.cram.structure.SliceBlocksReadStreams;
import htsjdk.samtools.cram.structure.SliceBlocksWriteStreams;

import java.io.ByteArrayOutputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.StringJoiner;

public class CramGolombDump {

    public static void main(final String[] args) {
        System.out.println("# CramGolombDump: the three codecs the specification is removing");

        // Golomb, over the values and the divisors where the two remainder widths both appear.
        for (final int m : new int[] {2, 3, 4, 5, 8, 10}) {
            for (final int value : new int[] {0, 1, 2, 3, 7, 8, 15}) {
                golomb(0, m, value);
            }
        }
        golomb(0, 2, 100);
        golomb(0, 10, 1000);
        golomb(5, 4, -5);
        golomb(5, 4, 0);
        golomb(0, 4, -1);
        golomb(0, 4, -4);

        // Golomb-Rice, whose parameter is a power rather than a divisor.
        for (final int log2m : new int[] {0, 1, 2, 3}) {
            for (final int value : new int[] {0, 1, 2, 3, 7, 8, 15, 16}) {
                rice(0, log2m, value);
            }
        }
        rice(0, 8, 255);
        rice(0, 8, 256);
        rice(4, 2, -4);

        // Golomb-Long, the same arithmetic on a long.
        for (final int m : new int[] {2, 4, 10}) {
            for (final long value : new long[] {0L, 1L, 7L, 100L}) {
                golombLong(0, m, value);
            }
        }
        golombLong(0, 1000, 1000000L);
        golombLong(10, 4, -10L);

        // Several values in a row, which is where the unary prefixes butt against each other.
        seqGolomb(0, 4, new int[] {0, 1, 2, 3});
        seqRice(0, 2, new int[] {0, 1, 2, 3});
        seqLong(0, 4, new long[] {0L, 1L, 2L, 3L});

        // The encoding parameters, which are two ITF8s whatever the codec.
        ser("golomb", 0, 4);
        ser("golomb", 10, 300);
        ser("rice", 0, 3);
        ser("long", 0, 4);

        // What each refuses, and where one refuses and another does not.
        errGolomb(0, 1, 0);
        errGolomb(0, 0, 0);
        errGolomb(0, -1, 0);
        errLong(0, 1, 0L);
        // Golomb-Rice takes the same parameter without a word.
        rice(0, 1, 0);
        errReadLength("golomb", GolombIntegerEncoding.fromSerializedEncodingParams(itf8(0, 4)));
        errReadLength("rice", new GolombRiceIntegerEncoding(0, 2));
        errReadLength("long", new GolombLongEncoding(0, 4));
    }

    static void golomb(final int offset, final int m, final int value) {
        final CRAMEncoding<Integer> encoding =
                GolombIntegerEncoding.fromSerializedEncodingParams(itf8(offset, m));
        final Result result = round(encoding, one(value));
        if (result.error != null) {
            System.out.printf("err\tgolomb\toffset=%d m=%d\t%d\t%s\t%s%n", offset, m, value,
                    result.error.getClass().getSimpleName(),
                    String.valueOf(result.error.getMessage()));
            return;
        }
        System.out.printf("golomb\t%d\t%d\t%d\t%s\t%s%n", offset, m, value, result.hex,
                result.values);
    }

    static void rice(final int offset, final int log2m, final int value) {
        final CRAMEncoding<Integer> encoding = new GolombRiceIntegerEncoding(offset, log2m);
        final Result result = round(encoding, one(value));
        if (result.error != null) {
            System.out.printf("err\trice\toffset=%d log2m=%d\t%d\t%s\t%s%n", offset, log2m, value,
                    result.error.getClass().getSimpleName(),
                    String.valueOf(result.error.getMessage()));
            return;
        }
        System.out.printf("rice\t%d\t%d\t%d\t%s\t%s%n", offset, log2m, value, result.hex,
                result.values);
    }

    static void golombLong(final int offset, final int m, final long value) {
        final CRAMEncoding<Long> encoding = new GolombLongEncoding(offset, m);
        final List<Long> values = new ArrayList<>();
        values.add(value);
        final Result result = round(encoding, values);
        if (result.error != null) {
            System.out.printf("err\tlong\toffset=%d m=%d\t%d\t%s\t%s%n", offset, m, value,
                    result.error.getClass().getSimpleName(),
                    String.valueOf(result.error.getMessage()));
            return;
        }
        System.out.printf("long\t%d\t%d\t%d\t%s\t%s%n", offset, m, value, result.hex,
                result.values);
    }

    static void seqGolomb(final int offset, final int m, final int[] values) {
        seq("golomb", String.format("offset=%d m=%d", offset, m),
                GolombIntegerEncoding.fromSerializedEncodingParams(itf8(offset, m)), boxed(values));
    }

    static void seqRice(final int offset, final int log2m, final int[] values) {
        seq("rice", String.format("offset=%d log2m=%d", offset, log2m),
                new GolombRiceIntegerEncoding(offset, log2m), boxed(values));
    }

    static void seqLong(final int offset, final int m, final long[] values) {
        final List<Long> boxed = new ArrayList<>();
        for (final long value : values) {
            boxed.add(value);
        }
        seq("long", String.format("offset=%d m=%d", offset, m),
                new GolombLongEncoding(offset, m), boxed);
    }

    static <T> void seq(final String name, final String params, final CRAMEncoding<T> encoding,
            final List<T> values) {
        final Result result = round(encoding, values);
        final StringJoiner shown = new StringJoiner(",");
        for (final T value : values) {
            shown.add(String.valueOf(value));
        }
        if (result.error != null) {
            System.out.printf("err\t%s\t%s\t%s\t%s\t%s%n", name, params, shown,
                    result.error.getClass().getSimpleName(),
                    String.valueOf(result.error.getMessage()));
            return;
        }
        System.out.printf("seq\t%s\t%s\t%s\t%s\t%s%n", name, params, shown, result.hex,
                result.values);
    }

    static void ser(final String name, final int offset, final int m) {
        final CRAMEncoding<?> encoding;
        final byte[] serialized;
        switch (name) {
            case "golomb":
                encoding = GolombIntegerEncoding.fromSerializedEncodingParams(itf8(offset, m));
                serialized = encoding.toSerializedEncodingParams();
                System.out.printf("ser\tgolomb\t%d\t%d\t%s\t%s%n", offset, m, hex(serialized),
                        hex(GolombIntegerEncoding.fromSerializedEncodingParams(serialized)
                                .toSerializedEncodingParams()));
                return;
            case "rice":
                serialized = new GolombRiceIntegerEncoding(offset, m).toSerializedEncodingParams();
                System.out.printf("ser\trice\t%d\t%d\t%s\t%s%n", offset, m, hex(serialized),
                        hex(GolombRiceIntegerEncoding.fromSerializedEncodingParams(serialized)
                                .toSerializedEncodingParams()));
                return;
            default:
                serialized = new GolombLongEncoding(offset, m).toSerializedEncodingParams();
                System.out.printf("ser\tlong\t%d\t%d\t%s\t%s%n", offset, m, hex(serialized),
                        hex(GolombLongEncoding.fromSerializedEncodingParams(serialized)
                                .toSerializedEncodingParams()));
        }
    }

    static void errGolomb(final int offset, final int m, final int value) {
        golomb(offset, m, value);
    }

    static void errLong(final int offset, final int m, final long value) {
        golombLong(offset, m, value);
    }

    static void errReadLength(final String name, final CRAMEncoding<?> encoding) {
        try {
            final CRAMCodec<?> codec = encoding.buildCodec(null, writeStreams());
            System.out.printf("err\t%s\tread-length\t4\t-\t%s%n", name,
                    String.valueOf(codec.read(4)));
        } catch (final Throwable t) {
            System.out.printf("err\t%s\tread-length\t4\t%s\t%s%n", name,
                    t.getClass().getSimpleName(), String.valueOf(t.getMessage()));
        }
    }

    static List<Integer> one(final int value) {
        final List<Integer> values = new ArrayList<>();
        values.add(value);
        return values;
    }

    static List<Integer> boxed(final int[] values) {
        final List<Integer> boxed = new ArrayList<>(values.length);
        for (final int value : values) {
            boxed.add(value);
        }
        return boxed;
    }

    static byte[] itf8(final int... values) {
        final ByteArrayOutputStream out = new ByteArrayOutputStream();
        for (final int value : values) {
            ITF8.writeUnsignedITF8(value, out);
        }
        return out.toByteArray();
    }

    static SliceBlocksWriteStreams writeStreams() {
        return new SliceBlocksWriteStreams(new CompressionHeader());
    }

    /** Write the values through the codec, take the raw core block, and read them back. */
    static <T> Result round(final CRAMEncoding<T> encoding, final List<T> values) {
        final Result result = new Result();
        try {
            final SliceBlocksWriteStreams streams = writeStreams();
            final CRAMCodec<T> writer = encoding.buildCodec(null, streams);
            for (final T value : values) {
                writer.write(value);
            }
            final SliceBlocks blocks = streams.flushStreamsToBlocks();
            result.hex = hex(blocks.getCoreBlock().getRawContent());

            final SliceBlocksReadStreams readStreams =
                    new SliceBlocksReadStreams(blocks, new CompressorCache());
            final CRAMCodec<T> reader = encoding.buildCodec(readStreams, null);
            final StringJoiner joiner = new StringJoiner(",");
            for (int i = 0; i < values.size(); i++) {
                joiner.add(String.valueOf(reader.read()));
            }
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
        final StringBuilder builder = new StringBuilder(bytes.length * 2);
        for (final byte value : bytes) {
            builder.append(String.format("%02x", value));
        }
        return builder.length() == 0 ? "-" : builder.toString();
    }
}
