/*
 * Canonical Huffman, the fourth codec written on the CRAM core bit stream.
 *
 * A CRAM file does not carry a Huffman tree. It carries an alphabet and one code word length per
 * symbol, and both writer and reader rebuild the same code words from that pair alone. So the
 * thing to pin is the rebuilding, not the tree that produced the lengths.
 *
 * Seven things here are decisions rather than layout.
 *
 *   - THE CODE WORDS ARE DERIVED, NOT STORED. Symbols are grouped by length, sorted inside each
 *     group by their own natural order, and handed consecutive integers, shifted left whenever the
 *     length grows. Two alphabets with the same lengths in a different order get the same codes;
 *   - BYTE SYMBOLS SORT SIGNED. The grouping is a TreeSet, so for the byte alphabet 0x80 sorts
 *     below 0x01 and takes the earlier code word;
 *   - A ONE-SYMBOL ALPHABET WRITES NOTHING. Its length is zero, so every write emits zero bits and
 *     every read consumes zero bits and returns the only symbol;
 *   - THE OVERFLOW CHECK COUNTS SET BITS, NOT BIT LENGTH. `Integer.bitCount(codeValue) >
 *     bitLength` is what refuses an impossible length table, so it fires later than a check on the
 *     value's width would;
 *   - READING WALKS THE LENGTHS IN ORDER, consuming only the difference between one length and the
 *     next, and matching against a table indexed by code word;
 *   - AN UNMATCHED CODE WORD RUNS OFF THE END OF THAT TABLE, which is sized to the largest code
 *     word and not to the largest the bits can hold. So a truncated or foreign stream comes out as
 *     an ArrayIndexOutOfBoundsException, and the codec's own "unable to map" message is reachable
 *     only with an empty alphabet;
 *   - A LENGTH TABLE SHORTER THAN THE ALPHABET IS NOT REFUSED. The scan is driven by the lengths,
 *     so the surplus symbols are silently dropped; the other way round it throws;
 *   - THE ENCODING PARAMETERS ARE ASYMMETRIC BETWEEN THE TWO FLAVOURS. Integer symbols are
 *     serialized as ITF8, byte symbols as raw bytes; both lengths are ITF8.
 *
 * Output:
 *
 *     canon\t<symbols>\t<lengths>\t<symbol:code/length,...>
 *     canonerr\t<symbols>\t<lengths>\t<class>\t<message>
 *     ser\tint|byte\t<symbols>\t<lengths>\t<params hex>\t<reparsed params hex>
 *     round\tint|byte\t<symbols>\t<lengths>\t<values>\t<core block hex>\t<values read back>
 *     cross\t<write symbols>\t<write lengths>\t<values>\t<read symbols>\t<read lengths>\t<result>
 *     err\t<what>\t<symbols>\t<lengths>\t<detail>\t<class>\t<message>
 *
 * Usage: CramHuffmanDump
 */

import htsjdk.samtools.cram.encoding.CRAMCodec;
import htsjdk.samtools.cram.encoding.CRAMEncoding;
import htsjdk.samtools.cram.encoding.core.CanonicalHuffmanByteEncoding;
import htsjdk.samtools.cram.encoding.core.CanonicalHuffmanIntegerEncoding;
import htsjdk.samtools.cram.encoding.core.huffmanUtils.HuffmanBitCode;
import htsjdk.samtools.cram.encoding.core.huffmanUtils.HuffmanCanoncialCodeGenerator;
import htsjdk.samtools.cram.encoding.core.huffmanUtils.HuffmanParams;
import htsjdk.samtools.cram.structure.CompressionHeader;
import htsjdk.samtools.cram.structure.CompressorCache;
import htsjdk.samtools.cram.structure.SliceBlocks;
import htsjdk.samtools.cram.structure.SliceBlocksReadStreams;
import htsjdk.samtools.cram.structure.SliceBlocksWriteStreams;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

public class CramHuffmanDump {

    public static void main(final String[] args) {
        System.out.println("# CramHuffmanDump: canonical Huffman on the core bit stream");

        // The code words, rebuilt from an alphabet and one length per symbol.
        canon(new int[] {1, 2, 3, 4}, new int[] {2, 2, 2, 2});
        canon(new int[] {4, 3, 2, 1}, new int[] {2, 2, 2, 2});
        canon(new int[] {1, 2, 3}, new int[] {1, 2, 2});
        canon(new int[] {3, 2, 1}, new int[] {2, 2, 1});
        canon(new int[] {5, 1, 3}, new int[] {2, 2, 2});
        canon(new int[] {1, 2, 3, 4, 5}, new int[] {1, 2, 3, 4, 4});
        canon(new int[] {42}, new int[] {0});
        canon(new int[] {7}, new int[] {1});
        canon(new int[] {0, 1}, new int[] {1, 1});
        canon(new int[] {1000000, 2, 3}, new int[] {1, 2, 2});
        canon(new int[] {-1, -2, 0}, new int[] {1, 2, 2});
        canon(new int[] {1, 2, 3, 4, 5, 6, 7, 8}, new int[] {3, 3, 3, 3, 3, 3, 3, 3});
        canon(new int[] {1, 2, 3, 4}, new int[] {1, 3, 3, 2});

        // Where the length table is refused, and where it is not.
        canonErr(new int[] {1, 2, 3, 4}, new int[] {1, 1, 1, 1});
        canonErr(new int[] {1, 2, 3}, new int[] {1, 1, 1});
        canonErr(new int[] {1, 2, 3, 4, 5, 6, 7, 8, 9}, new int[] {3, 3, 3, 3, 3, 3, 3, 3, 3});
        canonErr(new int[] {1, 2, 3}, new int[] {1, 2});
        canonErr(new int[] {1, 2}, new int[] {1, 2, 2});
        canonErr(new int[] {}, new int[] {});

        // The byte alphabet, whose symbols sort signed.
        canonByte(new byte[] {'A', 'C', 'G', 'T'}, new int[] {2, 2, 2, 2});
        canonByte(new byte[] {0x41, (byte) 0x80, 0x7f, 0x00}, new int[] {2, 2, 2, 2});
        canonByte(new byte[] {(byte) 0xff, 0x01}, new int[] {1, 1});
        canonByte(new byte[] {'N'}, new int[] {0});

        // The encoding parameters, which are ITF8 for one flavour and raw for the other.
        serInt(new int[] {1, 2, 3, 4}, new int[] {2, 2, 2, 2});
        serInt(new int[] {42}, new int[] {0});
        serInt(new int[] {1000000, 2, 3}, new int[] {1, 2, 2});
        serInt(new int[] {-1, 0}, new int[] {1, 1});
        serByte(new byte[] {'A', 'C', 'G', 'T'}, new int[] {2, 2, 2, 2});
        serByte(new byte[] {(byte) 0x80, 0x7f}, new int[] {1, 1});
        serByte(new byte[] {'N'}, new int[] {0});

        // What lands in the core block, and what comes back out of it.
        roundInt(new int[] {1, 2, 3, 4}, new int[] {2, 2, 2, 2}, new int[] {1, 2, 3, 4});
        roundInt(new int[] {1, 2, 3, 4}, new int[] {2, 2, 2, 2}, new int[] {4, 4, 4, 4, 4});
        roundInt(new int[] {1, 2, 3}, new int[] {1, 2, 2}, new int[] {1, 1, 2, 3});
        roundInt(new int[] {1, 2, 3, 4, 5}, new int[] {1, 2, 3, 4, 4}, new int[] {5, 4, 3, 2, 1});
        roundInt(new int[] {42}, new int[] {0}, new int[] {42, 42, 42});
        roundInt(new int[] {7}, new int[] {1}, new int[] {7, 7});
        roundInt(new int[] {1000000, 2, 3}, new int[] {1, 2, 2}, new int[] {1000000, 3, 2});
        roundInt(new int[] {-1, -2, 0}, new int[] {1, 2, 2}, new int[] {-1, -2, 0, -1});
        roundInt(new int[] {1, 2, 3, 4, 5, 6, 7, 8}, new int[] {3, 3, 3, 3, 3, 3, 3, 3},
                new int[] {1, 8, 1, 8, 1, 8, 1, 8});
        roundInt(new int[] {0, 1}, new int[] {1, 1}, new int[] {0});

        roundByte(new byte[] {'A', 'C', 'G', 'T'}, new int[] {2, 2, 2, 2},
                new byte[] {'A', 'C', 'G', 'T', 'A'});
        roundByte(new byte[] {0x41, (byte) 0x80, 0x7f, 0x00}, new int[] {2, 2, 2, 2},
                new byte[] {(byte) 0x80, 0x00, 0x7f, 0x41});
        roundByte(new byte[] {(byte) 0xff, 0x01}, new int[] {1, 1},
                new byte[] {0x01, (byte) 0xff, 0x01});
        roundByte(new byte[] {'N'}, new int[] {0}, new byte[] {'N', 'N'});

        // Written with one alphabet, read with another.
        cross(new int[] {1, 2, 3, 4}, new int[] {2, 2, 2, 2}, new int[] {4},
                new int[] {1, 2, 3}, new int[] {1, 2, 2});
        cross(new int[] {1, 2, 3, 4}, new int[] {2, 2, 2, 2}, new int[] {1},
                new int[] {5, 6, 7, 8}, new int[] {2, 2, 2, 2});
        cross(new int[] {1, 2}, new int[] {1, 1}, new int[] {1},
                new int[] {1, 2, 3, 4, 5}, new int[] {1, 2, 3, 4, 4});
        // A code word past the end of the code-word-to-symbol table.
        cross(new int[] {1, 2, 3, 4}, new int[] {2, 2, 2, 2}, new int[] {4},
                new int[] {1, 2, 3}, new int[] {2, 2, 2});
        cross(new int[] {1, 2}, new int[] {1, 1}, new int[] {2, 1, 1, 1},
                new int[] {1, 2, 3, 4, 5}, new int[] {1, 2, 3, 4, 4});
        // A code word inside the table but matching no length, which runs off the end of the
        // length scan rather than reaching the codec's own "unable to map" message.
        cross(new int[] {1, 2}, new int[] {1, 1}, new int[] {2, 2, 1, 2},
                new int[] {1, 2, 3}, new int[] {1, 2, 4});
        // An empty alphabet, the one way to reach the codec's own message.
        cross(new int[] {1, 2}, new int[] {1, 1}, new int[] {1}, new int[] {}, new int[] {});

        // What each refuses.
        errWriteUnknown(new int[] {1, 2, 3, 4}, new int[] {2, 2, 2, 2}, 9);
        errWriteUnknown(new int[] {42}, new int[] {0}, 43);
        errWriteUnknownByte(new byte[] {'A', 'C'}, new int[] {1, 1}, (byte) 'G');
        errReadEmpty(new int[] {1, 2, 3, 4}, new int[] {2, 2, 2, 2});
        errReadLength(new int[] {1, 2, 3, 4}, new int[] {2, 2, 2, 2});
        errReadLengthByte(new byte[] {'A', 'C'}, new int[] {1, 1});
    }

    static void canon(final int[] symbols, final int[] lengths) {
        final HuffmanCanoncialCodeGenerator<Integer> generator =
                new HuffmanCanoncialCodeGenerator<>(intParams(symbols, lengths));
        System.out.printf("canon\t%s\t%s\t%s%n", ints(symbols), ints(lengths),
                codes(generator.getCanonicalCodeWords()));
    }

    static void canonByte(final byte[] symbols, final int[] lengths) {
        final HuffmanCanoncialCodeGenerator<Byte> generator =
                new HuffmanCanoncialCodeGenerator<>(byteParams(symbols, lengths));
        System.out.printf("canon\tbytes:%s\t%s\t%s%n", bytes(symbols), ints(lengths),
                codes(generator.getCanonicalCodeWords()));
    }

    static void canonErr(final int[] symbols, final int[] lengths) {
        try {
            final HuffmanCanoncialCodeGenerator<Integer> generator =
                    new HuffmanCanoncialCodeGenerator<>(intParams(symbols, lengths));
            System.out.printf("canon\t%s\t%s\t%s%n", ints(symbols), ints(lengths),
                    codes(generator.getCanonicalCodeWords()));
        } catch (final Throwable t) {
            System.out.printf("canonerr\t%s\t%s\t%s\t%s%n", ints(symbols), ints(lengths),
                    t.getClass().getSimpleName(), String.valueOf(t.getMessage()));
        }
    }

    static void serInt(final int[] symbols, final int[] lengths) {
        final CanonicalHuffmanIntegerEncoding encoding =
                new CanonicalHuffmanIntegerEncoding(symbols, lengths);
        final byte[] params = encoding.toSerializedEncodingParams();
        final byte[] again = CanonicalHuffmanIntegerEncoding.fromSerializedEncodingParams(params)
                .toSerializedEncodingParams();
        System.out.printf("ser\tint\t%s\t%s\t%s\t%s%n", ints(symbols), ints(lengths), hex(params),
                hex(again));
    }

    static void serByte(final byte[] symbols, final int[] lengths) {
        final CanonicalHuffmanByteEncoding encoding =
                new CanonicalHuffmanByteEncoding(symbols, lengths);
        final byte[] params = encoding.toSerializedEncodingParams();
        final byte[] again = CanonicalHuffmanByteEncoding.fromSerializedEncodingParams(params)
                .toSerializedEncodingParams();
        System.out.printf("ser\tbyte\t%s\t%s\t%s\t%s%n", bytes(symbols), ints(lengths), hex(params),
                hex(again));
    }

    static void roundInt(final int[] symbols, final int[] lengths, final int[] values) {
        final CanonicalHuffmanIntegerEncoding encoding =
                new CanonicalHuffmanIntegerEncoding(symbols, lengths);
        final List<Integer> boxed = Arrays.stream(values).boxed().collect(Collectors.toList());
        final Result result = round(encoding, boxed);
        if (result.error != null) {
            System.out.printf("err\tround\t%s\t%s\t%s\t%s\t%s%n", ints(symbols), ints(lengths),
                    ints(values), result.error.getClass().getSimpleName(),
                    String.valueOf(result.error.getMessage()));
            return;
        }
        System.out.printf("round\tint\t%s\t%s\t%s\t%s\t%s%n", ints(symbols), ints(lengths),
                ints(values), result.hex, result.values);
    }

    static void roundByte(final byte[] symbols, final int[] lengths, final byte[] values) {
        final CanonicalHuffmanByteEncoding encoding =
                new CanonicalHuffmanByteEncoding(symbols, lengths);
        final List<Byte> boxed = new ArrayList<>(values.length);
        for (final byte value : values) {
            boxed.add(value);
        }
        final Result result = round(encoding, boxed);
        if (result.error != null) {
            System.out.printf("err\tround\tbytes:%s\t%s\t%s\t%s\t%s%n", bytes(symbols),
                    ints(lengths), bytes(values), result.error.getClass().getSimpleName(),
                    String.valueOf(result.error.getMessage()));
            return;
        }
        System.out.printf("round\tbyte\t%s\t%s\t%s\t%s\t%s%n", bytes(symbols), ints(lengths),
                bytes(values), result.hex, result.values);
    }

    /** Write with one alphabet, then read the same core block with another. */
    static void cross(final int[] writeSymbols, final int[] writeLengths, final int[] values,
            final int[] readSymbols, final int[] readLengths) {
        String outcome;
        try {
            final SliceBlocks blocks = write(
                    new CanonicalHuffmanIntegerEncoding(writeSymbols, writeLengths),
                    Arrays.stream(values).boxed().collect(Collectors.toList()));
            final SliceBlocksReadStreams readStreams =
                    new SliceBlocksReadStreams(blocks, new CompressorCache());
            final CRAMCodec<Integer> reader =
                    new CanonicalHuffmanIntegerEncoding(readSymbols, readLengths)
                            .buildCodec(readStreams, null);
            outcome = String.valueOf(reader.read());
        } catch (final Throwable t) {
            outcome = t.getClass().getSimpleName() + ": " + String.valueOf(t.getMessage());
        }
        System.out.printf("cross\t%s\t%s\t%s\t%s\t%s\t%s%n", ints(writeSymbols),
                ints(writeLengths), ints(values), ints(readSymbols), ints(readLengths), outcome);
    }

    static void errWriteUnknown(final int[] symbols, final int[] lengths, final int value) {
        try {
            write(new CanonicalHuffmanIntegerEncoding(symbols, lengths),
                    Arrays.asList(Integer.valueOf(value)));
            System.out.printf("err\twrite-unknown\t%s\t%s\t%d\t-\t-%n", ints(symbols),
                    ints(lengths), value);
        } catch (final Throwable t) {
            System.out.printf("err\twrite-unknown\t%s\t%s\t%d\t%s\t%s%n", ints(symbols),
                    ints(lengths), value, t.getClass().getSimpleName(),
                    String.valueOf(t.getMessage()));
        }
    }

    static void errWriteUnknownByte(final byte[] symbols, final int[] lengths, final byte value) {
        try {
            write(new CanonicalHuffmanByteEncoding(symbols, lengths),
                    Arrays.asList(Byte.valueOf(value)));
            System.out.printf("err\twrite-unknown\tbytes:%s\t%s\t%d\t-\t-%n", bytes(symbols),
                    ints(lengths), value);
        } catch (final Throwable t) {
            System.out.printf("err\twrite-unknown\tbytes:%s\t%s\t%d\t%s\t%s%n", bytes(symbols),
                    ints(lengths), value, t.getClass().getSimpleName(),
                    String.valueOf(t.getMessage()));
        }
    }

    /** Read from a core block that carries no bits at all. */
    static void errReadEmpty(final int[] symbols, final int[] lengths) {
        try {
            final CanonicalHuffmanIntegerEncoding encoding =
                    new CanonicalHuffmanIntegerEncoding(symbols, lengths);
            final SliceBlocks blocks = write(encoding, new ArrayList<>());
            final SliceBlocksReadStreams readStreams =
                    new SliceBlocksReadStreams(blocks, new CompressorCache());
            final CRAMCodec<Integer> reader = encoding.buildCodec(readStreams, null);
            System.out.printf("err\tread-empty\t%s\t%s\t-\t-\t%s%n", ints(symbols), ints(lengths),
                    String.valueOf(reader.read()));
        } catch (final Throwable t) {
            System.out.printf("err\tread-empty\t%s\t%s\t-\t%s\t%s%n", ints(symbols), ints(lengths),
                    t.getClass().getSimpleName(), String.valueOf(t.getMessage()));
        }
    }

    static void errReadLength(final int[] symbols, final int[] lengths) {
        try {
            final CRAMCodec<Integer> codec =
                    new CanonicalHuffmanIntegerEncoding(symbols, lengths)
                            .buildCodec(null, writeStreams());
            System.out.printf("err\tread-length\t%s\t%s\t4\t-\t%s%n", ints(symbols), ints(lengths),
                    String.valueOf(codec.read(4)));
        } catch (final Throwable t) {
            System.out.printf("err\tread-length\t%s\t%s\t4\t%s\t%s%n", ints(symbols), ints(lengths),
                    t.getClass().getSimpleName(), String.valueOf(t.getMessage()));
        }
    }

    static void errReadLengthByte(final byte[] symbols, final int[] lengths) {
        try {
            final CRAMCodec<Byte> codec = new CanonicalHuffmanByteEncoding(symbols, lengths)
                    .buildCodec(null, writeStreams());
            System.out.printf("err\tread-length\tbytes:%s\t%s\t4\t-\t%s%n", bytes(symbols),
                    ints(lengths), String.valueOf(codec.read(4)));
        } catch (final Throwable t) {
            System.out.printf("err\tread-length\tbytes:%s\t%s\t4\t%s\t%s%n", bytes(symbols),
                    ints(lengths), t.getClass().getSimpleName(), String.valueOf(t.getMessage()));
        }
    }

    static HuffmanParams<Integer> intParams(final int[] symbols, final int[] lengths) {
        return new HuffmanParams<>(Arrays.stream(symbols).boxed().collect(Collectors.toList()),
                Arrays.stream(lengths).boxed().collect(Collectors.toList()));
    }

    static HuffmanParams<Byte> byteParams(final byte[] symbols, final int[] lengths) {
        final List<Byte> boxed = new ArrayList<>(symbols.length);
        for (final byte symbol : symbols) {
            boxed.add(symbol);
        }
        return new HuffmanParams<>(boxed,
                Arrays.stream(lengths).boxed().collect(Collectors.toList()));
    }

    static SliceBlocksWriteStreams writeStreams() {
        return new SliceBlocksWriteStreams(new CompressionHeader());
    }

    static <T> SliceBlocks write(final CRAMEncoding<T> encoding, final List<T> values) {
        final SliceBlocksWriteStreams streams = writeStreams();
        final CRAMCodec<T> writer = encoding.buildCodec(null, streams);
        for (final T value : values) {
            writer.write(value);
        }
        return streams.flushStreamsToBlocks();
    }

    /** Write the values through the codec, take the raw core block, and read them back. */
    static <T> Result round(final CRAMEncoding<T> encoding, final List<T> values) {
        final Result result = new Result();
        try {
            final SliceBlocks blocks = write(encoding, values);
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

    static <T> String codes(final List<HuffmanBitCode<T>> bitCodes) {
        final StringJoiner joiner = new StringJoiner(",");
        for (final HuffmanBitCode<T> code : bitCodes) {
            joiner.add(String.format("%s:%d/%d", String.valueOf(code.getSymbol()),
                    code.getCodeWord(), code.getCodeWordBitLength()));
        }
        return joiner.length() == 0 ? "-" : joiner.toString();
    }

    static String ints(final int[] values) {
        final StringJoiner joiner = new StringJoiner(",");
        for (final int value : values) {
            joiner.add(Integer.toString(value));
        }
        return joiner.length() == 0 ? "-" : joiner.toString();
    }

    static String bytes(final byte[] values) {
        final StringJoiner joiner = new StringJoiner(",");
        for (final byte value : values) {
            joiner.add(Integer.toString(value));
        }
        return joiner.length() == 0 ? "-" : joiner.toString();
    }

    static String hex(final byte[] bytes) {
        final StringBuilder builder = new StringBuilder(bytes.length * 2);
        for (final byte value : bytes) {
            builder.append(String.format("%02x", value));
        }
        return builder.length() == 0 ? "-" : builder.toString();
    }
}
