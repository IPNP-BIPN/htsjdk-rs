/*
 * The codecs written on external blocks rather than on the core bit stream.
 *
 * Six of the encoding map's identifiers put their data in a block of its own, named by a content
 * id, and two more are built out of the first six. A block is bytes, not bits, so these codecs
 * carry no alignment problem and no offset; what they carry instead is a set of decisions about
 * where one value ends and the next begins.
 *
 * Seven things here are decisions rather than layout.
 *
 *   - EACH CODEC NAMES A BLOCK, and two codecs naming the same content id share one. That is what
 *     lets ByteArrayLen put its lengths in one block and its bytes in another, or in the same one;
 *   - INTEGER IS ITF8 AND LONG IS LTF8, straight onto the block, so an external integer costs
 *     between one and five bytes and nothing marks which;
 *   - EXTERNAL BYTE CANNOT SEE THE END OF ITS BLOCK. It returns `(byte) stream.read()`, and at the
 *     end that is `(byte) -1`, which is indistinguishable from a byte of 0xFF that is really
 *     there;
 *   - EXTERNAL BYTE ARRAY HAS NO LENGTH OF ITS OWN, so its no-argument read refuses outright and
 *     only a caller who already knows the length can use it;
 *   - BYTE ARRAY STOP TRUSTS THE DATA. It appends a stop byte after each array and reads until it
 *     sees one, so an array containing that byte is split in two and nothing reports it;
 *   - A STOPPED ARRAY THAT RUNS OUT OF BLOCK IS NOT AN ERROR EITHER. The read ends on the end of
 *     the stream exactly as it would on a stop byte;
 *   - BYTE ARRAY LEN IS A PAIR OF CODECS, so a length can live on the core bit stream while the
 *     bytes live in an external block. The two halves are independent encodings.
 *
 * Output:
 *
 *     ext\t<flavour>\t<id>\t<values>\t<blocks>\t<values read back>
 *     stop\t<stop>\t<id>\t<arrays>\t<blocks>\t<arrays read back>
 *     len\t<length encoding>\t<byte encoding>\t<arrays>\t<blocks>\t<arrays read back>
 *     ser\t<encoding>\t<params>\t<hex>\t<reparsed hex>
 *     err\t<what>\t<detail>\t<class>\t<message>
 *
 * A block column is `core=<hex>;<id>:<method>=<hex>` over every block the write produced, in
 * content id order with the core block first, an external block reported after its compression is
 * undone. An empty byte array prints as `.`, and so does an empty block.
 *
 * Usage: CramExternalCodecDump
 */

import htsjdk.samtools.cram.encoding.ByteArrayLenEncoding;
import htsjdk.samtools.cram.encoding.CRAMCodec;
import htsjdk.samtools.cram.encoding.CRAMEncoding;
import htsjdk.samtools.cram.encoding.core.CanonicalHuffmanIntegerEncoding;
import htsjdk.samtools.cram.encoding.external.ByteArrayStopEncoding;
import htsjdk.samtools.cram.encoding.external.ExternalByteArrayEncoding;
import htsjdk.samtools.cram.encoding.external.ExternalByteEncoding;
import htsjdk.samtools.cram.encoding.external.ExternalIntegerEncoding;
import htsjdk.samtools.cram.encoding.external.ExternalLongEncoding;
import htsjdk.samtools.cram.structure.CompressionHeader;
import htsjdk.samtools.cram.structure.CompressorCache;
import htsjdk.samtools.cram.structure.SliceBlocks;
import htsjdk.samtools.cram.structure.block.Block;
import htsjdk.samtools.cram.structure.SliceBlocksReadStreams;
import htsjdk.samtools.cram.structure.SliceBlocksWriteStreams;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.StringJoiner;

public class CramExternalCodecDump {

    public static void main(final String[] args) {
        System.out.println("# CramExternalCodecDump: the codecs written on external blocks");

        methods();

        // Integer: ITF8 onto the block, at every width it has.
        extInt(1, new int[] {0, 1, 127});
        extInt(1, new int[] {128, 16383});
        extInt(1, new int[] {16384, 2097151});
        extInt(1, new int[] {2097152, 268435455});
        extInt(1, new int[] {268435456, Integer.MAX_VALUE});
        extInt(1, new int[] {-1, -2, Integer.MIN_VALUE});
        extInt(7, new int[] {1, 2, 3});

        // Long: LTF8, which has three more widths than ITF8.
        extLong(1, new long[] {0L, 1L, 127L});
        extLong(1, new long[] {128L, 4294967295L});
        extLong(1, new long[] {4294967296L, Long.MAX_VALUE});
        extLong(1, new long[] {-1L, Long.MIN_VALUE});

        // Byte: one byte per value, and no way to see the end of the block.
        extByte(1, new byte[] {0, 1, 127});
        extByte(1, new byte[] {(byte) 0x80, (byte) 0xff});
        extByte(1, new byte[] {'A'});

        // Byte array: no length of its own, so the caller supplies one.
        extBytes(1, new byte[][] {{1, 2, 3}, {4, 5}});
        extBytes(1, new byte[][] {{}, {'A'}});
        extBytes(1, new byte[][] {{(byte) 0xff, 0x00}});

        // Byte array stop: a separator that the data is trusted not to contain.
        stop(1, (byte) 0x00, new byte[][] {{1, 2, 3}, {4, 5}});
        stop(1, (byte) 0x00, new byte[][] {{}, {}});
        stop(1, (byte) 0x00, new byte[][] {{'A', 'B'}, {}, {'C'}});
        stop(1, (byte) 0x09, new byte[][] {{'A', 'B'}, {'C'}});
        stop(1, (byte) 0xff, new byte[][] {{1, 2}, {3}});
        // The stop byte inside the data, which nothing reports.
        stop(1, (byte) 0x00, new byte[][] {{1, 0, 2}});
        stop(1, (byte) 0x09, new byte[][] {{'A', 0x09, 'B'}, {'C'}});

        // Byte array len: two encodings, and they need not share a block or even a kind of block.
        len("ext-int:1", "ext-bytes:2", new ExternalIntegerEncoding(1),
                new ExternalByteArrayEncoding(2), new byte[][] {{1, 2, 3}, {4, 5}});
        len("ext-int:1", "ext-bytes:1", new ExternalIntegerEncoding(1),
                new ExternalByteArrayEncoding(1), new byte[][] {{1, 2, 3}, {4, 5}});
        len("ext-int:1", "ext-bytes:2", new ExternalIntegerEncoding(1),
                new ExternalByteArrayEncoding(2), new byte[][] {{}, {'A'}});
        len("huffman:2,3", "ext-bytes:1",
                new CanonicalHuffmanIntegerEncoding(new int[] {2, 3}, new int[] {1, 1}),
                new ExternalByteArrayEncoding(1), new byte[][] {{1, 2}, {3, 4, 5}});
        len("huffman:3", "ext-bytes:1",
                new CanonicalHuffmanIntegerEncoding(new int[] {3}, new int[] {0}),
                new ExternalByteArrayEncoding(1), new byte[][] {{1, 2, 3}, {4, 5, 6}});
        len("ext-int:1", "stop:0/2", new ExternalIntegerEncoding(1),
                new ByteArrayStopEncoding((byte) 0x00, 2), new byte[][] {{1, 2}, {3}});

        // The encoding parameters, and what they carry back.
        serExt("ext-int", "id=1", new ExternalIntegerEncoding(1),
                params -> ExternalIntegerEncoding.fromSerializedEncodingParams(params));
        serExt("ext-int", "id=128", new ExternalIntegerEncoding(128),
                params -> ExternalIntegerEncoding.fromSerializedEncodingParams(params));
        serExt("ext-byte", "id=2", new ExternalByteEncoding(2),
                params -> ExternalByteEncoding.fromSerializedEncodingParams(params));
        serExt("ext-long", "id=3", new ExternalLongEncoding(3),
                params -> ExternalLongEncoding.fromSerializedEncodingParams(params));
        serExt("ext-bytes", "id=4", new ExternalByteArrayEncoding(4),
                params -> ExternalByteArrayEncoding.fromSerializedEncodingParams(params));
        serStop((byte) 0x00, 1);
        serStop((byte) 0xff, 300);
        serLen("ext-int:1", "ext-bytes:2",
                new ByteArrayLenEncoding(new ExternalIntegerEncoding(1),
                        new ExternalByteArrayEncoding(2)));

        // What each refuses, and the one that does not refuse when it might.
        errBytesUnknownLength();
        errReadLength("ext-int", new ExternalIntegerEncoding(1));
        errReadLength("ext-byte", new ExternalByteEncoding(1));
        errReadLength("ext-long", new ExternalLongEncoding(1));
        errStopReadLength();
        errLenReadLength();
        errBytesPastEnd();
        errBytePastEnd();
        errStopPastEnd();
    }

    static void extInt(final int id, final int[] values) {
        final ExternalIntegerEncoding encoding = new ExternalIntegerEncoding(id);
        final SliceBlocks blocks = write(encoding, boxInts(values));
        final StringJoiner back = new StringJoiner(",");
        final CRAMCodec<Integer> reader = encoding.buildCodec(readStreams(blocks), null);
        for (int i = 0; i < values.length; i++) {
            back.add(String.valueOf(reader.read()));
        }
        System.out.printf("ext\tint\t%d\t%s\t%s\t%s%n", id, ints(values), blocks(blocks), back);
    }

    static void extLong(final int id, final long[] values) {
        final ExternalLongEncoding encoding = new ExternalLongEncoding(id);
        final List<Long> boxed = new ArrayList<>();
        for (final long value : values) {
            boxed.add(value);
        }
        final SliceBlocks blocks = write(encoding, boxed);
        final StringJoiner back = new StringJoiner(",");
        final CRAMCodec<Long> reader = encoding.buildCodec(readStreams(blocks), null);
        for (int i = 0; i < values.length; i++) {
            back.add(String.valueOf(reader.read()));
        }
        final StringJoiner shown = new StringJoiner(",");
        for (final long value : values) {
            shown.add(Long.toString(value));
        }
        System.out.printf("ext\tlong\t%d\t%s\t%s\t%s%n", id, shown, blocks(blocks), back);
    }

    static void extByte(final int id, final byte[] values) {
        final ExternalByteEncoding encoding = new ExternalByteEncoding(id);
        final List<Byte> boxed = new ArrayList<>();
        for (final byte value : values) {
            boxed.add(value);
        }
        final SliceBlocks blocks = write(encoding, boxed);
        final StringJoiner back = new StringJoiner(",");
        final CRAMCodec<Byte> reader = encoding.buildCodec(readStreams(blocks), null);
        for (int i = 0; i < values.length; i++) {
            back.add(String.valueOf(reader.read()));
        }
        System.out.printf("ext\tbyte\t%d\t%s\t%s\t%s%n", id, bytes(values), blocks(blocks), back);
    }

    static void extBytes(final int id, final byte[][] values) {
        final ExternalByteArrayEncoding encoding = new ExternalByteArrayEncoding(id);
        final SliceBlocks blocks = write(encoding, Arrays.asList(values));
        final StringJoiner back = new StringJoiner(",");
        final CRAMCodec<byte[]> reader = encoding.buildCodec(readStreams(blocks), null);
        for (final byte[] value : values) {
            back.add(hex(reader.read(value.length)));
        }
        System.out.printf("ext\tbytes\t%d\t%s\t%s\t%s%n", id, arrays(values), blocks(blocks), back);
    }

    static void stop(final int id, final byte stopByte, final byte[][] values) {
        final ByteArrayStopEncoding encoding = new ByteArrayStopEncoding(stopByte, id);
        final SliceBlocks blocks = write(encoding, Arrays.asList(values));
        final StringJoiner back = new StringJoiner(",");
        final CRAMCodec<byte[]> reader = encoding.buildCodec(readStreams(blocks), null);
        for (int i = 0; i < values.length; i++) {
            back.add(hex(reader.read()));
        }
        System.out.printf("stop\t%d\t%d\t%s\t%s\t%s%n", stopByte & 0xff, id, arrays(values),
                blocks(blocks), back);
    }

    static void len(final String lenName, final String byteName,
            final CRAMEncoding<Integer> lenEncoding, final CRAMEncoding<byte[]> byteEncoding,
            final byte[][] values) {
        final ByteArrayLenEncoding encoding = new ByteArrayLenEncoding(lenEncoding, byteEncoding);
        final SliceBlocks blocks = write(encoding, Arrays.asList(values));
        final StringJoiner back = new StringJoiner(",");
        try {
            final CRAMCodec<byte[]> reader = encoding.buildCodec(readStreams(blocks), null);
            for (int i = 0; i < values.length; i++) {
                back.add(hex(reader.read()));
            }
        } catch (final Throwable t) {
            System.out.printf("err\tlen-read\t%s %s\t%s\t%s%n", lenName, byteName,
                    t.getClass().getSimpleName(), String.valueOf(t.getMessage()));
            return;
        }
        System.out.printf("len\t%s\t%s\t%s\t%s\t%s%n", lenName, byteName, arrays(values),
                blocks(blocks), back);
    }

    static void serExt(final String name, final String params, final CRAMEncoding<?> encoding,
            final java.util.function.Function<byte[], CRAMEncoding<?>> reparse) {
        final byte[] serialized = encoding.toSerializedEncodingParams();
        final byte[] again = reparse.apply(serialized).toSerializedEncodingParams();
        System.out.printf("ser\t%s\t%s\t%s\t%s%n", name, params, hex(serialized), hex(again));
    }

    static void serStop(final byte stopByte, final int id) {
        final ByteArrayStopEncoding encoding = new ByteArrayStopEncoding(stopByte, id);
        final byte[] serialized = encoding.toSerializedEncodingParams();
        final byte[] again = ByteArrayStopEncoding.fromSerializedEncodingParams(serialized)
                .toSerializedEncodingParams();
        System.out.printf("ser\tstop\tstop=%d id=%d\t%s\t%s%n", stopByte & 0xff, id,
                hex(serialized), hex(again));
    }

    static void serLen(final String lenName, final String byteName,
            final ByteArrayLenEncoding encoding) {
        final byte[] serialized = encoding.toSerializedEncodingParams();
        final byte[] again = ByteArrayLenEncoding.fromSerializedEncodingParams(serialized)
                .toSerializedEncodingParams();
        System.out.printf("ser\tlen\t%s %s\t%s\t%s%n", lenName, byteName, hex(serialized),
                hex(again));
    }

    /** The byte array codec has no length of its own, so its no-argument read refuses. */
    static void errBytesUnknownLength() {
        try {
            final ExternalByteArrayEncoding encoding = new ExternalByteArrayEncoding(1);
            final SliceBlocks blocks = write(encoding, Arrays.asList(new byte[] {1, 2}));
            final CRAMCodec<byte[]> reader = encoding.buildCodec(readStreams(blocks), null);
            System.out.printf("err\tbytes-unknown-length\tblock=0102\t-\t%s%n", hex(reader.read()));
        } catch (final Throwable t) {
            System.out.printf("err\tbytes-unknown-length\tblock=0102\t%s\t%s%n",
                    t.getClass().getSimpleName(), String.valueOf(t.getMessage()));
        }
    }

    static void errReadLength(final String name, final CRAMEncoding<?> encoding) {
        try {
            final CRAMCodec<?> codec = encoding.buildCodec(null, writeStreams());
            System.out.printf("err\tread-length\t%s\t-\t%s%n", name, String.valueOf(codec.read(4)));
        } catch (final Throwable t) {
            System.out.printf("err\tread-length\t%s\t%s\t%s%n", name,
                    t.getClass().getSimpleName(), String.valueOf(t.getMessage()));
        }
    }

    static void errStopReadLength() {
        errReadLength("stop", new ByteArrayStopEncoding((byte) 0x00, 1));
    }

    static void errLenReadLength() {
        errReadLength("len", new ByteArrayLenEncoding(new ExternalIntegerEncoding(1),
                new ExternalByteArrayEncoding(2)));
    }

    /** More bytes asked for than the block holds. */
    static void errBytesPastEnd() {
        try {
            final ExternalByteArrayEncoding encoding = new ExternalByteArrayEncoding(1);
            final SliceBlocks blocks = write(encoding, Arrays.asList(new byte[] {1, 2}));
            final CRAMCodec<byte[]> reader = encoding.buildCodec(readStreams(blocks), null);
            System.out.printf("err\tbytes-past-end\tblock=0102 length=4\t-\t%s%n",
                    hex(reader.read(4)));
        } catch (final Throwable t) {
            System.out.printf("err\tbytes-past-end\tblock=0102 length=4\t%s\t%s%n",
                    t.getClass().getSimpleName(),
                    String.valueOf(t.getMessage()));
        }
    }

    /** A byte read past the end of the block, which is not an error at all. */
    static void errBytePastEnd() {
        try {
            final ExternalByteEncoding encoding = new ExternalByteEncoding(1);
            final SliceBlocks blocks = write(encoding, Arrays.asList((byte) 0x41));
            final CRAMCodec<Byte> reader = encoding.buildCodec(readStreams(blocks), null);
            reader.read();
            System.out.printf("err\tbyte-past-end\tblock=41\t-\t%s%n", String.valueOf(reader.read()));
        } catch (final Throwable t) {
            System.out.printf("err\tbyte-past-end\tblock=41\t%s\t%s%n", t.getClass().getSimpleName(),
                    String.valueOf(t.getMessage()));
        }
    }

    /** A stopped array whose block ends before a stop byte does. */
    static void errStopPastEnd() {
        try {
            final ByteArrayStopEncoding encoding = new ByteArrayStopEncoding((byte) 0x00, 1);
            final SliceBlocks blocks = write(new ExternalByteArrayEncoding(1),
                    Arrays.asList(new byte[] {'A', 'B'}));
            final CRAMCodec<byte[]> reader = encoding.buildCodec(readStreams(blocks), null);
            final String first = hex(reader.read());
            System.out.printf("err\tstop-past-end\tblock=4142 stop=0\t-\t%s then %s%n", first,
                    hex(reader.read()));
        } catch (final Throwable t) {
            System.out.printf("err\tstop-past-end\tblock=4142 stop=0\t%s\t%s%n",
                    t.getClass().getSimpleName(),
                    String.valueOf(t.getMessage()));
        }
    }

    static List<Integer> boxInts(final int[] values) {
        final List<Integer> boxed = new ArrayList<>(values.length);
        for (final int value : values) {
            boxed.add(value);
        }
        return boxed;
    }

    static SliceBlocksWriteStreams writeStreams() {
        return new SliceBlocksWriteStreams(new CompressionHeader());
    }

    static SliceBlocksReadStreams readStreams(final SliceBlocks blocks) {
        return new SliceBlocksReadStreams(blocks, new CompressorCache());
    }

    static <T> SliceBlocks write(final CRAMEncoding<T> encoding, final List<T> values) {
        final SliceBlocksWriteStreams streams = writeStreams();
        final CRAMCodec<T> writer = encoding.buildCodec(null, streams);
        for (final T value : values) {
            writer.write(value);
        }
        return streams.flushStreamsToBlocks();
    }

    /** The blocks a write actually put bytes in, core first and then the externals in id order.
     *
     * A write stream is created for every data series the compression header knows, so most of the
     * blocks a single codec produces are empty and only the ones it wrote to say anything. The core
     * block is raw by definition; an external block is not, and what a row records is the content
     * after its compression is undone, which is what the codec wrote.
     */
    static String blocks(final SliceBlocks blocks) {
        final StringJoiner joiner = new StringJoiner(";");
        final byte[] core = blocks.getCoreBlock().getRawContent();
        if (core.length > 0) {
            joiner.add("core=" + hex(core));
        }
        final List<Integer> ids = new ArrayList<>(blocks.getExternalContentIDs());
        ids.sort(null);
        final CompressorCache cache = new CompressorCache();
        for (final Integer id : ids) {
            final Block block = blocks.getExternalBlock(id);
            final byte[] content = block.getUncompressedContent(cache);
            if (content.length > 0) {
                joiner.add(id + "=" + hex(content));
            }
        }
        return joiner.length() == 0 ? "-" : joiner.toString();
    }

    /** The compressor the default compression header picks for each external block. */
    static void methods() {
        final SliceBlocks blocks = write(new ExternalIntegerEncoding(1), boxInts(new int[] {1}));
        final List<Integer> ids = new ArrayList<>(blocks.getExternalContentIDs());
        ids.sort(null);
        final StringJoiner joiner = new StringJoiner(";");
        for (final Integer id : ids) {
            joiner.add(id + "=" + blocks.getExternalBlock(id).getCompressionMethod());
        }
        System.out.printf("methods\t%s%n", joiner);
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

    static String arrays(final byte[][] values) {
        final StringJoiner joiner = new StringJoiner(",");
        for (final byte[] value : values) {
            joiner.add(hex(value));
        }
        return joiner.length() == 0 ? "-" : joiner.toString();
    }

    /** An empty array is `.`, so a row can tell one apart from a missing one. */
    static String hex(final byte[] bytes) {
        if (bytes == null) {
            return "null";
        }
        final StringBuilder builder = new StringBuilder(bytes.length * 2);
        for (final byte value : bytes) {
            builder.append(String.format("%02x", value));
        }
        return builder.length() == 0 ? "." : builder.toString();
    }
}
