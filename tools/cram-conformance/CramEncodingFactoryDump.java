/*
 * EncodingFactory: which codec a data series type and an encoding identifier resolve to.
 *
 * A CRAM compression header names an encoding per data series, by identifier and parameters. The
 * factory is what turns that pair into a codec, and it is the last thing between a file's bytes and
 * the codecs. Forty lines of Java, and one missing keyword in them.
 *
 * Four things here are decisions rather than layout.
 *
 *   - THE SWITCH FALLS THROUGH. Only the BYTE arm ends in a break, so an INT that matches nothing
 *     falls into the LONG arm and then into the BYTE_ARRAY arm, and a LONG that matches nothing
 *     falls into BYTE_ARRAY. An INT data series named with BYTE_ARRAY_LEN therefore gets a byte
 *     array encoding rather than the refusal the last line promises;
 *   - THE REFUSAL NAMES BOTH HALVES, so a file that asks for something unreachable says which type
 *     and which identifier;
 *   - THE PARAMETERS ARE NOT VALIDATED AGAINST THE TYPE. Whatever encoding the switch lands on
 *     parses the bytes its own way, so the same parameters mean different things depending on the
 *     arm reached;
 *   - NULL IS AN IDENTIFIER LIKE ANY OTHER and matches nothing anywhere, so it always reaches the
 *     refusal.
 *
 * Every row carries the parameters it was given, in hex, so nothing here can be rebuilt from a
 * label.
 *
 * Output:
 *
 *     make\t<type>\t<encoding id>\t<params hex>\t<encoding class>\t<encoding toString>
 *     err\t<type>\t<encoding id>\t<params hex>\t<class>\t<message>
 *
 * Usage: CramEncodingFactoryDump
 */

import htsjdk.samtools.cram.encoding.CRAMEncoding;
import htsjdk.samtools.cram.encoding.EncodingFactory;
import htsjdk.samtools.cram.io.ITF8;
import htsjdk.samtools.cram.structure.DataSeriesType;
import htsjdk.samtools.cram.structure.EncodingDescriptor;
import htsjdk.samtools.cram.structure.EncodingID;

import java.io.ByteArrayOutputStream;

public class CramEncodingFactoryDump {

    public static void main(final String[] args) {
        System.out.println("# CramEncodingFactoryDump: what a type and an identifier resolve to");

        // Every pair of data series type and encoding identifier, with parameters that the
        // identifier's own encoding would accept.
        for (final DataSeriesType type : DataSeriesType.values()) {
            for (final EncodingID id : EncodingID.values()) {
                make(type, id, paramsFor(id));
            }
        }

        // The same parameters through every type, which is where the fall-through shows: one byte
        // means a content id in one arm and the first of a Huffman alphabet in another.
        for (final DataSeriesType type : DataSeriesType.values()) {
            make(type, EncodingID.EXTERNAL, itf8(7));
        }

        // More than one symbol, and a nested pair, to pin how each prints itself.
        make(DataSeriesType.INT, EncodingID.HUFFMAN, itf8(3, 1, 2, 3, 3, 1, 2, 2));
        make(DataSeriesType.BYTE, EncodingID.HUFFMAN, concat(itf8(2), new byte[] {0x41, 0x43},
                itf8(2, 1, 1)));
        make(DataSeriesType.BYTE_ARRAY, EncodingID.BYTE_ARRAY_LEN,
                concat(itf8(3, 4), itf8(1, 42, 1, 0), itf8(5, 2), new byte[] {0x00, 0x02}));

        // A descriptor rather than a loose pair, which is the way a compression header supplies it.
        makeFromDescriptor(DataSeriesType.INT, EncodingID.EXTERNAL, itf8(3));
        makeFromDescriptor(DataSeriesType.BYTE_ARRAY, EncodingID.BYTE_ARRAY_STOP,
                new byte[] {0x00, 0x01});
    }

    /** Parameters an encoding of this identifier parses without complaint. */
    static byte[] paramsFor(final EncodingID id) {
        switch (id) {
            case EXTERNAL:
                return itf8(1);
            case GOLOMB:
            case GOLOMB_RICE:
                return itf8(0, 4);
            case HUFFMAN:
                // One symbol, one length: the shortest alphabet there is.
                return itf8(1, 42, 1, 0);
            case BYTE_ARRAY_LEN:
                // An external integer for the length, an external byte array for the bytes.
                return concat(itf8(1, 1), itf8(1), itf8(1, 1), itf8(2));
            case BYTE_ARRAY_STOP:
                return new byte[] {0x00, 0x01};
            case BETA:
                return itf8(0, 8);
            case SUBEXPONENTIAL:
                return itf8(0, 2);
            case GAMMA:
                return itf8(0);
            default:
                return new byte[0];
        }
    }

    static void make(final DataSeriesType type, final EncodingID id, final byte[] params) {
        try {
            final CRAMEncoding<?> encoding = EncodingFactory.createCRAMEncoding(type, id, params);
            System.out.printf("make\t%s\t%s\t%s\t%s\t%s%n", type.name(), id.name(), hex(params),
                    encoding.getClass().getSimpleName(), String.valueOf(encoding));
        } catch (final Throwable t) {
            System.out.printf("err\t%s\t%s\t%s\t%s\t%s%n", type.name(), id.name(), hex(params),
                    t.getClass().getSimpleName(), String.valueOf(t.getMessage()));
        }
    }

    static void makeFromDescriptor(final DataSeriesType type, final EncodingID id,
            final byte[] params) {
        try {
            final CRAMEncoding<?> encoding = EncodingFactory.createCRAMEncoding(type,
                    new EncodingDescriptor(id, params));
            System.out.printf("make\t%s\t%s\t%s\t%s\t%s%n", type.name(), id.name(), hex(params),
                    encoding.getClass().getSimpleName(), String.valueOf(encoding));
        } catch (final Throwable t) {
            System.out.printf("err\t%s\t%s\t%s\t%s\t%s%n", type.name(), id.name(), hex(params),
                    t.getClass().getSimpleName(), String.valueOf(t.getMessage()));
        }
    }

    static byte[] itf8(final int... values) {
        final ByteArrayOutputStream out = new ByteArrayOutputStream();
        for (final int value : values) {
            ITF8.writeUnsignedITF8(value, out);
        }
        return out.toByteArray();
    }

    static byte[] concat(final byte[]... parts) {
        final ByteArrayOutputStream out = new ByteArrayOutputStream();
        for (final byte[] part : parts) {
            out.write(part, 0, part.length);
        }
        return out.toByteArray();
    }

    static String hex(final byte[] bytes) {
        final StringBuilder builder = new StringBuilder(bytes.length * 2);
        for (final byte value : bytes) {
            builder.append(String.format("%02x", value));
        }
        return builder.length() == 0 ? "-" : builder.toString();
    }
}
