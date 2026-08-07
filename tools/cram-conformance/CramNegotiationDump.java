/*
 * Codec negotiation: which encoding a writer picks for each data series, and which compressor.
 *
 * The reader takes what a file names. The writer chooses, and the choice is what makes one CRAM
 * of a set of records rather than another. Two things do the choosing: a fixed table for the
 * thirty data series, and a measurement for the tag series, which are not known until the records
 * are.
 *
 * Five things here are decisions rather than layout.
 *
 *   - THE DATA SERIES TABLE IS FIXED, not derived from the data. Every series gets an external
 *     encoding on a content id of its own, whatever the records hold;
 *   - THE COMPRESSOR IS CHOSEN BY TRYING ALL THREE. GZIP, rANS order 0 and rANS order 1 are each
 *     run over the data and the smallest output wins;
 *   - THE TIE-BREAK IS rANS 0, THEN rANS 1, THEN GZIP, which is the order of the comparisons and
 *     not the order of the compressions. Equal sizes therefore do not pick the first one tried, and
 *     every row carries all three lengths so the rule can be checked without a gzip;
 *   - A TAG'S ENCODING COMES FROM ITS TYPE, and for the two variable-length types from the RANGE
 *     of its values: one size and it is a length-prefixed array with a zero-bit Huffman length,
 *     several sizes and a Z becomes a stop-byte array while a B stays length-prefixed;
 *   - THE STOP BYTE IS A TAB, chosen and not searched for, so a Z tag whose text contains a tab is
 *     split by its own encoding.
 *
 * Output:
 *
 *     series\t<name>\t<type>\t<encoding>\t<params hex>
 *     compressor\t<label>\t<data hex or length>\t<gzip length>\t<rANS 0 length>\t<rANS 1 length>\t<chosen>
 *     tag\t<tag>\t<values>\t<encoding>\t<params hex>
 *     dictionary\t<label>\t<records' tags>\t<dictionary>
 *     err\t<what>\t<detail>\t<class>\t<message>
 *
 * Usage: CramNegotiationDump
 */

import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.cram.build.CompressionHeaderFactory;
import htsjdk.samtools.cram.structure.CRAMCompressionRecord;
import htsjdk.samtools.cram.structure.CRAMEncodingStrategy;
import htsjdk.samtools.cram.structure.CompressionHeader;
import htsjdk.samtools.cram.structure.CompressionHeaderEncodingMap;
import htsjdk.samtools.cram.structure.CompressorCache;
import htsjdk.samtools.cram.structure.block.BlockCompressionMethod;
import htsjdk.samtools.cram.structure.DataSeries;
import htsjdk.samtools.cram.structure.EncodingDescriptor;
import htsjdk.samtools.cram.structure.ReadTag;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

public class CramNegotiationDump {

    public static void main(final String[] args) {
        System.out.println("# CramNegotiationDump: what a writer picks, and why");

        // The fixed table: every data series the default strategy names.
        final CompressionHeaderEncodingMap map =
                new CompressionHeaderEncodingMap(new CRAMEncodingStrategy());
        for (final DataSeries series : DataSeries.values()) {
            final EncodingDescriptor descriptor = map.getEncodingDescriptorForDataSeries(series);
            System.out.printf("series\t%s\t%s\t%s\t%s%n", series.getCanonicalName(),
                    series.getType().name(),
                    descriptor == null ? "-" : descriptor.getEncodingID().name(),
                    descriptor == null ? "-" : hex(descriptor.getEncodingParameters()));
        }

        // The compressor, chosen by trying all three.
        compressor("empty", new byte[0]);
        compressor("one-byte", new byte[] {0x41});
        compressor("all-the-same", repeat((byte) 0x41, 1000));
        compressor("two-symbols", alternating(1000));
        compressor("ascending", ascending(256));
        compressor("ascending-x4", ascending(1024));
        compressor("text", "the quick brown fox jumps over the lazy dog".getBytes());

        // A tag's encoding, by its type and by the range of its values.
        tag("XAc", 'c', Arrays.asList(new byte[] {1}, new byte[] {2}));
        tag("XAC", 'C', Arrays.asList(new byte[] {1}));
        tag("XAs", 's', Arrays.asList(new byte[] {1, 2}));
        tag("XAi", 'i', Arrays.asList(new byte[] {1, 2, 3, 4}));
        tag("XAf", 'f', Arrays.asList(new byte[] {1, 2, 3, 4}));
        tag("XAA", 'A', Arrays.asList(new byte[] {0x41}));
        // Z of one size, then of two.
        tag("XAZ", 'Z', Arrays.asList("abc\0".getBytes(), "def\0".getBytes()));
        tag("XBZ", 'Z', Arrays.asList("abc\0".getBytes(), "defgh\0".getBytes()));
        // B of one size, then of two, then of two both above the hundred-byte threshold.
        tag("XAB", 'B', Arrays.asList(new byte[] {'c', 2, 0, 0, 0, 1, 2}));
        tag("XBB", 'B', Arrays.asList(new byte[] {'c', 1, 0, 0, 0, 1},
                new byte[] {'c', 2, 0, 0, 0, 1, 2}));
        tag("XCB", 'B', Arrays.asList(bigArray(101), bigArray(102)));

        // The dictionary, which is what a record's tag list index points into.
        dictionary("one-record-two-tags", Arrays.asList(Arrays.asList("XAc", "XBc")));
        dictionary("two-records-same-tags", Arrays.asList(Arrays.asList("XAc", "XBc"),
                Arrays.asList("XAc", "XBc")));
        dictionary("two-records-different-tags", Arrays.asList(Arrays.asList("XAc"),
                Arrays.asList("XBc")));
        dictionary("same-tags-opposite-order", Arrays.asList(Arrays.asList("XAc", "XBc"),
                Arrays.asList("XBc", "XAc")));
    }

    /** The three candidates' lengths as well as the winner, because the rule is the lengths. */
    static void compressor(final String label, final byte[] data) {
        final CRAMEncodingStrategy strategy = new CRAMEncodingStrategy();
        final CompressionHeaderFactory factory = new CompressionHeaderFactory(strategy);
        final CompressorCache cache = new CompressorCache();
        final int gzip = cache.getCompressorForMethod(BlockCompressionMethod.GZIP,
                strategy.getGZIPCompressionLevel()).compress(data).length;
        final int rans0 = cache.getCompressorForMethod(BlockCompressionMethod.RANS, 0)
                .compress(data).length;
        final int rans1 = cache.getCompressorForMethod(BlockCompressionMethod.RANS, 1)
                .compress(data).length;
        System.out.printf("compressor\t%s\t%s\t%d\t%d\t%d\t%s%n", label,
                data.length <= 32 ? hex(data) : data.length + " bytes", gzip, rans0, rans1,
                factory.getBestExternalCompressor(data).getMethod().name());
    }

    /** One record per value, each carrying the tag, through the whole header factory. */
    static void tag(final String name, final char type, final List<byte[]> values) {
        final List<CRAMCompressionRecord> records = new ArrayList<>();
        final int tagId = ReadTag.name3BytesToInt(new byte[] {(byte) name.charAt(0),
                (byte) name.charAt(1), (byte) type});
        long index = 0;
        for (final byte[] value : values) {
            records.add(recordWithTags(index++, Arrays.asList(
                    new ReadTag(tagId, value, ValidationStringency.SILENT))));
        }

        final CompressionHeader header = new CompressionHeaderFactory(new CRAMEncodingStrategy())
                .createCompressionHeader(records, true);
        final EncodingDescriptor descriptor = header.getTagEncodingMap().get(tagId);
        final StringJoiner shown = new StringJoiner(",");
        for (final byte[] value : values) {
            shown.add(hex(value));
        }
        System.out.printf("tag\t%s%c\t%s\t%s\t%s%n", name.substring(0, 2), type, shown,
                descriptor == null ? "-" : descriptor.getEncodingID().name(),
                descriptor == null ? "-" : hex(descriptor.getEncodingParameters()));
    }

    /** The tag id dictionary a set of records produces, and the index each record is given. */
    static void dictionary(final String label, final List<List<String>> recordTags) {
        final List<CRAMCompressionRecord> records = new ArrayList<>();
        long index = 0;
        final StringJoiner input = new StringJoiner(";");
        for (final List<String> names : recordTags) {
            final List<ReadTag> tags = new ArrayList<>();
            final StringJoiner shown = new StringJoiner(",");
            for (final String name : names) {
                final int tagId = ReadTag.name3BytesToInt(new byte[] {(byte) name.charAt(0),
                        (byte) name.charAt(1), (byte) name.charAt(2)});
                tags.add(new ReadTag(tagId, new byte[] {1}, ValidationStringency.SILENT));
                shown.add(name);
            }
            input.add(shown.toString());
            records.add(recordWithTags(index++, tags));
        }

        final CompressionHeader header = new CompressionHeaderFactory(new CRAMEncodingStrategy())
                .createCompressionHeader(records, true);
        final StringJoiner dictionary = new StringJoiner(";");
        for (final byte[][] group : header.getTagIDDictionary()) {
            final StringBuilder builder = new StringBuilder();
            for (final byte[] id : group) {
                builder.append(new String(id));
            }
            dictionary.add(builder.length() == 0 ? "." : builder.toString());
        }
        final StringJoiner indexes = new StringJoiner(",");
        for (final CRAMCompressionRecord record : records) {
            indexes.add(Integer.toString(record.getTagIdsIndex().value));
        }
        System.out.printf("dictionary\t%s\t%s\t%s\tindexes=%s%n", label, input,
                dictionary.length() == 0 ? "-" : dictionary.toString(), indexes);
    }

    static CRAMCompressionRecord recordWithTags(final long index, final List<ReadTag> tags) {
        final byte[] bases = new byte[10];
        Arrays.fill(bases, (byte) 'A');
        final byte[] scores = new byte[10];
        Arrays.fill(scores, (byte) 30);
        return new CRAMCompressionRecord(index, 0, CRAMCompressionRecord.CF_DETACHED,
                "r" + index, 10, 0, 100, 0, 40, scores, bases, tags, null, -1, 0, -1, 0, -1);
    }

    static byte[] repeat(final byte value, final int length) {
        final byte[] bytes = new byte[length];
        Arrays.fill(bytes, value);
        return bytes;
    }

    static byte[] alternating(final int length) {
        final byte[] bytes = new byte[length];
        for (int i = 0; i < length; i++) {
            bytes[i] = (byte) (i % 2 == 0 ? 'A' : 'C');
        }
        return bytes;
    }

    static byte[] ascending(final int length) {
        final byte[] bytes = new byte[length];
        for (int i = 0; i < length; i++) {
            bytes[i] = (byte) i;
        }
        return bytes;
    }

    /** A B array whose payload is long enough to be over the stop-encoding threshold. */
    static byte[] bigArray(final int payload) {
        final byte[] bytes = new byte[payload + 5];
        bytes[0] = 'c';
        bytes[1] = (byte) payload;
        for (int i = 5; i < bytes.length; i++) {
            bytes[i] = (byte) i;
        }
        return bytes;
    }

    static String hex(final byte[] bytes) {
        final StringBuilder builder = new StringBuilder(bytes.length * 2);
        for (final byte value : bytes) {
            builder.append(String.format("%02x", value));
        }
        return builder.length() == 0 ? "-" : builder.toString();
    }
}
