import htsjdk.samtools.SBIIndexWriter;
import java.io.ByteArrayOutputStream;
import java.security.MessageDigest;

/**
 * The bytes `SBIIndexWriter` writes, as a digest and a length, plus the header fields a reader
 * would take back out of them.
 *
 *   sbi <granularity> <records> <finalOffset> <fileLength> <bytes> <md5> <recordsField> <offsetCount>
 *
 * The corpus varies granularity against record counts that do and do not divide by it, because the
 * offset count is `ceil(records / granularity) + 1` and the header's record count is neither of
 * those two numbers. An empty index is included: it still carries the final offset.
 */
public class SbiDump {
  static String hex(byte[] b) {
    StringBuilder sb = new StringBuilder();
    for (byte x : b) sb.append(String.format("%02x", x));
    return sb.toString();
  }

  static long readLong(byte[] bytes, int offset) {
    long value = 0;
    for (int i = 7; i >= 0; i--) {
      value = (value << 8) | (bytes[offset + i] & 0xFFL);
    }
    return value;
  }

  public static void main(String[] args) throws Exception {
    long[] granularities = {1, 2, 3, 4096};
    int[] recordCounts = {0, 1, 2, 5, 10};
    for (long granularity : granularities) {
      for (int records : recordCounts) {
        ByteArrayOutputStream sink = new ByteArrayOutputStream();
        SBIIndexWriter writer = new SBIIndexWriter(sink, granularity);
        for (int i = 0; i < records; i++) {
          // A virtual offset with a block address and an offset inside it, increasing.
          writer.processRecord(((long) i << 16) | (i % 7));
        }
        long finalOffset = ((long) records << 16);
        writer.finish(finalOffset, 1_000_000L + records);
        byte[] bytes = sink.toByteArray();
        MessageDigest md = MessageDigest.getInstance("MD5");
        md.update(bytes);
        System.out.printf(
            "sbi\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%d%n",
            granularity,
            records,
            finalOffset,
            1_000_000L + records,
            bytes.length,
            hex(md.digest()),
            readLong(bytes, 44),
            readLong(bytes, 60));
      }
    }
  }
}
