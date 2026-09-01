import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.security.MessageDigest;

/**
 * The reference side of the I/O floor benchmark: htsjdk's BGZF writer and reader, timed on the
 * same payload the Rust side uses, in the same container.
 *
 * Prints one `name=value` line per measurement so the runner can read it without parsing prose,
 * and an md5 of every stream it writes so the comparison is byte equality rather than "both
 * finished". A speed number for a path whose bytes were never compared is not a measurement of
 * this port.
 */
public class BgzfBench {
  /** The same 64-bit LCG the conformance harnesses use, so the payload is one both sides agree on. */
  static byte[] lcg(int n, long seed, int shift) {
    byte[] b = new byte[n];
    long s = seed;
    for (int i = 0; i < n; i++) {
      s = s * 6364136223846793005L + 1442695040888963407L;
      b[i] = (byte) (s >>> shift);
    }
    return b;
  }

  /** Repetitive, SAM-shaped text: the compressible end of what a real BAM's blocks look like. */
  static byte[] text(int n) {
    byte[] p = "ACGTNacgtn\tSAMrecord\tRG:Z:rg1\n".getBytes();
    byte[] b = new byte[n];
    for (int i = 0; i < n; i++) b[i] = p[i % p.length];
    return b;
  }

  static String hex(byte[] b) {
    StringBuilder s = new StringBuilder();
    for (byte x : b) s.append(String.format("%02x", x));
    return s.toString();
  }

  static String md5(byte[] b, int len) throws Exception {
    MessageDigest md = MessageDigest.getInstance("MD5");
    md.update(b, 0, len);
    return hex(md.digest());
  }

  static byte[] deflate(byte[] input, int level) throws Exception {
    ByteArrayOutputStream sink = new ByteArrayOutputStream(input.length);
    BlockCompressedOutputStream out = new BlockCompressedOutputStream(sink, (File) null, level);
    out.write(input);
    out.close();
    return sink.toByteArray();
  }

  static byte[] inflate(byte[] framed) throws Exception {
    BlockCompressedInputStream in = new BlockCompressedInputStream(new ByteArrayInputStream(framed));
    ByteArrayOutputStream sink = new ByteArrayOutputStream();
    byte[] buf = new byte[65536];
    int n;
    while ((n = in.read(buf)) > 0) sink.write(buf, 0, n);
    in.close();
    return sink.toByteArray();
  }

  public static void main(String[] args) throws Exception {
    int megabytes = args.length > 0 ? Integer.parseInt(args[0]) : 64;
    int reps = args.length > 1 ? Integer.parseInt(args[1]) : 3;
    int size = megabytes * 1024 * 1024;

    byte[][] payloads = {text(size), lcg(size, 12345L, 58)};
    String[] names = {"text", "lcg"};

    for (int p = 0; p < payloads.length; p++) {
      byte[] input = payloads[p];
      System.out.printf("payload_%s_md5=%s%n", names[p], md5(input, input.length));
      for (int level : new int[] {1, 5, 6, 9}) {
        byte[] framed = null;
        // The first run is thrown away: it pays the JIT, and the port has no equivalent cost.
        deflate(input, level);
        for (int r = 0; r < reps; r++) {
          long t0 = System.nanoTime();
          framed = deflate(input, level);
          long ns = System.nanoTime() - t0;
          System.out.printf(
              "java_deflate_%s_level%d_run%d_mbps=%.2f%n",
              names[p], level, r, (megabytes * 1e9) / ns);
        }
        System.out.printf(
            "java_deflate_%s_level%d_bytes=%d md5=%s%n",
            names[p], level, framed.length, md5(framed, framed.length));

        inflate(framed);
        for (int r = 0; r < reps; r++) {
          long t0 = System.nanoTime();
          byte[] back = inflate(framed);
          long ns = System.nanoTime() - t0;
          if (back.length != input.length) throw new IllegalStateException("round trip lost bytes");
          System.out.printf(
              "java_inflate_%s_level%d_run%d_mbps=%.2f%n",
              names[p], level, r, (megabytes * 1e9) / ns);
        }
      }
    }
  }
}
