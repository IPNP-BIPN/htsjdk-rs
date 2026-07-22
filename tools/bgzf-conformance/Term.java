/*
 * Oracle for BlockCompressedInputStream.checkTermination. For a few payloads written by htsjdk's
 * BlockCompressedOutputStream at the default level (5, which BgzfWriter matches byte-for-byte), and
 * for a few truncation variants of each, emits the FileTermination htsjdk reports, as Rust tuples:
 *
 *   ("multi", "full", "HAS_TERMINATOR_BLOCK"),
 *
 * The Rust conformance test rebuilds the identical bytes with BgzfWriter, applies the same
 * truncation, and asserts check_termination returns the same variant.
 *
 *   java -cp htsjdk.jar:. Term
 */
import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import java.io.*;
import java.nio.file.*;

public class Term {
  static byte[] lcg(int n,long seed,int sh){ byte[] b=new byte[n]; long s=seed;
    for(int i=0;i<n;i++){ s=s*6364136223846793005L+1442695040888963407L; b[i]=(byte)(s>>>sh);} return b; }
  static byte[] runs(int n){ byte[] b=new byte[n]; for(int i=0;i<n;i++) b[i]=(byte)((i/64)%7); return b; }
  static byte[] text(int n){ byte[] p="ACGTNacgtn\tSAMrecord\n".getBytes(); byte[] b=new byte[n];
    for(int i=0;i<n;i++) b[i]=p[i%p.length]; return b; }

  static void report(String name, byte[] full, String variant, byte[] bytes) throws Exception {
    File f = File.createTempFile("term", ".bgzf"); f.deleteOnExit();
    try (FileOutputStream o = new FileOutputStream(f)) { o.write(bytes); }
    BlockCompressedInputStream.FileTermination t = BlockCompressedInputStream.checkTermination(f);
    System.out.printf("  (\"%s\", \"%s\", \"%s\"),%n", name, variant, t.name());
  }

  static byte[] trim(byte[] b, int n) { return java.util.Arrays.copyOf(b, Math.max(0, b.length - n)); }

  public static void main(String[] a) throws Exception {
    Object[][] cases = {
      {"tiny",   text(10)},
      {"exact1", runs(65498)},
      {"multi",  lcg(200000,12345L,58)},
      {"big",    text(500000)},
    };
    for (Object[] c : cases) {
      String name=(String)c[0]; byte[] in=(byte[])c[1];
      ByteArrayOutputStream bos=new ByteArrayOutputStream();
      BlockCompressedOutputStream os=new BlockCompressedOutputStream(bos,(File)null,5);
      os.write(in); os.close();
      byte[] full=bos.toByteArray();
      report(name, full, "full", full);
      report(name, full, "no_terminator", trim(full, 28));
      report(name, full, "truncated", trim(full, 33));
    }
  }
}
