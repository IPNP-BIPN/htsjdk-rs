import htsjdk.samtools.util.BlockCompressedOutputStream;
import java.io.*;

/** Writes BGZF fixture files with htsjdk itself, for the read-path conformance test. */
public class MakeFixtures {
  static byte[] lcg(int n,long seed,int sh){ byte[] b=new byte[n]; long s=seed;
    for(int i=0;i<n;i++){ s=s*6364136223846793005L+1442695040888963407L; b[i]=(byte)(s>>>sh);} return b; }
  static byte[] runs(int n){ byte[] b=new byte[n]; for(int i=0;i<n;i++) b[i]=(byte)((i/64)%7); return b; }
  static byte[] text(int n){ byte[] p="ACGTNacgtn\tSAMrecord\n".getBytes(); byte[] b=new byte[n];
    for(int i=0;i<n;i++) b[i]=p[i%p.length]; return b; }
  public static void main(String[] a) throws Exception {
    String outDir = a.length > 0 ? a[0] : "/out";
    new File(outDir).mkdirs();
    Object[][] cases = {
      {"empty", new byte[0]}, {"tiny", text(10)}, {"exact1", runs(65498)},
      {"over1", runs(65499)}, {"multi", lcg(200000,12345L,58)},
      {"incompr", lcg(200000,999L,56)}, {"big", text(500000)},
    };
    for (int lvl : new int[]{0,1,5,6,9}) {
      for (Object[] c : cases) {
        String name=(String)c[0]; byte[] in=(byte[])c[1];
        File f = new File(outDir, name + "_l" + lvl + ".bgzf");
        BlockCompressedOutputStream os = new BlockCompressedOutputStream(new FileOutputStream(f), (File)null, lvl);
        os.write(in); os.close();
        System.out.println(f.getName() + " " + f.length());
      }
    }
  }
}
