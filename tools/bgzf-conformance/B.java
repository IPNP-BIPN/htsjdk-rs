import htsjdk.samtools.util.BlockCompressedOutputStream;
import java.io.*; import java.security.MessageDigest;
public class B {
  static byte[] lcg(int n,long seed,int sh){ byte[] b=new byte[n]; long s=seed;
    for(int i=0;i<n;i++){ s=s*6364136223846793005L+1442695040888963407L; b[i]=(byte)(s>>>sh);} return b; }
  static byte[] runs(int n){ byte[] b=new byte[n]; for(int i=0;i<n;i++) b[i]=(byte)((i/64)%7); return b; }
  static byte[] text(int n){ byte[] p="ACGTNacgtn\tSAMrecord\n".getBytes(); byte[] b=new byte[n];
    for(int i=0;i<n;i++) b[i]=p[i%p.length]; return b; }
  public static void main(String[] a) throws Exception {
    Object[][] cases = {
      {"empty",   new byte[0]},
      {"tiny",    text(10)},
      {"exact1",  runs(65498)},          // exactly one full block
      {"over1",   runs(65499)},          // one full block + 1 byte
      {"multi",   lcg(200000,12345L,58)},// several blocks
      {"incompr", lcg(200000,999L,56)},  // incompressible, exercises fallback
      {"big",     text(500000)},
    };
    for (int lvl : new int[]{0,1,5,6,9}) {
      for (Object[] c : cases) {
        String name=(String)c[0]; byte[] in=(byte[])c[1];
        ByteArrayOutputStream bos=new ByteArrayOutputStream();
        BlockCompressedOutputStream os=new BlockCompressedOutputStream(bos,(File)null,lvl);
        os.write(in); os.close();
        byte[] out=bos.toByteArray();
        MessageDigest md=MessageDigest.getInstance("MD5"); md.update(out);
        System.out.printf("  (\"%s\", %d, %d, \"%s\"),%n", name, lvl, out.length, hex(md.digest()));
      }
    }
  }
  static String hex(byte[] b){ StringBuilder s=new StringBuilder(); for(byte x:b) s.append(String.format("%02x",x)); return s.toString(); }
}
