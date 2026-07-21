import java.util.zip.Deflater;
import java.security.MessageDigest;
public class Z2 {
  static byte[] lcg(int n, long seed, int shift){ byte[] b=new byte[n]; long s=seed;
    for(int i=0;i<n;i++){ s=s*6364136223846793005L+1442695040888963407L; b[i]=(byte)(s>>>shift);} return b; }
  static byte[] zeros(int n){ return new byte[n]; }
  static byte[] runs(int n){ byte[] b=new byte[n]; for(int i=0;i<n;i++) b[i]=(byte)((i/64)%7); return b; }
  static byte[] text(int n){ byte[] p="ACGTNacgtn\tSAMrecord\n".getBytes(); byte[] b=new byte[n];
    for(int i=0;i<n;i++) b[i]=p[i%p.length]; return b; }
  public static void main(String[] a) throws Exception {
    Object[][] payloads = {
      {"lcg64k",  lcg(65536,12345L,58)},
      {"rand64k", lcg(65536,999L,56)},
      {"zeros64k",zeros(65536)},
      {"runs64k", runs(65536)},
      {"text64k", text(65536)},
      {"empty",   new byte[0]},
      {"single",  new byte[]{0x42}},
    };
    for (Object[] p : payloads) {
      String name=(String)p[0]; byte[] in=(byte[])p[1];
      MessageDigest mi=MessageDigest.getInstance("MD5"); mi.update(in);
      System.out.printf("PAYLOAD %s len=%d md5=%s%n", name, in.length, hex(mi.digest()));
      for (int lvl=0; lvl<=9; lvl++) {
        Deflater d=new Deflater(lvl,true); d.setInput(in); d.finish();
        byte[] out=new byte[Math.max(in.length*2,256)]; int tot=0;
        while(!d.finished()){ int k=d.deflate(out,tot,out.length-tot); if(k==0) break; tot+=k; }
        d.end();
        MessageDigest md=MessageDigest.getInstance("MD5"); md.update(out,0,tot);
        System.out.printf("  (\"%s\", %d, %d, \"%s\"),%n", name, lvl, tot, hex(md.digest()));
      }
    }
  }
  static String hex(byte[] b){ StringBuilder s=new StringBuilder(); for(byte x:b) s.append(String.format("%02x",x)); return s.toString(); }
}
