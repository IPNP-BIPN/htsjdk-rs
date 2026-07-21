import java.util.zip.Deflater;
import java.security.MessageDigest;
public class Z {
  public static void main(String[] args) throws Exception {
    byte[] in = new byte[65536];
    long s = 12345L;
    for (int i = 0; i < in.length; i++) { s = s*6364136223846793005L + 1442695040888963407L; in[i] = (byte)(s >>> 58); }
    MessageDigest md0 = MessageDigest.getInstance("MD5"); md0.update(in);
    System.out.printf("input  len=%d md5=%s%n", in.length, hex(md0.digest()));
    for (int lvl : new int[]{1,5,6,9}) {
      Deflater d = new Deflater(lvl, true);   // nowrap=true, as BGZF uses
      d.setInput(in); d.finish();
      byte[] out = new byte[in.length*2];
      int n = d.deflate(out);
      d.end();
      MessageDigest md = MessageDigest.getInstance("MD5"); md.update(out, 0, n);
      System.out.printf("level=%d outlen=%d md5=%s%n", lvl, n, hex(md.digest()));
    }
  }
  static String hex(byte[] b){ StringBuilder sb=new StringBuilder(); for(byte x:b) sb.append(String.format("%02x",x)); return sb.toString(); }
}
