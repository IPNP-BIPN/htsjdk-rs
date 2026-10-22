import htsjdk.samtools.CRAMFileWriter;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.cram.compression.CompressionUtils;
import htsjdk.samtools.cram.compression.fqzcomp.FQZCompDecode;
import htsjdk.samtools.cram.compression.fqzcomp.FQZCompEncode;
import htsjdk.samtools.cram.compression.nametokenisation.NameTokenisationDecode;
import htsjdk.samtools.cram.compression.nametokenisation.NameTokenisationEncode;
import htsjdk.samtools.cram.compression.range.RangeDecode;
import htsjdk.samtools.cram.compression.range.RangeEncode;
import htsjdk.samtools.cram.compression.range.RangeParams;
import htsjdk.samtools.cram.compression.rans.RANS4x8Decode;
import htsjdk.samtools.cram.compression.rans.RANS4x8Encode;
import htsjdk.samtools.cram.compression.rans.RANS4x8Params;
import htsjdk.samtools.cram.compression.rans.RANSNx16Decode;
import htsjdk.samtools.cram.compression.rans.RANSNx16Encode;
import htsjdk.samtools.cram.compression.rans.RANSNx16Params;
import htsjdk.samtools.cram.compression.rans.RANSParams;
import htsjdk.samtools.cram.ref.ReferenceSource;
import java.io.BufferedOutputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.Locale;

/**
 * The reference side of `cram-bench`: htsjdk 5.0.0's CRAM 3.1 writer and reader on the same reads
 * and reference, and its CRAM 3.1 codecs on the same blocks. Run with `-Dsamjdk.use_libdeflate=true`,
 * 5.0.0's default, which the oracle image turns off for its goldens: the port's CRAM writer deflates
 * its GZIP blocks the way 5.0.0 does by default.
 *
 * Usage: CramBench dir reps
 */
public class CramBench {
  static File dir;
  static int reps;

  static <T> T time(String name, Callable<T> op) throws Exception {
    T out = op.call();
    for (int r = 0; r < reps; r++) {
      long t0 = System.nanoTime();
      op.call();
      System.out.printf(Locale.ROOT, "java_%s_run%d_secs=%.6f%n", name, r, (System.nanoTime() - t0) / 1e9);
    }
    return out;
  }

  static File port(String op, String suffix) {
    return new File(dir, "port_" + op + suffix);
  }

  static File java(String op, String suffix) {
    return new File(dir, "java_" + op + suffix);
  }

  static SAMFileHeader header;

  static List<SAMRecord> read(File f, File reference) throws Exception {
    List<SAMRecord> out = new ArrayList<>();
    SamReaderFactory factory = SamReaderFactory.makeDefault()
        .validationStringency(ValidationStringency.SILENT);
    if (reference != null) factory = factory.referenceSequence(reference);
    try (SamReader in = factory.open(f)) {
      header = in.getFileHeader();
      for (SAMRecord rec : in) out.add(rec);
    }
    return out;
  }

  static void writeSam(File f, SAMFileHeader h, List<SAMRecord> records) {
    try (SAMFileWriter w = new SAMFileWriterFactory().makeSAMWriter(h, true, f)) {
      for (SAMRecord rec : records) w.addAlignment(rec);
    }
  }

  static void cram() throws Exception {
    File fasta = new File(dir, "cram_ref.fasta");
    List<SAMRecord> records = read(new File(dir, "cram_input.bam"), null);
    SAMFileHeader h = header;
    File out = java("cram_write", ".cram");
    time("cram_write", () -> {
      try (OutputStream os = new BufferedOutputStream(new FileOutputStream(out))) {
        // The name is the file definition's 20-byte id; the port writes the same one.
        CRAMFileWriter w = new CRAMFileWriter(os, new ReferenceSource(fasta), h, "cram_write.cram");
        for (SAMRecord rec : records) w.addAlignment(rec);
        w.close();
      }
      return null;
    });

    List<SAMRecord> decoded = time("cram_read", () -> read(port("cram_write", ".cram"), fasta));
    writeSam(java("cram_read", ".sam"), header, decoded);
  }

  /** One codec direction: the block's index and the block in, one block out. */
  interface Codec {
    byte[] apply(int index, byte[] in) throws Exception;
  }

  static byte[] bytes(ByteBuffer b) {
    byte[] out = new byte[b.remaining()];
    b.duplicate().get(out);
    return out;
  }

  /** The `u32` length + bytes framing `cram-bench` writes. */
  static List<byte[]> unframe(File f) throws Exception {
    ByteBuffer in = ByteBuffer.wrap(Files.readAllBytes(f.toPath())).order(ByteOrder.LITTLE_ENDIAN);
    List<byte[]> out = new ArrayList<>();
    while (in.hasRemaining()) {
      byte[] part = new byte[in.getInt()];
      in.get(part);
      out.add(part);
    }
    return out;
  }

  static void frame(File f, List<byte[]> parts) throws Exception {
    ByteArrayOutputStream out = new ByteArrayOutputStream();
    ByteBuffer len = ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN);
    for (byte[] p : parts) {
      out.write(len.putInt(0, p.length).array());
      out.write(p);
    }
    Files.write(f.toPath(), out.toByteArray());
  }

  static List<byte[]> each(List<byte[]> in, Codec codec) throws Exception {
    List<byte[]> out = new ArrayList<>();
    for (int i = 0; i < in.size(); i++) out.add(codec.apply(i, in.get(i)));
    return out;
  }

  /**
   * One codec both ways. The blocks to compress are the port's decompressed ones, which the port
   * asserted equal its input; the ones to decompress are the port's compressed ones.
   */
  static void codec(String name, Codec compress, Codec uncompress) throws Exception {
    List<byte[]> packed = unframe(port(name + "_compress", ".bin"));
    List<byte[]> plain = unframe(port(name + "_uncompress", ".bin"));
    frame(java(name + "_compress", ".bin"), time(name + "_compress", () -> each(plain, compress)));
    frame(java(name + "_uncompress", ".bin"), time(name + "_uncompress", () -> each(packed, uncompress)));
  }

  /** The per-block `i32` arrays `cram-bench` writes beside the fqzcomp blocks. */
  static List<int[]> ints(String name) throws Exception {
    List<int[]> out = new ArrayList<>();
    for (byte[] part : unframe(new File(dir, "fqzcomp_" + name + ".bin"))) {
      int[] values = new int[part.length / 4];
      ByteBuffer.wrap(part).order(ByteOrder.LITTLE_ENDIAN).asIntBuffer().get(values);
      out.add(values);
    }
    return out;
  }

  public static void main(String[] args) throws Exception {
    dir = new File(args[0]);
    reps = args.length > 1 ? Integer.parseInt(args[1]) : 3;
    cram();

    codec("rans4x8_o0",
        (i, x) -> new RANS4x8Encode().compress(x, new RANS4x8Params(RANSParams.ORDER.ZERO)),
        (i, x) -> new RANS4x8Decode().uncompress(x));
    codec("rans4x8_o1",
        (i, x) -> new RANS4x8Encode().compress(x, new RANS4x8Params(RANSParams.ORDER.ONE)),
        (i, x) -> new RANS4x8Decode().uncompress(x));
    codec("ransnx16_o0",
        (i, x) -> new RANSNx16Encode().compress(x, new RANSNx16Params(0)),
        (i, x) -> new RANSNx16Decode().uncompress(x));
    codec("ransnx16_o1",
        (i, x) -> new RANSNx16Encode().compress(x, new RANSNx16Params(1)),
        (i, x) -> new RANSNx16Decode().uncompress(x));
    codec("range_o0",
        (i, x) -> bytes(new RangeEncode().compress(CompressionUtils.wrap(x), new RangeParams(0))),
        (i, x) -> bytes(new RangeDecode().uncompress(CompressionUtils.wrap(x))));
    codec("range_o1",
        (i, x) -> bytes(new RangeEncode().compress(CompressionUtils.wrap(x), new RangeParams(1))),
        (i, x) -> bytes(new RangeDecode().uncompress(CompressionUtils.wrap(x))));
    codec("name_tok",
        (i, x) -> bytes(new NameTokenisationEncode().compress(CompressionUtils.wrap(x), false, (byte) 0)),
        (i, x) -> new NameTokenisationDecode().uncompress(CompressionUtils.wrap(x), (byte) 0));
    List<int[]> lengths = ints("lengths");
    List<int[]> flags = ints("flags");
    codec("fqzcomp",
        (i, x) -> bytes(new FQZCompEncode().compress(CompressionUtils.wrap(x), lengths.get(i), flags.get(i))),
        (i, x) -> bytes(FQZCompDecode.uncompress(CompressionUtils.wrap(x))));
  }
}
