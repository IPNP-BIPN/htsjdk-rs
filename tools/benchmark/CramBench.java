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
import htsjdk.samtools.cram.compression.nametokenisation.NameTokenisationDecode;
import htsjdk.samtools.cram.compression.nametokenisation.NameTokenisationEncode;
import htsjdk.samtools.cram.compression.range.RangeDecode;
import htsjdk.samtools.cram.compression.range.RangeEncode;
import htsjdk.samtools.cram.compression.range.RangeParams;
import htsjdk.samtools.cram.compression.rans.RANSParams;
import htsjdk.samtools.cram.compression.rans.rans4x8.RANS4x8Decode;
import htsjdk.samtools.cram.compression.rans.rans4x8.RANS4x8Encode;
import htsjdk.samtools.cram.compression.rans.rans4x8.RANS4x8Params;
import htsjdk.samtools.cram.compression.rans.ransnx16.RANSNx16Decode;
import htsjdk.samtools.cram.compression.rans.ransnx16.RANSNx16Encode;
import htsjdk.samtools.cram.compression.rans.ransnx16.RANSNx16Params;
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

/**
 * The reference side of `cram-bench`: htsjdk 4.2.0's CRAM writer and reader on the same reads and
 * reference, and its CRAM 3.1 codecs on the same blocks.
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
      System.out.printf("java_%s_run%d_secs=%.6f%n", name, r, (System.nanoTime() - t0) / 1e9);
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
        CRAMFileWriter w = new CRAMFileWriter(os, new ReferenceSource(fasta), h, out.getName());
        for (SAMRecord rec : records) w.addAlignment(rec);
        w.close();
      }
      return null;
    });
    writeSam(java("cram_write", ".sam"), h, read(out, fasta));

    // htsjdk 4.2.0 cannot read the port's CRAM 3.1, so both sides read this CRAM 3.0 file.
    List<SAMRecord> decoded = time("cram_read", () -> read(out, fasta));
    writeSam(java("cram_read", ".sam"), header, decoded);
  }

  interface Codec {
    byte[] apply(byte[] in) throws Exception;
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
    for (byte[] x : in) out.add(codec.apply(x));
    return out;
  }

  /**
   * One codec both ways. The blocks to compress are the port's decompressed ones, which the port
   * asserted equal its input; the ones to decompress are the port's compressed ones.
   */
  static void codec(String name, Codec compress, Codec uncompress) throws Exception {
    List<byte[]> packed = unframe(port(name + "_compress", ".bin"));
    List<byte[]> plain = unframe(port(name + "_uncompress", ".bin"));
    if (compress != null) {
      frame(java(name + "_compress", ".bin"), time(name + "_compress", () -> each(plain, compress)));
      // Where the port follows a later encoder, the check is this htsjdk decoding its streams.
      if (port(name + "_compress", ".dec").exists()) {
        frame(java(name + "_compress", ".dec"), each(packed, uncompress));
      }
    }
    frame(java(name + "_uncompress", ".bin"), time(name + "_uncompress", () -> each(packed, uncompress)));
  }

  public static void main(String[] args) throws Exception {
    dir = new File(args[0]);
    reps = args.length > 1 ? Integer.parseInt(args[1]) : 3;
    cram();

    codec("rans4x8_o0",
        x -> bytes(new RANS4x8Encode().compress(CompressionUtils.wrap(x),
            new RANS4x8Params(RANSParams.ORDER.ZERO))),
        x -> bytes(new RANS4x8Decode().uncompress(CompressionUtils.wrap(x))));
    codec("rans4x8_o1",
        x -> bytes(new RANS4x8Encode().compress(CompressionUtils.wrap(x),
            new RANS4x8Params(RANSParams.ORDER.ONE))),
        x -> bytes(new RANS4x8Decode().uncompress(CompressionUtils.wrap(x))));
    codec("ransnx16_o0",
        x -> bytes(new RANSNx16Encode().compress(CompressionUtils.wrap(x), new RANSNx16Params(0))),
        x -> bytes(new RANSNx16Decode().uncompress(CompressionUtils.wrap(x))));
    codec("ransnx16_o1",
        x -> bytes(new RANSNx16Encode().compress(CompressionUtils.wrap(x), new RANSNx16Params(1))),
        x -> bytes(new RANSNx16Decode().uncompress(CompressionUtils.wrap(x))));
    codec("range_o0",
        x -> bytes(new RangeEncode().compress(CompressionUtils.wrap(x), new RangeParams(0))),
        x -> bytes(new RangeDecode().uncompress(CompressionUtils.wrap(x))));
    codec("range_o1",
        x -> bytes(new RangeEncode().compress(CompressionUtils.wrap(x), new RangeParams(1))),
        x -> bytes(new RangeDecode().uncompress(CompressionUtils.wrap(x))));
    codec("name_tok",
        x -> bytes(new NameTokenisationEncode().compress(CompressionUtils.wrap(x), false, (byte) 0)),
        x -> new NameTokenisationDecode().uncompress(CompressionUtils.wrap(x), (byte) 0));
    codec("fqzcomp", null, x -> bytes(FQZCompDecode.uncompress(CompressionUtils.wrap(x))));
  }
}
