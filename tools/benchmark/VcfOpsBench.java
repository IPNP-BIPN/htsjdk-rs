import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;
import java.util.concurrent.Callable;

/**
 * The reference side of `vcf-ops-bench`: htsjdk on the same inputs, writing `java_<op><suffix>`
 * beside each `port_<op><suffix>`.
 *
 * Usage: VcfOpsBench dir reps
 */
public class VcfOpsBench {
  static File dir;
  static int reps;

  interface Op<T> {
    T run(List<VariantContext> records) throws Exception;
  }

  /**
   * One untimed run, the check and the JIT's warm-up, then `reps` timed ones. Every run gets a
   * fresh read of the corpus outside the timing: decoding a record changes which branch the
   * encoders take, so a run must not see records an earlier one touched.
   */
  static <T> T time(String name, File input, Op<T> op) throws Exception {
    T out = op.run(input == null ? null : read(input));
    for (int r = 0; r < reps; r++) {
      List<VariantContext> records = input == null ? null : read(input);
      long t0 = System.nanoTime();
      op.run(records);
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

  static VCFHeader header;

  static List<VariantContext> read(File f) {
    List<VariantContext> out = new ArrayList<>();
    try (VCFFileReader in = new VCFFileReader(f, false)) {
      header = in.getFileHeader();
      for (VariantContext vc : in) out.add(vc);
    }
    return out;
  }

  static void write(File f, VCFHeader h, List<VariantContext> records, boolean index) {
    VariantContextWriterBuilder builder = new VariantContextWriterBuilder().setOutputFile(f);
    if (index) {
      builder
          .setReferenceDictionary(h.getSequenceDictionary())
          .setOptions(EnumSet.of(Options.INDEX_ON_THE_FLY));
    } else {
      builder.setOptions(EnumSet.noneOf(Options.class));
    }
    try (VariantContextWriter w = builder.build()) {
      w.writeHeader(h);
      for (VariantContext vc : records) w.add(vc);
    }
  }

  public static void main(String[] args) throws Exception {
    dir = new File(args[0]);
    reps = args.length > 1 ? Integer.parseInt(args[1]) : 3;
    File corpus = new File(dir, "corpus.vcf");

    time("bcf_write", corpus, records -> {
      write(java("bcf_write", ".bcf"), header, records, false);
      return null;
    });
    List<VariantContext> bcf = time("bcf_read", null, x -> read(port("bcf_write", ".bcf")));
    write(java("bcf_read", ".vcf"), header, bcf, false);

    time("vcf_idx_write", corpus, records -> {
      write(java("vcf_idx_write", ".vcf"), header, records, true);
      return null;
    });
    time("vcf_gz_write", corpus, records -> {
      write(java("vcf_gz_write", ".vcf.gz"), header, records, true);
      return null;
    });
    List<VariantContext> gz = time("vcf_gz_read", null, x -> read(port("vcf_gz_write", ".vcf.gz")));
    write(java("vcf_gz_read", ".vcf"), header, gz, false);

    File bed = new File(dir, "corpus.bed");
    List<BEDFeature> features = time("bed_decode", null, x -> {
      BEDCodec codec = new BEDCodec();
      List<BEDFeature> out = new ArrayList<>();
      try (BufferedReader in = new BufferedReader(new FileReader(bed))) {
        for (String line = in.readLine(); line != null; line = in.readLine()) {
          BEDFeature f = codec.decode(line);
          if (f != null) out.add(f);
        }
      }
      return out;
    });
    StringBuilder text = new StringBuilder();
    for (BEDFeature f : features) {
      text.append(f.getContig()).append('\t').append(f.getStart()).append('\t').append(f.getEnd())
          .append('\t').append(f.getName()).append('\t').append(Float.floatToIntBits(f.getScore()))
          .append('\t').append(f.getStrand().toString()).append('\n');
    }
    Files.write(java("bed_decode", ".txt").toPath(), text.toString().getBytes("ISO-8859-1"));
  }
}
