import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import java.io.File;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;
import java.util.Locale;

/**
 * The reference side of the BAM and VCF rows: htsjdk reading and writing the very files the port's
 * `record-bench` and `vcf-bench` generated, in the same container.
 *
 * It writes its output next to the port's, and the runner compares the two files byte for byte
 * before it prints a single duration. Durations are printed in seconds; the runner divides the
 * byte and record counts the port printed, so both columns are over the same denominator.
 *
 * Usage: CodecBench dir reps
 */
public class CodecBench {
  static int reps;

  interface Op {
    void run() throws Exception;
  }

  /** One discarded run for the JIT, as `BgzfBench` does, then `reps` timed ones. */
  static void time(String name, Op op) throws Exception {
    op.run();
    for (int r = 0; r < reps; r++) {
      long t0 = System.nanoTime();
      op.run();
      System.out.printf(Locale.ROOT, "java_%s_run%d_secs=%.6f%n", name, r, (System.nanoTime() - t0) / 1e9);
    }
  }

  /**
   * Eager decoding is the port's decoder: every field of every record, up front. SILENT because the
   * port's reader validates nothing, and STRICT would time `isValid` on top of the decode.
   */
  static SamReader openBam(File f) {
    return SamReaderFactory.makeDefault()
        .enable(SamReaderFactory.Option.EAGERLY_DECODE)
        .validationStringency(ValidationStringency.SILENT)
        .open(f);
  }

  static void bam(File dir) throws Exception {
    File corpus = new File(dir, "corpus.bam");
    time(
        "bam_decode",
        () -> {
          try (SamReader in = openBam(corpus)) {
            long n = 0;
            for (SAMRecord rec : in) n++;
            if (n == 0) throw new IllegalStateException("no records");
          }
        });

    List<SAMRecord> records = new ArrayList<>();
    SAMFileHeader header;
    try (SamReader in = openBam(corpus)) {
      header = in.getFileHeader();
      for (SAMRecord rec : in) records.add(rec);
    }
    System.out.printf(Locale.ROOT, "java_bam_records=%d%n", records.size());

    // htsjdk's own defaults: the JDK deflater at level 5, no index, no async writer.
    File out = new File(dir, "java_out.bam");
    time(
        "bam_encode",
        () -> {
          try (SAMFileWriter w = new SAMFileWriterFactory().makeBAMWriter(header, true, out)) {
            for (SAMRecord rec : records) w.addAlignment(rec);
          }
        });
  }

  static List<VariantContext> readVcf(File f, boolean decode) {
    List<VariantContext> out = new ArrayList<>();
    try (VCFFileReader in = new VCFFileReader(f, false)) {
      for (VariantContext vc : in) {
        // `getGenotype` reaches `LazyGenotypesContext.decode()`, the port's `genotypes.first()`.
        if (decode) vc.getGenotype(0);
        out.add(vc);
      }
    }
    return out;
  }

  static void writeVcf(File f, VCFHeader header, List<VariantContext> records) {
    try (VariantContextWriter w =
        new VariantContextWriterBuilder()
            .setOutputFile(f)
            .setOptions(EnumSet.noneOf(Options.class))
            .build()) {
      w.writeHeader(header);
      for (VariantContext vc : records) w.add(vc);
    }
  }

  static void vcf(File dir) throws Exception {
    File corpus = new File(dir, "corpus.vcf");
    VCFHeader header;
    try (VCFFileReader in = new VCFFileReader(corpus, false)) {
      header = in.getFileHeader();
    }
    time("vcf_read_lazy", () -> readVcf(corpus, false));
    time("vcf_read_decoded", () -> readVcf(corpus, true));

    System.out.printf(Locale.ROOT, "java_vcf_records=%d%n", readVcf(corpus, false).size());
    for (String mode : new String[] {"copy", "touched"}) {
      File out = new File(dir, "java_" + mode + ".vcf");
      boolean decode = mode.equals("touched");
      // Run -1 is the JIT's. A fresh read before every run, untimed: decoding is what makes the
      // encoder leave the verbatim branch, so a write must not see records an earlier one touched.
      for (int r = -1; r < reps; r++) {
        List<VariantContext> records = readVcf(corpus, decode);
        long t0 = System.nanoTime();
        writeVcf(out, header, records);
        if (r >= 0) {
          System.out.printf(
              "java_vcf_write_%s_run%d_secs=%.6f%n", mode, r, (System.nanoTime() - t0) / 1e9);
        }
      }
    }
  }

  public static void main(String[] args) throws Exception {
    File dir = new File(args[0]);
    reps = args.length > 1 ? Integer.parseInt(args[1]) : 3;
    bam(dir);
    vcf(dir);
  }
}
