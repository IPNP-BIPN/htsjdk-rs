import htsjdk.samtools.BAMIndexer;
import htsjdk.samtools.BamFileIoUtils;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.samtools.SAMRecordQueryNameComparator;
import htsjdk.samtools.SBIIndexWriter;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.fastq.BasicFastqWriter;
import htsjdk.samtools.fastq.FastqEncoder;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.reference.FastaReferenceWriter;
import htsjdk.samtools.reference.FastaReferenceWriterBuilder;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.samtools.util.SequenceUtil;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.Locale;

/**
 * The reference side of `bam-ops-bench`: htsjdk on the same inputs, writing `java_<op><suffix>`
 * beside each `port_<op><suffix>` for the runner to compare before it reads a duration.
 *
 * Usage: BamOpsBench dir reps
 */
public class BamOpsBench {
  static File dir;
  static int reps;

  /** One untimed run, whose answer is the check and the JIT's warm-up, then `reps` timed ones. */
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

  static SamReader open(File f) {
    return SamReaderFactory.makeDefault()
        .enable(SamReaderFactory.Option.EAGERLY_DECODE)
        .validationStringency(ValidationStringency.SILENT)
        .open(f);
  }

  static SAMFileHeader header;

  static List<SAMRecord> read(File f) throws Exception {
    List<SAMRecord> out = new ArrayList<>();
    try (SamReader in = open(f)) {
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

  static void writeBam(File f, SAMFileHeader h, List<SAMRecord> records) {
    try (SAMFileWriter w = new SAMFileWriterFactory().makeBAMWriter(h, true, f)) {
      for (SAMRecord rec : records) w.addAlignment(rec);
    }
  }

  static void write(File f, String text) throws Exception {
    Files.write(f.toPath(), text.getBytes("ISO-8859-1"));
  }

  static String hex(byte[] b) {
    StringBuilder sb = new StringBuilder();
    for (byte x : b) sb.append(String.format("%02x", x));
    return sb.toString();
  }

  static void sam(List<SAMRecord> records, SAMFileHeader h) throws Exception {
    time("sam_write", () -> {
      writeSam(java("sam_write", ".sam"), h, records);
      return null;
    });
    List<SAMRecord> read = time("sam_read", () -> read(port("sam_write", ".sam")));
    writeSam(java("sam_read", ".sam"), header, read);
  }

  /** Queryname order from the corpus's order, then coordinate order from queryname order. */
  static void sorting(List<SAMRecord> records, SAMFileHeader h) throws Exception {
    List<SAMRecord> input = records;
    for (String name : new String[] {"queryname", "coordinate"}) {
      List<SAMRecord> from = input;
      List<SAMRecord> sorted = time("sort_" + name, () -> {
        List<SAMRecord> copy = new ArrayList<>(from);
        copy.sort(name.equals("coordinate")
            ? new SAMRecordCoordinateComparator()
            : new SAMRecordQueryNameComparator());
        return copy;
      });
      SAMFileHeader sortedHeader = h.clone();
      sortedHeader.setSortOrder(SAMFileHeader.SortOrder.valueOf(name));
      writeBam(java("sort_" + name, ".bam"), sortedHeader, sorted);
      input = sorted;
    }
  }

  static long sbiOffset(long i) {
    return ((i / 40) * 9000) << 16 | ((i % 40) * 1500);
  }

  static void indexes() throws Exception {
    // As picard BuildBamIndex opens it: sources in records so the file pointers are there.
    File sorted = port("sort_coordinate", ".bam");
    time("bai_build", () -> {
      try (SamReader reader = SamReaderFactory.makeDefault()
          .disable(SamReaderFactory.Option.EAGERLY_DECODE)
          .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS)
          .open(sorted)) {
        BAMIndexer.createIndex(reader, java("bai_build", ".bai"));
      }
      return null;
    });

    long count = 10_000_000L;
    time("sbi_write", () -> {
      try (OutputStream out =
          new BufferedOutputStream(new FileOutputStream(java("sbi_write", ".sbi")))) {
        SBIIndexWriter writer = new SBIIndexWriter(out, 4096);
        for (long i = 0; i < count; i++) writer.processRecord(sbiOffset(i));
        writer.finish(sbiOffset(count), 1L << 40);
      }
      return null;
    });
  }

  static void blockCopy() throws Exception {
    Path sorted = port("sort_coordinate", ".bam").toPath();
    read(sorted.toFile());
    SAMFileHeader h = header.clone();
    h.addComment("benchmark");
    time("bam_reheader", () -> {
      BamFileIoUtils.reheaderBamFile(h, sorted, java("bam_reheader", ".bam").toPath(), false, false);
      return null;
    });
    List<File> inputs = List.of(
        new File(dir, "corpus.bam"), port("sort_coordinate", ".bam"), port("sort_queryname", ".bam"));
    time("bam_gather", () -> {
      BamFileIoUtils.gatherWithBlockCopying(inputs, java("bam_gather", ".bam"), false, false);
      return null;
    });
  }

  static void fastq(List<SAMRecord> records) throws Exception {
    List<FastqRecord> fastq = new ArrayList<>();
    for (SAMRecord rec : records) fastq.add(FastqEncoder.asFastqRecord(rec));
    time("fastq_write", () -> {
      try (BasicFastqWriter w = new BasicFastqWriter(java("fastq_write", ".fastq"))) {
        for (FastqRecord r : fastq) w.write(r);
      }
      return null;
    });
    List<FastqRecord> read = time("fastq_read", () -> {
      List<FastqRecord> out = new ArrayList<>();
      try (FastqReader in = new FastqReader(port("fastq_write", ".fastq"))) {
        for (FastqRecord r : in) out.add(r);
      }
      return out;
    });
    try (BasicFastqWriter w = new BasicFastqWriter(java("fastq_read", ".fastq"))) {
      for (FastqRecord r : read) w.write(r);
    }
  }

  static List<ReferenceSequence> readFasta(File f) {
    List<ReferenceSequence> out = new ArrayList<>();
    try (FastaSequenceFile in = new FastaSequenceFile(f, true)) {
      for (ReferenceSequence s = in.nextSequence(); s != null; s = in.nextSequence()) out.add(s);
    }
    return out;
  }

  static List<ReferenceSequence> fasta() throws Exception {
    File portFasta = port("fasta_write", ".fasta");
    List<ReferenceSequence> genome = readFasta(portFasta);
    time("fasta_write", () -> {
      try (FastaReferenceWriter w = new FastaReferenceWriterBuilder()
          .setFastaOutput(Files.newOutputStream(java("fasta_write", ".fasta").toPath()))
          .setIndexOutput(Files.newOutputStream(java("fasta_write", ".fasta.fai").toPath()))
          .setDictOutput(Files.newOutputStream(java("fasta_write", ".dict").toPath()))
          .setBasesPerLine(60)
          .setEmitMd5(true)
          .build()) {
        for (ReferenceSequence s : genome) w.startSequence(s.getName()).appendBases(s.getBases());
      }
      return null;
    });

    List<ReferenceSequence> read = time("fasta_read", () -> readFasta(portFasta));
    StringBuilder text = new StringBuilder();
    for (ReferenceSequence s : read) {
      text.append(s.getName()).append('\t').append(s.length()).append('\t')
          .append(hex(MessageDigest.getInstance("MD5").digest(s.getBases()))).append('\n');
    }
    write(java("fasta_read", ".txt"), text.toString());

    List<String[]> queries = lines(new File(dir, "fasta_queries.txt"));
    byte[] digest = time("fasta_query", () -> {
      MessageDigest md5 = MessageDigest.getInstance("MD5");
      try (IndexedFastaSequenceFile f = new IndexedFastaSequenceFile(portFasta.toPath())) {
        for (String[] q : queries) {
          md5.update(f.getSubsequenceAt(q[0], Long.parseLong(q[1]), Long.parseLong(q[2])).getBases());
        }
      }
      return md5.digest();
    });
    write(java("fasta_query", ".txt"), hex(digest) + "\n");
    return genome;
  }

  static List<String[]> lines(File f) throws Exception {
    List<String[]> out = new ArrayList<>();
    for (String line : Files.readAllLines(f.toPath())) out.add(line.split("\t"));
    return out;
  }

  static void mdNm(List<SAMRecord> records, byte[] chr1) throws Exception {
    time("md_nm", () -> {
      for (SAMRecord rec : records) SequenceUtil.calculateMdAndNmTags(rec, chr1, true, true);
      return null;
    });
    StringBuilder text = new StringBuilder();
    for (SAMRecord rec : records) {
      text.append(rec.getStringAttribute("MD")).append('\t')
          .append(rec.getIntegerAttribute("NM")).append('\n');
    }
    write(java("md_nm", ".txt"), text.toString());
  }

  static void intervals() throws Exception {
    File input = port("interval_list_write", ".interval_list");
    IntervalList list = IntervalList.fromFile(input);
    time("interval_list_write", () -> {
      list.write(java("interval_list_write", ".interval_list"));
      return null;
    });
    IntervalList read = time("interval_list_read", () -> IntervalList.fromFile(input));
    read.write(java("interval_list_read", ".interval_list"));
    time("interval_list_sort", list::sorted).write(java("interval_list_sort", ".interval_list"));
    time("interval_list_uniq", () -> list.uniqued())
        .write(java("interval_list_uniq", ".interval_list"));

    List<Interval> all = list.getIntervals();
    List<Interval> targets = all.subList(0, 50_000);
    List<Interval> queries = all.subList(0, 100_000);
    int[] counts = time("overlap_detector", () -> {
      OverlapDetector<Integer> detector = new OverlapDetector<>(0, 0);
      for (int i = 0; i < targets.size(); i++) detector.addLhs(i, targets.get(i));
      int[] out = new int[queries.size()];
      for (int i = 0; i < out.length; i++) out[i] = detector.getOverlaps(queries.get(i)).size();
      return out;
    });
    StringBuilder text = new StringBuilder();
    for (int c : counts) text.append(c).append('\n');
    write(java("overlap_detector", ".txt"), text.toString());

    List<String[]> lifts = lines(new File(dir, "liftover_intervals.txt"));
    File chain = new File(dir, "liftover.chain");
    List<Interval> lifted = time("liftover", () -> {
      LiftOver lift = new LiftOver(chain);
      List<Interval> out = new ArrayList<>();
      for (String[] q : lifts) {
        out.add(lift.liftOver(new Interval(q[0], Integer.parseInt(q[1]), Integer.parseInt(q[2]))));
      }
      return out;
    });
    text = new StringBuilder();
    for (Interval i : lifted) {
      if (i == null) {
        text.append("null\n");
      } else {
        text.append(i.getContig()).append('\t').append(i.getStart()).append('\t')
            .append(i.getEnd()).append('\t').append(i.isNegativeStrand() ? '-' : '+').append('\n');
      }
    }
    write(java("liftover", ".txt"), text.toString());
  }

  public static void main(String[] args) throws Exception {
    // LiftOver logs every interval it fails to lift at INFO, which would time the logger.
    Log.setGlobalLogLevel(Log.LogLevel.WARNING);
    dir = new File(args[0]);
    reps = args.length > 1 ? Integer.parseInt(args[1]) : 3;
    List<SAMRecord> records = read(new File(dir, "corpus.bam"));
    SAMFileHeader h = header;
    sam(records, h);
    sorting(records, h);
    indexes();
    blockCopy();
    fastq(records);
    List<ReferenceSequence> genome = fasta();
    mdNm(records, genome.get(0).getBases());
    intervals();
  }
}
