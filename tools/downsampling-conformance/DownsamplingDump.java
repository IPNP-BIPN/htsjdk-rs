import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.DownsamplingIteratorFactory;
import htsjdk.samtools.DownsamplingIterator;
import java.util.ArrayList;
import java.util.List;

/**
 * What `ConstantMemoryDownsamplingIterator` keeps, and the statistics it reports for keeping it.
 *
 *   kept    <proportion> <seed> <record>  <true|false>
 *   stats   <proportion> <seed> seen=<n> accepted=<n> discarded=<n>
 *
 * The iterator is package-private, so it is reached through `DownsamplingIteratorFactory` with
 * `ConstantMemory` and an accuracy that keeps that strategy: the factory is what Picard's
 * DownsampleSam calls, and the strategy is its default.
 *
 * The proportions are chosen for the arithmetic rather than for realism. `1.0` is the one that
 * catches a port using 64-bit arithmetic: the threshold is
 * `Integer.MIN_VALUE + (int) Math.round(range * proportion)`, and at 1 that cast wraps back to
 * `Integer.MAX_VALUE`. `0.0` is its mirror. The rest bracket the middle.
 */
public class DownsamplingDump {
  static SAMFileHeader header() {
    return new SAMFileHeader();
  }

  static List<SAMRecord> records(SAMFileHeader h, int n) {
    List<SAMRecord> out = new ArrayList<>();
    for (int i = 0; i < n; i++) {
      SAMRecord r = new SAMRecord(h);
      r.setReadName("read" + i);
      r.setReadUnmappedFlag(true);
      r.setReadBases("ACGT".getBytes());
      r.setBaseQualities(new byte[] {30, 30, 30, 30});
      out.add(r);
    }
    return out;
  }

  public static void main(String[] args) {
    SAMFileHeader h = header();
    int n = 40;
    double[] proportions = {0.0, 0.01, 0.25, 0.5, 0.75, 0.99, 1.0};
    int[] seeds = {1, 42};

    for (double proportion : proportions) {
      for (int seed : seeds) {
        List<SAMRecord> corpus = records(h, n);
        DownsamplingIterator it =
            DownsamplingIteratorFactory.make(
                corpus.iterator(),
                DownsamplingIteratorFactory.Strategy.ConstantMemory,
                proportion,
                0.0001,
                seed);
        boolean[] kept = new boolean[n];
        while (it.hasNext()) {
          SAMRecord r = it.next();
          kept[Integer.parseInt(r.getReadName().substring(4))] = true;
        }
        for (int i = 0; i < n; i++) {
          System.out.printf("kept\t%s\t%d\t%d\t%b%n", Double.toString(proportion), seed, i, kept[i]);
        }
        System.out.printf(
            "stats\t%s\t%d\tseen=%d accepted=%d discarded=%d%n",
            Double.toString(proportion), seed, it.getSeenCount(), it.getAcceptedCount(), it.getDiscardedCount());
      }
    }
  }
}
