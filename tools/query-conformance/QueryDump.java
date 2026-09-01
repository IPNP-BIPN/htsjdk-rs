import htsjdk.samtools.BAMQueryMultipleIntervalsIteratorFilter;
import htsjdk.samtools.Chunk;
import htsjdk.samtools.GenomicIndexUtil;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.QueryInterval;
import java.util.BitSet;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * `QueryInterval`'s ordering and merging, and `Chunk`'s overlap, adjacency and coalescing.
 *
 *   compare   <a> <b> <int>              QueryInterval.compareTo, which is a difference
 *   overlaps  <a> <b> <true|false>
 *   abuts     <a> <b> <true|false>
 *   optimize  <inputs> <outputs>         QueryInterval.optimizeIntervals
 *   chunkcmp  <a> <b> <int>
 *   chunkovl  <a> <b> <true|false>
 *   chunkadj  <a> <b> <true|false>
 *   chunkopt  <minimumOffset> <inputs> <outputs>
 *
 * An interval is written `ref:start-end` and a chunk `blockAddress:blockOffset-blockAddress:blockOffset`,
 * both as htsjdk's own toString does, so a mismatch reads as the objects rather than as numbers.
 *
 * `Chunk.overlaps`, `isAdjacentTo` and `optimizeChunkList` are package-private, so they are reached
 * by reflection: the alternative is to test them through a real index, which would measure the
 * index rather than the arithmetic.
 */
public class QueryDump {
  static QueryInterval qi(int ref, int start, int end) {
    return new QueryInterval(ref, start, end);
  }

  static String show(QueryInterval i) {
    return i.toString();
  }

  static String show(List<QueryInterval> list) {
    if (list.isEmpty()) return "[]";
    StringBuilder sb = new StringBuilder();
    for (QueryInterval i : list) {
      if (sb.length() > 0) sb.append(",");
      sb.append(show(i));
    }
    return sb.toString();
  }

  static String showChunk(Chunk c) {
    return c.toString();
  }

  static String showChunks(List<Chunk> list) {
    if (list.isEmpty()) return "[]";
    StringBuilder sb = new StringBuilder();
    for (Chunk c : list) {
      if (sb.length() > 0) sb.append(",");
      sb.append(showChunk(c));
    }
    return sb.toString();
  }

  static long vfp(long block, int offset) {
    return (block << 16) | offset;
  }

  /** `GenomicIndexUtil.regionToBins`, as a comma-separated list or `null`. */
  static String bins(int start, int end) {
    BitSet set = GenomicIndexUtil.regionToBins(start, end);
    if (set == null) return "null";
    StringBuilder sb = new StringBuilder();
    for (int i = set.nextSetBit(0); i >= 0; i = set.nextSetBit(i + 1)) {
      if (sb.length() > 0) sb.append(",");
      sb.append(i);
    }
    return sb.length() == 0 ? "[]" : sb.toString();
  }

  static SAMFileHeader header() {
    SAMFileHeader h = new SAMFileHeader();
    SAMSequenceDictionary d = new SAMSequenceDictionary();
    d.addSequence(new SAMSequenceRecord("chr1", 100000));
    d.addSequence(new SAMSequenceRecord("chr2", 100000));
    h.setSequenceDictionary(d);
    return h;
  }

  /**
   * Record `i` of the filter corpus: reference `i % 2`, start `100 * (i + 1)`, and every third one
   * unmapped but PLACED at that start, which is the case the whole filter turns on.
   */
  static SAMRecord filterRecord(SAMFileHeader h, int i) {
    SAMRecord r = new SAMRecord(h);
    r.setReadName("read" + i);
    r.setReferenceIndex(i % 2);
    r.setAlignmentStart(100 * (i + 1));
    r.setReadBases("ACGTACGTAC".getBytes());
    r.setBaseQualities(new byte[] {30, 30, 30, 30, 30, 30, 30, 30, 30, 30});
    if (i % 3 == 2) {
      r.setReadUnmappedFlag(true);
    } else {
      r.setCigarString("10M");
      r.setMappingQuality(60);
    }
    return r;
  }

  public static void main(String[] args) throws Exception {
    QueryInterval[] intervals = {
      qi(0, 100, 200), qi(0, 150, 250), qi(0, 201, 300), qi(0, 100, 0),
      qi(0, 300, 0), qi(1, 100, 200), qi(1, 100, 100), qi(0, 1, 1),
      qi(2, 500, 400),
    };
    for (int a = 0; a < intervals.length; a++) {
      for (int b = 0; b < intervals.length; b++) {
        System.out.printf(
            "compare\t%s\t%s\t%d%n", show(intervals[a]), show(intervals[b]),
            intervals[a].compareTo(intervals[b]));
        System.out.printf(
            "overlaps\t%s\t%s\t%b%n", show(intervals[a]), show(intervals[b]),
            intervals[a].overlaps(intervals[b]));
        System.out.printf(
            "abuts\t%s\t%s\t%b%n", show(intervals[a]), show(intervals[b]),
            intervals[a].endsAtStartOf(intervals[b]));
      }
    }

    int[][] subsets = {{0, 1}, {0, 2}, {0, 3}, {3, 4}, {5, 6}, {0, 1, 2, 3, 4, 5, 6, 7}, {8}, {}};
    for (int[] subset : subsets) {
      QueryInterval[] input = new QueryInterval[subset.length];
      for (int k = 0; k < subset.length; k++) input[k] = intervals[subset[k]];
      List<QueryInterval> before = new ArrayList<>(Arrays.asList(input));
      QueryInterval[] optimized = QueryInterval.optimizeIntervals(input.clone());
      System.out.printf(
          "optimize\t%s\t%s%n", show(before), show(Arrays.asList(optimized)));
    }

    int[][] regions = {
      {1, 100}, {1, 0}, {0, 0}, {100, 50}, {1, 16384}, {16385, 32768},
      {1, 536870912}, {536870911, 536870912}, {-5, 100}, {1, -1}, {1, 1},
    };
    for (int[] region : regions) {
      System.out.printf("bins\t%d\t%d\t%s%n", region[0], region[1], bins(region[0], region[1]));
    }

    SAMFileHeader fh = header();
    QueryInterval[] filterIntervals = {qi(0, 100, 200), qi(0, 500, 900), qi(1, 100, 100000)};
    for (int f = 0; f < filterIntervals.length; f++) {
      for (int i = 0; i < 12; i++) {
        System.out.printf(
            "cmprec\t%s\t%d\t%s%n",
            show(filterIntervals[f]), i,
            BAMQueryMultipleIntervalsIteratorFilter.compareIntervalToRecord(
                filterIntervals[f], filterRecord(fh, i)));
      }
    }
    for (boolean contained : new boolean[] {false, true}) {
      BAMQueryMultipleIntervalsIteratorFilter filter =
          new BAMQueryMultipleIntervalsIteratorFilter(filterIntervals, contained);
      for (int i = 0; i < 12; i++) {
        System.out.printf(
            "filter\t%b\t%d\t%s%n", contained, i, filter.compareToFilter(filterRecord(fh, i)));
      }
    }

    Chunk[] chunks = {
      new Chunk(vfp(10, 0), vfp(20, 40)),
      new Chunk(vfp(20, 40), vfp(30, 0)),
      new Chunk(vfp(20, 41), vfp(30, 0)),
      new Chunk(vfp(1, 0), vfp(2, 0)),
      new Chunk(vfp(50, 0), vfp(60, 0)),
      new Chunk(vfp(10, 0), vfp(20, 40)),
    };
    Method overlaps = Chunk.class.getDeclaredMethod("overlaps", Chunk.class);
    Method adjacent = Chunk.class.getDeclaredMethod("isAdjacentTo", Chunk.class);
    overlaps.setAccessible(true);
    adjacent.setAccessible(true);
    for (int a = 0; a < chunks.length; a++) {
      for (int b = 0; b < chunks.length; b++) {
        System.out.printf(
            "chunkcmp\t%s\t%s\t%d%n", showChunk(chunks[a]), showChunk(chunks[b]),
            chunks[a].compareTo(chunks[b]));
        System.out.printf(
            "chunkovl\t%s\t%s\t%b%n", showChunk(chunks[a]), showChunk(chunks[b]),
            (Boolean) overlaps.invoke(chunks[a], chunks[b]));
        System.out.printf(
            "chunkadj\t%s\t%s\t%b%n", showChunk(chunks[a]), showChunk(chunks[b]),
            (Boolean) adjacent.invoke(chunks[a], chunks[b]));
      }
    }

    Method optimizeChunks = Chunk.class.getDeclaredMethod("optimizeChunkList", List.class, long.class);
    optimizeChunks.setAccessible(true);
    long[] minimums = {0, vfp(2, 0), vfp(10, 0), vfp(25, 0)};
    int[][] chunkSubsets = {{0, 1}, {0, 2}, {3, 4}, {0, 1, 2, 3, 4, 5}, {}};
    for (long minimum : minimums) {
      for (int[] subset : chunkSubsets) {
        List<Chunk> input = new ArrayList<>();
        for (int k : subset) input.add(new Chunk(chunks[k].getChunkStart(), chunks[k].getChunkEnd()));
        List<Chunk> before = new ArrayList<>(input);
        @SuppressWarnings("unchecked")
        List<Chunk> result = (List<Chunk>) optimizeChunks.invoke(null, input, minimum);
        System.out.printf(
            "chunkopt\t%d\t%s\t%s%n", minimum, showChunks(before), showChunks(result));
      }
    }
  }
}
