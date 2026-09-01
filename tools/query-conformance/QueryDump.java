import htsjdk.samtools.Chunk;
import htsjdk.samtools.QueryInterval;
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
