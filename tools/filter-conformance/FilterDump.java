import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.filter.AlignedFilter;
import htsjdk.samtools.filter.ReadNameFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.filter.TagFilter;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

/**
 * Every answer the three record filters give, over a corpus built by index so the port can rebuild
 * the same one without a fixture file.
 *
 * Both forms of `filterOut` are dumped, because the pair form is not the single form applied twice:
 * each filter is asymmetric in its own way, and those asymmetries are the whole reason a port of
 * these classes can be wrong while every single-record row agrees.
 *
 *   single <filter> <record> <true|false>
 *   pair   <filter> <first> <second> <true|false>
 *
 * `true` means the record is DROPPED, which is htsjdk's sense of the word and the easiest thing to
 * invert by accident.
 */
public class FilterDump {
  static SAMFileHeader header() {
    SAMFileHeader h = new SAMFileHeader();
    SAMSequenceDictionary d = new SAMSequenceDictionary();
    d.addSequence(new SAMSequenceRecord("chr1", 1000));
    h.setSequenceDictionary(d);
    return h;
  }

  /**
   * Record `i` of the corpus: name `read<i>`, unmapped when `i` is odd, and an `RG` of `rg1` when
   * `i % 3 == 0`, `rg2` when `i % 3 == 1`, and no tag at all when `i % 3 == 2`. Three residues
   * against two mapping states covers every combination the filters can see.
   */
  static SAMRecord record(SAMFileHeader h, int i) {
    SAMRecord r = new SAMRecord(h);
    r.setReadName("read" + i);
    r.setReadBases("ACGT".getBytes());
    r.setBaseQualities(new byte[] {30, 30, 30, 30});
    if (i % 2 == 1) {
      r.setReadUnmappedFlag(true);
    } else {
      r.setReferenceIndex(0);
      r.setAlignmentStart(100 + i);
      r.setCigarString("4M");
      r.setMappingQuality(60);
    }
    if (i % 3 == 0) r.setAttribute("RG", "rg1");
    else if (i % 3 == 1) r.setAttribute("RG", "rg2");
    return r;
  }

  public static void main(String[] args) {
    SAMFileHeader h = header();
    int n = 12;

    Set<String> names = new LinkedHashSet<>(Arrays.asList("read0", "read3", "read4", "read11"));
    List<Object> values = Collections.singletonList((Object) "rg1");

    String[] labels = {
      "aligned_include", "aligned_exclude",
      "name_include", "name_exclude",
      "tag_include", "tag_exclude",
    };
    SamRecordFilter[] filters = {
      new AlignedFilter(true), new AlignedFilter(false),
      new ReadNameFilter(names, true), new ReadNameFilter(names, false),
      new TagFilter("RG", values, true), new TagFilter("RG", values, false),
    };

    for (int f = 0; f < filters.length; f++) {
      for (int i = 0; i < n; i++) {
        System.out.printf("single\t%s\t%d\t%b%n", labels[f], i, filters[f].filterOut(record(h, i)));
      }
      for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        System.out.printf(
            "pair\t%s\t%d\t%d\t%b%n",
            labels[f], i, j, filters[f].filterOut(record(h, i), record(h, j)));
      }
    }
  }
}
