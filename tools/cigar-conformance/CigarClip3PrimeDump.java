import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CigarUtil;

/**
 * `CigarUtil.softClip3PrimeEndOfRead`, which is a record transform rather than a cigar function.
 *
 *   clip3 <cigar> <strand> <start> <clipFrom> <newCigar> <newStart> <unmapped> <nmDropped>
 *
 * `strand` is `+` or `-`, and it decides which end of the stored bases is the read's three-prime
 * end: for a negative-strand record the cigar is reversed, clipped, and reversed back, and the
 * alignment start then moves by however much the reference span shrank.
 *
 * `unmapped` is what happens when nothing aligned survives the clip: the record loses its cigar,
 * its reference, its start and its mapping quality. `nmDropped` is the separate rule that NM, MD
 * and UQ are invalidated whenever the reference length changed at all, which is true far more often
 * than the record is unmapped.
 */
public class CigarClip3PrimeDump {
  static SAMFileHeader header() {
    SAMFileHeader h = new SAMFileHeader();
    SAMSequenceDictionary d = new SAMSequenceDictionary();
    d.addSequence(new SAMSequenceRecord("chr1", 100000));
    h.setSequenceDictionary(d);
    return h;
  }

  public static void main(String[] args) {
    SAMFileHeader h = header();
    String[] cigars = {
      "50M", "10S40M", "40M10S", "20M2I28M", "20M3D30M", "10M5N35M", "5H45M", "45M5H",
      "25M25S", "50S",
    };
    boolean[] strands = {false, true};
    int[] clipFroms = {1, 2, 10, 25, 40, 50};
    int start = 1000;

    for (String cigarText : cigars) {
      for (boolean negative : strands) {
        for (int clipFrom : clipFroms) {
          SAMRecord r = new SAMRecord(h);
          r.setReadName("read");
          r.setReferenceIndex(0);
          r.setAlignmentStart(start);
          r.setCigarString(cigarText);
          r.setMappingQuality(60);
          r.setReadNegativeStrandFlag(negative);
          byte[] bases = new byte[r.getCigar().getReadLength()];
          byte[] quals = new byte[bases.length];
          for (int i = 0; i < bases.length; i++) {
            bases[i] = (byte) "ACGT".charAt(i % 4);
            quals[i] = 30;
          }
          r.setReadBases(bases);
          r.setBaseQualities(quals);
          r.setAttribute("NM", 3);
          r.setAttribute("MD", "50");
          r.setAttribute("UQ", 42);

          String result;
          try {
            CigarUtil.softClip3PrimeEndOfRead(r, clipFrom);
            result =
                String.format(
                    "%s\t%d\t%b\t%b",
                    r.getCigarString(),
                    r.getAlignmentStart(),
                    r.getReadUnmappedFlag(),
                    r.getAttribute("NM") == null);
          } catch (Exception e) {
            result = "EXCEPTION\t" + e.getClass().getName();
          }
          System.out.printf(
              "clip3\t%s\t%s\t%d\t%d\t%s%n",
              cigarText, negative ? "-" : "+", start, clipFrom, result);
        }
      }
    }
  }
}
