import htsjdk.samtools.DuplicateScoringStrategy;
import htsjdk.samtools.DuplicateScoringStrategy.ScoringStrategy;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

/**
 * Every score and every comparison `DuplicateScoringStrategy` gives, over a corpus built by index.
 *
 *   score   <strategy> <record> <short>
 *   compare <strategy> <first> <second> <int>
 *
 * The corpus is chosen for the arithmetic rather than for realism: qualities on both sides of the
 * Q15 threshold, a read long enough to reach the Short.MAX_VALUE/2 clamp, vendor-failed records
 * that take the Short.MIN_VALUE/2 discount, paired and unpaired ends so `compare`'s first branch
 * fires, and names that differ so its tie-break does.
 *
 * `assumeMateCigar` is false throughout: the true branch reads the MC tag through
 * SAMUtils.getMateCigar, which throws when there is none, and a corpus that carried one would be
 * measuring the tag rather than the score.
 */
public class DuplicateScoringDump {
  static SAMFileHeader header() {
    SAMFileHeader h = new SAMFileHeader();
    SAMSequenceDictionary d = new SAMSequenceDictionary();
    d.addSequence(new SAMSequenceRecord("chr1", 100000));
    h.setSequenceDictionary(d);
    return h;
  }

  /** Record `i`: see the class comment for what each residue is for. */
  static SAMRecord record(SAMFileHeader h, int i) {
    SAMRecord r = new SAMRecord(h);
    r.setReadName("read" + i);
    int length = (i == 7) ? 5000 : 8;
    byte[] bases = new byte[length];
    byte[] quals = new byte[length];
    for (int b = 0; b < length; b++) {
      bases[b] = (byte) "ACGT".charAt(b % 4);
      // Below the Q15 threshold for even b when i is even, above it otherwise: both sides of the
      // filter appear in the same record.
      quals[b] = (byte) ((i % 2 == 0 && b % 2 == 0) ? 10 : 20 + (i % 5));
    }
    r.setReadBases(bases);
    r.setBaseQualities(quals);
    if (i % 4 == 3) {
      r.setReadUnmappedFlag(true);
    } else {
      r.setReferenceIndex(0);
      r.setAlignmentStart(100 + i);
      r.setCigarString(length + "M");
      r.setMappingQuality(60);
    }
    if (i % 3 != 2) {
      r.setReadPairedFlag(true);
      if (i % 3 == 0) r.setFirstOfPairFlag(true);
      else r.setSecondOfPairFlag(true);
      r.setMateUnmappedFlag(true);
    }
    if (i % 5 == 4) r.setReadFailsVendorQualityCheckFlag(true);
    return r;
  }

  public static void main(String[] args) {
    SAMFileHeader h = header();
    int n = 10;
    for (ScoringStrategy strategy : ScoringStrategy.values()) {
      for (int i = 0; i < n; i++) {
        System.out.printf(
            "score\t%s\t%d\t%d%n",
            strategy.name(), i, DuplicateScoringStrategy.computeDuplicateScore(record(h, i), strategy));
      }
      for (int i = 0; i < n; i++) {
        int j = (i + 3) % n;
        System.out.printf(
            "compare\t%s\t%d\t%d\t%d%n",
            strategy.name(),
            i,
            j,
            DuplicateScoringStrategy.compare(record(h, i), record(h, j), strategy));
      }
    }
  }
}
