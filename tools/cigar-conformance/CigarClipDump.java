/*
 * Oracle dump harness for CigarUtil.softClipEndOfRead conformance in htsjdk-rs.
 *
 * For each (clipFrom, cigar) case it prints `clip<TAB>clipFrom:inputCigar<TAB>resultCigar`. The
 * committed corpus is `java ... CigarClipDump | gzip > tests/data/cigar_clip.txt.gz`, regenerated
 * and compared in CI. The Rust port parses the input cigar, soft-clips from clipFrom, and must
 * reproduce the result cigar's text.
 *
 *   java -cp htsjdk.jar:. CigarClipDump
 */
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.util.CigarUtil;
import htsjdk.samtools.TextCigarCodec;

import java.util.List;

public class CigarClipDump {
    public static void main(final String[] args) {
        // (clipFrom, cigar): a spread of straddle positions, indels at the boundary, and
        // leading/trailing clips already present.
        final Object[][] cases = {
            {30, "36M"}, {25, "36M"}, {1, "36M"}, {36, "36M"},
            {20, "10M5I21M"}, {12, "10M5I21M"}, {11, "10M5I21M"},
            {15, "10M2D26M"}, {20, "19M1I16M"},
            {10, "5S31M"}, {30, "30M6S"}, {28, "30M6S"},
            {30, "31M5H"}, {20, "5H31M"}, {25, "5H26M5H"},
            {8, "3M2I3M2I3M"},
        };
        for (final Object[] c : cases) {
            final int clipFrom = (Integer) c[0];
            final String cigarStr = (String) c[1];
            final Cigar cigar = TextCigarCodec.decode(cigarStr);
            final List<CigarElement> result = CigarUtil.softClipEndOfRead(clipFrom, cigar.getCigarElements());
            System.out.println("clip\t" + clipFrom + ":" + cigarStr + "\t" + new Cigar(result).toString());
        }
    }
}
