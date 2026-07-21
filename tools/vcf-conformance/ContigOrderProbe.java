/*
 * Is VCFHeaderLine's ordering a total order?
 *
 * VCFContigHeaderLine.compareTo sorts contigs by contigIndex, and everything else by rendered
 * string. Mixing the two can make the comparator inconsistent: contig A may sort before contig
 * B by index while an INFO line sorts between them by string. A TreeSet with an inconsistent
 * comparator produces an order that depends on insertion sequence.
 *
 * This inserts the same lines in two different orders and prints both results.
 */

import htsjdk.variant.vcf.*;
import htsjdk.variant.variantcontext.writer.*;
import java.io.ByteArrayOutputStream;
import java.util.*;

public class ContigOrderProbe {

    static VCFContigHeaderLine contig(final String id, final int index) {
        return new VCFContigHeaderLine(
                new LinkedHashMap<String, String>() {{ put("ID", id); put("length", "1000"); }},
                index);
    }

    static String render(List<VCFHeaderLine> lines) {
        Set<VCFHeaderLine> set = new LinkedHashSet<>(lines);
        VCFHeader header = new VCFHeader(set);
        ByteArrayOutputStream out = new ByteArrayOutputStream();
        VariantContextWriter w = new VariantContextWriterBuilder()
                .setOutputVCFStream(out).unsetOption(Options.INDEX_ON_THE_FLY).build();
        w.writeHeader(header);
        w.close();
        StringBuilder sb = new StringBuilder();
        for (String l : out.toString().split("\n")) {
            if (l.startsWith("##") && !l.contains("fileformat")) sb.append(l).append(" | ");
        }
        return sb.toString();
    }

    public static void main(String[] args) {
        // Contig index order and string order disagree: zzz has the LOW index, aaa the high one.
        List<VCFHeaderLine> lines = new ArrayList<>();
        lines.add(contig("zzz", 0));
        lines.add(contig("aaa", 1));
        // A line whose KEY is also "contig" but which is not a VCFContigHeaderLine, so it
        // compares by string against the real ones. Its rendered form sorts between them.
        lines.add(new VCFHeaderLine("contig", "<ID=mmm,length=1000>"));

        List<VCFHeaderLine> forward = new ArrayList<>(lines);
        List<VCFHeaderLine> reversed = new ArrayList<>(lines);
        Collections.reverse(reversed);

        String a = render(forward);
        String b = render(reversed);
        System.out.println("inserted forward : " + a);
        System.out.println("inserted reversed: " + b);
        System.out.println("SAME_OUTPUT=" + a.equals(b));

        // And the pairwise comparisons that make it inconsistent.
        VCFHeaderLine zzz = contig("zzz", 0), aaa = contig("aaa", 1);
        VCFHeaderLine mmm = new VCFHeaderLine("contig", "<ID=mmm,length=1000>");
        System.out.println("zzz vs aaa (by index)  = " + Integer.signum(zzz.compareTo(aaa)));
        System.out.println("zzz vs mmm (by string) = " + Integer.signum(zzz.compareTo(mmm)));
        System.out.println("mmm vs aaa (by string) = " + Integer.signum(mmm.compareTo(aaa)));
    }
}
