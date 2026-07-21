/*
 * Probe: does re-sorting an already-sorted VCF header change it?
 *
 * VCFWriter.setHeader, when doNotWriteGenotypes is set, rebuilds the header:
 *
 *     this.mHeader = new VCFHeader(header.getMetaDataInSortedOrder());
 *
 * getMetaDataInSortedOrder() hands back a TreeSet's iteration order, and the VCFHeader
 * constructor pours it into a *new* TreeSet with the same comparator. For a consistent
 * comparator that is a no-op. Decision 0016 established that VCFHeaderLine.compareTo is not
 * consistent - VCFContigHeaderLine compares to other contigs by index and to everything else by
 * rendered string, which admits a cycle - and re-inserting a sequence into a TreeSet under an
 * inconsistent comparator is not guaranteed to reproduce it.
 *
 * So this asks the oracle directly: write the same header both ways and compare the bytes.
 *
 * Prints RESORT_STABLE=true or false, and the two headers when they differ.
 */

import htsjdk.variant.vcf.*;
import htsjdk.variant.variantcontext.writer.*;
import java.io.ByteArrayOutputStream;
import java.util.*;

public class ResortProbe {

    static String write(VCFHeader header, boolean doNotWriteGenotypes) {
        ByteArrayOutputStream out = new ByteArrayOutputStream();
        VariantContextWriterBuilder b = new VariantContextWriterBuilder()
                .setOutputVCFStream(out)
                .unsetOption(Options.INDEX_ON_THE_FLY);
        if (doNotWriteGenotypes) b.setOption(Options.DO_NOT_WRITE_GENOTYPES);
        VariantContextWriter w = b.build();
        w.writeHeader(header);
        w.close();
        return out.toString();
    }

    /** The three contigs from decision 0016, whose comparison cycles. */
    static Set<VCFHeaderLine> cyclic() {
        Set<VCFHeaderLine> lines = new LinkedHashSet<>();
        String[] ids = {"mmm", "zzz", "aaa"};
        for (int i = 0; i < ids.length; i++) {
            final String id = ids[i];
            lines.add(new VCFContigHeaderLine(
                    new LinkedHashMap<String, String>() {{ put("ID", id); }}, i));
        }
        // Something that is not a contig, so the string-comparison branch is reachable.
        lines.add(new VCFInfoHeaderLine("DP", 1, VCFHeaderLineType.Integer, "Depth"));
        return lines;
    }

    static void compare(String label, Set<VCFHeaderLine> lines) {
        VCFHeader h = new VCFHeader(lines, Arrays.asList("s1"));

        // The ordinary path.
        String direct = write(h, false);
        // The path that rebuilds the header from its own sorted order. The genotype columns
        // differ by construction, so only the metadata lines are compared.
        String resorted = write(h, true);

        String a = metadataOnly(direct);
        String b = metadataOnly(resorted);
        System.out.println(label + "_RESORT_STABLE=" + a.equals(b));
        if (!a.equals(b)) {
            System.out.println("  direct  : " + a.replace("\n", " | "));
            System.out.println("  resorted: " + b.replace("\n", " | "));
        }
    }

    static String metadataOnly(String vcf) {
        StringBuilder sb = new StringBuilder();
        for (String line : vcf.split("\n")) {
            if (line.startsWith("##")) sb.append(line).append("\n");
        }
        return sb.toString();
    }

    public static void main(String[] args) {
        compare("CYCLIC", cyclic());

        // A control: contigs whose index order and string order agree, so the comparator is
        // consistent over this set and re-sorting must be a no-op.
        Set<VCFHeaderLine> ordered = new LinkedHashSet<>();
        String[] ids = {"aaa", "mmm", "zzz"};
        for (int i = 0; i < ids.length; i++) {
            final String id = ids[i];
            ordered.add(new VCFContigHeaderLine(
                    new LinkedHashMap<String, String>() {{ put("ID", id); }}, i));
        }
        ordered.add(new VCFInfoHeaderLine("DP", 1, VCFHeaderLineType.Integer, "Depth"));
        compare("CONSISTENT", ordered);

        // A realistic dictionary, where the index order is the natural genomic one and the
        // string order is not: chr10 sorts before chr2.
        Set<VCFHeaderLine> realistic = new LinkedHashSet<>();
        String[] chroms = {"chr1", "chr2", "chr3", "chr10", "chr11", "chrX"};
        for (int i = 0; i < chroms.length; i++) {
            final String id = chroms[i];
            realistic.add(new VCFContigHeaderLine(
                    new LinkedHashMap<String, String>() {{ put("ID", id); put("length", "1000"); }}, i));
        }
        realistic.add(new VCFInfoHeaderLine("DP", 1, VCFHeaderLineType.Integer, "Depth"));
        compare("REALISTIC", realistic);
    }
}
