/*
 * What a record's genotype columns are written from, and what changes it.
 *
 * `VCFCodec` does not parse the genotype columns of a line: it hands the substring to a
 * `LazyGenotypesContext` and returns. `VCFEncoder` then has two paths, and the first copies the
 * FORMAT column too:
 *
 *     if (gc.isLazyWithData() && ... instanceof String) { append(FIELD_SEPARATOR); append(text); }
 *     else { calcVCFGenotypeKeys, GT first, the rest sorted, trailing missing fields trimmed }
 *
 * So a record read and written untouched keeps the INPUT's key order, and a record whose genotypes
 * were looked at gets the recomputed one. This dump measures three things per input:
 *
 *   - COPY: read, write, nothing else. The verbatim branch.
 *   - TOUCH: read, call getGenotypes().get(0) -- a READ, not a modification -- then write.
 *   - REBUILD: read, rebuild the genotypes through a builder, then write.
 *
 * The inputs are chosen so the two orders differ: a FORMAT of GT:GQ:DP sorts to GT:DP:GQ, and a
 * record with a trailing missing field keeps it verbatim and loses it when rebuilt.
 *
 * Output:
 *     line\t<case>\t<mode>\t<the data line, escaped>
 *     lazy\t<case>\t<whether the context was lazy with data before the write>
 */

import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.SynchronousLineReader;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.LazyGenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFEncoder;
import htsjdk.variant.vcf.VCFHeader;

import java.io.StringReader;
import java.util.ArrayList;
import java.util.List;

public class VcfLazyGenotypesDump {

    static String escape(final String s) {
        return s.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
    }

    static String header(final String format) {
        return "##fileformat=VCFv4.2\n"
                + "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
                + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                + "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Quality\">\n"
                + "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
                + "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Depths\">\n"
                + "##contig=<ID=chr1,length=100000>\n"
                + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\n";
    }

    public static void main(final String[] args) {
        // The FORMAT order the file carries is NOT the order the encoder computes: GT is moved to
        // the front and the rest sort, so GQ:DP becomes DP:GQ.
        run("unsorted-format", "chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\tGT:GQ:DP\t0/1:60:10\t1/1:70:20");
        // Already in the computed order: the two paths agree, which is what makes the case above a
        // measurement of the branch rather than of the sort.
        run("sorted-format", "chr1\t200\t.\tC\tG\t50\tPASS\tDP=10\tGT:DP:GQ\t0/1:10:60\t1/1:20:70");
        // A trailing field that is missing in one sample: the verbatim branch keeps whatever the
        // file had, the recomputed branch applies the trailing-field trim per sample.
        run("trailing-missing", "chr1\t300\t.\tG\tA\t50\tPASS\tDP=10\tGT:GQ:AD\t0/1:60:5,6\t1/1:.:.");
        // A field declared in the header but absent from this line's FORMAT.
        run("short-format", "chr1\t400\t.\tT\tC\t50\tPASS\tDP=10\tGT\t0/1\t1/1");
    }

    static void run(final String label, final String line) {
        final String text = header(null) + line + "\n";
        emit(label, "copy", text, Mode.COPY);
        emit(label, "touch", text, Mode.TOUCH);
        emit(label, "rebuild", text, Mode.REBUILD);
    }

    enum Mode { COPY, TOUCH, REBUILD }

    static void emit(final String label, final String mode, final String text, final Mode what) {
        final VCFCodec codec = new VCFCodec();
        final LineIterator it =
                new LineIteratorImpl(new SynchronousLineReader(new StringReader(text)));
        final VCFHeader header = (VCFHeader) codec.readActualHeader(it);
        // `outputTrailingFormatFields = false`, which is what `VariantContextWriter` builds:
        // passing true here disables the trailing-missing-field trim and measures an encoder no
        // writer constructs.
        final VCFEncoder encoder = new VCFEncoder(header, true, false);
        while (it.hasNext()) {
            VariantContext vc = codec.decode(it.next());
            if (vc == null) {
                continue;
            }
            switch (what) {
                case COPY:
                    break;
                case TOUCH:
                    // A READ. `getGenotypes().get(0)` runs `ensureSampleOrdering`, which decodes,
                    // and `decode()` ends with `unparsedGenotypeData = null`.
                    vc.getGenotypes().get(0);
                    break;
                case REBUILD:
                    final List<Genotype> rebuilt = new ArrayList<>();
                    for (final Genotype g : vc.getGenotypes()) {
                        rebuilt.add(new GenotypeBuilder(g).make());
                    }
                    vc = new VariantContextBuilder(vc)
                            .genotypes(GenotypesContext.create(new ArrayList<>(rebuilt)))
                            .make();
                    break;
            }
            final GenotypesContext gc = vc.getGenotypes();
            final boolean lazy = gc.isLazyWithData()
                    && ((LazyGenotypesContext) gc).getUnparsedGenotypeData() instanceof String;
            System.out.printf("lazy\t%s\t%s\t%s%n", label, mode, lazy);
            System.out.printf("line\t%s\t%s\t%s%n", label, mode, escape(encoder.encode(vc)));
        }
    }
}
