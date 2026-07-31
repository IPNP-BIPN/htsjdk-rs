/*
 * VCFCodec.canDecode and VCF3Codec.canDecode, taken from the reference.
 *
 * This is the third codec GATK's -L reaches, and the odd one out: it OPENS THE FILE. BEDCodec and
 * IntervalListCodec answer on the path alone, so a `.list` holding a BED body is not a Feature
 * file and dies in the interval reader; a `.list` holding a VCF body IS one, because this codec
 * reads the first eighteen bytes and finds the magic there. Two files with the same extension and
 * different contents take different branches of -L, and only this codec makes that true.
 *
 * Three attempts in order: plain, then GZIPInputStream, then BlockCompressedInputStream. A BGZF
 * file is a gzip file, so the second attempt already answers for it and the third is unreachable
 * for well-formed input.
 *
 * `nread` is computed and never used, so a file shorter than the magic leaves zero bytes in the
 * tail of the buffer, the comparison fails, and the answer is false rather than an exception. An
 * absent path is the same: canDecodeFile catches IOException and answers false.
 *
 * The magic carries no minor version (`##fileformat=VCFv4`), so 4.0 through 4.3 all match and 3.3
 * does not; VCF3Codec carries the other magic, and GATK registers both.
 *
 * Output:
 *
 *     candecode\t<label>\t<extension>\t<VCFCodec>\t<VCF3Codec>
 *
 * Usage: VcfCanDecodeDump
 */

import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.variant.vcf.VCF3Codec;
import htsjdk.variant.vcf.VCFCodec;

import java.io.ByteArrayOutputStream;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.zip.GZIPOutputStream;

public class VcfCanDecodeDump {

    static final String VCF4 = "##fileformat=VCFv4.2\n"
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            + "chr1\t100\t.\tA\tC\t.\t.\t.\n";
    static final String VCF3 = "##fileformat=VCFv3.3\n"
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    public static void main(final String[] args) throws Exception {
        System.out.println("# VcfCanDecodeDump: VCFCodec.canDecode and VCF3Codec.canDecode");

        final Path dir = Path.of("vcf-candecode-dump");
        deleteRecursively(dir);
        Files.createDirectories(dir);

        // The extension is not consulted, so the same body answers the same under any of them.
        probe(dir, "vcf4-plain", "a.vcf", VCF4.getBytes());
        probe(dir, "vcf4-list-extension", "a.list", VCF4.getBytes());
        probe(dir, "vcf4-bed-extension", "a.bed", VCF4.getBytes());
        probe(dir, "vcf4-no-extension", "plain", VCF4.getBytes());

        // Nor is the contents' plausibility: only the first eighteen bytes decide.
        probe(dir, "magic-only", "b.vcf", "##fileformat=VCFv4".getBytes());
        probe(dir, "magic-then-junk", "c.vcf", "##fileformat=VCFv4NONSENSE".getBytes());
        probe(dir, "leading-space", "d.vcf", (" " + VCF4).getBytes());
        probe(dir, "leading-newline", "e.vcf", ("\n" + VCF4).getBytes());

        // Versions. The magic has no minor version, and VCF 3 has its own codec.
        probe(dir, "vcf40", "f.vcf", "##fileformat=VCFv4.0\n".getBytes());
        probe(dir, "vcf43", "g.vcf", "##fileformat=VCFv4.3\n".getBytes());
        probe(dir, "vcf3", "h.vcf", VCF3.getBytes());
        probe(dir, "vcf5", "i.vcf", "##fileformat=VCFv5.0\n".getBytes());

        // Short and empty files, which are false rather than errors.
        probe(dir, "truncated-magic", "j.vcf", "##fileformat=V".getBytes());
        probe(dir, "one-byte", "k.vcf", "#".getBytes());
        probe(dir, "empty", "l.vcf", new byte[0]);

        // Compression: gzip and BGZF, which the second attempt already covers.
        probe(dir, "vcf4-gzip", "m.vcf.gz", gzip(VCF4.getBytes()));
        probe(dir, "vcf4-bgzf", "n.vcf.gz", bgzf(VCF4.getBytes()));
        probe(dir, "vcf4-bgzf-bed-extension", "o.bed", bgzf(VCF4.getBytes()));
        probe(dir, "gzip-not-vcf", "p.vcf.gz", gzip("hello there\n".getBytes()));
        // Gzip magic followed by nothing that inflates.
        probe(dir, "broken-gzip", "q.vcf.gz", new byte[] {0x1f, (byte) 0x8b, 0x08, 0, 0, 0, 0, 0});

        // Bodies that are not VCF at all.
        probe(dir, "bed-body", "r.bed", "chr1\t0\t10\n".getBytes());
        probe(dir, "interval-list-body", "s.interval_list",
                "@HD\tVN:1.6\nchr1\t1\t10\t+\t.\n".getBytes());

        // A path that does not exist: caught, and false.
        System.out.printf("candecode\t%s\t%s\t%b\t%b%n", "absent", ".vcf",
                new VCFCodec().canDecode(dir.resolve("absent.vcf").toString()),
                new VCF3Codec().canDecode(dir.resolve("absent.vcf").toString()));
        // And a directory, which opens and fails to read.
        System.out.printf("candecode\t%s\t%s\t%b\t%b%n", "directory", "",
                new VCFCodec().canDecode(dir.toString()),
                new VCF3Codec().canDecode(dir.toString()));
    }

    static void probe(final Path dir, final String label, final String name, final byte[] body)
            throws Exception {
        final Path path = dir.resolve(name);
        Files.write(path, body);
        final String extension = name.contains(".") ? name.substring(name.indexOf('.')) : "";
        System.out.printf("candecode\t%s\t%s\t%b\t%b%n", label, extension,
                new VCFCodec().canDecode(path.toString()),
                new VCF3Codec().canDecode(path.toString()));
    }

    static byte[] gzip(final byte[] body) throws Exception {
        final ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        try (OutputStream out = new GZIPOutputStream(bytes)) {
            out.write(body);
        }
        return bytes.toByteArray();
    }

    static byte[] bgzf(final byte[] body) throws Exception {
        final ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        try (OutputStream out = new BlockCompressedOutputStream(bytes, (Path) null)) {
            out.write(body);
        }
        return bytes.toByteArray();
    }

    static void deleteRecursively(final Path path) throws Exception {
        if (!Files.exists(path)) {
            return;
        }
        try (java.util.stream.Stream<Path> walk = Files.walk(path)) {
            walk.sorted(java.util.Comparator.reverseOrder()).forEach(p -> {
                try {
                    Files.delete(p);
                } catch (final Exception e) {
                    throw new IllegalStateException(e);
                }
            });
        }
    }
}
