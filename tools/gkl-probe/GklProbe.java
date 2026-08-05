/*
 * Which deflater actually produces GATK's default BAM bytes, and is it reproducible without
 * porting a compression algorithm?
 *
 * `libgkl_compression.so` bundles BOTH ISA-L's igzip (igzip_base.c, encode_deflate_icf,
 * IGZIP_DIST_TABLE_SIZE) and zlib (deflate_fast, deflate_medium, deflate_slow). Which one runs
 * depends on the compression level, and htsjdk's BGZF default level decides which one a BAM is
 * actually written with.
 *
 * So the decisive experiment is not to read the source but to compress the same bytes both ways
 * and compare:
 *
 *   - IntelDeflater at each level, against
 *   - java.util.zip.Deflater at the same level (the JDK's zlib).
 *
 * If they agree at the level htsjdk uses, then "GKL-exact deflate" is not a port of ISA-L at all:
 * it is zlib with the right settings, and the whole milestone changes size.
 */
import com.intel.gkl.compression.IntelDeflater;
import com.intel.gkl.compression.IntelDeflaterFactory;
import htsjdk.samtools.util.BlockCompressedStreamConstants;

import java.util.Arrays;
import java.util.Random;
import java.util.zip.Deflater;

public class GklProbe {

    public static void main(final String[] args) throws Exception {
        System.out.println("htsjdk default BGZF compression level = "
                + BlockCompressedStreamConstants.DEFAULT_COMPRESSION_LEVEL);

        // Data with the shape a BAM block has: repetitive, but not trivially so.
        final byte[] data = new byte[60000];
        final Random random = new Random(7);
        for (int i = 0; i < data.length; i++) {
            data[i] = (byte) ("ACGT".charAt(random.nextInt(4)));
        }

        final IntelDeflaterFactory factory = new IntelDeflaterFactory();
        System.out.println("GKL native library loaded = " + factory.usingIntelDeflater());

        for (int level = 1; level <= 9; level++) {
            final byte[] gkl = deflate(factory.makeDeflater(level, true), data);
            final byte[] jdk = deflate(new Deflater(level, true), data);
            System.out.printf("level %d: gkl=%d bytes jdk=%d bytes  %s%n", level, gkl.length,
                    jdk.length, Arrays.equals(gkl, jdk) ? "IDENTICAL" : "different");
        }
    }

    static byte[] deflate(final Deflater deflater, final byte[] data) {
        deflater.reset();
        deflater.setInput(data);
        deflater.finish();
        final byte[] out = new byte[data.length * 2 + 1024];
        final int n = deflater.deflate(out, 0, out.length);
        deflater.end();
        return Arrays.copyOf(out, n);
    }
}
