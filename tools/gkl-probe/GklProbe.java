/*
 * Which deflater produces GATK's default BAM bytes, and are those bytes a property of the
 * algorithm or of the CPU that ran it?
 *
 * `libgkl_compression.so` bundles BOTH ISA-L's igzip (igzip_base.c, encode_deflate_icf,
 * IGZIP_DIST_TABLE_SIZE) and zlib (deflate_fast, deflate_medium, deflate_slow, and the string
 * "deflate 1.2.13 Copyright 1995-2022 Jean-loup Gailly and Mark Adler"). Which one runs depends on
 * the compression level, and htsjdk's BGZF default level decides which one a BAM is written with.
 *
 * Decision 0028 measured the split by compressing the same bytes both ways: levels 1 to 6 differ
 * from the JDK's zlib, levels 7 to 9 are byte-identical to it. Decision 0029 then read the level
 * branch out of the library rather than inferring it, and the split is not where 0028 said: only
 * levels 1 and 2 reach igzip, and 3 to 9 reach a zlib 1.2.13 carrying Intel's `deflate_medium`
 * patch, which the JDK's zlib 1.3.2 disagrees with below level 7. The default is 5, so the
 * default path is that patched zlib.
 *
 * THIS PROBE ANSWERS THE SECOND HALF, which 0028 named and deliberately left open. igzip ships
 * hand-written AVX2 and AVX512 kernels and dispatches on CPU features at load time. If those
 * kernels emit different bytes from the base C path, then "igzip output" is not one thing and a
 * port has no fixed target: the same file compressed on two machines would differ, exactly as
 * `Math.pow` does in decision 0007. That would move H.4 out of the reproducible column entirely,
 * which is a larger result than any amount of porting.
 *
 * It is answered the way 0007 answered pow: run here, run on a real x86-64 CI host, diff. This
 * laptop translates linux/amd64 through Rosetta, which does not implement AVX2, so the local run
 * exercises igzip's SSE path while the CI run exercises its AVX2 path. The two columns are the
 * experiment.
 *
 * Output is a hash per (fixture, level, backend), not a length, because two deflate streams of
 * equal length are not equal streams and a length-only comparison would report an agreement it
 * never checked.
 *
 * ONE ASSERTION GUARDS THE WHOLE THING. If the native library fails to load, GKL falls back to the
 * JDK deflater silently, every gkl line would equal its jdk line, and the probe would report a
 * beautifully stable result about nothing. `usingIntelDeflater` is checked, and a false answer is
 * a failure rather than a footnote.
 *
 * Output:
 *
 *     env\t<gkl native>\t<java version>\t<os.arch>
 *     cpu\t<model name, or "unavailable">
 *     default-bgzf-level\t<level>
 *     fixture\t<name>\t<bytes>\t<sha256 of the input>
 *     deflate\t<fixture>\t<level>\t<gkl|jdk>\t<bytes>\t<sha256 of the output>
 *     bgzf\t<fixture>\t<level>\t<gkl|jdk>\t<bytes>\t<sha256 of the whole BGZF stream>
 *
 * Usage: GklProbe
 */
import com.intel.gkl.compression.IntelDeflaterFactory;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.BlockCompressedStreamConstants;
import htsjdk.samtools.util.zip.DeflaterFactory;

import java.io.ByteArrayOutputStream;
import java.io.OutputStream;
import java.math.BigInteger;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Random;
import java.util.zip.Deflater;

public class GklProbe {

    public static void main(final String[] args) throws Exception {
        final IntelDeflaterFactory gkl = new IntelDeflaterFactory();
        if (!gkl.usingIntelDeflater()) {
            // Not a warning. Without the native library every comparison below is the JDK against
            // itself, which agrees trivially and means nothing.
            System.err.println("FATAL: GKL native deflater did not load; this probe would compare "
                    + "the JDK against itself.");
            System.exit(1);
        }
        System.out.printf("env\tgkl-native=true\tjava=%s\tarch=%s%n",
                System.getProperty("java.version"), System.getProperty("os.arch"));
        System.out.printf("cpu\t%s%n", cpuModel());
        System.out.printf("default-bgzf-level\t%d%n",
                BlockCompressedStreamConstants.DEFAULT_COMPRESSION_LEVEL);

        final Map<String, byte[]> fixtures = new LinkedHashMap<>();
        // The shape a BAM's sequence column has: four symbols, no structure beyond that.
        fixtures.put("acgt", bases(60000, 7));
        // Long runs, which is where a match-finder's choice among equally long matches matters
        // most, and so where two implementations are likeliest to disagree while producing streams
        // of the same length.
        fixtures.put("runs", runs(60000));
        // Incompressible, so the stored path is exercised and not only the huffman one.
        fixtures.put("random", noise(60000, 11));
        // Larger than one BGZF block's 65280-byte payload, so the multi-block framing and the
        // deflater reuse between blocks are under test rather than a single shot.
        fixtures.put("acgt-2blocks", bases(200000, 13));

        for (final Map.Entry<String, byte[]> fixture : fixtures.entrySet()) {
            System.out.printf("fixture\t%s\t%d\t%s%n", fixture.getKey(), fixture.getValue().length,
                    sha256(fixture.getValue()));
        }

        for (final Map.Entry<String, byte[]> fixture : fixtures.entrySet()) {
            for (int level = 1; level <= 9; level++) {
                emit("deflate", fixture.getKey(), level, "gkl",
                        deflate(gkl.makeDeflater(level, true), fixture.getValue()));
                emit("deflate", fixture.getKey(), level, "jdk",
                        deflate(new Deflater(level, true), fixture.getValue()));
            }
        }

        // The real path. A BAM is not a deflate stream; it is BGZF, a sequence of gzip members
        // each carrying a CRC and an extra field holding the block length. Hashing that covers the
        // framing as well as the compressed bytes, and the framing is where a block-size decision
        // would show up that the raw deflate comparison cannot see.
        for (final Map.Entry<String, byte[]> fixture : fixtures.entrySet()) {
            for (final int level : new int[] {
                    1, BlockCompressedStreamConstants.DEFAULT_COMPRESSION_LEVEL, 9 }) {
                emit("bgzf", fixture.getKey(), level, "gkl", bgzf(fixture.getValue(), level, gkl));
                emit("bgzf", fixture.getKey(), level, "jdk",
                        bgzf(fixture.getValue(), level, new DeflaterFactory()));
            }
        }
    }

    static void emit(final String kind, final String fixture, final int level, final String backend,
                     final byte[] out) throws Exception {
        System.out.printf("%s\t%s\t%d\t%s\t%d\t%s%n", kind, fixture, level, backend, out.length,
                sha256(out));
    }

    static byte[] deflate(final Deflater deflater, final byte[] data) {
        deflater.reset();
        deflater.setInput(data);
        deflater.finish();
        final byte[] out = new byte[data.length * 2 + 1024];
        int n = 0;
        while (!deflater.finished() && n < out.length) {
            final int written = deflater.deflate(out, n, out.length - n);
            if (written == 0) {
                break;
            }
            n += written;
        }
        deflater.end();
        return Arrays.copyOf(out, n);
    }

    static byte[] bgzf(final byte[] data, final int level, final DeflaterFactory factory)
            throws Exception {
        final ByteArrayOutputStream sink = new ByteArrayOutputStream();
        try (final OutputStream out =
                new BlockCompressedOutputStream(sink, (Path) null, level, factory)) {
            out.write(data);
        }
        return sink.toByteArray();
    }

    static byte[] bases(final int length, final long seed) {
        final byte[] data = new byte[length];
        final Random random = new Random(seed);
        for (int i = 0; i < length; i++) {
            data[i] = (byte) "ACGT".charAt(random.nextInt(4));
        }
        return data;
    }

    static byte[] runs(final int length) {
        final byte[] data = new byte[length];
        for (int i = 0; i < length; i++) {
            data[i] = (byte) "ACGT".charAt((i / 300) % 4);
        }
        return data;
    }

    static byte[] noise(final int length, final long seed) {
        final byte[] data = new byte[length];
        new Random(seed).nextBytes(data);
        return data;
    }

    static String sha256(final byte[] data) throws Exception {
        final byte[] digest = MessageDigest.getInstance("SHA-256").digest(data);
        final String hex = new BigInteger(1, digest).toString(16);
        return "0".repeat(64 - hex.length()) + hex;
    }

    static String cpuModel() {
        // Self-describing output: the CPU is this experiment's independent variable, so a run that
        // does not say which one it was cannot be compared against another.
        try {
            for (final String line : Files.readAllLines(Path.of("/proc/cpuinfo"),
                    StandardCharsets.UTF_8)) {
                if (line.startsWith("model name")) {
                    return line.substring(line.indexOf(':') + 1).trim();
                }
            }
        } catch (final Exception ignored) {
            // Not Linux, or no procfs. The line still prints, saying so.
        }
        return "unavailable";
    }
}
