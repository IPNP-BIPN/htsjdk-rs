/*
 * The reference side of a differential fuzzer over the byte-level parsers.
 *
 * Reads one hex string per line on stdin, runs each through the named parser, and prints what the
 * reference did with it: a value, or the exception it threw. The port's side of the same corpus is
 * `crates/htsjdk-cram/examples/differential_fuzz.rs`, and the two are diffed line by line.
 *
 * The parsers here are the ones a hostile file reaches first, and each is a pure function of its
 * bytes: no state, no environment, nothing to seed. That is what makes a divergence a bug rather
 * than a difference of setup.
 *
 * Output, one line per input:
 *
 *     <hex>\t<parser>\tok:<value>
 *     <hex>\t<parser>\terr:<exception class>
 *
 * Usage: FuzzDriver <itf8|ltf8|crai|block> < corpus.hex
 */

import htsjdk.samtools.cram.CRAIEntry;
import htsjdk.samtools.cram.io.ITF8;
import htsjdk.samtools.cram.io.LTF8;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;

public class FuzzDriver {

    public static void main(final String[] args) throws Exception {
        final String parser = args.length > 0 ? args[0] : "itf8";
        try (final BufferedReader reader =
                new BufferedReader(new InputStreamReader(System.in, StandardCharsets.UTF_8))) {
            String line;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) {
                    continue;
                }
                System.out.printf("%s\t%s\t%s%n", line, parser, outcome(parser, unhex(line)));
            }
        }
    }

    static String outcome(final String parser, final byte[] bytes) {
        try {
            switch (parser) {
                case "itf8":
                    // The InputStream overload, which returns -1 past the end rather than
                    // throwing: that is the one the port models, and the one a CRAM reader uses.
                    return "ok:" + ITF8.readUnsignedITF8(new ByteArrayInputStream(bytes));
                case "ltf8":
                    return "ok:" + LTF8.readUnsignedLTF8(new ByteArrayInputStream(bytes));
                case "crai":
                    return "ok:" + new CRAIEntry(new String(bytes, StandardCharsets.UTF_8))
                            .toString().replace('\t', ' ');
                default:
                    return "err:UnknownParser";
            }
        } catch (final Throwable t) {
            return "err:" + t.getClass().getSimpleName();
        }
    }

    static byte[] unhex(final String hex) {
        final byte[] bytes = new byte[hex.length() / 2];
        for (int i = 0; i < bytes.length; i++) {
            bytes[i] = (byte) Integer.parseInt(hex.substring(i * 2, i * 2 + 2), 16);
        }
        return bytes;
    }
}
