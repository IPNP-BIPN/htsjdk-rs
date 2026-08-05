/*
 * The CRAM substitution matrix: five bytes inside the preservation map that decide which
 * substitution gets the shortest code.
 *
 * cram-preservation-map measured the SM key as five bytes and left them opaque. This is what they
 * mean. Each byte is a packed vector of four two-bit codes, one per possible substitute of that
 * reference base, in the spec's base order A, C, G, T, N minus the base itself. The codes are
 * RANKS by observed frequency, so the commonest substitution gets code 0 and the shortest ITF8.
 *
 * Six things are decisions rather than layout.
 *
 *   - THE SORT IS RUN TWICE, AND THE SECOND RUN IS A SORT BY ORDINAL WEARING THE FREQUENCY
 *     COMPARATOR'S CLOTHES. substitutionCodeVector sorts by frequency, writes each entry's rank,
 *     then ZEROES EVERY FREQUENCY and sorts again with the same comparator, which now falls
 *     through to the ordinal tie-break every time. A port that sorts once and packs in sorted
 *     order emits the codes in the wrong slots;
 *   - THE COMPARATOR SUBTRACTS TWO LONGS AND CASTS TO INT. `(int) (o2.freq - o1.freq)`, so two
 *     frequencies whose difference is a multiple of 2^32 compare EQUAL and fall through to the
 *     ordinal tie-break. A substitution seen 4294967296 times can rank as if it were tied with one
 *     seen never;
 *   - A LOWER CASE REFERENCE BASE CAN BE DECODED BUT NOT ENCODED. `base` accepts it, because the
 *     reading constructor copies each row to its lower case twin; `code` throws on it by name;
 *   - `base`'s FAILURE MESSAGE BLAMES THE WRONG ARGUMENT. When the lookup lands on NO_BASE it
 *     formats the REFERENCE base into "Attempt to retrieve a substitution base for invalid base",
 *     even though the reference base is the one thing that was valid enough to index with;
 *   - THE DEFAULT MATRIX IS NOT ZEROES. With no substitutions observed, every frequency ties, the
 *     ordinal order wins, and each byte becomes 0x1b, which is the ranks 0, 1, 2, 3 packed two bits
 *     apiece;
 *   - THE MATRIX IS A SQUARE OF 128 BY 128 FOR FIVE BASES. Both lookup tables cover the whole
 *     symbol space so a base can index them directly, which is why an invalid base is caught by a
 *     sign test rather than by a bounds failure.
 *
 * Output:
 *
 *     sizes\t<basesSize>\t<codesPerBase>\t<symbolSpaceSize>
 *     vector\t<label>\t<refBase>\t<frequencies A,C,G,T,N>\t<code vector, unsigned>\t<code per substitute>
 *     decode\t<matrix hex>\t<refBase>\t<the four bases by code, in code order>
 *     lowercase\t<matrix hex>\t<refBase>\t<the four bases by code for the lower case twin>
 *     code\t<matrix hex>\t<refBase>\t<readBase>\t<code, or ERR class: message>
 *     base\t<matrix hex>\t<refBase>\t<code>\t<base, or ERR class: message>
 *     tostring\t<matrix hex>\t<toString with tabs replaced by spaces>
 *     err\t<label>\t<class>\t<message>
 *
 * Usage: CramSubstitutionMatrixDump
 */

import htsjdk.samtools.cram.structure.SubstitutionMatrix;

import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.util.Arrays;
import java.util.StringJoiner;

public class CramSubstitutionMatrixDump {

    private static final byte[] BASES = {'A', 'C', 'G', 'T', 'N'};

    public static void main(final String[] args) throws Exception {
        System.out.println("# CramSubstitutionMatrixDump: the five bytes behind the SM key");

        System.out.printf("sizes\t%d\t%d\t%d%n", SubstitutionMatrix.BASES_SIZE,
                SubstitutionMatrix.BASES_SIZE - 1, symbolSpaceSize());

        // No substitutions at all: every frequency ties and the ordinal order decides.
        vector("all-zero", frequencies(0, 0, 0, 0, 0));
        // One dominant substitute, then a strict ordering, then a partial tie.
        vector("one-dominant", frequencies(0, 100, 0, 0, 0));
        vector("strictly-ordered", frequencies(40, 30, 20, 10, 5));
        vector("reversed", frequencies(5, 10, 20, 30, 40));
        vector("tied-pair", frequencies(10, 10, 1, 1, 0));

        // The comparator subtracts two longs and casts to int, so a difference that is a multiple
        // of 2^32 reads as zero and the ordinal tie-break decides instead of the frequency.
        vector("difference-is-two-to-the-32", frequencies(0, 4294967296L, 0, 0, 0));
        vector("difference-is-two-to-the-32-plus-one", frequencies(0, 4294967297L, 0, 0, 0));

        // What a matrix decodes to, including the lower case rows the reading constructor fills.
        for (final byte[] matrix : new byte[][] {
                {0x1b, 0x1b, 0x1b, 0x1b, 0x1b},
                {0x00, 0x00, 0x00, 0x00, 0x00},
                {(byte) 0xe4, (byte) 0xe4, (byte) 0xe4, (byte) 0xe4, (byte) 0xe4},
                {0x1b, 0x27, 0x4e, (byte) 0x93, 0x1b}}) {
            decode(matrix);
        }

        // The two accessors, and where they disagree.
        final byte[] standard = {0x1b, 0x1b, 0x1b, 0x1b, 0x1b};
        final SubstitutionMatrix matrix = new SubstitutionMatrix(standard);
        for (final byte refBase : new byte[] {'A', 'C', 'G', 'T', 'N', 'a', 'c', 0, -1}) {
            for (final byte readBase : new byte[] {'C', 'G', 0}) {
                code(standard, matrix, refBase, readBase);
            }
            for (final byte code : new byte[] {0, 3}) {
                base(standard, matrix, refBase, code);
            }
        }

        // Escaped rather than trimmed: toString() ends in a tab and, for `n`, in four NULs, and
        // String.trim() would silently drop both because it strips everything at or below a space.
        System.out.printf("tostring\t%s\t%s%n", hex(standard), escape(matrix.toString()));
    }

    /** A frequency array indexed by base value, as `substitutionCodeVector` expects it. */
    static long[] frequencies(final long a, final long c, final long g, final long t, final long n)
            throws Exception {
        final long[] out = new long[symbolSpaceSize()];
        out['A'] = a;
        out['C'] = c;
        out['G'] = g;
        out['T'] = t;
        out['N'] = n;
        return out;
    }

    /**
     * `substitutionCodeVector` for each reference base, with the same frequencies each time. The
     * method is private and has a documented side effect on the code lookup, so both are recorded.
     */
    static void vector(final String label, final long[] frequencies) throws Exception {
        for (final byte refBase : BASES) {
            final SubstitutionMatrix matrix = new SubstitutionMatrix(new byte[] {0, 0, 0, 0, 0});
            final Method method = SubstitutionMatrix.class
                    .getDeclaredMethod("substitutionCodeVector", byte.class, long[].class);
            method.setAccessible(true);
            final byte codeVector = (byte) method.invoke(matrix, refBase, frequencies);

            final StringJoiner counts = new StringJoiner(",");
            final StringJoiner codes = new StringJoiner(",");
            for (final byte base : BASES) {
                counts.add(Long.toString(frequencies[base]));
                if (base != refBase) {
                    codes.add((char) base + "=" + codeByBase(matrix, refBase, base));
                }
            }
            System.out.printf("vector\t%s\t%c\t%s\t%d\t%s%n", label, (char) refBase, counts,
                    0xFF & codeVector, codes);
        }
    }

    /** What a serialized matrix decodes to, upper case and lower. */
    static void decode(final byte[] encoded) throws Exception {
        final SubstitutionMatrix matrix = new SubstitutionMatrix(encoded);
        for (final byte refBase : BASES) {
            final StringJoiner upper = new StringJoiner(",");
            final StringJoiner lower = new StringJoiner(",");
            for (byte code = 0; code < SubstitutionMatrix.BASES_SIZE - 1; code++) {
                upper.add(describe(baseByCode(matrix, refBase, code)));
                lower.add(describe(baseByCode(matrix,
                        (byte) Character.toLowerCase((char) refBase), code)));
            }
            System.out.printf("decode\t%s\t%c\t%s%n", hex(encoded), (char) refBase, upper);
            System.out.printf("lowercase\t%s\t%c\t%s%n", hex(encoded), (char) refBase, lower);
        }
    }

    static void code(final byte[] encoded, final SubstitutionMatrix matrix, final byte refBase,
            final byte readBase) {
        String answer;
        try {
            answer = Byte.toString(matrix.code(refBase, readBase));
        } catch (final Exception e) {
            answer = "ERR " + e.getClass().getName() + ": " + escape(e.getMessage());
        }
        System.out.printf("code\t%s\t%d\t%d\t%s%n", hex(encoded), refBase, readBase, answer);
    }

    static void base(final byte[] encoded, final SubstitutionMatrix matrix, final byte refBase,
            final byte code) {
        String answer;
        try {
            answer = describe(matrix.base(refBase, code));
        } catch (final Exception e) {
            answer = "ERR " + e.getClass().getName() + ": " + escape(e.getMessage());
        }
        System.out.printf("base\t%s\t%d\t%d\t%s%n", hex(encoded), refBase, code, answer);
    }

    /**
     * A message with every character outside printable ASCII written as \\uXXXX. Java formats an
     * invalid base with %c, so these messages carry a NUL or a U+FFFF, and a golden that held them
     * raw would be a text file with control bytes in it.
     */
    static String escape(final String message) {
        if (message == null) {
            return "<null>";
        }
        final StringBuilder out = new StringBuilder(message.length());
        for (int i = 0; i < message.length(); i++) {
            final char c = message.charAt(i);
            if (c >= 0x20 && c <= 0x7e) {
                out.append(c);
            } else {
                out.append(String.format("\\u%04x", (int) c));
            }
        }
        return out.toString();
    }

    static String describe(final byte base) {
        return base == 0 ? "0" : Character.toString((char) base);
    }

    static int symbolSpaceSize() throws Exception {
        final Field field = SubstitutionMatrix.class.getDeclaredField("SYMBOL_SPACE_SIZE");
        field.setAccessible(true);
        return field.getInt(null);
    }

    static int codeByBase(final SubstitutionMatrix matrix, final byte refBase, final byte readBase)
            throws Exception {
        final Field field = SubstitutionMatrix.class.getDeclaredField("codeByBase");
        field.setAccessible(true);
        return ((byte[][]) field.get(matrix))[refBase][readBase];
    }

    static byte baseByCode(final SubstitutionMatrix matrix, final byte refBase, final byte code)
            throws Exception {
        final Field field = SubstitutionMatrix.class.getDeclaredField("baseByCode");
        field.setAccessible(true);
        return ((byte[][]) field.get(matrix))[refBase][code];
    }

    static String hex(final byte[] bytes) {
        final StringBuilder b = new StringBuilder(bytes.length * 2);
        for (final byte value : bytes) {
            b.append(String.format("%02x", value));
        }
        return b.toString();
    }
}
