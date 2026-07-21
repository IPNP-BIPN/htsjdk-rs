/*
 * Is htsjdk's metrics number formatting locale-dependent?
 *
 * FormatUtil builds its formatters from NumberFormat.getNumberInstance(), which takes the
 * default locale. If that is so, then the bytes of every Picard metrics file depend on the
 * JVM's locale, and the oracle contract has to pin it.
 *
 * This prints the same values under several locales. It asserts nothing; the answer is
 * whatever comes out.
 */

import htsjdk.samtools.util.FormatUtil;
import java.util.Locale;

public class FormatProbe {

    static final double[] VALUES = {
            0.0, 1.0, -1.0, 0.5, 1.0 / 3.0, 2.0 / 3.0, 1234567.891234567,
            0.0000001, 0.5000005, 0.4999995, 1.5, 2.5, -1.5, -2.5,
            Double.NaN, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY,
            Double.MIN_VALUE, Double.MAX_VALUE};

    public static void main(final String[] args) {
        System.out.println("default_locale\t" + Locale.getDefault());
        for (final String tag : new String[]{"en-US", "fr-FR", "de-DE", "ar-EG", "hi-IN"}) {
            Locale.setDefault(Locale.forLanguageTag(tag));
            final FormatUtil f = new FormatUtil();
            final StringBuilder sb = new StringBuilder();
            for (final double v : VALUES) {
                sb.append(f.format(v)).append('|');
            }
            sb.append("INT:").append(f.format(1234567L)).append('|');
            sb.append("BOOL:").append(f.format(true)).append(f.format(false));
            System.out.println(tag + "\t" + sb);
        }
    }
}
