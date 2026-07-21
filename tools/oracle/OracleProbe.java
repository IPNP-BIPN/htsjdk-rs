import com.intel.gkl.compression.IntelDeflaterFactory;
import com.intel.gkl.compression.IntelInflaterFactory;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.TreeSet;

/**
 * Asserts and records the oracle environment before any golden is produced.
 *
 * <p>The point of this class is that it <em>fails</em>. Intel GKL logs a WARNING and falls back
 * to the JDK deflater when its native library cannot be extracted or loaded; it does not throw.
 * An oracle without this check will emit goldens from a degraded configuration that look
 * completely normal. See docs/decisions/0004-oracle-platform.md, where exactly that happened
 * because commons-io was missing from a classpath.
 *
 * <p>Exit code 0 means the environment matches the declared contract. Any other exit code means
 * no golden produced in this environment may be trusted.
 */
public final class OracleProbe {

    /** The contract every golden in this project is produced under. */
    private static final String EXPECTED_ARCH = "amd64";
    private static final String EXPECTED_JAVA_MAJOR = "17";
    private static final boolean EXPECTED_GKL_AVAILABLE = true;
    private static final List<String> REQUIRED_CPU_FLAGS = Arrays.asList("avx", "avx2", "sse4_2");

    public static void main(String[] args) throws Exception {
        final List<String> failures = new ArrayList<>();

        final String arch = System.getProperty("os.arch");
        final String javaVersion = System.getProperty("java.version");
        final String javaVendor = System.getProperty("java.vendor");
        final String javaMajor = javaVersion.split("\\.")[0];

        final boolean intelDeflater = new IntelDeflaterFactory().usingIntelDeflater();
        final boolean intelInflater = new IntelInflaterFactory().usingIntelInflater();
        final TreeSet<String> cpuFlags = readCpuFlags();

        if (!EXPECTED_ARCH.equals(arch)) {
            failures.add("os.arch is '" + arch + "', expected '" + EXPECTED_ARCH + "'");
        }
        if (!EXPECTED_JAVA_MAJOR.equals(javaMajor)) {
            failures.add("java major is '" + javaMajor + "', expected '" + EXPECTED_JAVA_MAJOR + "'");
        }
        if (intelDeflater != EXPECTED_GKL_AVAILABLE) {
            failures.add("usingIntelDeflater is " + intelDeflater + ", expected "
                    + EXPECTED_GKL_AVAILABLE
                    + ". GKL degrades silently; check commons-io and commons-logging are on the"
                    + " classpath and that the native library extracted.");
        }
        if (intelInflater != EXPECTED_GKL_AVAILABLE) {
            failures.add("usingIntelInflater is " + intelInflater + ", expected "
                    + EXPECTED_GKL_AVAILABLE);
        }
        for (final String flag : REQUIRED_CPU_FLAGS) {
            if (!cpuFlags.isEmpty() && !cpuFlags.contains(flag)) {
                failures.add("CPU flag '" + flag + "' is absent");
            }
        }

        // Provenance record, emitted on stdout so a caller can capture it verbatim.
        final StringBuilder json = new StringBuilder();
        json.append("{\n");
        json.append("  \"os_arch\": \"").append(arch).append("\",\n");
        json.append("  \"java_version\": \"").append(javaVersion).append("\",\n");
        json.append("  \"java_vendor\": \"").append(javaVendor).append("\",\n");
        json.append("  \"using_intel_deflater\": ").append(intelDeflater).append(",\n");
        json.append("  \"using_intel_inflater\": ").append(intelInflater).append(",\n");
        json.append("  \"cpu_flags_checked\": \"").append(String.join(",", REQUIRED_CPU_FLAGS)).append("\",\n");
        json.append("  \"avx\": ").append(cpuFlags.contains("avx")).append(",\n");
        json.append("  \"avx2\": ").append(cpuFlags.contains("avx2")).append(",\n");
        json.append("  \"avx512f\": ").append(cpuFlags.contains("avx512f")).append(",\n");
        json.append("  \"contract_satisfied\": ").append(failures.isEmpty()).append("\n");
        json.append("}");
        System.out.println(json);

        if (!failures.isEmpty()) {
            System.err.println("ORACLE CONTRACT VIOLATED. No golden produced here may be trusted.");
            for (final String f : failures) {
                System.err.println("  - " + f);
            }
            System.exit(2);
        }
    }

    /** Returns the flags of the first CPU, or an empty set where /proc is unavailable. */
    private static TreeSet<String> readCpuFlags() {
        final TreeSet<String> flags = new TreeSet<>();
        try (BufferedReader r = new BufferedReader(new FileReader("/proc/cpuinfo"))) {
            String line;
            while ((line = r.readLine()) != null) {
                if (line.startsWith("flags")) {
                    final int colon = line.indexOf(':');
                    if (colon >= 0) {
                        flags.addAll(Arrays.asList(line.substring(colon + 1).trim().split("\\s+")));
                    }
                    break;
                }
            }
        } catch (final Exception ignored) {
            // Non-Linux or restricted /proc: the flag assertions are then skipped, and the
            // recorded provenance says so by reporting every flag as false.
        }
        return flags;
    }
}
