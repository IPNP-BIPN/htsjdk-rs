import htsjdk.samtools.SamFileValidator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * What `SamFileValidator` says about each SAM in `cases/`, one line per error.
 *
 *   error <case> <SAMValidationError.toString>
 *   clean <case>
 *
 * The corpus is read from files rather than embedded, so this harness and the Rust test measure the
 * same bytes rather than two transcriptions of them. Each case is one file broken in one way: a
 * version the header may not carry, no read groups at all, a missing and an invalid platform, a
 * record with no read group and a record with no qualities, records out of coordinate order, a pair
 * whose mate fields disagree, a paired read whose mate never arrives, and flags a single-end or
 * unmapped read may not set.
 *
 * VERBOSE is what is dumped, because it is the mode that prints the errors themselves; the SUMMARY
 * histogram is Picard's rendering of the same list and belongs to the tool.
 */
public class ValidatorDump {
  public static void main(String[] args) throws Exception {
    // The harness directory is mounted at /harness, and only <class>.java is copied into /work.
    Path dir = Paths.get(args.length > 0 ? args[0] : "/harness/cases");
    List<Path> cases = new ArrayList<>();
    Files.list(dir).filter(p -> p.toString().endsWith(".sam")).forEach(cases::add);
    Collections.sort(cases);

    for (Path path : cases) {
      String name = path.getFileName().toString().replace(".sam", "");
      SamReaderFactory factory =
          SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
      try (SamReader reader = factory.open(path.toFile())) {
        StringWriter captured = new StringWriter();
        PrintWriter out = new PrintWriter(captured);
        SamFileValidator v = new SamFileValidator(out, 8000);
        v.setVerbose(true, Integer.MAX_VALUE);
        v.validateSamFileVerbose(reader, null);
        out.flush();
        boolean any = false;
        for (String line : captured.toString().split("\n")) {
          if (line.isEmpty() || line.startsWith("No errors found")) {
            continue;
          }
          System.out.printf("error\t%s\t%s%n", name, line);
          any = true;
        }
        if (!any) {
          System.out.printf("clean\t%s%n", name);
        }
      }
    }
  }
}
