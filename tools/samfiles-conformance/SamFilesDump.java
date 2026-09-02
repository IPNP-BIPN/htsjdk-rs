import htsjdk.samtools.SamFiles;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

/**
 * Dumps what `htsjdk.samtools.SamFiles.findIndex` answers for a directory laid out around one
 * data file.
 *
 * Each case is a layout and a query, and the layout travels with the answer so the port rebuilds
 * the same directory rather than keeping a second copy of the table:
 *
 *   case  <id>  <query>  <entry;entry;...>  <answer relative to the root, or null>
 *
 * An entry is `f:<name>` for a regular file, `d:<name>` for a directory, and
 * `l:<name>-&gt;<target>` for a symbolic link, where a name may contain `/`.
 *
 * The cases are chosen for the four places the search is not what its summary suggests: the
 * extension is REPLACED before it is APPENDED, so a `.csi` beside a `.bam` beats a `.bam.bai`
 * beside the same file; a CRAM has no replaced `.csi` at all; `Files.isRegularFile` refuses a
 * directory, so a directory named like an index is skipped rather than returned; and a failure to
 * find anything resolves the path's symbolic links and searches again at the real location.
 */
public class SamFilesDump {

  /** One layout, its query and the answer, printed as a row. */
  static void run(String id, String query, String... entries) throws IOException {
    Path root = Files.createTempDirectory("samfiles");
    for (String entry : entries) {
      String kind = entry.substring(0, 1);
      String body = entry.substring(2);
      if (kind.equals("f")) {
        Path file = root.resolve(body);
        Files.createDirectories(file.getParent());
        Files.write(file, new byte[] {0});
      } else if (kind.equals("d")) {
        Files.createDirectories(root.resolve(body));
      } else if (kind.equals("l")) {
        int arrow = body.indexOf("->");
        Path link = root.resolve(body.substring(0, arrow));
        Files.createDirectories(link.getParent());
        Files.createSymbolicLink(link, root.resolve(body.substring(arrow + 2)));
      } else {
        throw new IllegalArgumentException(entry);
      }
    }
    Path answer = SamFiles.findIndex(root.resolve(query));
    System.out.printf("case\t%s\t%s\t%s\t%s%n", id, query, String.join(";", entries), rel(root, answer));
  }

  /** The answer as a name under the root, so the row says nothing about where the root was. */
  static String rel(Path root, Path answer) throws IOException {
    if (answer == null) {
      return "null";
    }
    // The symlink retry answers a canonical path, so both sides are canonicalized before they are
    // compared: on a machine whose temporary directory is itself a link, the raw paths differ by a
    // prefix that has nothing to do with the search.
    return root.toRealPath().relativize(answer.toAbsolutePath().toRealPath()).toString();
  }

  public static void main(String[] args) throws IOException {
    // A BAM: the extension is replaced first, and `.bai` before `.csi` within that.
    run("bam-replaced-bai", "reads.bam", "f:reads.bam", "f:reads.bai");
    run("bam-appended-bai", "reads.bam", "f:reads.bam", "f:reads.bam.bai");
    run("bam-replaced-beats-appended", "reads.bam", "f:reads.bam", "f:reads.bai", "f:reads.bam.bai");
    run("bam-replaced-csi", "reads.bam", "f:reads.bam", "f:reads.csi");
    run("bam-replaced-bai-beats-replaced-csi", "reads.bam", "f:reads.bam", "f:reads.bai", "f:reads.csi");
    // The one a summary of the method gets wrong: the replaced pair is exhausted before anything
    // is appended, so a `.csi` beside the BAM wins over a `.bam.bai` beside it.
    run("bam-replaced-csi-beats-appended-bai", "reads.bam", "f:reads.bam", "f:reads.csi", "f:reads.bam.bai");
    run("bam-appended-csi", "reads.bam", "f:reads.bam", "f:reads.bam.csi");
    run("bam-appended-bai-beats-appended-csi", "reads.bam", "f:reads.bam", "f:reads.bam.bai", "f:reads.bam.csi");
    run("bam-nothing", "reads.bam", "f:reads.bam");
    // `Files.isRegularFile` is false for a directory, so the search steps over it.
    run("bam-index-is-a-directory", "reads.bam", "f:reads.bam", "d:reads.bai");
    run("bam-directory-then-appended", "reads.bam", "f:reads.bam", "d:reads.bai", "f:reads.bam.bai");
    // The data file is never opened, so it need not exist for its name to imply an index.
    run("bam-data-file-absent", "reads.bam", "f:reads.bai");
    // The extension test is case-sensitive: `.BAM` takes the fallthrough only.
    run("bam-uppercase-extension", "reads.BAM", "f:reads.BAM", "f:reads.bai");
    run("bam-uppercase-appended", "reads.BAM", "f:reads.BAM", "f:reads.BAM.bai");

    // A CRAM: `.crai` replaced, then `.crai` APPENDED, and only then the shared fallthrough.
    run("cram-replaced-crai", "reads.cram", "f:reads.cram", "f:reads.crai");
    run("cram-appended-crai", "reads.cram", "f:reads.cram", "f:reads.cram.crai");
    run("cram-replaced-beats-appended-crai", "reads.cram", "f:reads.cram", "f:reads.crai", "f:reads.cram.crai");
    run("cram-appended-bai", "reads.cram", "f:reads.cram", "f:reads.cram.bai");
    run("cram-appended-csi", "reads.cram", "f:reads.cram", "f:reads.cram.csi");
    run("cram-appended-crai-beats-appended-bai", "reads.cram", "f:reads.cram", "f:reads.cram.crai", "f:reads.cram.bai");
    // There is no replaced `.csi` outside the BAM branch, so this one is not found at all.
    run("cram-no-replaced-csi", "reads.cram", "f:reads.cram", "f:reads.csi");
    run("cram-no-replaced-bai", "reads.cram", "f:reads.cram", "f:reads.bai");

    // Every other name reaches the fallthrough and nothing else.
    run("sam-appended-bai", "reads.sam", "f:reads.sam", "f:reads.sam.bai");
    run("sam-no-replaced-bai", "reads.sam", "f:reads.sam", "f:reads.bai");
    run("sam-appended-csi", "reads.sam", "f:reads.sam", "f:reads.sam.csi");
    run("no-extension-appended", "reads", "f:reads", "f:reads.bai");

    // The symlink retry, which runs only when the first search found nothing.
    run("symlink-index-beside-target", "link.bam", "f:real/reads.bam", "f:real/reads.bai", "l:link.bam->real/reads.bam");
    run("symlink-index-beside-link", "link.bam", "f:real/reads.bam", "f:link.bai", "l:link.bam->real/reads.bam");
    run("symlink-link-beats-target", "link.bam", "f:real/reads.bam", "f:real/reads.bai", "f:link.bai", "l:link.bam->real/reads.bam");
    run("symlink-appended-beside-target", "link.bam", "f:real/reads.bam", "f:real/reads.bam.bai", "l:link.bam->real/reads.bam");
    // A link to nothing cannot be resolved, and the exception is caught rather than thrown.
    run("symlink-broken", "link.bam", "l:link.bam->missing/reads.bam");
  }
}
