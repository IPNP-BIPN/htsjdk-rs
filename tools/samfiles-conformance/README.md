# SamFiles.findIndex conformance harness

Dumps the index `htsjdk.samtools.SamFiles.findIndex` picks for a directory laid out around one data
file. Each row carries the layout as well as the answer, so the port rebuilds the same directory
from the dump rather than keeping a second copy of the case table.

The cases are chosen for the four places the search is not what its one-line summary suggests:

* the extension is **replaced** before it is **appended**, so a `.csi` beside `reads.bam` beats a
  `reads.bam.bai` beside the same file;
* a CRAM has no replaced `.csi` and no replaced `.bai`: its replacement is `.crai` alone, followed
  by an appended `.crai` before the shared fallthrough;
* `Files.isRegularFile` refuses a directory, so a directory named `reads.bai` is stepped over
  rather than returned;
* finding nothing resolves the path's symbolic links and searches again at the real location, so
  an index beside the link wins and an index beside the target is still found.

## Run

```sh
python3 ../conformance/run_suite.py --suites samfiles
```
