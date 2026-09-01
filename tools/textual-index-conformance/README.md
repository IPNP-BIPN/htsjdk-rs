# TextualBAMIndexWriter conformance harness

Dumps each `.bai` twice: as bytes, and as the text `TextualBAMIndexWriter` prints for it. The BAMs
and their indexes are the ones `tools/build-index-conformance` builds, so the two suites measure
the same files from two sides: one that the index bytes are right, one that the text form of those
bytes is right.

The writer is package-private, so the harness reaches it the way GATK does, through
`BAMIndexer.createAndWriteIndex(input, output, true)`.

A blank line is part of the format (an empty bin prints one) and is emitted as `<blank>` so it
survives being split into lines on the other side.

## Run

```sh
python3 ../conformance/run_suite.py --suites textual-index
```
