# TabixIndexCreator conformance harness

Dumps the `.tbi` `htsjdk.tribble.index.tabix.TabixIndexCreator` builds for a stream of features,
twice: the little-endian body the port composes, and the same body inside the BGZF stream that
lands beside the feature file.

The creator is fed features and their file positions **directly** rather than through a codec. What
is measured is the index; a VCF reader in front of it would measure the reader too, and the reader
is measured by its own suites.

The cases are chosen for where the layout is decided rather than described: a feature is indexed
one feature late, because a chunk needs both ends and the end of one feature is the start of the
next; the bin comes from a zero-based half-open region, because `getIndexingBin` returns null; a
feature with no end is one base for the bin **and** a shifted window for the linear index; the
linear index fills its gaps with the last non-empty offset; and the name block counts its null
terminators.

Four refusals are here too, all of them `IllegalArgumentException`, including the one that says
equal starts are **not** out of order because `compareTo` never looks at the end.

## Run

```sh
python3 ../conformance/run_suite.py --suites tabix-index
```
