# htsjdk.samtools.filter conformance harness

Dumps `AlignedFilter`, `ReadNameFilter` and `TagFilter`'s answers over a corpus built by index, in
both forms of `filterOut`: the single-record one and the pair one. The pair form is not the single
form applied twice, and each filter is asymmetric differently, so a port can agree on every
single-record row and still be wrong.

`true` in the dump means the record is **dropped**.

## Run

```sh
python3 ../conformance/run_suite.py --suites filter
```
