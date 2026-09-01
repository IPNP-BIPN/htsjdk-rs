# SBIIndexWriter conformance harness

Dumps the bytes `SBIIndexWriter` writes, as a length and an MD5, plus the two header counters a
reader takes back out: the record count and the offset count.

They are different numbers and neither is the other. An offset is written every `granularity`
records **starting with the first**, and `finish` writes the final offset as well, so an index over
`n` records holds `ceil(n / granularity) + 1` offsets while the header records `n`. The corpus
varies granularity against record counts that do and do not divide by it, and includes an empty
index, which still carries the final offset.

## Run

```sh
python3 ../conformance/run_suite.py --suites sbi
```
