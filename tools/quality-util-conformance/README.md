# QualityUtil conformance harness

Dumps `htsjdk.samtools.util.QualityUtil`'s answers as bit patterns: the 101-entry error-probability
table, `getPhredScoreFromErrorProbability` over a sweep, and `getPhredScoreFromObsAndErrors` over
pairs a metrics file actually produces.

The sweep is chosen for the places a port goes wrong rather than for coverage: a probability above
one makes the argument to `Math.round` negative, where Java rounds half **up** and Rust's
`f64::round` rounds half **away from zero**; zero makes it infinite, where `(int) Math.round` wraps
rather than saturating; and the table itself is built with `Math.pow`, which decision 0007 deferred.

## Run

```sh
python3 ../conformance/run_suite.py --suites quality-util
```
