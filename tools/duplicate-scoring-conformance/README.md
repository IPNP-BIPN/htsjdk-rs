# DuplicateScoringStrategy conformance harness

Dumps every score and every comparison the three strategies give, over a corpus built by index and
chosen for the arithmetic: qualities on both sides of the Q15 threshold, a read long enough to
reach the `Short.MAX_VALUE / 2` clamp, vendor-failed records that take the `Short.MIN_VALUE / 2`
discount, and paired and unpaired ends so `compare`'s first branch and its name tie-break both fire.

`assumeMateCigar` is false throughout: its true branch reads `MC` through `SAMUtils.getMateCigar`,
which throws when there is none.

## Run

```sh
python3 ../conformance/run_suite.py --suites duplicate-scoring
```
