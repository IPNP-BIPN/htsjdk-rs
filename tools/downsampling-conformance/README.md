# ConstantMemoryDownsamplingIterator conformance harness

Dumps which records the constant-memory downsampler keeps, and the seen/accepted/discarded counts
it reports, over seven proportions and two seeds.

The proportions are chosen for the arithmetic: the threshold is
`Integer.MIN_VALUE + (int) Math.round(range * proportion)`, and at a proportion of 1 that `(int)`
cast wraps back to `Integer.MAX_VALUE`, which is the only value that keeps everything. A port doing
the same sum in 64 bits keeps everything at proportions where the reference does not.

## Run

```sh
python3 ../conformance/run_suite.py --suites downsampling
```
