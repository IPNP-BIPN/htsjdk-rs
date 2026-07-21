# Java math conformance

Generates the corpus that `crates/jmath` is validated against.

## Run

```sh
docker run --rm --platform linux/amd64 \
  -v "$PWD":/harness:ro -v <jars>:/j:ro -v <out>:/out -w /work htsjdk-rs-oracle:4.2.0 \
  'cp /harness/JMathDump.java . && javac -cp /j/commons-math3-3.6.1.jar -d . JMathDump.java \
   && java -cp .:/j/commons-math3-3.6.1.jar JMathDump /out/jmath.csv'
gzip -9 jmath.csv
```

Requires `commons-math3-3.6.1.jar` (sha256
`1e56d7b058d28b65abd256b8458e3885b674c1d588fa43cd7d1cbb9c7ef2b308`) for `FastMath`.

## Format

`function,input_bits,math_bits,strictmath_bits,fastmath_bits`, all hex raw bit patterns.

Bits and not decimal, deliberately: decimal rendering discards exactly the last-place
difference the corpus exists to measure.

## Why three columns of output

`java.lang.Math`, `java.lang.StrictMath` and commons-math3 `FastMath` are three different
functions that disagree pairwise. GATK calls all three depending on the call site. See
`../../docs/decisions/0005-java-math-has-three-implementations.md`.

## Sampling

Stratified rather than uniform, because the failures live at the edges: special values, the
`(0,1]` probability range that `log10` sees in every genotype likelihood, phred-scale
magnitudes for `exp`, a dense sweep around 1.0, subnormals, and 40,000 pseudo-random bit
patterns. 44,996 inputs, 809,930 rows.
