# zlib conformance harness

Checks that the Rust deflate backend produces byte-identical output to the JDK's
`java.util.zip.Deflater` in `nowrap` mode, which is what htsjdk's
`BlockCompressedOutputStream` uses for BGZF blocks.

See `../../docs/decisions/0001-deflate-backend.md` for the result and the decision.

## Run

The comparison CI runs is the `zlib` suite, which regenerates the Java side in the pinned
container and compares its 70 rows against the table in
`crates/htsjdk-bgzf/tests/zlib_conformance.rs`:

```sh
python3 ../conformance/run_suite.py --suites zlib
```

By hand, both halves print the same lines and `diff` is the whole comparison:

```sh
# Reference side, in the pinned container (Z2: 7 payloads x 10 levels)
docker run --rm --platform linux/amd64 -v "$PWD":/harness:ro -w /work htsjdk-rs-oracle:4.2.0 \
  'cp /harness/Z2.java . && javac -d . Z2.java && java -cp . Z2' > /tmp/z_java.txt

# Rust side (zlib backend, as pinned in the decision record)
(cd rust && cargo run --release) > /tmp/z_rust.txt
diff /tmp/z_java.txt /tmp/z_rust.txt
```

`Z.java` is the original four-vector smoke test and is kept because decision 0001's first
experiment is written against it.

## Status

**Re-derived on every push** by the `zlib` suite, in the pinned `linux/amd64` container: 70 of 70
vectors, keyed on payload and level. It was previously confirmed there once, by hand, which is
what decision 0001's third follow-up records.
