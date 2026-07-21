# zlib conformance harness

Checks that the Rust deflate backend produces byte-identical output to the JDK's
`java.util.zip.Deflater` in `nowrap` mode, which is what htsjdk's
`BlockCompressedOutputStream` uses for BGZF blocks.

See `../../docs/decisions/0001-deflate-backend.md` for the result and the decision.

## Run

```sh
# Reference side
javac -d . Z.java && java -cp . Z

# Rust side (zlib backend, as pinned in the decision record)
cd rust && cargo run --release
```

Every `level=N ... md5=` line must match between the two.

## Status

Validated locally against OpenJDK 17.0.19 on macOS arm64.
**Not yet confirmed inside the pinned `linux/amd64` container**, which is the authoritative
oracle. See the "remaining work" section of the decision record.
