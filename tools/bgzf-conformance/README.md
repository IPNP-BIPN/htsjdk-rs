# BGZF conformance harness

Generates golden vectors from the official htsjdk jar for
`crates/htsjdk-bgzf/tests/bgzf_conformance.rs`.

## Run

```sh
curl -sSLO https://repo1.maven.org/maven2/com/github/samtools/htsjdk/4.2.0/htsjdk-4.2.0.jar
# sha256 52c9eb1a568d8261767ebf888d6ebafa60911bd44e0c9242d413fbff1d1e2398
javac -cp htsjdk-4.2.0.jar -d . B.java && java -cp .:htsjdk-4.2.0.jar B
```

Paste the emitted rows into the `GOLDEN` table.

## Corpus design

The payloads are chosen to sit on boundaries, not to be representative:

| Case | Purpose |
|---|---|
| `empty` | output must be the 28-byte terminator alone |
| `tiny` | single short block |
| `exact1` | exactly one full uncompressed block (65498) |
| `over1` | one full block plus one byte |
| `multi` | several blocks |
| `incompr` | incompressible, forces the no-compression fallback |
| `big` | many blocks, highly compressible |

`incompr` is the one that matters most: it caught the boundary bug recorded in
`../../docs/decisions/0003-deflate-fallback-is-a-status-not-a-length.md`, which the other six
cases were happy to miss.
