# Oracle

The environment every golden in this repository is produced under, and the assertion that
protects it.

## Why the assertion exists

Intel GKL does not throw when its native library fails to load. It logs a `WARNING` and falls
back to the JDK deflater. An oracle without an explicit check will therefore emit goldens from a
degraded configuration that are indistinguishable from good ones.

This is not hypothetical: it happened while writing
[decision 0004](../../docs/decisions/0004-oracle-platform.md). A missing `commons-io` on the
classpath produced `usingIntelDeflater = false`, which would have "confirmed" a wrong assumption
in the program plan. `OracleProbe` exits non-zero rather than logging, and the Docker build runs
it, so a degraded image cannot be built at all.

## Contract

| Property | Required |
|---|---|
| `os.arch` | `amd64` |
| Java major | 17 |
| `usingIntelDeflater` / `usingIntelInflater` | true |
| CPU flags | `avx`, `avx2`, `sse4_2` |

Everything is pinned: the base image by digest, every jar by sha256.

## Use

```sh
docker build --platform linux/amd64 -t htsjdk-rs-oracle:4.2.0 .
./run.sh ../bgzf-conformance/B.java goldens/bgzf.txt
```

`run.sh` checks the contract first and refuses to run the harness if it fails. Provenance
(image id, platform, JVM, GKL state, CPU flags, timestamp) goes to stderr; harness output goes
to stdout or the named file.

## Note on the platform

`linux/amd64` runs under Rosetta translation on Apple Silicon and that is fine: AVX and AVX2 are
translated and GKL's x86 natives load. AVX-512 is **not** available, so any work depending on it
(PairHMM implementation selection) must pin its choice explicitly rather than rely on
`FASTEST_AVAILABLE`.

## Not yet here

GATK 4.6.2.0 and Picard 3.4.0 themselves. This image covers the htsjdk layer, which is what
`htsjdk-rs` needs. The tool-level oracle extends this image and belongs with `picard-rs` and
`gatk-rs`.
