# igzip oracle

ISA-L 2.31.0 built for x86-64. It writes raw deflate exactly as GKL configures
`isal_deflate_stateless` at levels 1 and 2. `crates/gkl-deflate/src/igzip.rs` was checked against it
(decision 0050). On Apple Silicon it runs under Rosetta, which reports SSE4.2 and AVX2, so ISA-L takes
the kernel path GKL takes.

```sh
docker build -t igzip-oracle tools/igzip-oracle
docker run --rm -v "$PWD/corpus":/c igzip-oracle sh -c 'igzip-oracle /c/*'   # writes <file>.l1, <file>.l2
cargo run -p gkl-deflate --release --example igzip_corpus -- corpus/*        # writes <file>.r1, <file>.r2
for f in corpus/*.l1; do cmp -s "$f" "${f%.l1}.r1" || echo "$f"; done
```
