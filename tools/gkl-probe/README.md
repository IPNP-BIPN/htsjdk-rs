# GKL probe

Two questions about the deflater that writes GATK's BAM files, both answered by compressing the
same bytes and hashing the result.

1. **Which deflater is it?** `libgkl_compression.so` bundles ISA-L's igzip *and* zlib 1.2.13.
   Levels 1 to 6 differ from `java.util.zip.Deflater` and are igzip; levels 7 to 9 are
   byte-identical to it and are zlib. htsjdk's BGZF default is level 5, so the default path is
   igzip.
2. **Is igzip's output a property of the algorithm or of the CPU?** igzip dispatches on CPU
   features and ships AVX2 and AVX512 kernels. If those emit different bytes from the base path,
   there is no single byte sequence for a port to target.

See `../../docs/decisions/0028-gkl-is-igzip-below-seven-and-zlib-above.md`.

## Run

```sh
docker build --platform linux/amd64 -t htsjdk-rs-oracle:4.2.0 ../oracle
docker run --rm --platform linux/amd64 -v "$PWD":/harness:ro -w /work htsjdk-rs-oracle:4.2.0 \
  'cp /harness/GklProbe.java . && javac -cp "$ORACLE_CP" -d . GklProbe.java \
   && java -cp "$ORACLE_CP":. GklProbe'
```

## `real-x86-64.txt`

The committed column, and a golden: `crates/gkl-deflate` compares against it. Its `cpu` line reads
`AMD EPYC 7763 64-Core Processor` and its `flags` line `avx,avx2,pclmulqdq,sse4_1,sse4_2`, because
it was derived by the `igzip-portability` job on a GitHub runner and downloaded from that run's
`gkl-probe-real-x86-64` artefact. Decision 0008: a golden comes from the pinned container on a real
x86-64 runner and never from a developer machine.

It did not start out that way, and the way it went wrong is worth keeping. The file was first
written on Apple Silicon under Rosetta as a *reference* for the portability diff, and this README
said so: "nothing in the Rust tests reads it." Then `gkl-deflate` started reading it, which made it
a golden and made that sentence false in the same commit. Every measured row was identical between
the two columns, so nothing was wrong with the numbers; what was wrong was a claim about them that
had quietly stopped being true.

The job still re-derives every row on each run and fails on any difference, which is what the
portability question needs. What changed is that the file it compares against now comes from the
same class of host it is compared on.

## Which backend a level reaches

The probe measures *that* a level differs from the JDK, not *why*. The why is in the library, and
decision 0029 corrects the inference 0028 drew from the numbers alone. Reproduce it with:

```sh
unzip -o gkl.jar -d jar
S=jar/com/intel/gkl/native/libgkl_compression.so
# the level branch: (level - 1) <= 1 goes to isal_deflate_stateless_init, else deflateInit2_
llvm-objdump -d --start-address=0x53e0 --stop-address=0x5628 "$S"
# Intel's zlib patches, absent from stock zlib
llvm-nm "$S" | grep -E 'deflate_medium|slide_hash_sse|longest_match'
strings "$S" | grep -m1 'deflate 1\.'
```
