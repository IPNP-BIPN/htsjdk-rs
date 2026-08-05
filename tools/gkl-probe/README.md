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

## `emulated.txt`

The committed reference column, produced on Apple Silicon where Docker translates `linux/amd64`
through Rosetta. Its `cpu` line reads `VirtualApple @ 2.50GHz`, and Rosetta implements no AVX, so
this column is igzip's SSE path.

It is a **reference, not a golden**: nothing in the Rust tests reads it. Its only consumer is the
`igzip-portability` CI job, which reruns the probe on a real x86-64 host and diffs. Question 2 is
that diff.

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
