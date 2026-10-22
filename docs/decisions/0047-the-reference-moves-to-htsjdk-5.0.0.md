# 0047. The reference moves to htsjdk 5.0.0

**Status:** accepted
**Date:** 2026-10-22
**Supersedes:** the pin of [0042](0042-the-reference-version-is-pinned-by-the-consumers-not-by-this-repository.md)
**Resolves:** [#182](https://github.com/IPNP-BIPN/htsjdk-rs/issues/182)
**Result:** 87 of 87 suites oracle-backed against htsjdk 5.0.0

## The decision

The reference is htsjdk `5.0.0`, the version GATK 4.7.0.0 and Picard 3.5.0 pin. Decision 0042 held
the pin at 4.2.0 until all three repositories could move together. picard-rs moved to Picard 3.5.0
on its own (its decision 0013), so the three repositories already name two sets of pins, and
picard-rs reaches htsjdk-rs for the behaviours 5.0.0 changed, CRAM 3.1 above all. gatk-rs stays at
GATK 4.6.2.0: it pins htsjdk-rs by revision, so it keeps 4.2.0's behaviour until it moves that pin,
and its own goldens decide when.

## How it was measured

The same way picard-rs measured its move: the oracle image moved first, alone, so the first CI
run compared htsjdk 5.0.0 with every 4.2.0 golden before any port code changed.

The image takes htsjdk 5.0.0 and its runtime dependencies at its own `build.gradle`'s versions,
each pinned by sha256: jlibdeflate 0.1.0, snappy-java 1.1.10.5, xz 1.9, json 20231013 and
commons-logging 1.3.0. 5.0.0 deflates with libdeflate by default, which `OracleProbe` refused to
build the image under; `JAVA_TOOL_OPTIONS=-Dsamjdk.use_libdeflate=false` turns it off for every
JVM in the image, so the deflater a golden records is still the one its dump chose.
`Deflater::Libdeflate` stays the port of what 5.0.0 does when it is on.

That run failed 18 of 87 suites, the same 18 on this laptop and on the x86-64 runner. Seven of
them were dumps that no longer compiled against 5.0.0's API; once ported to it, the two rANS 4x8
suites and the core, Huffman and Golomb codec suites came back identical to their 4.2.0 goldens.

## What moved, and why

Thirteen goldens were re-frozen from the runner's candidates, each byte-identical to the local
dump:

| golden | the difference | the cause |
|---|---|---|
| `bam_tag_read` | an unterminated string tag is `SAMFormatException: Null-terminated string tag value is not null terminated.` | `readNullTerminatedString` scans the array instead of running off the buffer |
| `cram_encoding_map`, `cram_compression_header`, `cram_record_read`, `cram_record_write`, `cram_container`, `cram_block`, `cram_preservation_map`, `cram_tag_encoding_map` | every data series is on another content id | `DataSeries` takes htslib's `cram_DS_ID` numbering, 11 to 46 |
| `cram_container`, `cram_block` | the files are CRAM 3.1, with rANS Nx16, fqzcomp and the name tokeniser | 5.0.0 writes 3.1 by default |
| `cram_negotiation` | a tag's compressor is GZIP or rANS Nx16, and every `Z` is tab-stopped | a per-tag `TrialCompressor`; `Z` no longer takes the fixed-size encoding |
| `cram_external_codecs` | a short `readFully` is `IllegalStateException`, a missing stop byte a `CRAMException` | `CRAMByteReader` replaces the `InputStream`s |
| `cram_record_flags` | a pair at 100 and 200 has TLEN 150, not 101 | `computeInsertSize` is htslib's span of both ends |
| `cram_slice_header` | a slice header has no tags | 5.0.0 no longer writes `BD`, `SD`, `B1`, `S1`, `B5`, `S5` |

## What the suites do not reach

A suite measures what its dump calls. To find what 5.0.0 changed that no suite calls, every Java
class the port cites was compared between the two tags with comments, imports, modifier order and
line wrapping removed. Of 142 classes, 112 are identical, 2 moved package (the rANS encoders), and
28 changed. Each change was read:

* **behaviour, now ported:** the unterminated string tag above, and `calculateMdAndNm`, which stops
  a deletion at the end of the reference where 4.2.0 indexed past it and threw. The port already
  stopped there.
* **behaviour, not in this repository:** `SamPairUtil.getPairOrientation` and `SAMRecord.toString`,
  which picard-rs ports beside the Picard code that reaches them (its decision 0013).
* **performance, same output:** the nibble lookup in `SAMUtils`, the thread-local MD5 in
  `SequenceUtil`, the tag list append in `BinaryTagCodec`, the header in `BAMRecordCodec`.
* **CRAM 3.1:** the remaining CRAM classes, whose 5.0.0 behaviour `reader.rs` and `writer.rs` were
  written from (decision 0039's extension surface), and which picard-rs measures byte for byte
  against Picard 3.5.0's CRAM output.

## What this does not change

Every module still says which htsjdk it was transcribed from, and most say 4.2.0: that is where the
code came from, and the suites are what say it matches 5.0.0. A module that was changed for 5.0.0
says so.
