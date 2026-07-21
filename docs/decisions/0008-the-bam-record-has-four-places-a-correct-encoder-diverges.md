# 0008. The BAM record has four places where a correct encoder still diverges

**Status:** accepted, measured against htsjdk
**Date:** 2026-07-21
**Follows:** [0003](0003-deflate-fallback-is-a-status-not-a-length.md)

## Finding

The BAM record layout is fully specified and short enough to implement from the specification
in an afternoon. Doing so produces a file that `samtools view` reads, that carries every read
correctly, and that is **not** htsjdk's file.

Four fields carry the divergence. Each was confirmed against htsjdk 4.2.0 in the pinned oracle
container, and each is covered by a golden in `crates/htsjdk-bam/tests/data/bam_codec.txt.gz`.

### 1. Integer tag width is chosen from the value, by a ladder that is not monotone

`BinaryTagCodec.getIntegerType` picks the narrowest representation, testing in a fixed order.
Read as ranges:

| value | type | width |
|---|---|---|
| `-2^31 .. -32769` | `i` | 4 |
| `-32768 .. -129` | `s` | 2 |
| `-128 .. 127` | `c` | 1 |
| `128 .. 255` | `C` | 1, unsigned |
| **`256 .. 32767`** | **`s`** | **2, signed** |
| `32768 .. 65535` | `S` | 2, unsigned |
| `65536 .. 2^31-1` | `i` | 4 |
| `2^31 .. 2^32-1` | `I` | 4, unsigned |

The bolded row is the trap. A value of 300 fits an unsigned short, and htsjdk writes a
**signed** one. Any encoder implementing "smallest type that fits, preferring unsigned" gets
`S` there. Confirmed from the oracle: `int_300` ends `...58 49 73 2c 01`, type `0x73 = 's'`,
while `int_200` ends `...58 49 43 c8`, type `0x43 = 'C'`.

The declared Java type has no influence at all: `getTagValueType` reaches every integral box
type through `((Number) value).longValue()`, so `Byte(100)`, `Short(100)`, `Integer(100)` and
`Long(100)` produce identical bytes. Arrays follow the **opposite** rule: their element width
comes from the array's class and is never narrowed.

### 2. Tags are ordered by the packed short, which weights the second character

`SAMTag.makeBinaryTag` is `(char[1] << 8) | char[0]`. The two characters land on disk in
reading order, so the packing is invisible in a hex dump, but the **numeric** value used for
ordering puts the second character in the high byte.

So `ZA` sorts before `AZ`, the reverse of the string order. Confirmed from the oracle: a record
given `ZA AZ NM MD AS XS SA` in that order is written `SA ZA MD NM AS XS`.

### 3. The trailing nibble of an odd-length read is `=`, not `N`

`bytesToCompressedBases` allocates `(len + 1) / 2` zeroed bytes and writes only the high nibble
of the last one. Zero decodes as `=`. Padding with `0x0F` (`N`), which looks like the safer
choice, gives a file that reads back identically and hashes differently.

### 4. The indexing `bin` is computed, and readers never check it

`computeIndexingBin` converts the start to 0-based and passes the still-1-based end through as
the half-open exclusive end. Converting *both* to 0-based, the symmetric-looking thing to do,
shifts every bin by one boundary. A wrong bin is invisible to any reader that scans linearly,
which is most of them.

## Why these four are one finding and not four

They share a shape, and it is the same shape as decisions 0001, 0002 and 0003: **the wrong
answer is valid.** Not "valid-looking": actually valid, parseable, semantically correct, and
accepted by every other tool in the ecosystem. Three of the four are what a careful engineer
implementing the published specification would naturally write.

That is the argument for porting from source rather than from the specification, and it is
worth restating because the specification is *right there* and porting is slower. The
specification defines what a BAM file may contain. It does not define what htsjdk writes.

## How it was verified

`tools/bam-conformance/BamCodecDump.java` runs `BAMRecordCodec.encode` inside the pinned
oracle image over 136 records chosen to hit each boundary: the whole promotion ladder including
every transition value, both box-type variants, tag sets whose packed order differs from their
string order, odd and even read lengths, every CIGAR operator, alignment starts on every bin
boundary at all six levels, and CIGARs of 65,535 / 65,536 / 65,537 operators around the `CG`
displacement threshold. Bin boundaries are reached with deletions rather than with matches, so
a 64 Mb reference span costs four read bases instead of 64 million; the first version of the
harness produced a 222 MB corpus that GitHub refused outright. It records; it asserts nothing. The comparison lives on the Rust side
so a bug in the harness cannot define the expected answer.

All 136 encode byte-identically, and all 136 survive a decode/re-encode round trip.

The harness was then checked against itself by sabotage, because a conformance suite that
passes on the first run has not yet demonstrated it can fail:

| sabotage | cases that diverged | first differing byte |
|---|---:|---|
| `region_to_bin` off by one at level 5 | 120 of 136 | 14, the `bin` field |
| tag order by string instead of packed short | 2 of 136 | 52, the start of the tag block |

Both point at exactly the right field. The second is the more instructive number: only two
records in the corpus carry more than one tag, so that trap is detectable by a **1.5%** slice
of the suite. Coverage of the parameter space is not the same as coverage of the failure modes,
and a corpus assembled without the failure modes in mind would plausibly have missed it
entirely.

## Consequence

`getReadNameLength` counts UTF-16 units and `StringUtil.stringToBytes` truncates each to one
byte, so `Z` tags and read names are neither UTF-8 nor ASCII. `café` occupies 4 bytes plus a
terminator. The port reproduces this rather than correcting it; a corrected encoder would
disagree with htsjdk on the length field, not merely on the content.

One item is deferred and recorded here so it is not lost: `BAMRecordCodec` forces the bin to 0
for reference sequences longer than `BIN_GENOMIC_SPAN`, after warning once. That needs the
sequence dictionary, so it belongs to the writer rather than to the record, and it lands with
the BAM file writer.
