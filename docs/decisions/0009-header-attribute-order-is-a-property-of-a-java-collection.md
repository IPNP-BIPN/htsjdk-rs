# 0009. The SAM header's byte order is decided by a Java collection choice

**Status:** accepted, measured against htsjdk
**Date:** 2026-07-21
**Follows:** [0008](0008-the-bam-record-has-four-places-a-correct-encoder-diverges.md)

## Finding

`AbstractSAMHeaderRecord` holds a header line's attributes in a **`LinkedHashMap`**:

```java
private final Map<String,String> mAttributes = new LinkedHashMap<String, String>();
```

`SAMTextHeaderCodec.encodeTags` iterates that map directly. So the byte order of every `@HD`,
`@SQ`, `@RG` and `@PG` line is whatever order the attributes were *inserted* in, and nothing
in the SAM specification, the codec, or the header classes says so.

Two consequences, and the second is the one that bites:

1. **Attributes are not sorted.** A `BTreeMap`, the natural Rust choice for a small string
   map, gives alphabetical order. Every multi-attribute line comes out different.
2. **Overwriting a key does not move it.** `LinkedHashMap.put` on an existing key replaces the
   value and leaves the entry in its original position. A `Vec` with remove-then-append, which
   is the obvious way to get "insertion order" in Rust, moves it to the end.

Confirmed from the oracle. Setting `SM`, then `LB`, then `SM` again, then `PL` gives:

```
@RG	ID:rgx	SM:second	LB:lib	PL:ILLUMINA
```

`SM` carries its second value and keeps its first position.

This matters in practice rather than only in principle: reading a header and writing it back
is what most tools do, and `setAttribute` on an already-present key is how any tool that
normalises a field reaches this path.

## Why this is worth its own record

Decisions 0001, 0002, 0003 and 0008 all describe places where the wrong answer is a valid
file. This one is different in kind, and worth naming separately: the correct behaviour is not
written down anywhere in htsjdk. It is not in a method, a comment, or a constant. It is an
emergent property of a data-structure choice made once, in a field declaration, and every
header htsjdk has ever written depends on it.

Porting from documentation cannot recover this. Porting from the source only recovers it if
the port reads the *declaration* and not merely the method that uses it. That is the concrete
argument for the rule that every feature must name the symbol it ports.

## Verification

`tools/bam-conformance/BamFileDump.java` emits 9 header shapes and 4 complete BAM files from
htsjdk in the pinned oracle container. All 13 are reproduced exactly.

The file cases are the first end-to-end comparison in the port: they compose the BGZF writer
(decisions 0001, 0003), the record codec (0008) and the header encoder, and framing errors
live in the seams rather than in the parts. The largest is 20,000 records over several BGZF
blocks, 128,733 bytes, byte-identical.

Sabotage check, sorting attributes alphabetically instead of by insertion:

| | diverged |
|---|---:|
| headers | 6 of 9 |
| whole files | 4 of 4 |

Every file diverges, because every file carries an `@RG`- or `@SQ`-bearing header, which is to
say: this trap is not an edge case, it is the common path.

## Consequence

`Attributes` in `crates/htsjdk-bam/src/header.rs` is a `Vec<(String, String)>` with explicit
update-in-place semantics, and the test `overwriting_an_attribute_keeps_its_original_position`
pins it. It is deliberately not a `HashMap`, a `BTreeMap`, or an `IndexMap`: the first two are
wrong, and the third would be right by accident of its `insert` semantics rather than by a
choice this port made on purpose.
