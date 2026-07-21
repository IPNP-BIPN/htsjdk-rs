# 0015. htsjdk's SAM writer emits records its own reader rejects

**Status:** accepted; stringency modelled, both behaviours reproduced
**Date:** 2026-07-21
**Follows:** [0008](0008-the-bam-record-has-four-places-a-correct-encoder-diverges.md)

## Finding

The whole-file SAM conformance suite reads htsjdk's own SAM output back. One case failed:

```
unmapped: Inconsistent("MAPQ must be zero if RNAME is not specified")
```

The record has `RNAME = *` and `MAPQ = 60`. `SAMFileWriter` wrote it without complaint.
`SAMLineParser.parseLine` checks exactly that combination:

```java
} else {
    if (pos != 0)   reportErrorParsingLine("POS must be zero if RNAME is not specified");
    if (mapq != 0)  reportErrorParsingLine("MAPQ must be zero if RNAME is not specified");
    if (!cigar.equals("*")) reportErrorParsingLine("CIGAR must be '*' if RNAME is not specified");
}
```

and `reportErrorParsingLine` throws when the stringency is `STRICT`. `ValidationStringency`
declares:

```java
public static final ValidationStringency DEFAULT_STRINGENCY = STRICT;
```

So **htsjdk, at its default settings, writes a file it will not read.**

## Why this is a port problem and not a curiosity

A port that is unconditionally strict cannot read every file htsjdk produces. A port that is
unconditionally lenient cannot reproduce htsjdk's default behaviour, which is to refuse. Both
are wrong, and neither is discoverable from the writer alone: the writer performs no such check,
so reading the write path gives no hint that the read path has an opinion.

This generalises past this one field. htsjdk's validation lives on the **read** side and on
`SAMRecord.isValid()`, not in the encoder. The port had, until now, only encoders and readers
that reproduced encoder behaviour. Every consistency rule in `SAMLineParser` is a rule the
writer does not enforce.

## Decision

`ValidationStringency` is ported, with `Strict` as the default to match
`DEFAULT_STRINGENCY`, and every consistency check routed through one place so the stringency
governs them all together, exactly as `reportErrorParsingLine` does.

One distinction is reproduced carefully. htsjdk has **two** error paths:

| method | governed by stringency? |
|---|---|
| `reportErrorParsingLine` | yes — throws only at `STRICT` |
| `reportFatalErrorParsingLine` | **no** — throws always |

A malformed integer or an unparseable tag takes the second path. So a lenient parser still
refuses `POS = NOTANUMBER`, and a test pins that at all three stringencies. Collapsing the two
into one switch would make a lenient reader accept garbage.

## Consequence for the conformance suite

The whole-file test reads htsjdk's output at `SILENT`, and the reason is written at the call
site rather than left as a convenience. That is the honest shape: the suite is not relaxing its
standard, it is reading a file that cannot be read at the default standard, and saying so.

The finding also means the port's `write_sam` is not obliged to refuse what htsjdk writes. It
does not, which is why the BAM → SAM direction passed while the SAM → records direction failed.
Those two tests disagreeing is what surfaced this at all.
