# 0042: The reference version is pinned by the consumers, not by this repository

**Status.** Accepted. Records why the target stays at htsjdk 4.2.0, and what the move costs when it
comes. Closes #182.

## The pin

This port targets htsjdk **4.2.0**, the version GATK 4.6.2.0's `build.gradle` names. GATK 4.7.0.0
(2026-08-18) pins htsjdk 5.0.0 and Picard 3.5.0.

The three repositories are coherent only because all three name **one** set of pins:

| Repo | Ports | Pins |
|---|---|---|
| `htsjdk-rs` | htsjdk 4.2.0 | (the floor) |
| `picard-rs` | Picard 3.4.0 | htsjdk 4.2.0 |
| `gatk-rs` | GATK 4.6.2.0 | Picard 3.4.0, htsjdk 4.2.0 |

Moving this repository alone would leave the other two reproducing a GATK that no longer exists:
every golden above here was taken against tools built on 4.2.0, so a 5.0.0 htsjdk under a 4.6.2.0
GATK is a combination nobody ships and nobody can check against. The decision is therefore not
this repository's to make, and it is recorded once, in IPNP-BIPN/gatk-rs#810: **the target moves
when the tool ports are done, and all three move together.**

## Why the issue is closed rather than left open

An issue is a thing to do. This is a thing not to do yet, decided elsewhere, and an open issue
tracking someone else's decision is a queue entry that can only be closed by that decision
changing. The content that mattered -- what will be waiting here on the move -- is below, where it
survives without anyone triaging it.

## What is waiting

**`jlibdeflate` is the default DEFLATE engine in 5.0.0.** Every BGZF and BAM byte here depends on
which deflater ran, and the dumps already print a `deflater\t<class>` line into the golden, which
is the guard that makes the bump survivable. On the move that pinning is re-verified rather than
assumed: the first question is whether libdeflate at a given level writes the bytes
`java.util.zip.Deflater` writes at that level. It almost certainly does not, and then a golden
taken with the 5.0.0 default is a *different artefact* from one taken with the JDK deflater, and
the manifest has to say which one it holds.

**CRAM 3.1 on the write path, defaulting to 3.1.** rANS Nx16, the adaptive arithmetic Range coder,
FQZComp, Name Tokenisation and STRIPE, over four profiles, with a `TrialCompressor` that picks the
codec per data series *by trial compression*. Reproducing the container means reproducing the
trial, because the choice is data-dependent and lands in the bytes. Decision 0039 already sized
this as the extension surface.

**Changes that are behaviour rather than defects here.** The LTF-8 nine-byte write used `>> 28`
where it needed `>> 24`, corrupting CRAM offsets past 256 GB; `SamLocusIterator` no longer offsets
the read position; `SamPairUtil.getPairOrientation` no longer answers asymmetrically on dovetail
pairs. This port reproduces the 4.2.0 behaviour, correctly, and each of those is a *target change*
on the move rather than a bug to fix now.

**The rest.** NM and MD stripped on encode and regenerated on decode; `TLEN` computed as htslib
computes it; slice headers without BD/SD/B5/S5/B1/S1; slices bounded by bases as well as reads;
`SAMRecord.toString()` returning the whole SAM line; SRA support removed.

## The rule

A version pin in a port is a property of the thing being reproduced, not a dependency to keep
current. Bumping it is a re-measurement of every golden in three repositories, and the only
question that decides when is whether the consumers are ready.
