# 0035. Reading a VCF does not give you the file's header, and the version is the part it loses

**Status:** accepted; scopes H.2's read side and constrains what a byte claim over a VCF header means
**Date:** 2026-08-05
**Follows:** [0016](0016-the-vcf-header-comparator-is-not-a-total-order.md)

## The question

H.2 has a suite for the header frame, one for a data line and one for a genotype block. What it did
not have was the loop that joins them, and the entry sized that loop as thin: read the frame, then
call the line decoder once per line.

Before writing it, one thing was worth measuring rather than assuming, because every downstream
claim rests on it: **is the header a reader hands back the header the file contains?**

## It is not, in three separate ways

Measured in the pinned container, `tools/vcf-conformance/VcfReadDump.java`:

| the file says | the reader hands back |
|---|---|
| `##fileformat=VCFv4.0` | `fileformat=VCFv4.2` |
| `##fileformat=VCFv4.2` then `##fileformat=VCFv4.1` | `fileformat=VCFv4.2`, and the codec's version is **4.1** |
| `##INFO=<ID=DP,Number=1,Type=Float,Description="Depth">` | `##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">` |
| `##INFO=<ID=DP,Number=1,Type=Integer,Description="my own depth">` | unchanged |

**The fileformat line is deleted and a different one is put back.** The `VCFHeader` constructor
calls `removeVCFVersionLines`, so the stored metadata never holds it; `getMetaDataInInputOrder`
then *prepends* a synthesized line, which says `VCFv4.3` when the header's own version is 4.3 or
later and `VCFv4.2` in every other case. A v4.0 file therefore reads back declaring v4.2, and the
last of two declarations is the one the codec believes while neither is the one it reports.

**Standard INFO and FORMAT lines are replaced.** `doOnTheFlyModifications` defaults to true, so
every read runs `repairStandardHeaderLines`, and for the eighteen IDs htsjdk holds a definition for,
a disagreement on count or type replaces the **whole line**, description included. A disagreement on
the description alone does not, because `REPAIR_BAD_DESCRIPTIONS` is a private `false`. That is what
makes the rewrite hard to notice: the common case looks like a passthrough.

**The header forgets which version it is.** The repair rebuilds through the `VCFHeader` constructor
that takes no version and re-attaches it **only from 4.3 up**, so below 4.3 `getVCFHeaderVersion()`
returns null while `codec.getVersion()` still answers. Two different answers to the same question,
and a port with one field cannot give both.

## Why that last one is not a curiosity

It is load-bearing twice over.

`VCFWriter.rejectVCFV43Headers` tests `getVCFHeaderVersion()`, so **a 4.3 file can be read and
cannot be written back**: `IllegalArgumentException: Writing VCF version VCF4_3 is not implemented`,
where the message interpolates the enum constant and not the version string. The writer can refuse
4.3 precisely because a 4.2 header no longer knows what it is.

And `VCFUtils.smartMergeHeaders`, ported in htsjdk-rs #82, reads the same field. That suite already
recorded that the merge's version comes from a **field** set at parse time rather than from a
`##fileformat` line, and this is the other half of it: for anything below 4.3 the field is **empty
on every header a reader produced**, so the merge's version policy is unreachable through the read
path. It was reachable in that suite only because the headers there were assembled in memory.

## Decision

**The reader carries both versions, as two fields with two names.** `VcfFile::codec_version` and
`VcfFile::header_version`, the second an `Option`, and the doc comment on each says which upstream
call it answers. The repair table is its own module rather than a branch inside the reader, because
it is a data table with its own rule and it is consulted from the header path only.

**A byte claim about a VCF header is a claim about a header that has been through this.** Not the
file's bytes: htsjdk's rewrite of them. A round trip is measured in the same suite and is not the
identity, at any version, for two reasons that are both in this record.

## The other thing the same experiment settled

The line counter, which is the one piece of codec state a single-line suite cannot see. `lineNo` is
incremented in `readActualHeader` per header line, and again in `parseVCFLine` on entry. The column
check in `decodeLine` runs **before** that second increment and `generateException` runs after it,
so the same position reports two different numbers depending on which check refuses it: measured,
`Line 12` for a short line and `line number 13` for a bad INFO field, on two files with identical
twelve-line headers and one data line each. A `#` line in the body increments nothing at all,
because `decodeLine` returns null before reaching the counter, which also means such a line is a
**silently dropped record** rather than a refusal.

And one message is wrong upstream, reproduced as it is: the count in "there aren't enough columns"
is `header == null ? 8 : 9` while the check that produced it is `hasGenotypingData() ? 9 : 8`. A
sites-only file is checked against 8 and told it was expecting 9. A port that formats the checked
value into the message agrees with the check and disagrees with htsjdk.
