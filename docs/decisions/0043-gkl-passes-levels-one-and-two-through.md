# 0043. GKL passes levels 1 and 2 through to ISA-L

**Status:** accepted
**Date:** 2026-09-23
**Corrects:** [0031](0031-gkl-levels-one-and-two-are-isal-level-one.md)

## What 0031 got wrong

0031 found that Java levels 1 and 2 produced identical bytes on its four fixtures, and it explained
that by saying GKL does not pass the level through to ISA-L. The configuration it recorded (ISA-L
level 1, `ISAL_DEF_LVL1_DEFAULT`, stateless, end of stream) is right for Java level 1. The
explanation is wrong, and so is using that configuration for Java level 2.

## How it showed

gatk-rs measured `FilterMutectCalls` writing a `.vcf.gz`. GATKConfig sets the level to 2, and the
BGZF block holding the 15 KB header came out 18 bytes longer than the reference's, even though the
decompressed text was identical. Every `.vcf.gz` measured before that was under 1.4 KB. Run in the
pinned container on that block's text, `IntelDeflater` gives:

| Java level | bytes | equals |
|---|---|---|
| 1 | 5321 | ISA-L level 1, stateless |
| 2 | 5303 | ISA-L level 2, stateless, and not level 1 at any buffer size |

The four fixtures of 0031 happen to compress identically at ISA-L levels 1 and 2. That coincidence
is what made the two Java levels look the same.

## The answer

`deflate_gkl` passes the Java level to ISA-L, with `ISAL_DEF_LVL1_DEFAULT` as the token buffer at
level 1 and `ISAL_DEF_LVL2_DEFAULT` at level 2. The eight (fixture, level) pairs of
`levels_one_and_two_go_through_isal_and_match_gkl` still match GKL's recorded column. The canary is
level 1 and is unchanged.

## What is still open

The recorded fixtures cannot tell level 1 from level 2, so they cannot catch this regression. A
fixture that separates the two levels belongs in the pinned-container column the next time it is
regenerated. Until then, gatk-rs's `FilterMutectCalls` covering array is the measurement that does.
