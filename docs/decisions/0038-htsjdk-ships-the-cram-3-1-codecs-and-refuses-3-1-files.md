# 0038. htsjdk ships the CRAM 3.1 codecs and refuses 3.1 files, which sizes H.3

**Status:** accepted; scopes Milestone H.3
**Date:** 2026-08-05
**Follows:** [0037](0037-noodles-reads-htsjdks-cram-and-0036-said-otherwise.md)

## The question

`htsjdk/samtools/cram/compression/` contains directories named `fqzcomp`, `nametokenisation`,
`range` and `rans/ransnx16`. Those are the **CRAM 3.1** codecs: the specification that added them is
CRAMcodecs 3.1, and none of them can appear in a 3.0 file. Their presence reads as "htsjdk 4.2.0
supports CRAM 3.1", which would put 21 files of entropy coding and quality-score modelling inside
this milestone.

It does not. Measured three ways in the pinned container.

## The measurement

```text
supported  2.1  true
supported  3.0  true
supported  3.1  false
supported  4.0  false
default    3.0

class  ...rans.ransnx16.RANSNx16Encode              present
class  ...range.RangeEncode                         present
class  ...fqzcomp.FQZCompEncode                     present
class  ...nametokenisation.NameTokenisationEncode   present

open   3.0  htsjdk.samtools.util.RuntimeEOFException  null
open   3.1  java.lang.RuntimeException  CRAM version 3.1 is not supported
```

`CramVersions.isSupportedVersion` answers **false** for 3.1. The default is 3.0. And a file whose
definition declares 3.1 is refused at the version check with `CRAM version 3.1 is not supported`.

The 3.0 row is the control that makes the 3.1 row mean something: a 3.0 definition gets *past* the
version check and dies later, on truncation, because the probe's file is 26 bytes of header and
nothing else. So the 3.1 failure is the version gate and not the truncation.

**And all four codec classes load.** htsjdk ships the 3.1 codecs and will not read a file that could
contain them. They are reachable from their own unit tests and from nowhere else in a read path.

## What it scopes

**The port needs CRAM 2.1 and 3.0.** rANS **4x8** is 3.0's entropy codec and is required. rANS
Nx16, the range coder, fqzcomp and the name tokeniser are **not required for byte-identity with
htsjdk 4.2.0**, because no file htsjdk will open can contain a block encoded with them.

| directory | files | needed |
|---|---|---|
| `rans/rans4x8` | 3 | **yes** |
| `rans/ransnx16` | 3 | no |
| `range` | 6 | no |
| `fqzcomp` | 8 | no |
| `nametokenisation` | 4 | no |

Twenty-one of the 169 files in `cram/` are out of scope on this evidence, and they are the hardest
twenty-one: adaptive arithmetic coding and a quality-score model.

## What it does not license

"Not required" is not "will never be required". It is required the moment the oracle moves to an
htsjdk that supports 3.1, and this record is what says why the directory was skipped rather than
missed. The suite that would notice is a version probe: assert that
`isSupportedVersion(new CRAMVersion(3, 1))` is still false, so a jar bump that quietly adds 3.1
fails the run rather than silently invalidating this scoping.

## The other side of it

[noodles](https://github.com/zaeleus/noodles) supports CRAM 3.1 and htsjdk does not, so on that
axis the third-party crate is **ahead** of the oracle. That is harmless for the read side under
0037's rule: a 3.1 file is one htsjdk refuses, so no GATK run ever sees one, and there is nothing to
be identical to. It is worth stating because it is the first place where the reference
implementation is the narrower of the two.
