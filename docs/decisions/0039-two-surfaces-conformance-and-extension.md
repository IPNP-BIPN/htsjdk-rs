# 0039. Two surfaces: what is byte-identical to htsjdk, and what is deliberately more

**Status:** accepted; sets the architecture for everything this port does beyond htsjdk
**Date:** 2026-08-05
**Follows:** [0004](0004-oracle-platform.md), [0037](0037-noodles-reads-htsjdks-cram-and-0036-said-otherwise.md), [0038](0038-htsjdk-ships-the-cram-3-1-codecs-and-refuses-3-1-files.md)

## The question

Every record before this one answers the same question: does this port produce the bytes htsjdk
produces? That question has a mechanical answer, an oracle, and thirty-nine suites.

The new requirement is different: support formats htsjdk **does not**. CRAM 3.1, which 0038 measured
htsjdk refusing outright. VCF 4.3 on the *write* side, which decision 0035 measured htsjdk refusing
with `Writing VCF version VCF4_3 is not implemented`. There is no oracle for either, because the
oracle's answer is an exception.

Both requirements are legitimate and they are not the same requirement. The danger is not doing the
second one; it is doing it in a way that makes the first one unverifiable, by leaving a reader unable
to tell which claim a given function carries.

## Decision

**The port has two surfaces, and every public entry point belongs to exactly one.**

### The conformance surface

Byte-identity with htsjdk 4.2.0 is the contract. Validated by the oracle: a Java dump in the pinned
container on a real x86-64 runner, a golden, and a Rust test that compares. Everything shipped so
far is here, and nothing on this surface may change because of anything on the other one.

The rule that enforces it: **an extension may add entry points and may not alter existing ones.** If
supporting VCF 4.3 output required `write_vcf` to behave differently for a 4.2 header, the extension
would be refused rather than the suite adjusted.

### The extension surface

Deliberately beyond htsjdk. Byte-identity is **not** the contract, because there are no bytes to be
identical to. What replaces it is stated per feature, and it is never nothing:

 * **the specification**, cited by section, because for these features the specification is the only
   normative text;
 * **an independent reader**, under 0037's rule. This is where noodles earns its keep for a second
   time: 0038 found it supports CRAM 3.1 where htsjdk does not, so on the extension surface the
   crate that was the second opinion becomes the *reference*;
 * **our own oracle-backed reader**, where one exists. Writing VCF 4.3 is checked by reading it back
   with the reader of decision 0035, which htsjdk validated. That is a real constraint: it says the
   extension emits what the conformance surface already agrees is a 4.3 file.

### Why writing cannot be delegated, measured rather than argued

`crates/htsjdk-vcf/tests/why_noodles_differs.rs` writes the simplest VCF there is with both this
port and `noodles-vcf`, and asserts the difference. Three lines differ, and **two of the three are
forcible through noodles' API**:

| | forcible |
|---|---|
| noodles writes the newest version it knows, `VCFv4.5`; htsjdk always writes `VCFv4.2` | **yes**, `set_file_format` |
| noodles keeps INFO fields in insertion order; htsjdk sorts them | **yes**, insert them sorted |
| noodles writes header lines grouped by kind; htsjdk sorts across kinds | **no** |

The third is hard-coded in `write_header`, which iterates `infos()`, then `filters()`, then
`formats()`, then `contigs()`. Insertion order decides the order within a kind and nothing decides
the order between kinds, so no sequence of `add_*` calls produces htsjdk's `FILTER` before `INFO`.

And a fork would not be enough. Decision 0016 measured htsjdk's comparator and found it **is not a
total order**: it has a cycle, so on some headers the output depends on a Java `TreeSet`'s insertion
order. A correct implementation cannot reproduce that; only a port of the same broken comparator
can. That is what `VcfHeader::write` is, and it is why it exists.

None of this is a defect in noodles. The specification fixes none of the three, so two correct
implementations differ — and VCF is the *text* format, the mildest case. A BAM adds a deflate stream
on top, where decision 0001 already measured two zlibs disagreeing.

### How a caller tells them apart

Three things, none of them optional:

 * **the module**. Extension code lives under an `ext` module in its crate, never beside the ported
   code it extends;
 * **the doc comment**. Conformance modules say "Ported from `<class>` at htsjdk 4.2.0". Extension
   modules say what htsjdk does instead, and cite the specification. `tools/audit/provenance.py`
   already enforces the first half; the second is why an extension module must not carry a
   `Ported from` line it cannot honour;
 * **the name**. An extension entry point names its version or its capability, so
   `write_vcf_v43` sits beside `write_vcf` and neither can be reached by accident.

## What this costs, said plainly

A feature on the extension surface has weaker evidence behind it than one on the conformance
surface, and no amount of testing changes that: an oracle says "this is what the reference does" and
a specification says "this is what the text can be read to mean". The two are not the same kind of
claim, and this record exists so that no future reader has to guess which one a function is making.

It also costs a real risk that is worth naming: if htsjdk later *adds* a feature this port already
extended, the extension becomes a conformance obligation and will almost certainly disagree, because
it was written from the text and htsjdk will have been written by someone else. That is not a reason
to refuse the extension. It is a reason for the version probe 0038 asks for, so the day the oracle
moves is a failing run rather than a silent divergence.

## Sequencing, which is not negotiable

An extension extends something. CRAM 3.1 support cannot precede CRAM 3.0 in this port, because
3.1 is 3.0's container model with four more codecs, and what exists today is the ITF8 and LTF8
integers underneath both. The order is therefore: the CRAM container model on the conformance
surface, then 3.1 on the extension surface.

VCF 4.3 output has no such dependency. The reader is oracle-backed, the writer exists, and the only
thing between them is htsjdk's refusal, so it is the first feature on this surface.
