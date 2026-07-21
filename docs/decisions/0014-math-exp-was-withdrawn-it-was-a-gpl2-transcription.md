# 0014. `Math.exp` is withdrawn: it was a transcription of GPL2-only source

**Status:** accepted; code removed, function reclassified as unported
**Date:** 2026-07-21
**Follows:** [0013](0013-the-last-divergences-are-blocked-by-a-licence-not-by-difficulty.md)

## What happened

Decision 0013 established that `FloatingDecimal` cannot be ported because `java.base` is GPL2
and the OpenJDK Assembly Exception grants permission to *link*, not to translate. It closed by
observing that the gap was general: "every place the port needs JVM behaviour rather than
htsjdk behaviour hits the same wall", and it named `Math` intrinsics as one of those places.

An audit of every symbol the port claims to have taken from a source ran immediately afterwards.
It found this, in already-merged, already-published code:

```rust
//! `java.lang.Math.exp`, reproducing HotSpot's x86 intrinsic.
//!
//! Ported from openjdk/jdk `src/hotspot/cpu/x86/macroAssembler_x86_exp.cpp` at `jdk-17-ga`
```

and, in the companion table:

```rust
//! Lookup table for `exp`, transcribed from HotSpot's `_Tbl_addr`.
```

The header of that source file, fetched from `openjdk/jdk17u` at the tag rather than recalled:

```
* Copyright (c) 2016, Intel Corporation.
* Intel Math Library (LIBM) Source Code
*
* This code is free software; you can redistribute it and/or modify it
* under the terms of the GNU General Public License version 2 only, as
* published by the Free Software Foundation.
```

**GPL version 2 only.** No Classpath Exception — that exception covers the class libraries, and
HotSpot is not the class libraries. So `crates/jmath/src/exp.rs` and
`crates/jmath/src/exp_table.rs` were a line-by-line translation of GPL2-only code, published in
an MIT-licensed repository. That is a licence violation, and it was mine.

## What was done

Both files are removed. `jmath::math::exp` no longer exists. `Math.exp` now has the same status
as `Math.pow`: unported, with the reason recorded rather than the gap left unexplained.

The conformance corpus keeps its `exp` points, routed to the system libm, so the function stays
in the reported table instead of quietly vanishing:

```
agreement with java.lang.Math: exp=99.9711%  pow=99.9378%  log1p=99.5889% …
```

**That 0.0289% is the measured price of the licence.** It is not an estimate: it is the fraction
of the 44,987-point `exp` corpus where the system libm and `java.lang.Math` disagree, and it is
what the port now gets wrong on every `exp` call site.

## Why the neighbours are clean, and why that was luck

The audit checked the whole crate.

| module | provenance | status |
|---|---|---|
| `sqrt` | IEEE-754 mandates the rounding | clean, no source consulted |
| `log`, `log10` | correctly-rounded double-double, constants generated at 400-bit precision | clean, independent implementation |
| `dd` | Dekker/Knuth two-sum and two-product, classical | clean |
| `exp` | **transcription of GPL2-only HotSpot** | **removed** |
| `pow` | never written, deferred by decision 0007 | clean by accident of scheduling |

`log` and `log10` are clean for a reason that had nothing to do with licensing. Decision 0006
asked "is the intrinsic correctly rounded?", found that it was, and concluded that *rounding the
true result* was a smaller job than porting the algorithm. That decision was made on effort
grounds. It happens also to be the only reason those two functions are not in the same position
as `exp`.

And `pow` is clean only because decision 0007 deferred it for a *different* reason — that its
intrinsic depends on an approximate hardware instruction — before anyone got round to
transcribing 2,220 lines of the same GPL2 file.

So two of the four escaped by luck and one by unrelated caution. That is not a system, and the
consequence is below.

## The rule this establishes

**Before porting any symbol, check the licence of the file it lives in.** Not the licence of the
project: the licence of the file. htsjdk is MIT, Picard is MIT, GATK is Apache 2.0, and none of
that says anything about the JDK sources those three run on.

Concretely, a symbol is portable into this program only if it comes from the pinned htsjdk,
Picard or GATK clones. Anything reached through `java.lang`, `java.util`, `java.text` or HotSpot
is GPL2 and is **not** portable. When such behaviour is observable in an output byte, the
options are the ones decision 0013 lists: reproduce it from first principles as decision 0006
did, obtain permission, change the oracle, or quarantine the values as bio-identical.

`tools/audit/provenance.py` runs that check mechanically over every `Ported from` claim in the
tree, and CI runs it, so the next one is caught before it is published rather than eleven
commits later.

## On the disclosure

The infringing code was public for a matter of hours, in a repository with no users. It is
removed from `main`. It remains in the git history of `IPNP-BIPN/htsjdk-rs`, which a history
rewrite could excise; that is a judgement call for the repository owner, and this record exists
so the call can be made with the facts rather than discovered later.
