# 0052. htscodecs and HTSlib are sources for decoders, never for encoders

**Status:** accepted
**Date:** 2026-11-25

## Context

The CRAM 3.1 codecs (rANS 4x8 and Nx16, the adaptive arithmetic coder, fqzcomp, the name
tokeniser), BGZF and the FASTA indexes all have one correct decoded output: the format fixes it.
How fast a decoder reaches it is free. htscodecs 1.6 and HTSlib 1.24 are the fastest decoders of
these formats, and pasta (a MIT fork of noodles) has Rust ports of them.

An encoder is different. The bytes it writes are a choice, and this library must make htsjdk
5.0.0's and GKL's choices (igzip level 2, htsjdk's CRAM writer, `VCFEncoder`). htscodecs makes
other ones.

## Decision

- htscodecs and HTSlib join `ALLOWED` in `tools/audit/provenance.py`. htscodecs is under the
  Modified (3-clause) BSD licence, with its rANS and range coder primitives in the public domain;
  HTSlib is under the MIT/Expat licence, with its CRAM code under the 3-clause BSD licence. Both
  permit redistribution of altered versions with the notice kept, which `NOTICE` does.
- A port from them may only **decode**: the same answer as htsjdk's decoder, reached faster. What
  htsjdk refuses must still be refused, with htsjdk's error. A fast path that cannot fail runs
  only where the reference could not fail either, and the reference's loop stays as the fallback.
- Encoders stay htsjdk's and GKL's. From htscodecs, HTSlib or pasta they may take buffering and
  threading, never a choice of bytes.

## Consequences

A decoder taken from htscodecs is checked against the first, symbol-by-symbol port of htsjdk's on
conforming, truncated and damaged streams, as `rans_nx16`'s tests do.
