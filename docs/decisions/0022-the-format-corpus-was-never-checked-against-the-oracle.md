# 0022. The format corpus was never checked against the oracle

**Status:** accepted
**Date:** 2026-07-29

## What was found

Moving the conformance suites into `tools/conformance/manifest.json` required naming, for each
committed golden, the CI step that regenerates it. Eighteen goldens had one. One did not:

| golden | rows | regenerated in CI |
|---|---:|---|
| `crates/htsjdk-metrics/tests/data/format.txt.gz` | 41,678 | **no** |

Its harness exists (`tools/metrics-conformance/FormatDump.java`) and its Rust test reads it
(`crates/htsjdk-metrics/tests/format_conformance.rs`), but no workflow step has ever run the
harness in the pinned container and compared the result. The corpus was produced once, committed,
and has been compared only against itself since.

That matters more here than it would for most corpora, because this is the corpus the README's
**99.73%** number-formatting figure rests on, and because decision 0011 established that metrics
number formatting is **locale-dependent**. The oracle image pins `en_US` and its probe refuses to
build otherwise; a developer machine pins nothing. So the one corpus whose correctness depends on
a property the container exists to guarantee is the one the container never saw.

`tools/zlib-conformance` is in the same position by its own admission: its README says the
comparison was "validated locally against OpenJDK 17.0.19 on macOS arm64" and is "not yet confirmed
inside the pinned `linux/amd64` container", which decision 0001 lists as remaining work. It is not
a golden-versus-oracle suite (it compares a Java run to a Rust run), so the manifest does not cover
it, and it stays open rather than being quietly counted as done.

## What was done

1. `format` is declared as a suite with `status: unchecked`, and generated into its own CI job
   whose name carries that status. Its next run is the first comparison it has ever had.
2. `tools/conformance/audit_goldens.py` fails the `guard` job if a committed golden belongs to no
   suite. Corpora checked by analysis rather than byte equality (the jmath corpus, whose
   per-function conditions come from decision 0007) are declared under `goldens_handled_elsewhere`
   with the job that checks them, so "not in a suite" cannot be used as a loophole.
3. The same audit and the same manifest shape now exist in picard-rs, where the equivalent search
   found sixteen unchecked goldens; see that repository's decision 0008.

## What this does not settle

Until the `format` suite has run green in CI, the 99.73% figure describes a measurement made on an
unpinned machine. If the first run diverges, the figure changes rather than the method: the number
that survives is the one the oracle produces.
