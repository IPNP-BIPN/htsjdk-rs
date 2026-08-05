# Integrating noodles, one pull request per format

A working prompt. Hand this file, whole, to a session that is to carry out any part of it. It is
written to be self-contained: every constraint it relies on is stated here, with the measurement
that established it, so nothing has to be rediscovered.

---

## 1. What this repository is, in three sentences

`htsjdk-rs` is a Rust port of htsjdk 4.2.0 whose contract is **byte-identity**: the bytes it writes
are the bytes htsjdk writes, verified by an oracle. `gatk-rs` and `picard-rs` sit on top with the
same contract against GATK 4.6.2.0 and Picard 3.4.0. Forty-one conformance suites enforce it, each
one a Java dump run in a pinned container against a golden committed from a real x86-64 CI runner.

**Nothing in this plan may weaken that.** Everything here is additive.

## 2. The two surfaces, and which one a change belongs to

Decision 0039. Every public entry point belongs to exactly one:

**The conformance surface.** Byte-identity with htsjdk is the contract, the oracle enforces it,
and an extension may add entry points and may never alter one. All 41 suites are here.

**The extension surface.** Deliberately beyond htsjdk, so there are no bytes to be identical to.
Validated by the specification, by an independent reader, and by our own oracle-backed reader
where one exists. Lives under an `ext` module, never beside the ported code it extends, and never
carries a `Ported from` line.

## 3. What has already been measured about noodles

Do not re-derive these. Cite them.

| finding | decision | evidence |
|---|---|---|
| noodles' CRAM **primitives** are unreachable | 0036 | `io::reader::num` and `codecs::rans_4x8::{encode,decode}` are `pub(crate)`; calling them is `error[E0603]` |
| noodles' ITF8 reader **errors on truncation**; htsjdk returns `-1` | 0036 | `read_exact` versus `InputStream.read()`; the `cram-varint` golden pins htsjdk's |
| noodles **reads the CRAM htsjdk writes** | 0037 | 1646-byte file, four unmapped reads, 4/4 records identical in name, flags and sequence |
| the axis is **reading against writing**, not crate against port | 0037 | a reader's output is records, which the format defines; a writer's output is bytes, which it does not |
| htsjdk **ships the CRAM 3.1 codecs and refuses 3.1 files** | 0038 | `isSupportedVersion(3.1)` is false; opening a 3.1 definition throws; all four codec classes load |
| noodles **cannot be forced** to write htsjdk's VCF | 0039 | 3 differences, 2 forcible; `write_header` iterates by kind and htsjdk sorts across kinds |
| htsjdk's header comparator **is not a total order** | 0016 | it has a cycle, so the output depends on a Java `TreeSet`'s insertion order |
| the oracle **cannot produce bzip2 or LZMA CRAM blocks** | 0036 | Commons Compress is not on its classpath: `NoClassDefFoundError` |
| our VCF writer is **at parity** with noodles' | this branch | 0.97–0.99×, 131 MiB/s, after removing four allocations per record |

noodles ships every one to two weeks and its **0.x minors carry semantic changes** — 0.90 changed
`Filters::is_pass` and fixed CRLF stripping, which is the area the `vcf-read` golden pins. Every
noodles dependency is pinned to a single minor, and must stay so.

## 4. Invariants every pull request in this plan must satisfy

1. **All 41 suites pass, unmodified.** If a suite needs changing, the change is wrong. Say so and
   stop rather than adjust the suite.
2. **noodles is a dev dependency on the conformance surface and may be a dependency on the
   extension surface.** Never the reverse.
3. **A golden may only come from the pinned container on a real x86-64 CI runner.** Two commits:
   first the dump, the port and a `golden-pending` manifest entry with `expect_rows` and
   `expect_contains`; push; then the artefact as the golden, the Rust test, and the flip to
   `oracle-backed`. Never generate a golden locally — this machine is Apple Silicon.
4. **Run `python3 tools/conformance/generate_ci.py` after every manifest change**, and
   `--check` guards `ci.yml` against drift.
5. **`python3 tools/audit/provenance.py` must stay clean.** A file citing a source not in its
   `ALLOWED` map fails.
6. **Never chain a verification behind `grep -c`.** It exits 1 on zero matches, so
   `cargo clippy | grep -c warning && git commit` skips the commit exactly when the lint is clean.
   Run each check as its own command and read the output.
7. **Use absolute paths, and check `pwd` before any `git` command.** Three repositories sit side by
   side and a commit has already landed in the wrong one this way.
8. Every commit ends with
   `Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>`, no history is rewritten,
   and the project board is updated in the same pass rather than afterwards.
9. **No em dashes** anywhere in prose, commit messages or documentation.

## 5. The sequence

One pull request per row. Do them in this order: each depends on what the previous established.

---

### PR 1. VCF: measure noodles' reader against the `vcf-read` golden

**Why first.** VCF is the format where both sides are most complete here, and the golden that
would judge noodles has 102 rows including every htsjdk quirk this port has found.

**The measurement, before any code.** For each of the 24 whole-file cases in
`crates/htsjdk-vcf/tests/data/vcf_read.txt.gz`, feed the same text to `noodles_vcf::io::Reader` and
compare, row by row, against the golden. Report three counts: rows where noodles agrees with
htsjdk, rows where it differs, and rows it cannot express.

**The expected shape of the answer**, which must be confirmed rather than assumed: noodles will
return the same **records** and not the same **quirks**. The quirks the golden pins are the
fileformat line the reader substitutes, the standard header lines it rewrites, the version the
header forgets below 4.3, the two different line numbers for one malformed line, and the `#` line
in the body that decodes to null.

**What the PR contains.** `crates/htsjdk-vcf/tests/noodles_reader_agreement.rs`, a test that asserts
the counts, so the day noodles changes its reader the number moves and the file says so. A decision
record stating the rule that comes out of it, in the form "noodles may be used for X and may not be
used for Y", with X and Y named by the measurement rather than by category.

**Acceptance.** The counts are asserted, not printed. The decision record cites row numbers.

---

### PR 2. VCF 4.3 and 4.4 output, on the extension surface

**Why.** htsjdk reads 4.3 and refuses to write it, measured in decision 0035:
`IllegalArgumentException: Writing VCF version VCF4_3 is not implemented`. A caller who reads a 4.3
file with htsjdk cannot write it back. noodles writes 4.3 and 4.4.

**What the PR contains.** `crates/htsjdk-vcf/src/ext/vcf43.rs`, an **adapter**, not an
implementation: it converts this crate's `VcfHeader` and `VariantContext` into noodles' types and
calls its writer. `noodles-vcf` becomes a real dependency of this crate, gated to the `ext` module.

**Do not write the percent escaping by hand.** A first attempt at this module did, and reading
noodles' writer found two bugs in it within minutes: the escape set is **per column**, since `;`
and `=` delimit in INFO while `:` delimits in a sample field, and a value that is exactly `.` must
be written `%2E` or it reads back as *missing*. Both are in noodles already.

**Acceptance.** A round trip through **this crate's own reader**, which decision 0035 validated
against htsjdk: what noodles writes, the ported decoder reads back unchanged. That is the strongest
available check because one end of it is oracle-backed. Plus the specification, cited by section,
for each escape. 4.4 is refused rather than written as 4.3, because emitting a file claiming a
version it was never checked against is the one thing an extension must not do.

---

### PR 3. CRAM: the file definition and the container header

**Why.** `cram-varint` is oracle-backed and is the only thing in `htsjdk-cram`. The next structure
up is the file definition, `CRAM` plus two version bytes plus a 20-byte file id, then the container
header, which is a run of the ITF8s that already exist.

**This is conformance surface, hand-ported.** noodles' container reader is public but its
primitives are not, and the target is htsjdk's reading of a container, not the specification's.

**The dump.** `tools/cram-conformance/CramContainerDump.java`. Build small CRAMs with htsjdk in the
pinned container, exactly as `CramRoundTrip.java` in the scratch probes did, and emit the file
definition fields and every container header field. Include a file with more than one container so
the landmark offsets are exercised.

**Acceptance.** Two commits, golden from CI. Plus a noodles cross-check under decision 0037: the
containers this port reads and the containers noodles reads are the same containers.

---

### PR 4. CRAM: the compression header and the encodings

**Why.** This is the largest single piece of H.3 and it is where the port either works or does not.

**Scope, from decision 0038.** rANS **4x8** is required. rANS Nx16, the range coder, fqzcomp and
the name tokeniser are **not**, because no file htsjdk will open can contain them: 21 of the 169
files in `cram/` are out of scope on that evidence, and they are the hardest 21.

**Add the version probe decision 0038 asks for**: assert that
`isSupportedVersion(new CRAMVersion(3, 1))` is still false, so a jar bump that quietly adds 3.1
fails the run rather than silently invalidating the scoping.

**The bzip2 and LZMA problem.** The oracle cannot produce those blocks: Commons Compress is not on
its classpath. Either add the jars to `tools/oracle` in this PR and say so, or declare those two
external compressors out of scope with this record as the reason. Do not leave it implicit.

---

### PR 5. CRAM 3.1, on the extension surface

**Strictly after PR 4.** An extension extends something, and 3.1 is 3.0's container model with four
more codecs.

**What the PR contains.** An `ext` module that reads and writes 3.1 through noodles, which supports
it where htsjdk does not. This is the first place where the third-party crate is **ahead** of the
oracle, and it is the cleanest case in this whole plan: there is nothing to be identical to, and a
maintained implementation exists.

**Acceptance.** Round trip through noodles, plus: a 3.1 file this port writes must be refused by
htsjdk with `CRAM version 3.1 is not supported`, asserted in the container, because that is what
says the file really is 3.1.

---

### PR 6. BAM and SAM: extend the noodles cross-check

**Why.** `crates/htsjdk-bam/tests/noodles_cross_check.rs` exists and checks three records and three
tags. The BAM suites are the oldest and largest in the repository and the cross-check is thin.

**What the PR contains.** Extend it to every case the `bam-file` and `bam-codec` suites already
build, so that every BAM this port writes is read by a second implementation. Every tag type,
including `B` arrays and the empty one. The unmapped and placed-unmapped shapes. A record with a
CIGAR long enough to trigger the sentinel.

**Do not** replace the BAM reader with noodles: PR 1's measurement decides that for VCF and the
same measurement has to be made here before any such change. If PR 1's answer is "records agree,
quirks do not", the same is almost certainly true here and must still be shown.

---

### PR 7. BGZF: the one place a dependency is already the right answer

**Why.** Decision 0001 pins `flate2` to a **vendored zlib**, the same one the JDK links, which is
why it is a dependency and not a port. `noodles-bgzf` is a different question: it is a framing
layer over deflate, and the framing is what `bgzf` and `bgzf-termination` pin.

**The measurement.** Does `noodles-bgzf` produce the same blocks as this crate for the same input at
the same level? Almost certainly not, because block sizing is an implementation choice. Measure it
and record the answer; if it agrees, that is a genuine finding and the crate can be used.

---

### PR 8. BCF, tabix and CSI

**Why.** htsjdk supports BCF 2.1 and 2.2 and this port has none of it. GATK reads BCF.

**Order.** BCF depends on the VCF header model, which is complete. Tabix and CSI are index formats
that sit beside BGZF, which is complete.

**Surface.** Conformance, hand-ported, because htsjdk supports them: the oracle can judge. Use
noodles as the second reader as in PR 6, not as the implementation.

---

### PR 9. The formats GATK does not have: FASTA, FASTQ, GFF, GTF, BED, htsget, refget

**All extension surface.** htsjdk has some of these and GATK reaches few; noodles has all of them,
maintained.

**One PR per format, or one PR for the group** if the adapters are thin, which they will be: these
are the formats with no htsjdk quirks to reproduce because there is no htsjdk implementation to be
identical to.

**Acceptance is the weakest in this plan and must say so.** The specification and noodles, with no
oracle. Every module says that in its first paragraph.

---

## 6. Two things to measure once and cite everywhere after

**Speed.** Milestone S.1 asks for a harness measuring both sides on the same inputs.
`crates/htsjdk-vcf/tests/write_speed.rs` is the pattern: 50,000 records, release build, best of five
with the first discarded, `#[ignore]`d because a timing in CI is a flaky test. Repeat it per format
as each lands. The VCF writer went from 1.44× slower than noodles to 0.98× by removing four
allocations a record, with all 41 suites still green, so parity is reachable and byte-identity is
not what costs the time.

**Maintenance.** Every noodles dependency pinned to one minor, and a note in the PR of what changed
in noodles since the last one. Its changelog is the source; it ships every one to two weeks.

## 7. What would make this plan wrong

Say so and stop, rather than working around it:

 * a conformance suite needing modification to accommodate a noodles change;
 * a measurement showing noodles' reader disagreeing with htsjdk on **records** rather than on
   quirks, which would retire decision 0037's rule;
 * htsjdk 4.3 or later appearing as the oracle, which would turn every extension into a conformance
   obligation and is what decision 0038's version probe is for.
