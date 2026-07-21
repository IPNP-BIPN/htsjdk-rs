# 0018. Four ordering rules in one library, none of them written down

**Status:** accepted; all four reproduced
**Date:** 2026-07-21
**Follows:** [0009](0009-the-sam-header-keeps-insertion-order.md),
[0016](0016-the-vcf-header-comparator-is-not-a-total-order.md)

## Finding

Porting `VCFEncoder` completes a set. htsjdk orders four different sequences of key/value data,
and uses a different rule for each:

| what | rule | where it is decided |
|---|---|---|
| SAM header attributes | **insertion order** | `LinkedHashMap` in `AbstractSAMHeaderRecord` |
| VCF header lines | **sorted by the rendered text of the whole line** | `TreeSet` + `VCFHeaderLine.compareTo` |
| VCF INFO fields | **sorted by key alone** | `new TreeMap<>()` in `VCFEncoder.write` |
| VCF FORMAT keys | **sorted, then `GT` prepended** | `calcVCFGenotypeKeys` |

No two are alike. None is stated in the SAM or VCF specifications, in a Javadoc comment, or
anywhere else in the source. Each is decided entirely by which collection class the author
reached for, in a single line that reads as an implementation detail:

```java
final Map<String, String> infoFields = new TreeMap<>();
```

## Why this matters more than it looks

The failure mode is uniform across all four. A port that picks one convention and applies it
throughout produces a file that is **valid**, that every other tool reads correctly, and that
differs from htsjdk's bytes. Round-tripping does not catch it. Schema validation does not catch
it. Only byte comparison does.

And the convention most likely to be picked is insertion order, because that is what an
`IndexMap`, a `Vec<(K, V)>`, or a struct field list gives for free in Rust. That choice is right
for the SAM header and wrong for the other three.

The FORMAT rule deserves its own note because it looks like domain knowledge and is not:

```java
List<String> sortedList = ParsingUtils.sortList(new ArrayList<>(keys));
if (sawGoodGT) {
    final List<String> newList = new ArrayList<>(sortedList.size() + 1);
    newList.add(VCFConstants.GENOTYPE_KEY);
    newList.addAll(sortedList);
    sortedList = newList;
}
```

`GT` is prepended **after** the sort, so it never participates in it. The familiar
`GT:AD:DP:GQ:PL` is not an ordering by significance, as it is often read; it is `GT` followed by
plain ASCII order. `GT:AD:DP:GQ:PL` and the alphabet agree by coincidence.

## Two more rules in the same shape

Two further behaviours found while porting the encoder have the same character, in that the
right answer is decided somewhere other than where a reader would look for it.

**Allele case.** `SimpleAllele`'s constructor uppercases the bases, except for symbolic alleles,
which it leaves alone. And "symbolic" is `wouldBeSymbolicAllele`, which is true for anything
starting with `<`, **ending** with `>`, containing `[` or `]`, or starting or ending with `.`.
So `<del>` keeps its case, `at>` keeps its case, and `at` becomes `AT`. Uppercasing
unconditionally fails 5 of the 54 corpus cases; not uppercasing at all fails 1.

**`formatVCFDouble`'s signed branches.**

```java
if (d < 1) {
    if (d < 0.01) { ... "%.3e" ... }
    else          { ... "%.3f" ... }
} else            { ... "%.2f" ... }
```

Both comparisons are signed, so the sign decides as much as the magnitude does. `-1e10` takes
the same branch as `1e-10` and prints `-1.000e+10`, while `1e10` prints `10000000000.00`. Only
the `Math.abs(d) >= 1e-20` guard looks at magnitude alone, and it sends everything nearer zero
than 1e-20, in either direction, to the literal string `0.00`.

## Decision

All four orderings and both formatting rules are reproduced exactly, each with the collection or
comparison that produces it, and each carrying a comment naming the htsjdk line that decides it.
Where the Rust-natural choice differs from htsjdk's, the comment says which and why, because the
next person to read the code will otherwise "fix" it.

The one place htsjdk cannot be reproduced is the VCF header's comparator when its inconsistency
is triggered, which is decision 0016 and unchanged by this.

## Verification

`crates/htsjdk-vcf/tests/record_conformance.rs`, 54 record cases against the pinned oracle, all
byte-equal, with case names asserted so a case cannot silently leave the list.

Sabotage-checked, since a conformance suite that cannot fail proves nothing:

| deliberate break | cases failed |
|---|---|
| INFO in insertion order | 1 |
| symbolic alleles uppercased | 5 |
| Rust's native `{:.2}` in place of the JVM model | 2 records, 3 sweep |
