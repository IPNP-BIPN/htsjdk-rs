# 0016. The VCF header comparator is not a total order

**Status:** accepted; reproduced where it is well-defined, diverged where htsjdk is not
**Date:** 2026-07-21
**Follows:** [0009](0009-header-attribute-order-is-a-property-of-a-java-collection.md)

## Two header formats, opposite rules

Decision 0009 found that a SAM header's attribute order is **insertion order**, from a
`LinkedHashMap`. The same library orders VCF header lines the other way:

```java
public Set<VCFHeaderLine> getMetaDataInSortedOrder() {
    return makeGetMetaDataSet(new TreeSet<VCFHeaderLine>(mMetaData));
}

public int compareTo(Object other) {
    return toString().compareTo(other.toString());
}
```

A `TreeSet` over the **rendered string of the whole line**. So `##FILTER` before `##FORMAT`
before `##INFO` before `##alsoUnstructured`, because `I` (0x49) sorts before `a` (0x61) in
plain ASCII.

Neither format states its rule anywhere except in the choice of collection. A port that guessed
"sorted" for SAM or "insertion" for VCF would be wrong in both directions, and would produce
valid files either way.

## And the VCF rule is not consistent

`VCFContigHeaderLine` overrides the comparison:

```java
public int compareTo(final Object other) {
    if (other instanceof VCFContigHeaderLine)
        return contigIndex.compareTo(((VCFContigHeaderLine) other).contigIndex);
    else
        return super.compareTo(other);   // by rendered string
}
```

Contigs compare to each other **by index** and to everything else **by string**. Those two
orders can disagree, and when they do the comparator has a cycle.

Constructed and measured in the pinned oracle. Three lines, with a non-contig line whose key is
also `contig` so it compares by string:

```
zzz vs aaa (by index)  = -1
zzz vs mmm (by string) = +1
mmm vs aaa (by string) = +1
```

which is `aaa < mmm < zzz < aaa`. Feeding the same three lines to `VCFHeader` in two different
orders:

```
inserted forward : ##contig=<ID=mmm> | ##contig=<ID=zzz> | ##contig=<ID=aaa>
inserted reversed: ##contig=<ID=aaa> | ##contig=<ID=mmm> | ##contig=<ID=zzz>
SAME_OUTPUT=false
```

**The same logical header serialises two ways depending on the order its lines were added.**

## Why this is not merely academic

It bears directly on the bit-identity claim. A VCF header is not a set of lines with a
canonical rendering; it is a set of lines plus the history of how they were assembled. Any tool
that builds a header by merging sources, which is most of them, can produce either output.

And there is a second consequence that lands on decision 0013. Reproducing htsjdk exactly in
the inconsistent case means reproducing what `TreeSet` does with an inconsistent comparator,
which is a property of `java.util.TreeMap`'s red-black balancing. That is `java.base`, GPL2, and
therefore **not portable** by the rule decision 0014 established. The pathological case is not
just hard to reproduce; it is legally out of reach.

## Decision

**Reproduce the rule where it is well-defined; use a genuine total order where htsjdk has none.**

`SortKey` gives contigs a synthetic string of `contig=` followed by their zero-padded index, and
every other line its rendered string. That:

- orders contigs among themselves by index, as htsjdk does;
- orders contigs against any line not beginning `contig=` identically to htsjdk, because the
  `contig=` prefix is preserved and the comparison resolves before the padding;
- is a total order, so the output does not depend on insertion sequence.

It diverges from htsjdk only for a header containing a non-contig line whose rendered string
begins `contig=` — which is exactly the header on which htsjdk itself is order-dependent, so
there is no single answer to match.

Real headers never hit this. No line other than a contig renders a string starting with
`contig=`, so the cross comparisons never interleave. That is a property of the key namespace,
not of the comparator, and **nothing enforces it**: `new VCFHeaderLine("contig", "<...>")` is a
public constructor.

## Note on how this was found

Not by reading `compareTo`. The conformance suite had a case with twelve contigs, `chr1` to
`chr12`, and it failed on `chr10` versus `chr2` — the port sorted them lexicographically and
htsjdk did not. Chasing that one mismatch led to the override, and the override suggested the
cycle, which was then constructed deliberately and measured.

The corpus contained twelve contigs because twelve is what a real header has, not because
anyone expected them to be interesting.
