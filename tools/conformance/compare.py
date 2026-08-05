#!/usr/bin/env python3
"""The one comparator for every conformance suite.

Before this file existed, each suite carried its own copy of the same thirty lines of parsing
and diffing inside `.github/workflows/ci.yml`: 783 lines here, 1556 in picard-rs, and the copies
had already drifted (some keyed on two columns, some on one, some read tuples out of a Rust
source with a regex). A comparison that quietly ignores part of a file is worthless, so what each
suite compares, and what it ignores, now lives in `manifest.json` and this file is the only thing
that applies it. picard-rs carries the same file; the two are kept in step by hand, because the
repositories do not share a build.

The dump format, produced by the `*Dump.java` harnesses:

    <kind>\\t<case>\\t<payload>

one row per line, where the payload joins the record lines of a file with a *literal* backslash-n
(two characters), not a newline. That is why `SEGMENT_SEP` below is `'\\\\n'` and not `'\\n'`: a
comparison that split on real newlines would see one segment and silently canonicalize nothing.

Canonicalization rules are a closed set. Adding one is a reviewable event, because canonicalizing
is how a bit-identity claim quietly weakens.
"""

import gzip
import json
import sys
from pathlib import Path

# The payload separator emitted by the Java harnesses: a literal backslash followed by 'n'.
SEGMENT_SEP = "\\n"


# --------------------------------------------------------------------------------------
# Canonicalization rules. Each takes the payload and its rule spec, and returns the payload
# that will be compared. Every rule must be declared in the manifest with a `why`.
# --------------------------------------------------------------------------------------


def _strip_line_prefixes(payload, spec):
    """Drop whole segments starting with any declared prefix.

    Used for metrics files, whose `# <Tool> <command line>` and `# Started on: <timestamp>`
    headers carry JVM temp paths and a clock reading. The prefix list is per suite and holds
    *every* tool the harness runs, not just the current one: a multi-tool harness emits more
    than one header, and stripping only the current tool's left the others' temp paths in
    place. That bug made CI report 14 mismatches that were not divergences at all.
    """
    prefixes = tuple(spec["prefixes"])
    return SEGMENT_SEP.join(
        seg for seg in payload.split(SEGMENT_SEP) if not seg.startswith(prefixes)
    )


def _strip_pg(payload, spec):
    """Drop `@PG` segments: the provenance record whose `CL:` is the command line."""
    return SEGMENT_SEP.join(
        seg for seg in payload.split(SEGMENT_SEP) if not seg.startswith("@PG")
    )


def _strip_ur(payload, spec):
    """Drop `UR:` tab-fields: the reference's `file:` URI, which is path-dependent."""
    out = []
    for seg in payload.split(SEGMENT_SEP):
        out.append("\t".join(f for f in seg.split("\t") if not f.startswith("UR:")))
    return SEGMENT_SEP.join(out)


def _strip_banner(payload, spec):
    """Drop the metrics banner: `## htsjdk...` and `# ...` segments.

    This is the stripping `CmpCorpus.java` already applies on the Java side; declaring it here
    keeps the two sides describing the same operation.
    """
    return SEGMENT_SEP.join(
        seg
        for seg in payload.split(SEGMENT_SEP)
        if not (seg.startswith("## htsjdk") or seg.startswith("# "))
    )


RULES = {
    "strip_line_prefixes": _strip_line_prefixes,
    "strip_pg": _strip_pg,
    "strip_ur": _strip_ur,
    "strip_banner": _strip_banner,
}


def canonicalize(kind, payload, rules):
    """Apply the suite's declared rules, in order, to one row's payload."""
    for spec in rules:
        name = spec["rule"]
        if name not in RULES:
            raise SystemExit(f"unknown canonicalization rule: {name}")
        # A rule may be restricted to certain row kinds, e.g. `strip_ur` applies to the `dict`
        # row and not to the inputs beside it.
        kinds = spec.get("kinds")
        if kinds is not None and kind not in kinds:
            continue
        payload = RULES[name](payload, spec)
    return payload


# --------------------------------------------------------------------------------------
# Readers
# --------------------------------------------------------------------------------------


def _open(path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def parse_keyed(path, compare_spec):
    """Read a dump into {(kind, case): canonical payload}."""
    rules = compare_spec.get("rules", [])
    skip_kinds = set(compare_spec.get("skip_kinds", []))
    skip_comments = compare_spec.get("skip_comment_lines", True)
    rows = {}
    with _open(path) as fh:
        for line in fh:
            if skip_comments and line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t", 2)
            if len(parts) < 3 or parts[0] in skip_kinds:
                continue
            kind, case, payload = parts
            rows[(kind, case)] = canonicalize(kind, payload, rules)
    return rows


def parse_named(path, compare_spec):
    """Read a dump into {name: rest-of-line}.

    The single-key shape: one case per line, everything after the first tab is the value. Used
    where the case name is already unique (a BAM record's name, a VCF file's case).
    """
    rules = compare_spec.get("rules", [])
    skip_comments = compare_spec.get("skip_comment_lines", True)
    allow_empty = compare_spec.get("allow_missing_value", False)
    rows = {}
    with _open(path) as fh:
        for line in fh:
            if skip_comments and line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t", 1)
            if len(parts) < 2:
                # `interval_list` emits cases whose value is empty (an empty list is a result),
                # so dropping valueless lines there would drop real cases.
                if not allow_empty or not parts[0]:
                    continue
                parts.append("")
            rows[parts[0]] = canonicalize(parts[0], parts[1], rules)
    return rows


def parse_regex(path, compare_spec):
    """Read tuples out of a Rust test source with a declared regex.

    Two BGZF corpora are committed as literal arrays inside the test file rather than as a
    `.txt.gz`, so the oracle comparison has to read the source. The pattern, and which groups are
    the key, are declared in the manifest instead of being buried in a YAML heredoc.
    """
    import re

    pattern = re.compile(compare_spec["pattern"])
    key_groups = compare_spec["key_groups"]
    rows = {}
    with open(path) as fh:
        text = fh.read()
    for match in pattern.findall(text):
        key = tuple(match[i] for i in key_groups)
        value = tuple(v for i, v in enumerate(match) if i not in key_groups)
        rows[key] = value
    return rows


# The one bit pattern a comparison may not assert: the sign of a NaN.
#
# `Double.doubleToRawLongBits` records what the JVM held, which is the right thing for a dump to
# do. It is the wrong thing to compare, because the sign of a NaN is not a property of the
# arithmetic that produced it. commons-math3's `FastMathCalc` divides an overflowed exponential by
# itself, and `inf / inf` gives x86's negative default quiet NaN when the FPU computes it and
# Java's positive canonical NaN when the JIT folds it, on the same architecture, in the same
# container. One CI host produced each.
#
# So a suite may declare `canonicalise_nan: true` and have every NaN-shaped field collapsed on both
# sides. Only NaN patterns are touched: the guard is on the exponent and mantissa, so a value that
# is an index, a count or any finite double travels unchanged. No index in any dump comes within
# nine quintillion of the smallest NaN pattern.
NAN_EXPONENT_MASK = 0x7FF0000000000000
CANONICAL_NAN = 0x7FF8000000000000


def _canonicalise_nan_fields(line):
    """Replace every NaN-shaped integer field of a tab-separated line with the canonical one."""
    fields = line.split("\t")
    for i, field in enumerate(fields):
        try:
            value = int(field)
        except ValueError:
            continue
        bits = value & 0xFFFFFFFFFFFFFFFF
        if bits & NAN_EXPONENT_MASK == NAN_EXPONENT_MASK and bits & 0x000FFFFFFFFFFFFF:
            fields[i] = str(CANONICAL_NAN)
    return "\t".join(fields)


def parse_lines(path, compare_spec):
    """Read a dump as a line sequence.

    Suites whose corpus has several rows of one kind per case (the shards of a split, the
    inputs of a merge) cannot use a (kind, case) key: it would collapse those rows and compare
    only the last one.
    """
    skip_comments = compare_spec.get("skip_comment_lines", True)
    canonicalise = compare_spec.get("canonicalise_nan", False)
    with _open(path) as fh:
        lines = [
            line.rstrip("\n")
            for line in fh
            if line.strip() and not (skip_comments and line.startswith("#"))
        ]
    return [_canonicalise_nan_fields(line) for line in lines] if canonicalise else lines


# --------------------------------------------------------------------------------------
# Comparison
# --------------------------------------------------------------------------------------


def compare_case(real_path, golden_path, compare_spec):
    """Compare one regenerated dump against its committed golden.

    Returns (ok, compared_count, message lines).
    """
    mode = compare_spec.get("mode", "keyed")
    out = []

    if mode == "lines":
        real = parse_lines(real_path, compare_spec)
        committed = parse_lines(golden_path, compare_spec)
        if real != committed:
            out.append(f"lines differ: real={len(real)} committed={len(committed)}")
            for i, (r, c) in enumerate(zip(real, committed)):
                if r != c:
                    out.append(f"  first diff at line {i}")
                    out.append(f"    real     ={r[:200]}")
                    out.append(f"    committed={c[:200]}")
                    break
            return False, len(real), out
        return True, len(real), out

    if mode == "regex":
        # The committed side is a Rust source file, and the two sides legitimately hold different
        # case sets: the test file carries cases the dump does not emit and vice versa. So only
        # the intersection is compared, and an empty intersection is a failure rather than a pass,
        # because "nothing compared" must never read as "nothing wrong".
        real = parse_regex(real_path, compare_spec)
        committed = parse_regex(golden_path, compare_spec)
        shared = set(real) & set(committed)
        bad = [k for k in shared if real[k] != committed[k]]
        for k in list(bad)[:10]:
            out.append(f"  {k}: real={real[k]} committed={committed[k]}")
        if not shared:
            out.append(f"nothing compared: real={len(real)} committed={len(committed)} cases")
            return False, 0, out
        return not bad, len(shared), out

    if mode == "named":
        real = parse_named(real_path, compare_spec)
        committed = parse_named(golden_path, compare_spec)
    elif mode == "keyed":
        real = parse_keyed(real_path, compare_spec)
        committed = parse_keyed(golden_path, compare_spec)
    else:
        raise SystemExit(f"unknown compare mode: {mode}")

    if set(real) != set(committed):
        out.append(f"row sets differ: {sorted(set(real) ^ set(committed))}")
        return False, len(real), out
    bad = [k for k in real if real[k] != committed[k]]
    for k in bad[:5]:
        out.append(f"  {k}")
        out.append(f"    real     ={real[k][:200]}")
        out.append(f"    committed={committed[k][:200]}")
    return not bad, len(real), out


def load_manifest(path=None):
    path = Path(path or Path(__file__).with_name("manifest.json"))
    with open(path) as fh:
        return json.load(fh)


def suite_by_id(manifest, suite_id):
    for suite in manifest["suites"]:
        if suite["id"] == suite_id:
            return suite
    raise SystemExit(f"no suite {suite_id!r} in the manifest")


def main(argv):
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest")
    ap.add_argument("--suite", required=True)
    ap.add_argument(
        "--real",
        action="append",
        required=True,
        metavar="DUMP=PATH",
        help="the regenerated dump for one case, e.g. QualityYieldDump=/tmp/qy.txt",
    )
    args = ap.parse_args(argv)

    manifest = load_manifest(args.manifest)
    suite = suite_by_id(manifest, args.suite)
    reals = dict(pair.split("=", 1) for pair in args.real)

    failed, total = 0, 0
    for case in suite["cases"]:
        dump = case["dump"]
        if dump not in reals:
            print(f"FAIL {suite['id']}/{dump}: no regenerated dump supplied")
            failed += 1
            continue
        ok, compared, messages = compare_case(reals[dump], case["golden"], suite["compare"])
        total += compared
        status = "ok  " if ok else "FAIL"
        print(f"{status} {suite['id']}/{dump}: compared={compared}")
        for line in messages:
            print(line)
        failed += 0 if ok else 1

    print(f"suite={suite['id']} cases={len(suite['cases'])} compared={total} failed={failed}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
