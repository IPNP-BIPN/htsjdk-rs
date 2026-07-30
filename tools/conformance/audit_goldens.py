#!/usr/bin/env python3
"""Check that every committed golden belongs to a declared conformance suite.

A golden that no suite regenerates has never been compared to anything but itself: the Rust test
reading it is asserting that the port reproduces a file of unknown provenance, which is weaker
than the claim the README makes. `format.txt.gz`, the 41,678-row corpus behind the 99.73% figure,
was in exactly that state until the manifest existed.

Corpora that are checked by an analysis job rather than by byte equality (the jmath corpus, whose
per-function pass conditions come from decision 0007) are declared in the manifest's
`goldens_handled_elsewhere`, with the job that checks them. Only committed files count
(`git ls-files`), so an in-flight slice's golden does not fail the audit before its suite lands.
"""

import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import compare as comparator  # noqa: E402

REPO = Path(__file__).resolve().parents[2]


def committed_goldens():
    out = subprocess.run(
        ["git", "ls-files", "crates/*/tests/data/*.gz"],
        cwd=REPO,
        capture_output=True,
        text=True,
        check=True,
    )
    return sorted(line for line in out.stdout.splitlines() if line)


def main():
    manifest = comparator.load_manifest()
    declared = {case["golden"] for suite in manifest["suites"] for case in suite["cases"]}
    elsewhere = {
        entry["path"]: entry["job"] for entry in manifest.get("goldens_handled_elsewhere", [])
    }
    committed = committed_goldens()

    undeclared = [g for g in committed if g not in declared and g not in elsewhere]
    missing = sorted(g for g in declared if not (REPO / g).exists())

    for golden in undeclared:
        print(f"FAIL undeclared golden: {golden}")
        print("     add a suite to tools/conformance/manifest.json that regenerates it, or")
        print("     declare it under goldens_handled_elsewhere with the job that checks it")
    for golden in missing:
        print(f"FAIL declared golden does not exist: {golden}")

    by_status = {}
    for suite in manifest["suites"]:
        by_status[suite["status"]] = by_status.get(suite["status"], 0) + len(suite["cases"])
    summary = " ".join(f"{status}={count}" for status, count in sorted(by_status.items()))
    print(
        f"goldens committed={len(committed)} declared={len(declared)} ({summary}) "
        f"handled elsewhere={len(elsewhere)}"
    )

    return 1 if (undeclared or missing) else 0


if __name__ == "__main__":
    sys.exit(main())
