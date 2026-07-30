#!/usr/bin/env python3
"""Run one conformance suite (or one probe) against the oracle image.

The CI job and a local run go through this same script, so a suite that passes locally and fails
in CI is a real difference in the environment rather than a difference in how it was invoked.

    python3 tools/conformance/run_suite.py --suites metrics
    python3 tools/conformance/run_suite.py --suites "metrics rnaseq snvq"
    python3 tools/conformance/run_suite.py --probe rnaseq-overlap-order
    python3 tools/conformance/run_suite.py --list

The oracle image must exist; build it with

    docker build --platform linux/amd64 -t picard-rs-oracle:3.4.0 tools/oracle

Goldens are only valid when produced on real x86-64 (docs/decisions/0004 and the README's
bit-identity contract), so a local run on Apple Silicon is a smoke test, not a source of goldens.
"""

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import compare as comparator  # noqa: E402

REPO = Path(__file__).resolve().parents[2]


def container_command(cls, props):
    """The in-container script: compile the harness against the pinned jar, then run it.

    `2>/dev/null` drops htsjdk's chatter on stderr; the dump itself goes to stdout.
    """
    prop_str = (" ".join(props) + " ") if props else ""
    return (
        f'cp /harness/{cls}.java . && javac -cp "$ORACLE_CP" -d . {cls}.java '
        f'&& java {prop_str}-cp ".:$ORACLE_CP" {cls} 2>/dev/null'
    )


def docker_run(manifest, harness, cls, props, stdout):
    oracle = manifest["oracle"]
    cmd = [
        "docker",
        "run",
        "--rm",
        "--platform",
        oracle["platform"],
        "-v",
        f"{REPO}/{harness}:/harness:ro",
        "-w",
        "/work",
        oracle["image"],
        container_command(cls, props),
    ]
    print("+ " + " ".join(cmd[:-1]) + f" '{cmd[-1][:60]}...'", flush=True)
    with open(stdout, "w") as fh:
        result = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, text=True)
    # A compile failure used to report itself as "0 rows" and nothing else, because javac's
    # diagnostics went to the container's stderr and the run only kept stdout. The dump's own
    # chatter is already dropped inside the container by `2>/dev/null`, so anything arriving here
    # is the harness failing to build or to start, and it is the only thing that says why.
    if result.returncode != 0 and result.stderr.strip():
        print(f"--- {cls}: the container wrote to stderr:", flush=True)
        for line in result.stderr.rstrip("\n").split("\n"):
            print("   ", line, flush=True)
    return result.returncode


PENDING_DIR = REPO / "tools" / "conformance" / "pending"


def run_pending(manifest, suite, workdir):
    """A suite with no golden yet: run it, check the shape, and leave the dump for CI to publish.

    The alternative was to generate the golden here and commit it. That is exactly what produced
    the goldens of decision 0008, whose provenance turned out to be a laptop rather than the pinned
    container, so it is refused: this prints what the reference did, asserts only what the suite
    declares it must contain, and says plainly that nothing was compared.
    """
    props = suite.get("java_props", manifest.get("default_java_props", []))
    PENDING_DIR.mkdir(parents=True, exist_ok=True)
    failed = 0
    for case in suite["cases"]:
        dump = case["dump"]
        out = PENDING_DIR / f"{suite['id']}.{dump}.txt"
        rc = docker_run(manifest, suite["harness"], dump, props, out)
        rows = [line for line in open(out) if line.strip() and not line.startswith("#")]
        print(f"--- {suite['id']}/{dump}: {len(rows)} rows, nothing compared (no golden yet)")
        for line in rows:
            print("   ", line.rstrip("\n")[:200])
        expected = suite.get("expect_rows")
        if rc != 0 or (expected is not None and len(rows) != expected):
            print(f"FAIL {suite['id']}/{dump}: exit {rc}, {len(rows)} rows, expected {expected}")
            failed += 1
            continue
        # A row count is not evidence: a dump whose every row says "cannot read file" has the
        # right count. The behaviours the suite exists for are named, and each must appear.
        body = "".join(rows)
        for phrase in suite.get("expect_contains", []):
            if phrase not in body:
                print(f"FAIL {suite['id']}/{dump}: expected a row containing {phrase!r}")
                failed += 1
    print(
        f"suite={suite['id']} status=golden-pending cases={len(suite['cases'])} failed={failed}; "
        f"dumps in {PENDING_DIR} are the candidate goldens, valid only from a real x86-64 run"
    )
    return failed


def run_suite(manifest, suite, workdir):
    if suite["status"] == "golden-pending":
        return run_pending(manifest, suite, workdir)
    props = suite.get("java_props", manifest.get("default_java_props", []))
    failed = 0
    reals = {}
    for case in suite["cases"]:
        dump = case["dump"]
        out = Path(workdir) / f"{suite['id']}.{dump}.txt"
        rc = docker_run(manifest, suite["harness"], dump, props, out)
        lines = sum(1 for _ in open(out))
        print(f"regenerated {dump}: {lines} lines (docker exit {rc})")
        if rc != 0 or lines == 0:
            print(f"FAIL {suite['id']}/{dump}: the oracle produced no dump")
            failed += 1
            continue
        reals[dump] = out

    total = 0
    for case in suite["cases"]:
        dump = case["dump"]
        if dump not in reals:
            continue
        ok, compared, messages = comparator.compare_case(
            reals[dump], REPO / case["golden"], suite["compare"]
        )
        total += compared
        print(f"{'ok  ' if ok else 'FAIL'} {suite['id']}/{dump}: compared={compared}")
        for line in messages:
            print(line)
        failed += 0 if ok else 1

    print(
        f"suite={suite['id']} status={suite['status']} cases={len(suite['cases'])} "
        f"compared={total} failed={failed}"
    )
    return failed


def run_probe(manifest, probe, workdir):
    out = Path(workdir) / f"probe.{probe['id']}.txt"
    props = probe.get("java_props", manifest.get("default_java_props", []))
    docker_run(manifest, probe["harness"], probe["class"], props, out)
    text = open(out).read()
    print(text)
    # A probe may assert more than one line: decision 0019 measures the cyclic header and the
    # realistic dictionary, and either regressing invalidates the decision.
    expected = probe.get("expect_all") or [probe["expect"]]
    missing = [e for e in expected if e not in text]
    if missing:
        print(f"FAIL probe {probe['id']}: expected {missing}")
        print(probe["on_failure"])
        return 1
    print(f"ok   probe {probe['id']}: {', '.join(expected)}")
    return 0


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest")
    ap.add_argument(
        "--suites",
        help="one or more suite ids, space separated. Several suites share a CI job because each "
        "job pays the oracle image restore once.",
    )
    ap.add_argument("--probe")
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args(argv)

    manifest = comparator.load_manifest(args.manifest)

    if args.list:
        for suite in manifest["suites"]:
            print(f"{suite['id']:28} {suite['status']:14} {len(suite['cases'])} case(s)")
        for probe in manifest.get("probes", []):
            print(f"{probe['id']:28} {'probe':14} {probe['class']}")
        return 0

    with tempfile.TemporaryDirectory() as workdir:
        if args.suites:
            ids = args.suites.split()
            # Every suite runs even after one fails: the run exists to say which suites diverge,
            # not that at least one does.
            failed = sum(
                run_suite(manifest, comparator.suite_by_id(manifest, suite_id), workdir)
                for suite_id in ids
            )
            print(f"suites={len(ids)} failing={failed}")
            return 1 if failed else 0
        if args.probe:
            for probe in manifest.get("probes", []):
                if probe["id"] == args.probe:
                    return run_probe(manifest, probe, workdir)
            raise SystemExit(f"no probe {args.probe!r} in the manifest")

    ap.error("pass --suites, --probe or --list")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
