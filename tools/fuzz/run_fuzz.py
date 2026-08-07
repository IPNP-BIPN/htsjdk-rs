#!/usr/bin/env python3
"""Drive a differential fuzzer over the byte-level parsers, and report what diverges.

Three steps.

1. **Generate** a corpus: a fixed set of interesting inputs, then mutations of them under a seed
   given on the command line. The seed is an argument rather than the clock, so a finding can be
   reproduced by rerunning the same command.
2. **Run both sides** over the same corpus. The reference runs in the pinned container through
   `FuzzDriver.java`; the port runs through `crates/htsjdk-cram/examples/differential_fuzz.rs`.
   Each prints one line per input, so the comparison is a diff rather than a report.
3. **Minimise** every divergence: bytes are dropped from the end for as long as the divergence
   survives, so what gets written to `findings/` is the shortest input that still shows it.

    python3 tools/fuzz/run_fuzz.py --parser itf8 --iterations 400 --seed 1

Findings are evidence, not goldens. They say what the reference did and where the port differs,
and one produced anywhere other than real x86-64 CI cannot be quoted (decision 0008).
"""

import argparse
import json
import random
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
IMAGE = "htsjdk-rs-oracle:4.2.0"

# The inputs worth having whatever the seed does: every ITF8 and LTF8 width boundary, the empty
# input, and the shapes a CRAI line can take.
SEEDS = {
    "itf8": ["", "00", "7f", "80", "8080", "bfff", "c04000", "dfffff", "e0200000",
             "efffffff", "f100000000", "ffffffffff", "ff", "ffff", "ffffff", "ffffffff"],
    "ltf8": ["", "00", "7f", "80", "8080", "c0", "e0", "f0", "f8", "fc", "fe", "ff",
             "ff7fffffffffffffff", "ffffffffffffffffff", "f100000000"],
    "crai": [
        "302039093530", "30093130300935300931303030093230093330300a",
        "3009310931093009300931",          # 0 1 1 0 0 1
        "2d3109300930093530300931300932300a",
        "78093130300935300931303030093230093330300a",  # a letter in the first column
        "",
        "09090909090909",
    ],
}


def corpus(parser, iterations, seed):
    """The seeds, then mutations of them: a byte flipped, a byte appended, a byte removed."""
    random = _random(seed)
    cases = list(SEEDS[parser])
    seen = set(cases)
    while len(cases) < iterations:
        base = random.choice(SEEDS[parser])
        bytes_ = bytearray.fromhex(base)
        for _ in range(random.randint(1, 3)):
            choice = random.randint(0, 2)
            if choice == 0 and bytes_:
                bytes_[random.randrange(len(bytes_))] = random.randrange(256)
            elif choice == 1 and len(bytes_) < 24:
                bytes_.append(random.randrange(256))
            elif bytes_:
                del bytes_[random.randrange(len(bytes_))]
        case = bytes_.hex()
        if case not in seen:
            seen.add(case)
            cases.append(case)
    return cases


def _random(seed):
    generator = random.Random()
    generator.seed(seed)
    return generator


def reference(parser, cases):
    command = (
        f'cp /harness/FuzzDriver.java . && javac -cp "$ORACLE_CP" -d . FuzzDriver.java '
        f'&& java -cp ".:$ORACLE_CP" FuzzDriver {parser} 2>/dev/null'
    )
    result = subprocess.run(
        ["docker", "run", "--rm", "-i", "--platform", "linux/amd64",
         "-v", f"{REPO}/tools/fuzz:/harness:ro", "-w", "/work", IMAGE, command],
        input="\n".join(cases) + "\n", capture_output=True, text=True,
    )
    if result.returncode != 0:
        print(result.stderr[-2000:], file=sys.stderr)
        raise SystemExit("the reference side did not run")
    return parse(result.stdout)


def port(parser, cases):
    result = subprocess.run(
        ["cargo", "run", "--quiet", "--release", "-p", "htsjdk-cram",
         "--example", "differential_fuzz", "--", parser],
        input="\n".join(cases) + "\n", capture_output=True, text=True, cwd=REPO,
    )
    if result.returncode != 0:
        print(result.stderr[-2000:], file=sys.stderr)
        raise SystemExit("the port side did not run")
    return parse(result.stdout)


def parse(text):
    outcomes = {}
    for line in text.splitlines():
        parts = line.split("\t")
        if len(parts) == 3:
            outcomes[parts[0]] = parts[2]
    return outcomes


def minimise(parser, case, expected, actual):
    """Drop bytes from the end for as long as the divergence survives."""
    shortest = case
    while len(shortest) >= 2:
        shorter = shortest[:-2]
        theirs = reference(parser, [shorter]).get(shorter)
        ours = port(parser, [shorter]).get(shorter)
        if theirs is None or ours is None or theirs == ours:
            break
        shortest, expected, actual = shorter, theirs, ours
    return shortest, expected, actual


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--parser", default="itf8", choices=sorted(SEEDS))
    ap.add_argument("--iterations", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args(argv)

    cases = corpus(args.parser, args.iterations, args.seed)
    print(f"{args.parser}: {len(cases)} inputs, seed {args.seed}")

    theirs = reference(args.parser, cases)
    ours = port(args.parser, cases)

    divergent = [case for case in cases
                 if theirs.get(case) != ours.get(case)]
    print(f"{len(cases) - len(divergent)} agree, {len(divergent)} differ")

    if not divergent:
        return 0

    findings = REPO / "tools" / "fuzz" / "findings"
    findings.mkdir(parents=True, exist_ok=True)
    report = []
    for case in divergent[:20]:
        shortest, expected, actual = minimise(
            args.parser, case, theirs.get(case), ours.get(case))
        report.append({"parser": args.parser, "input": shortest,
                       "reference": expected, "port": actual, "found_from": case})
        print(f"  {shortest}: reference={expected} port={actual}")
    (findings / f"{args.parser}-seed{args.seed}.json").write_text(
        json.dumps(report, indent=2) + "\n")
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
