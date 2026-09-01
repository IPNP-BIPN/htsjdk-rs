#!/usr/bin/env bash
# Times the port's BGZF against htsjdk's, on the same payloads, inside the pinned oracle container,
# and asserts the two write the same bytes in the same run.
#
# Both sides run in the container for the same reason picard-rs's benchmark does: timing a native
# binary against an emulated JVM measures the emulator. The port is cross-built for linux/amd64
# first, so what runs there is a real amd64 binary rather than an interpreted one.
#
# A number produced on a developer machine is a harness check, not a measurement. Issue #78's
# figures come from the `benchmark` job, which runs this on a real x86-64 runner.
#
# Usage: tools/benchmark/run.sh [megabytes] [reps]   (default 64 3)
set -euo pipefail
cd "$(dirname "$0")/../.."
MB="${1:-64}"
REPS="${2:-3}"
OUT="$PWD/target/benchmark"
mkdir -p "$OUT"

echo "== cross-building the port for linux/amd64"
docker run --rm --platform linux/amd64 -v "$PWD":/src -v "$OUT":/out -w /src rust:1.97 \
  bash -c 'cargo build --release --bin bgzf-bench --target-dir /out/amd64 2>&1 | tail -1'

echo "== htsjdk, in the pinned container"
docker run --rm --platform linux/amd64 \
  -v "$PWD/tools/benchmark":/harness:ro -w /work htsjdk-rs-oracle:4.2.0 \
  "cp /harness/BgzfBench.java . && javac -cp \$ORACLE_CP -d . BgzfBench.java \
   && java -cp .:\$ORACLE_CP BgzfBench $MB $REPS" | tee "$OUT/java.txt"

echo "== the port, in the same container"
docker run --rm --platform linux/amd64 -v "$OUT":/out -w /out htsjdk-rs-oracle:4.2.0 \
  "/out/amd64/release/bgzf-bench $MB $REPS" | tee "$OUT/rust.txt"

echo "== bytes first"
python3 - "$OUT" <<'PY'
import sys, re, pathlib

out = pathlib.Path(sys.argv[1])
def digests(path, side):
    found = {}
    for line in path.read_text().splitlines():
        m = re.match(rf"{side}_deflate_(\w+)_level(\d+)_bytes=(\d+) md5=([0-9a-f]{{32}})", line)
        if m:
            found[(m.group(1), int(m.group(2)))] = (int(m.group(3)), m.group(4))
        m = re.match(r"payload_(\w+)_md5=([0-9a-f]{32})", line)
        if m:
            found[("payload", m.group(1))] = m.group(2)
    return found

java = digests(out / "java.txt", "java")
rust = digests(out / "rust.txt", "rust")
if not java or not rust:
    raise SystemExit("one side produced no digests at all")
if java != rust:
    for key in sorted(set(java) | set(rust), key=str):
        if java.get(key) != rust.get(key):
            print(f"  {key}: htsjdk={java.get(key)} port={rust.get(key)}")
    raise SystemExit("the two sides wrote different bytes; the timing below means nothing")
print(f"byte-identical on {len(java)} keys (payload digests and framed streams)")
PY

echo "== throughput, MB/s, median of $REPS"
python3 - "$OUT" <<'PY'
import sys, re, pathlib, statistics

out = pathlib.Path(sys.argv[1])
rates = {}
for side, name in (("java", "htsjdk"), ("rust", "port")):
    for line in (out / f"{side}.txt").read_text().splitlines():
        m = re.match(rf"{side}_(deflate|inflate)_(\w+)_level(\d+)_run\d+_mbps=([\d.]+)", line)
        if m:
            key = (m.group(1), m.group(2), int(m.group(3)))
            rates.setdefault(key, {}).setdefault(name, []).append(float(m.group(4)))

print(f"{'operation':10} {'payload':8} {'level':>5} {'htsjdk':>9} {'port':>9} {'ratio':>7}")
for key in sorted(rates):
    row = rates[key]
    if "htsjdk" not in row or "port" not in row:
        continue
    j = statistics.median(row["htsjdk"])
    r = statistics.median(row["port"])
    print(f"{key[0]:10} {key[1]:8} {key[2]:>5} {j:>9.1f} {r:>9.1f} {r / j:>6.2f}x")
PY
