#!/usr/bin/env bash
# Times the port's BGZF, BAM and VCF codecs against htsjdk's, on the same inputs, inside the
# pinned oracle container, and asserts the two write the same bytes in the same run.
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
# One `-p` per binary, and not one build of all three: selecting binaries across the workspace
# unifies FEATURES across it too, which pulled `gkl-deflate/isal` into a build with no assembler
# and failed in `isal-sys`'s configure. Each package builds with the features it declares.
docker run --rm --platform linux/amd64 -v "$PWD":/src -v "$OUT":/out -w /src rust:1.97 \
  bash -c 'set -e
    cargo build --release -p htsjdk-bgzf --bin bgzf-bench --target-dir /out/amd64 2>&1 | tail -1
    cargo build --release -p htsjdk-bam --bin record-bench --target-dir /out/amd64 2>&1 | tail -1
    cargo build --release -p htsjdk-vcf --bin vcf-bench --target-dir /out/amd64 2>&1 | tail -1'

echo "== htsjdk, in the pinned container"
docker run --rm --platform linux/amd64 \
  -v "$PWD/tools/benchmark":/harness:ro -w /work htsjdk-rs-oracle:4.2.0 \
  "cp /harness/BgzfBench.java . && javac -cp \$ORACLE_CP -d . BgzfBench.java \
   && java -cp .:\$ORACLE_CP BgzfBench $MB $REPS" | tee "$OUT/java.txt"

echo "== the port, in the same container"
docker run --rm --platform linux/amd64 -v "$OUT":/out -w /out htsjdk-rs-oracle:4.2.0 \
  "/out/amd64/release/bgzf-bench $MB $REPS" | tee "$OUT/rust.txt"

# The record codec, which #78 names as the second path under everything: every read filter and
# every walker pays it per record. Its in-memory rows have no htsjdk column -- htsjdk's codec writes
# into a stream rather than answering with bytes -- and they print the digest of the corpus they
# encoded, which is what makes a speed change reviewable. With a directory it also writes that
# corpus as a BAM file and times reading and rewriting the whole file, which htsjdk can do too.
echo "== the record codec, in the same container"
docker run --rm --platform linux/amd64 -v "$OUT":/out -w /out htsjdk-rs-oracle:4.2.0 \
  "/out/amd64/release/record-bench 400000 $REPS /out" | tee "$OUT/records.txt"

# The VCF codec, #78's third path. The two write rows are the two the library has rather than two
# implementations of one: a record nothing looked at is written from the file's own text, and a
# record whose genotypes were read is rebuilt column by column. Every tool that copies records pays
# the first and every tool that inspects them pays the second.
echo "== the VCF codec, in the same container"
docker run --rm --platform linux/amd64 -v "$OUT":/out -w /out htsjdk-rs-oracle:4.2.0 \
  "/out/amd64/release/vcf-bench 200000 $REPS 4 /out" | tee "$OUT/vcf.txt"

# htsjdk on the files the two programs above just wrote, writing its own output next to theirs.
echo "== htsjdk's BAM and VCF codecs, on the same files"
docker run --rm --platform linux/amd64 -v "$PWD/tools/benchmark":/harness:ro -v "$OUT":/out \
  -w /work htsjdk-rs-oracle:4.2.0 \
  "cp /harness/CodecBench.java . && javac -cp \$ORACLE_CP -d . CodecBench.java \
   && java -cp .:\$ORACLE_CP CodecBench /out $REPS" | tee "$OUT/java-codec.txt"

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

python3 - "$OUT" <<'FILES'
import sys, re, pathlib, hashlib, gzip

out = pathlib.Path(sys.argv[1])
text = "".join((out / f).read_text() for f in ("records.txt", "vcf.txt", "java-codec.txt"))
counts = dict(re.findall(r"(\w+_records)=(\d+)", text))
for kind in ("bam", "vcf"):
    if counts.get(f"{kind}_records") != counts.get(f"java_{kind}_records"):
        raise SystemExit(f"{kind}: the port read {counts.get(f'{kind}_records')} records, "
                         f"htsjdk {counts.get(f'java_{kind}_records')}")
bad = False
for port, java in (("port_out.bam", "java_out.bam"), ("port_copy.vcf", "java_copy.vcf"),
                   ("port_touched.vcf", "java_touched.vcf")):
    a, b = (out / port).read_bytes(), (out / java).read_bytes()
    if a == b:
        print(f"{port} == {java}: {len(a)} bytes, md5={hashlib.md5(a).hexdigest()}")
        continue
    bad = True
    print(f"{port} != {java}: {len(a)} vs {len(b)} bytes")
    if port.endswith(".bam"):
        same = gzip.decompress(a) == gzip.decompress(b)
        print(f"  inflated: {'identical, so the deflater differs' if same else 'different too'}")
if bad:
    raise SystemExit("the two sides wrote different files; the timing below means nothing")
FILES

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

echo "== the record codec in memory, MB/s and records/s, median of $REPS"
python3 - "$OUT" <<'REC'
import sys, re, pathlib, statistics

out = pathlib.Path(sys.argv[1])
rows = {}
for line in (out / "records.txt").read_text().splitlines():
    m = re.match(r"rust_record_(decode|encode)_run\d+_mbps=([\d.]+) recs_per_sec=(\d+)", line)
    if m:
        rows.setdefault(m.group(1), []).append((float(m.group(2)), int(m.group(3))))
    m = re.match(r"records=(\d+) encoded_bytes=(\d+) md5=([0-9a-f]{32})", line)
    if m:
        print(f"corpus: {m.group(1)} records, {m.group(2)} bytes, md5={m.group(3)}")
for name in ("decode", "encode"):
    if name in rows:
        mbps = statistics.median(r[0] for r in rows[name])
        recs = statistics.median(r[1] for r in rows[name])
        print(f"{name:8} {mbps:>9.1f} MB/s {recs:>12,.0f} records/s")
REC

echo "== BAM and VCF files, htsjdk against the port, median of $REPS"
python3 - "$OUT" <<'CODEC'
import sys, re, pathlib, statistics

out = pathlib.Path(sys.argv[1])
text = "".join((out / f).read_text() for f in ("records.txt", "vcf.txt", "java-codec.txt"))
n = {k: int(v) for k, v in
     re.findall(r"\b(bam_records|bam_inflated_bytes|vcf_records|vcf_text_bytes)=(\d+)", text)}
secs = {}
for side, kind, op, value in re.findall(r"(java|rust)_(bam|vcf)_(\w+?)_run\d+_secs=([\d.]+)", text):
    secs.setdefault((kind, op), {}).setdefault(side, []).append(float(value))
# Every rate is over one denominator for both sides: the inflated BAM stream, the VCF both read, or
# the VCF both wrote (identical, checked above).
size = {
    ("bam", "decode"): n["bam_inflated_bytes"],
    ("bam", "encode"): n["bam_inflated_bytes"],
    ("vcf", "read_lazy"): n["vcf_text_bytes"],
    ("vcf", "read_decoded"): n["vcf_text_bytes"],
    ("vcf", "write_copy"): (out / "port_copy.vcf").stat().st_size,
    ("vcf", "write_touched"): (out / "port_touched.vcf").stat().st_size,
}
print(f"{'path':18} {'htsjdk MB/s':>12} {'port MB/s':>10} {'htsjdk rec/s':>13} {'port rec/s':>12} {'ratio':>7}")
for key, nbytes in size.items():
    row = secs.get(key, {})
    if "java" not in row or "rust" not in row:
        continue
    recs = n[f"{key[0]}_records"]
    j, r = statistics.median(row["java"]), statistics.median(row["rust"])
    mb = nbytes / (1024 * 1024)
    print(f"{key[0] + ' ' + key[1]:18} {mb / j:>12.1f} {mb / r:>10.1f} "
          f"{recs / j:>13,.0f} {recs / r:>12,.0f} {j / r:>6.2f}x")
CODEC
