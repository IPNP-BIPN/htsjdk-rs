#!/usr/bin/env bash
# Run a Java harness inside the pinned oracle image and capture provenance alongside its output.
#
# The contract check runs first and its failure is fatal. A golden produced in an unverified
# environment is worse than no golden: it looks exactly like a good one.
#
# Usage:
#   ./run.sh <harness.java> [output-file]
#
# Example:
#   ./run.sh ../zlib-conformance/Z2.java goldens/zlib.txt

set -euo pipefail

IMAGE="${ORACLE_IMAGE:-htsjdk-rs-oracle:4.2.0}"
PLATFORM=linux/amd64

if [ $# -lt 1 ]; then
    sed -n '2,12p' "$0" | sed 's/^# \{0,1\}//'
    exit 64
fi

harness_path=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
harness_dir=$(dirname "$harness_path")
harness_file=$(basename "$harness_path")
harness_class="${harness_file%.java}"
output="${2:-}"

if ! docker image inspect "$IMAGE" >/dev/null 2>&1; then
    echo "oracle image '$IMAGE' not found; build it with:" >&2
    echo "  docker build --platform $PLATFORM -t $IMAGE $(dirname "$0")" >&2
    exit 69
fi

# 1. Contract check. Provenance goes to stderr so stdout stays clean for the harness output.
provenance=$(docker run --rm --platform "$PLATFORM" "$IMAGE" \
    'java -cp "$ORACLE_CP" OracleProbe' 2>/dev/null) || {
    echo "ORACLE CONTRACT VIOLATED; refusing to generate goldens." >&2
    docker run --rm --platform "$PLATFORM" "$IMAGE" 'java -cp "$ORACLE_CP" OracleProbe' >/dev/null || true
    exit 3
}

image_digest=$(docker image inspect "$IMAGE" --format '{{.Id}}')

{
    echo "# oracle provenance"
    echo "# image: $IMAGE"
    echo "# image_id: $image_digest"
    echo "# platform: $PLATFORM"
    echo "# harness: $harness_file"
    echo "# generated: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "$provenance" | sed 's/^/# /'
} >&2

# 2. Run the harness.
result=$(docker run --rm --platform "$PLATFORM" \
    -v "$harness_dir":/harness:ro -w /work "$IMAGE" \
    "cp /harness/$harness_file . && javac -cp \"\$ORACLE_CP\" -d . $harness_file \
     && java -cp \".:\$ORACLE_CP\" $harness_class" 2>/dev/null)

if [ -n "$output" ]; then
    mkdir -p "$(dirname "$output")"
    printf '%s\n' "$result" > "$output"
    echo "wrote $(printf '%s\n' "$result" | wc -l | tr -d ' ') lines to $output" >&2
else
    printf '%s\n' "$result"
fi
