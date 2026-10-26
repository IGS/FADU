#!/usr/bin/env bash
# Regression test: run FADU on the tiny dataset and compare against expected output.
# Override the FADU command with FADU_CMD (e.g. FADU_CMD="fadu.jl" for an installed copy).
set -euo pipefail

DATA_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/data"
ROOT_DIR="$(cd "${DATA_DIR}/../.." && pwd)"
FADU_CMD="${FADU_CMD:-julia ${ROOT_DIR}/fadu.jl}"

OUT_DIR="$(mktemp -d)"
trap 'rm -rf "${OUT_DIR}"' EXIT

${FADU_CMD} -g "${DATA_DIR}/tiny.gff3" -b "${DATA_DIR}/tiny.bam" -o "${OUT_DIR}" -s "reverse" -f "CDS" -a "ID"

if diff "${DATA_DIR}/expected.counts.txt" "${OUT_DIR}/tiny.counts.txt"; then
    echo "PASS: output matches expected.counts.txt"
else
    echo "FAIL: output differs from expected.counts.txt" >&2
    exit 1
fi
