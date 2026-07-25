#!/usr/bin/env bash
# Rebuild (inq-study) + run the two identical vacuum WP runs with the CORRECTED
# 80-Bohr geometry (WP launches 10 sigma clear of the one-sided +z CAP / wrapped
# -z wall). no-CAP then CAP. Superseded old LZ=60/launch=-26 results are removed.
set -euo pipefail
cd "$(dirname "$0")"
export INQ_SOURCE=/local/data/public/skcb2/tddft/inq-study
export CUDA_VISIBLE_DEVICES=0
PATH="/local/data/public/skcb2/tddft/shared/bin:$PATH"

echo "== removing superseded results =="
rm -rf results/nocap results/cap results/comparison

echo "== build + run no-CAP (WP_ETA=0) =="
WP_OUT=nocap WP_ETA=0 inq-run

echo "== run CAP (WP_ETA=-0.7) with the same binary =="
env INQ_SHARE_PATH="${INQ_SHARE_PATH:-/local/data/public/skcb2/tddft/inq/install/share}" \
    PSEUDOPOD_SHARE_PATH="${PSEUDOPOD_SHARE_PATH:-/local/data/public/skcb2/tddft/inq/install/share/pseudopod}" \
    CUDA_VISIBLE_DEVICES=0 WP_OUT=cap WP_ETA=-0.7 ./run

echo "== both runs done =="
grep -H "launch_z\|cap_z\|clearance\|n_steps\|cell_bohr" results/nocap/run_summary.txt results/cap/run_summary.txt
