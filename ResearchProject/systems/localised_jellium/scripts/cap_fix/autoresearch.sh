#!/usr/bin/env bash
# LOCKED EVAL HARNESS — cap_fix campaign (autoresearch skill contract).
# OFF LIMITS to experiments: this script and run_metrics.py define the
# benchmark and the metric. Experiments vary ONLY the CAP config env vars.
#
# usage:  [EM_CAP_MODE=..] [EM_CAP_ETA=..] [EM_CAP_CENTER_BOHR=..]
#         [EM_CAP_WIDTH_BOHR=..] [EM_WRAP_WIDTH_BOHR=..] [EM_N_STEPS=..]
#         ./autoresearch.sh <run_name> <gpu_id>
#
# The benchmark is the FIXED witness configuration of the diagnosis campaign:
# slab_n52 GS + sigma1 WP (mass fork), ETRS, dt=0.04, 700 steps (t=28 a.u.),
# write_every=5 — the exact setup in which the artifact was reproduced
# (capon_weak_partial). Deterministic (no RNG in the propagator), so the noise
# floor is the propagator conservation floor, not seed variance.
set -euo pipefail

RUN_NAME="${1:?usage: autoresearch.sh <run_name> <gpu_id>}"
GPU_ID="${2:?usage: autoresearch.sh <run_name> <gpu_id>}"

cd "$(dirname "$0")"

export INQ_SOURCE=/local/data/public/skcb2/tddft/inq-study
export CUDA_VISIBLE_DEVICES="${GPU_ID}"
export EM_OUT="${RUN_NAME}"
export EM_N_STEPS="${EM_N_STEPS:-700}"
export EM_WRITE_EVERY="${EM_WRITE_EVERY:-5}"

# The binary is PRE-BUILT (inq-run) and run.cpp is frozen for the campaign —
# execute it directly. Rationale (2026-07-13): two concurrent inq-run calls
# race in the shared build/ dir and corrupt the device link; experiments are
# env-var configs, so no per-run rebuild is ever needed. If run.cpp ever
# changes (deliberate re-baseline), rebuild once with inq-run BEFORE looping.
[[ ./run -nt ./run.cpp ]] || { echo "FATAL: ./run stale vs run.cpp — rebuild with inq-run first" >&2; exit 3; }

# per-experiment wall-clock cap: 2 h (a 700-step run takes ~55 min)
timeout 7200 env \
    INQ_SHARE_PATH="${INQ_SHARE_PATH:-$INQ_SOURCE/install/share}" \
    PSEUDOPOD_SHARE_PATH="${PSEUDOPOD_SHARE_PATH:-$INQ_SOURCE/install/share/pseudopod}" \
    ./run

/local/data/public/skcb2/tddft/venv/bin/python3 run_metrics.py "results/${RUN_NAME}"
