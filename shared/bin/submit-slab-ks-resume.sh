#!/bin/bash
#
# submit-slab-ks-resume.sh — resume the 11 slab-KS runs that were killed when
# /rds hit 100 % on 2026-07-31.
# Plan: docs/plans/slab-ks-orbital-stopping-wrap.md
# Handover: docs/handovers/slab-ks-orbital-stopping-wrap.md
#
# WHAT HAPPENED. 16 concurrent runs writing five VTI streams each filled the
# 1099.5 GB project quota. 11 died mid-flight with
#   VTIImageDataWriter: failed while writing file: ...density_delta_tNNNNNN.vti
# (SIGABRT, exit 6) after 1-2 h. No physics or code fault — pure disk.
#
# WHAT CHANGED BEFORE THIS RESUME (user instruction 2026-07-31):
#   * density frames  301 -> ~119 per run   (LJ_SAVE_EVERY 15/12/10/9 -> 38/30/25/22)
#   * wavefunctions   100 -> ~19 per run    (LJ_WF_EVERY   45/36/30/27 -> 228/180/150/132)
#   * retained checkpoints 5 -> 3           (LJ_CKPT_EVERY = N/3, LJ_MAX_CKPT = 3,
#                                            oldest pruned as new ones are written)
#   * the FINAL step is now ALWAYS written as a step-stamped ckpt_step<N_STEPS>
#     in addition to the rolling `checkpoint`, and sorts last so pruning can
#     never remove it.
#   * freed 347 GB: deleted density_delta (exactly reconstructible from
#     density_total) and the interior checkpoints of the 5 COMPLETED runs.
#
# THROTTLING. The arrays carry %4, and the classical array waits on the WP one,
# so at most ~4 runs write concurrently instead of 16 — roughly a quarter of the
# peak footprint that blew the quota. Slower, but the failure cannot repeat.
#
# THE 11 RUNS (the other 5 are COMPLETE and are not touched):
#
#   half  idx  run          reached   target   resumes from
#   wp     0   n100_v2p0     2774      4529      2715
#   wp     1   n100_v2p5     2760      3623      2172
#   wp     2   n100_v3p0     2719      3020      2416
#   wp     4   n40_v2p0      4365      4529      3620
#   cl     0   n100_v2p0     2040      4529      1810
#   cl     1   n100_v2p5     1732      3623      1448
#   cl     2   n100_v3p0     1730      3020      1208
#   cl     3   n100_v3p5     1710      2588      1551
#   cl     4   n40_v2p0      2760      4529      2715
#   cl     5   n40_v2p5      2748      3623      2172
#   cl     6   n40_v3p0      2730      3020      2416
#
# Both smokes run FIRST because run.cpp changed (checkpoint policy) and the
# array tasks exec ./run directly — without a rebuild they would silently run
# the old binary.
#
#   Usage:  cd <repo root> && ./shared/bin/submit-slab-ks-resume.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"
[ -f shared/bin/csd3-env.sh ] || { echo "FATAL: run from the repo root." >&2; exit 2; }

# Refuse to launch into a nearly-full quota — that is what caused the failure.
avail=$(quota -s 2>/dev/null | awk '/user\/skcb2/ {gsub(/[!*]/,"",$2); print $3-$2; exit}')
echo "quota headroom: ${avail:-unknown} GB"
if [ -n "${avail:-}" ] && [ "$(printf '%.0f' "$avail")" -lt 150 ]; then
  echo "FATAL: under 150 GB free — free space before resuming." >&2
  exit 3
fi
echo

# Rebuild both binaries (run.cpp changed), then resume, throttled to 4 at a time.
j_wpsm=$(sbatch --parsable shared/bin/run-slab-ks-wp.slurm smoke)
echo "1  slabks-wp smoke (rebuild)    $j_wpsm"

j_clsm=$(sbatch --parsable shared/bin/run-slab-ks-classical.slurm smoke)
echo "2  slabks-cl smoke (rebuild)    $j_clsm"

j_wp=$(sbatch --parsable --dependency=afterok:$j_wpsm \
       --export=ALL,LJ_RESUME=1 --array=0,1,2,4%4 shared/bin/run-slab-ks-wp.slurm)
echo "3  slabks-wp resume 0,1,2,4     $j_wp   (afterok $j_wpsm, max 4 concurrent)"

j_cl=$(sbatch --parsable --dependency=afterok:$j_clsm,afterany:$j_wp \
       --export=ALL,LJ_RESUME=1 --array=0-6%4 shared/bin/run-slab-ks-classical.slurm)
echo "4  slabks-cl resume 0-6         $j_cl   (afterok $j_clsm + afterany $j_wp, max 4 concurrent)"

echo
echo "watch:  squeue -u \$USER -o '%.12i %.12j %.9T %.8M %.14R'"
echo "space:  quota -s | grep user/skcb2"
