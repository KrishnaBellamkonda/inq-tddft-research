#!/bin/bash
# submit-lz-bulk-sweep.sh — the slab->bulk L_slab sweep, submitted as one
# autonomous, dependency-chained, PILOT-GATED set of SLURM jobs.
#
#     ./shared/bin/submit-lz-bulk-sweep.sh
#     LZB_EXCLUDE=none ./shared/bin/submit-lz-bulk-sweep.sh   # no node excludes
#
# Plan: docs/plans/jellium-slab-extend-Lz.md
# Handover: docs/handovers/jellium-slab-extend-Lz.md
#
# WHY A PILOT (user decision 2026-08-05): "do one or two runs of velocity for
# all the Lz. Check if everything is alright before committing the massive
# number of GPU hours." One velocity (v = 3.0) runs at ALL four boxes on both
# halves (~20 GPU-h); an automated gate then checks completeness, ledger
# closure and S sanity against the L = 25 anchors. Only on PASS does SLURM
# release the remaining ~70 GPU-h (v = 2.0/2.5/3.5). A gate failure leaves the
# production arrays DependencyNeverSatisfied and writes PILOT_REPORT.md + email
# saying why — autonomous either way, silent in neither.
#
# THE CHAIN
#   1     gs s0p5_L15   builds the gs binary, then runs (smallest box first)
#   2-4   gs x3         the other boxes, exec mode          [afterok 1]
#   5     smoke         builds wp+cl binaries, 8 x t=0 gates [afterok 1-4]
#   6-13  pilot         v = 3.0: {4 boxes} x {wp, cl}        [afterok 5]
#   14    pilot gate    python pilot_gate.py                 [afterany 6-13]
#   15    wp production array 0-11%4  (v = 2.0/2.5/3.5)      [afterok 14]
#   16    cl production array 0-11%4                         [afterok 14]
#   17-20 vac x4        CAP-only baselines per box           [afterok 14]
#   21    finalize (1)  repair short runs, figures, report   [afterany 15-20]
#   22    finalize (2)  second attempt                       [afterany 21]
#
# COST (WARN, not a gate): pilot ~20 GPU-h, production ~70 GPU-h, projected
# from states x grid scaling of the measured 3.15 s/step at the 105-box; the
# pilot gate re-projects from MEASURED s/step and emails the number before
# production starts. Every run checkpoints every N/5 steps and supports
# LJ_RESUME=1; killing one costs at most one interval:
#     scancel <jobid>
#     sbatch --export=ALL,LZB_CFG=<preset>,LJ_RESUME=1 shared/bin/run-lzb-wp.slurm <idx>
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

LJ_BASE="$REPO_ROOT/ResearchProject/systems/localised_jellium"
CFGS=(s0p5_L15 s0p5_L35 s5p0_L15 s5p0_L35)
GS_TAGS=(slab_n60_L35x35x75_dx0p4_per2 slab_n140_L35x35x95_dx0p4_per2 \
         slab_n60_L35x35x95_dx0p4_per2 slab_n140_L35x35x115_dx0p4_per2)

# Node excludes: gpu-q-2 hung large reads and gpu-q-25 hung a checkpoint load
# (2026-08-03, sigma56). Override with LZB_EXCLUDE=none once they are known good.
EXCLUDE="${LZB_EXCLUDE:-gpu-q-2,gpu-q-25}"
if [ "$EXCLUDE" = "none" ]; then
  sub() { sbatch --parsable "$@"; }
else
  sub() { sbatch --parsable --exclude="$EXCLUDE" "$@"; }
fi

echo "=== lz_bulk_sweep — submitting the pilot-gated chain ==="
echo "repo: $REPO_ROOT   excludes: $EXCLUDE"
echo

ALL_JOBS=()

# ---- 1-4. ground states (build-once, then exec) ----------------------------
GS_DEPS=""
BUILDER=""
for i in 0 1 2 3; do
  CFG="${CFGS[$i]}"; TAG="${GS_TAGS[$i]}"
  if [ -d "$LJ_BASE/shared_gs/$TAG" ]; then
    echo "  gs $CFG        SKIPPED (already at shared_gs/$TAG)"
    continue
  fi
  if [ -z "$BUILDER" ]; then
    J=$(sub shared/bin/run-lzb-gs.slurm "$CFG")
    BUILDER="$J"
    echo "  gs $CFG        $J  (builds the binary)"
  else
    J=$(sub --dependency=afterok:$BUILDER shared/bin/run-lzb-gs.slurm "$CFG" exec)
    echo "  gs $CFG        $J  [afterok $BUILDER]"
  fi
  GS_DEPS="$GS_DEPS:$J"
  ALL_JOBS+=("$J")
done

# ---- 5. smoke: builds wp+cl binaries + 8 x t=0 gates -----------------------
if [ -n "$GS_DEPS" ]; then
  J_SM=$(sub --dependency=afterok$GS_DEPS shared/bin/run-lzb-smoke.slurm)
else
  J_SM=$(sub shared/bin/run-lzb-smoke.slurm)
fi
echo "  smoke          $J_SM  (wp+cl builds, 8 smokes)"
ALL_JOBS+=("$J_SM")

# ---- 6-13. pilot: v = 3.0 at every box, both halves ------------------------
PILOT_DEPS=""
for CFG in "${CFGS[@]}"; do
  JW=$(sub --dependency=afterok:$J_SM --export=ALL,LZB_CFG="$CFG" shared/bin/run-lzb-wp.slurm 2)
  echo "  pilot wp $CFG  $JW  (v = 3.0)"
  JC=$(sub --dependency=afterok:$J_SM --export=ALL,LZB_CFG="$CFG" shared/bin/run-lzb-cl.slurm 2)
  echo "  pilot cl $CFG  $JC  (v = 3.0)"
  PILOT_DEPS="$PILOT_DEPS:$JW:$JC"
  ALL_JOBS+=("$JW" "$JC")
done

# ---- 14. the pilot gate ----------------------------------------------------
# afterany, not afterok: a crashed pilot run must still be JUDGED (the gate
# fails on it and blocks production with a report), not silently skipped.
J_GATE=$(sub --dependency=afterany$PILOT_DEPS shared/bin/run-lzb-pilotgate.slurm)
echo "  pilot gate     $J_GATE  (releases production on PASS)"
ALL_JOBS+=("$J_GATE")

# ---- 15-16. production: the remaining three velocities ---------------------
J_WP=$(sub --dependency=afterok:$J_GATE --array=0-11%4 shared/bin/run-lzb-wp.slurm)
echo "  wp production  $J_WP  (array 0-11%4, v = 2.0/2.5/3.5)"
J_CL=$(sub --dependency=afterok:$J_GATE --array=0-11%4 shared/bin/run-lzb-cl.slurm)
echo "  cl production  $J_CL  (array 0-11%4)"
ALL_JOBS+=("$J_WP" "$J_CL")

# ---- 17-20. vacuum baselines ----------------------------------------------
VAC_DEPS=""
for CFG in "${CFGS[@]}"; do
  JV=$(sub --dependency=afterok:$J_GATE --export=ALL,LZB_CFG="$CFG" shared/bin/run-lzb-vac.slurm)
  echo "  vac $CFG       $JV"
  VAC_DEPS="$VAC_DEPS:$JV"
  ALL_JOBS+=("$JV")
done

# ---- 21-22. finalize: repair, build, report --------------------------------
J_F1=$(sub --dependency=afterany:$J_WP:$J_CL$VAC_DEPS shared/bin/run-lzb-finalize.slurm 1 2)
echo "  finalize (1)   $J_F1"
J_F2=$(sub --dependency=afterany:$J_F1 shared/bin/run-lzb-finalize.slurm 2 2)
echo "  finalize (2)   $J_F2"
ALL_JOBS+=("$J_F1" "$J_F2")

echo
echo "Kill the whole campaign with:"
echo "  scancel ${ALL_JOBS[*]}"
echo
echo "Check status at any time WITHOUT changing anything:"
echo "  cd ResearchProject/systems/localised_jellium/hypotheses/lz_bulk_sweep"
echo "  python finalize.py --status-only"
echo
echo "Watch with:  squeue -u \$USER -o '%.10i %.12j %.8T %.10M %R'"
