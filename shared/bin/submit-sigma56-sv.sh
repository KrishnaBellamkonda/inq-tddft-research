#!/bin/bash
# submit-sigma56-sv.sh — the sigma_WP = 5 and 6 Bohr classical+wavepacket TWIN
# campaign, submitted as one autonomous, dependency-chained set of SLURM jobs.
#
#     ./shared/bin/submit-sigma56-sv.sh
#     SKIP_GS=1 ./shared/bin/submit-sigma56-sv.sh    # GS already exists
#
# WHY THIS CAMPAIGN EXISTS (user instruction 2026-08-02). The existing sigma =
# 0.5/2/3 wavepacket runs have no width-matched classical twin, and their packets
# disperse so much in transit (x3.2 at sigma = 0.5, x2.7 at sigma = 2) that a
# single sigma label does not describe the object doing the physics. At sigma =
# 5/6 the packet grows by only x1.23 / x1.12 over the in-slab transit, so it is
# effectively constant-width and its LABEL agrees with its TIME-AVERAGE — which is
# the condition under which a classical projectile at fixed sigma_pot is a fair
# comparison at all. The question the campaign answers: at what width do the
# classical and quantum projectiles stop being distinguishable?
#
# WHAT CHANGES vs the 85-Bohr campaigns: L_z 85 -> 105 (pure vacuum; slab, N and
# r_s untouched) and launch z -24 -> -27.5, both forced by geometry — a sigma = 6
# packet has a density std of 4.243 Bohr and the old box could not hold it clear
# of both the absorber and the slab. Everything else is held.
#
# BOTH HALVES CARRY THE CAP (user decision 2026-08-02) so E_absorbed/L_slab is the
# same estimator on both, plus one CAP-free classical control per sigma at v = 3.0
# to measure what the absorber costs.
#
# THE OBSERVABLE: S = (E_total(t_final) - E_GS) / L_slab_z, L_slab_z = 25 Bohr,
# norm-corrected on the WP half.
#
# THE CHAIN
#   1  gs           ONE ground state serves all 18 runs (sigma does not enter the
#                   bath, and both halves load the same checkpoint). afterok-gates
#                   everything: a bad GS must not spawn 60 GPU-hours.
#   2  wp smoke s6  BUILDS the WP binary + runs the t=0 analytic gates.
#   3  cl smoke s6  BUILDS the classical binary (20 steps).
#   4  wp sweep s6  array 0-3, v = 2.0/2.5/3.0/3.5 in PARALLEL   [afterok 2]
#   5  cl sweep s6  array 0-3, the width-matched twins           [afterok 3]
#   6  vac s6       4 CAP-only baselines                          [afterany 4]
#   7  cl nocap s6  v = 3.0 CAP-free control                      [afterany 5]
#   8  wp smoke s5  t=0 gates at the other width (binary already built)
#   9  wp sweep s5  array 0-3                                     [afterok 8]
#  10  cl sweep s5  array 0-3                                     [afterany 5]
#  11  vac s5       4 CAP-only baselines                          [afterany 6]
#  12  cl nocap s5  v = 3.0                                       [afterany 7]
#  13  finalize     repair any short/missing run, then build figures + notebooks
#                   + CAMPAIGN_REPORT.md + email          [afterany 9,10,11,12]
#  14  finalize     second attempt, in case 13 hit its walltime mid-repair
#
# COST — WARNING, NOT A GATE (.claude/rules/checkpoint-dont-block.md).
# MEASURED on the smoke stages (2026-08-02): 3.15 s/step (WP) and 3.00 s/step
# (classical) on an A100, against 2.75 s/step at L_z = 85 — the grid grew
# 88x88x213 -> 88x88x264 (x1.24). Per sigma that is 13246 steps per half, but the
# four velocities run as a PARALLEL array, so wall clock is the longest single
# run: 4360 x 3.15 s = 3.8 h. Whole campaign ~10-14 h wall, ~60 GPU-h.
# EVERY run checkpoints every N/5 steps and supports LJ_RESUME=1, so killing one
# costs at most one interval:
#     scancel <jobid>
#     sbatch --export=ALL,LJ_SIGMA=6.0,LJ_RESUME=1 shared/bin/run-s56-wp.slurm 0
#
# Plan: docs/plans/sigma56-sv-twin.md
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

GS_DIR="$REPO_ROOT/ResearchProject/systems/localised_jellium/shared_gs/slab_n100_L35x35x105_dx0p4_per2"

sub() { sbatch --parsable "$@"; }

echo "=== sigma56_sv twin campaign — submitting the chain ==="
echo "repo: $REPO_ROOT"
echo

# ---- 1. ground state -------------------------------------------------------
if [ "${SKIP_GS:-0}" = "1" ] || [ -d "$GS_DIR" ]; then
  if [ ! -d "$GS_DIR" ]; then
    echo "FATAL: SKIP_GS=1 but no ground state at $GS_DIR" >&2; exit 2
  fi
  echo "  1  gs           SKIPPED (already at $GS_DIR)"
  GS_DEP=""
else
  J_GS=$(sub shared/bin/run-s56-gs.slurm 0.4)
  echo "  1  gs           $J_GS"
  GS_DEP="--dependency=afterok:$J_GS"
fi

# ---- sigma = 6 -------------------------------------------------------------
J_WS6=$(sub $GS_DEP --export=ALL,LJ_SIGMA=6.0 shared/bin/run-s56-wp.slurm smoke)
echo "  2  wp smoke s6   $J_WS6"
J_CS6=$(sub $GS_DEP --export=ALL,LJ_SIGMA=6.0 shared/bin/run-s56-cl.slurm smoke)
echo "  3  cl smoke s6   $J_CS6"

J_W6=$(sub --dependency=afterok:$J_WS6 --export=ALL,LJ_SIGMA=6.0 \
           --array=0-3 shared/bin/run-s56-wp.slurm)
echo "  4  wp sweep s6   $J_W6  (array 0-3)"
J_C6=$(sub --dependency=afterok:$J_CS6 --export=ALL,LJ_SIGMA=6.0 \
           --array=0-3 shared/bin/run-s56-cl.slurm)
echo "  5  cl sweep s6   $J_C6  (array 0-3)"

J_V6=$(sub --dependency=afterany:$J_W6 --export=ALL,LJ_SIGMA=6.0 \
           shared/bin/run-s56-vac.slurm)
echo "  6  vac s6        $J_V6"
J_N6=$(sub --dependency=afterany:$J_C6 --export=ALL,LJ_SIGMA=6.0,LJ_CAP_ETA=0 \
           shared/bin/run-s56-cl.slurm 2)
echo "  7  cl nocap s6   $J_N6  (v = 3.0)"

# ---- sigma = 5 -------------------------------------------------------------
# Queued behind the sigma = 6 sweeps rather than alongside them, so the campaign
# never asks for more than ~8 GPUs at once.
J_WS5=$(sub --dependency=afterany:$J_W6 --export=ALL,LJ_SIGMA=5.0 \
            shared/bin/run-s56-wp.slurm smoke)
echo "  8  wp smoke s5   $J_WS5"
J_W5=$(sub --dependency=afterok:$J_WS5 --export=ALL,LJ_SIGMA=5.0 \
           --array=0-3 shared/bin/run-s56-wp.slurm)
echo "  9  wp sweep s5   $J_W5  (array 0-3)"
J_C5=$(sub --dependency=afterany:$J_C6 --export=ALL,LJ_SIGMA=5.0 \
           --array=0-3 shared/bin/run-s56-cl.slurm)
echo " 10  cl sweep s5   $J_C5  (array 0-3)"
J_V5=$(sub --dependency=afterany:$J_V6 --export=ALL,LJ_SIGMA=5.0 \
           shared/bin/run-s56-vac.slurm)
echo " 11  vac s5        $J_V5"
J_N5=$(sub --dependency=afterany:$J_N6 --export=ALL,LJ_SIGMA=5.0,LJ_CAP_ETA=0 \
           shared/bin/run-s56-cl.slurm 2)
echo " 12  cl nocap s5   $J_N5  (v = 3.0)"

# ---- finalize: repair, build, report ---------------------------------------
# NOT a one-shot notebook stage. A SLURM chain gets runs LAUNCHED, not FINISHED:
# a walltime kill, a preempted node or a dependency that never fired all leave a
# short run, and a one-shot post-processor would happily build a figure out of
# whatever happened to be on disk. run-s56-finalize.slurm checks every expected
# run against its step target, RESUMES the short ones in place (by invoking these
# same dispatchers as plain bash — no sbatch from inside a job), and only then
# builds the figures and notebooks. It writes CAMPAIGN_REPORT.md and attempts an
# email on EVERY path, including total failure.
#
# TWO are chained: if the first hits its 36 h walltime mid-repair, the second
# picks up from the checkpoints it left behind. Both are bounded, so neither can
# idle a GPU waiting for work that is never coming.
J_F1=$(sub --dependency=afterany:$J_W5:$J_C5:$J_V5:$J_N5 \
           shared/bin/run-s56-finalize.slurm 1 2)
echo " 13  finalize (1)  $J_F1"
J_F2=$(sub --dependency=afterany:$J_F1 shared/bin/run-s56-finalize.slurm 2 2)
echo " 14  finalize (2)  $J_F2"

echo
echo "Kill the whole campaign with:"
echo "  scancel ${J_GS:-} $J_WS6 $J_CS6 $J_W6 $J_C6 $J_V6 $J_N6 $J_WS5 $J_W5 $J_C5 $J_V5 $J_N5 $J_F1 $J_F2"
echo
echo "Check status at any time WITHOUT changing anything:"
echo "  cd ResearchProject/systems/localised_jellium/hypotheses/sigma56_sv"
echo "  python finalize.py --status-only"
echo
echo "Watch with:  squeue -u \$USER -o '%.10i %.12j %.8T %.10M %R'"
