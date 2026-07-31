#!/bin/bash
#
# submit-slab-ks-wrap.sh — the full autonomous chain for the CAP-free
# wrap-around slab KS-orbital stopping study.
# Plan: docs/plans/slab-ks-orbital-stopping-wrap.md
#
#   Usage:  cd <repo root> && ./shared/bin/submit-slab-ks-wrap.sh
#
# Chain (each stage waits on the one it depends on; the two arrays run their
# eight points CONCURRENTLY, which is the "run in parallel" the user asked for):
#
#   1  slabks-tests        library correctness gate (slab occupancy + minimum-image
#                          Gaussian + projectile wrap). afterok: nothing runs if
#                          the kernels the study rests on are wrong.
#   2  slabks-gs N=40      new r_s = 5.67 ground state         afterok 1
#   3  slabks-gs N=100     r_s = 4.18 — SKIPS, already exists  afterok 1
#   4  slabks-wp smoke     builds the WP binary + t=0 gates    afterok 1
#   5  slabks-cl smoke     builds the classical binary         afterok 1
#   6  slabks-wp  0-7      8 WP points, concurrent             afterok 4,2,3
#   7  slabks-cl  0-7      8 classical twins, concurrent       afterok 5,2,3
#
# Stages 6 and 7 are independent of each other, so both arrays queue at once:
# 16 GPU tasks total, ~29 GPU-hours, ~3.5 h wall if the partition gives them all
# a slot.
#
# The matrix (density x velocity), identical for both halves:
#
#   idx  N    r_s    v      N_steps   t (a.u.)   path (Bohr)   plasma periods
#     0  100  4.18   2.0     4529      181.2        362.3          5.8
#     1  100  4.18   2.5     3623      144.9        362.3          4.7
#     2  100  4.18   3.0     3020      120.8        362.4          3.9
#     3  100  4.18   3.5     2588      103.5        362.3          3.3
#     4   40  5.67   2.0     4529      181.2        362.3          3.7
#     5   40  5.67   2.5     3623      144.9        362.3          3.0
#     6   40  5.67   3.0     3020      120.8        362.4          2.5
#     7   40  5.67   3.5     2588      103.5        362.3          2.1
#
# Resume after a kill: re-submit the single index with LJ_RESUME=1, e.g.
#   LJ_RESUME=1 sbatch shared/bin/run-slab-ks-wp.slurm 0
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

[ -f shared/bin/csd3-env.sh ] || { echo "FATAL: run from the repo root." >&2; exit 2; }

echo "repo: $REPO_ROOT"
echo

j_tests=$(sbatch --parsable shared/bin/run-slab-ks-tests.slurm)
echo "1  slabks-tests           $j_tests"

j_gs40=$(sbatch --parsable --dependency=afterok:$j_tests shared/bin/run-slab-ks-gs.slurm 40)
echo "2  slabks-gs N=40         $j_gs40   (afterok $j_tests)"

j_gs100=$(sbatch --parsable --dependency=afterok:$j_tests shared/bin/run-slab-ks-gs.slurm 100)
echo "3  slabks-gs N=100        $j_gs100   (afterok $j_tests; skips if present)"

# The smokes BUILD the binaries and run 20 steps against the N=100 ground state,
# which already exists — so they depend on the test gate only, NOT on the new
# N=40 ground state. Compile errors therefore surface immediately instead of
# after a ground-state run.
j_wpsm=$(sbatch --parsable --dependency=afterok:$j_tests shared/bin/run-slab-ks-wp.slurm smoke)
echo "4  slabks-wp smoke        $j_wpsm   (afterok $j_tests)"

j_clsm=$(sbatch --parsable --dependency=afterok:$j_tests shared/bin/run-slab-ks-classical.slurm smoke)
echo "5  slabks-cl smoke        $j_clsm   (afterok $j_tests)"

j_wp=$(sbatch --parsable --dependency=afterok:$j_wpsm:$j_gs40:$j_gs100 --array=0-7 shared/bin/run-slab-ks-wp.slurm)
echo "6  slabks-wp array 0-7    $j_wp   (afterok $j_wpsm,$j_gs40,$j_gs100)"

j_cl=$(sbatch --parsable --dependency=afterok:$j_clsm:$j_gs40:$j_gs100 --array=0-7 shared/bin/run-slab-ks-classical.slurm)
echo "7  slabks-cl array 0-7    $j_cl   (afterok $j_clsm,$j_gs40,$j_gs100)"

echo
echo "submitted. watch with:  squeue -u \$USER -o '%.10i %.14j %.8T %.10M %R'"
