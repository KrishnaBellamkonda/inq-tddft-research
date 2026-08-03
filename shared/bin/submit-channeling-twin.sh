#!/bin/bash
#
# submit-channeling-twin.sh — the full autonomous chain for the annular-tube
# CHANNELING TWIN: a matched classical + wavepacket pair shot down the bore of a
# periodic r_s = 3 jellium tube, to validate a KS-orbital definition of stopping
# power against the classical deltaE/ds one.
#
# Plan: docs/plans/cylindrical-channeling-ks-stopping.md
#
#   Usage:  cd <repo root> && ./shared/bin/submit-channeling-twin.sh
#           ./shared/bin/submit-channeling-twin.sh --no-tests   # skip the gate
#
# Chain (each stage waits on what it depends on; the two production halves run
# CONCURRENTLY, which is the whole point of a twin):
#
#   1  chan-tests        library gate: radial occupancy, minimum-image charge
#                        AND force, projectile wrap. Nothing runs if these fail.
#   2  chan-gs           the SHARED r_s = 3 tube ground state    afterok 1
#                        (skips instantly if the checkpoint exists)
#   3  chan-twin wp smoke        builds the WP binary + t=0 gates    afterok 1
#   4  chan-twin classical smoke builds the classical binary         afterok 1
#   5  chan-twin wp              1500-step production               afterok 3,2
#   6  chan-twin classical       1500-step production               afterok 4,2
#   7  chan-nb            twin parity gate -> 2 run notebooks -> the
#                         COMPARISON notebook                       afterany 5,6
#
# WHY THE SMOKES DEPEND ONLY ON THE TESTS. They build the binaries and would
# surface a compile error immediately; they need the ground state only to run, so
# they are ordered after stage 1 and the production halves wait on stage 2 as
# well. If the ground state does not exist yet the smokes exit 2 with a clear
# message and the chain stops there rather than after an 8-hour SCF.
#
# COST. ~2-4 GPU-hours per production half (768k grid, 104 states, 1500 steps),
# both concurrent, plus up to 8 h for the ground state. Everything is
# checkpointed every 500 steps and resumable
# (.claude/rules/checkpoint-dont-block.md), so a kill costs at most one interval:
#
#   kill:    scancel <jobid>
#   resume:  CH_RESUME=1 sbatch shared/bin/run-chan-twin.slurm wp
#   extend:  CH_RESUME=1 CH_N_STEPS=3000 sbatch shared/bin/run-chan-twin.slurm wp
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

[ -f shared/bin/csd3-env.sh ] || { echo "FATAL: run from the repo root." >&2; exit 2; }

RUN_TESTS=1
[ "${1:-}" = "--no-tests" ] && RUN_TESTS=0

echo "repo: $REPO_ROOT"
echo

if [ "$RUN_TESTS" -eq 1 ]; then
  j_tests=$(sbatch --parsable shared/bin/run-chan-tests.slurm)
  echo "1  chan-tests             $j_tests"
  DEP_TESTS="--dependency=afterok:$j_tests"
else
  echo "1  chan-tests             SKIPPED (--no-tests)"
  DEP_TESTS=""
fi

# shellcheck disable=SC2086
j_gs=$(sbatch --parsable $DEP_TESTS shared/bin/run-chan-gs.slurm)
echo "2  chan-gs                 $j_gs"

# --job-name PER HALF, not the script's default. Both halves come from ONE
# script with the half as an argument, and neither squeue nor `scontrol show job`
# displays a job's arguments — so without this override every twin job appears
# in the queue as an indistinguishable "chan-twin" and there is no way to tell
# whether the classical half was submitted at all. (It always was; you just
# could not see it. User raised exactly this, 2026-08-01.) The name also feeds
# %x in the script's --output pattern, so the log files separate too.

# BUILD runs CONCURRENTLY WITH THE GROUND STATE, not after it. It needs no
# checkpoint, so a compile error surfaces in minutes instead of after an
# 8-hour SCF. The first version merged build into the smoke and put it before the
# ground state existed — the script's own GS guard then killed it in 5 seconds
# WITHOUT COMPILING, and its afterok cancelled both production halves.
# shellcheck disable=SC2086
j_wpb=$(sbatch --parsable --job-name=chan-wp-build $DEP_TESTS shared/bin/run-chan-twin.slurm wp build)
echo "3  chan-wp-build           $j_wpb   (compiles; no GS needed)"

# shellcheck disable=SC2086
j_clb=$(sbatch --parsable --job-name=chan-cl-build $DEP_TESTS shared/bin/run-chan-twin.slurm classical build)
echo "4  chan-cl-build           $j_clb   (compiles; no GS needed)"

j_wpsm=$(sbatch --parsable --job-name=chan-wp-smoke \
         --dependency=afterok:$j_wpb:$j_gs shared/bin/run-chan-twin.slurm wp smoke)
echo "5  chan-wp-smoke           $j_wpsm   (afterok $j_wpb,$j_gs)"

j_clsm=$(sbatch --parsable --job-name=chan-cl-smoke \
         --dependency=afterok:$j_clb:$j_gs shared/bin/run-chan-twin.slurm classical smoke)
echo "6  chan-cl-smoke           $j_clsm   (afterok $j_clb,$j_gs)"

j_wp=$(sbatch --parsable --job-name=chan-wp \
       --dependency=afterok:$j_wpsm shared/bin/run-chan-twin.slurm wp)
echo "7  chan-wp  (wavepacket)   $j_wp   (afterok $j_wpsm)"

j_cl=$(sbatch --parsable --job-name=chan-cl \
       --dependency=afterok:$j_clsm shared/bin/run-chan-twin.slurm classical)
echo "8  chan-cl  (classical)    $j_cl   (afterok $j_clsm)"

j_nb=$(sbatch --parsable --dependency=afterany:$j_wp:$j_cl shared/bin/run-chan-notebooks.slurm)
echo "7  chan-nb                 $j_nb   (afterany $j_wp,$j_cl)"

echo
echo "submitted. watch with:  squeue -u \$USER -o '%.10i %.14j %.8T %.10M %R'"
