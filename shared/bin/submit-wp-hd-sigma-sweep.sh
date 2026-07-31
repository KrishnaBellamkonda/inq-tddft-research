#!/bin/bash
# submit-wp-hd-sigma-sweep.sh — the sigma_WP = 2 and 3 Bohr companion campaigns to
# the completed sigma_WP = 0.5 wavepacket sweep, submitted as one autonomous,
# dependency-chained set of SLURM jobs. Run once from anywhere:
#
#     ./shared/bin/submit-wp-hd-sigma-sweep.sh
#
# WHY THESE CAMPAIGNS EXIST (user instruction 2026-07-31). At sigma_WP = 0.5 the
# packet's free-dispersion rate is 1/(sqrt2 sigma) = 1.414 Bohr/a.u., so by the
# time it reaches the slab face it has grown from a density std of 0.35 Bohr to
# 4.7-8.1 Bohr depending on velocity, and by the exit face to 15-26 Bohr — wider
# than the 25 Bohr slab itself. It is no longer a localised projectile when it
# does the physics. sigma = 2 and 3 cut the dispersion rate to 0.354 and 0.236
# Bohr/a.u., so the packet meets the jellium at 1.8-2.5 Bohr and leaves at
# 3.2-6.6 Bohr: an object that stays recognisably itself across the transit.
#
# WHAT CHANGES: sigma_WP only. Same ground state (sigma does not enter the bath),
# same dx = 0.40, same launch z = -24, same CAPs (12.5 Bohr/face, eta = -1), same
# dt = 0.04 and step counts, same cadences, same four velocities. That is what
# makes the three sigma traces comparable point for point.
#
# THE OBSERVABLE: S = (E_total(t_final) - E_GS) / L_slab_z, the localised-jellium
# deposit definition, plotted for sigma = 0.5, 2.0 and 3.0 on one axis.
#
# THE CHAIN — one campaign's four velocities run in PARALLEL, then the next
# campaign's four, per the user's requested shape:
#
#   1. smoke2   20 steps at sigma = 2      rebuilds the binary and runs the t=0
#                                          analytic gates. afterok-gates stage 2,
#                                          so a bad packet costs 2 minutes, not
#                                          4 GPUs x 3 hours.
#   2. sweep2   array 0-3, sigma = 2       v = 2.0/2.5/3.0/3.5 in PARALLEL.
#   3. vac2     4 CAP-only controls        sigma = 2 baseline. Each campaign needs
#                                          its own: CAP attrition depends on the
#                                          spreading rate, which is what sigma sets.
#   4. smoke3   20 steps at sigma = 3      as (1), for the second campaign.
#   5. sweep3   array 0-3, sigma = 3       v = 2.0/2.5/3.0/3.5 in PARALLEL.
#   6. vac3     4 CAP-only controls        sigma = 3 baseline.
#   7. nb       build + EXECUTE notebooks  per-run notebooks for all 8 new runs,
#                                          per-campaign synthesis, and the cross-
#                                          sigma S_deposit(v) figure carrying all
#                                          three traces.
#
# Stages 2->3->4 are chained afterANY rather than afterok: one velocity failing
# should not cancel the other campaign. Only the smoke->sweep edges are afterok,
# because those are correctness gates (.claude/rules/checkpoint-dont-block.md:
# gate on correctness, never on cost).
#
# THE -z CAP CAVEAT AT sigma = 3. Launch z = -24, CAP inner edge z = -30, so 6
# Bohr of clearance = 2.8 density-std at sigma = 3 (vs 4.2 at sigma = 2). About
# 0.23 % of the packet starts inside the absorbing band and is removed over the
# first few a.u. Small, one-off, reproduced by vac3, and reported in the
# notebooks rather than folded silently into S.
#
# Plan: docs/plans/wavepacket-highdensity-sv-twin.md
set -euo pipefail

cd "$(dirname "$0")/../.."
REPO_ROOT="$(pwd)"

GS_DIR="$REPO_ROOT/ResearchProject/systems/localised_jellium/shared_gs/slab_n100_L35x35x85_dx0p4_per2"
if [ ! -d "$GS_DIR" ]; then
  echo "FATAL: production ground state missing at $GS_DIR" >&2
  echo "       It is sigma-independent and was built by the 0.5 campaign; if it is" >&2
  echo "       gone, run: sbatch shared/bin/run-wp-hd-gs.slurm 0.4" >&2
  exit 2
fi
echo "Ground state present (sigma-independent, reused by both campaigns):"
echo "  $GS_DIR"
echo ""
echo "Submitting the sigma = 2 and sigma = 3 wavepacket campaigns..."

SMOKE2=$(sbatch --parsable --export=ALL,LJ_SIGMA=2.0 \
          shared/bin/run-wp-hd-wp.slurm smoke)
echo "  1. smoke2  (sigma=2, build + t=0 gates)   job $SMOKE2"

SWEEP2=$(sbatch --parsable --dependency=afterok:"$SMOKE2" --export=ALL,LJ_SIGMA=2.0 \
          --array=0-3 shared/bin/run-wp-hd-wp.slurm)
echo "  2. sweep2  (sigma=2, 4 velocities || )    job $SWEEP2   [afterok $SMOKE2]"

VAC2=$(sbatch --parsable --dependency=afterany:"$SWEEP2" --export=ALL,LJ_SIGMA=2.0 \
        shared/bin/run-wp-hd-vac.slurm)
echo "  3. vac2    (sigma=2, CAP baselines)       job $VAC2   [afterany $SWEEP2]"

SMOKE3=$(sbatch --parsable --dependency=afterany:"$VAC2" --export=ALL,LJ_SIGMA=3.0 \
          shared/bin/run-wp-hd-wp.slurm smoke)
echo "  4. smoke3  (sigma=3, t=0 gates)           job $SMOKE3   [afterany $VAC2]"

SWEEP3=$(sbatch --parsable --dependency=afterok:"$SMOKE3" --export=ALL,LJ_SIGMA=3.0 \
          --array=0-3 shared/bin/run-wp-hd-wp.slurm)
echo "  5. sweep3  (sigma=3, 4 velocities || )    job $SWEEP3   [afterok $SMOKE3]"

VAC3=$(sbatch --parsable --dependency=afterany:"$SWEEP3" --export=ALL,LJ_SIGMA=3.0 \
        shared/bin/run-wp-hd-vac.slurm)
echo "  6. vac3    (sigma=3, CAP baselines)       job $VAC3   [afterany $SWEEP3]"

NB=$(sbatch --parsable --dependency=afterany:"$VAC3" \
      shared/bin/run-wp-hd-notebooks.slurm all)
echo "  7. nb      (all notebooks + sigma sweep)  job $NB   [afterany $VAC3]"

echo ""
echo "Chain submitted. Monitor with:"
echo "    squeue -u \$USER"
echo "    tail -f wp-hd-*.out"
echo ""
echo "Expected wall clock: each campaign ~3.2 h (the v=2.0 point is the long pole,"
echo "3623 steps; the four run concurrently), so ~7 h end to end including the"
echo "cheap vacuum controls and the notebook stage."
echo ""
echo "Each production run keeps 5 retained checkpoints. To extend a point:"
echo "    LJ_RESUME=1 LJ_SIGMA=2.0 sbatch shared/bin/run-wp-hd-wp.slurm <idx>"
