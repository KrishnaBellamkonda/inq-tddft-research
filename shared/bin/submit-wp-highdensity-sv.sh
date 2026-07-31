#!/bin/bash
# submit-wp-highdensity-sv.sh — submit the WHOLE wavepacket twin campaign as one
# dependency-chained, autonomous set of SLURM jobs. Run once from the repo root:
#
#     ./shared/bin/submit-wp-highdensity-sv.sh
#
# Nothing needs babysitting after this: each stage starts only when its
# predecessor succeeds, and a correctness failure anywhere stops the chain
# instead of silently seeding downstream runs with a bad input.
#
# THE CHAIN
#   1. gs05    ground state at dx = 0.50            FIDELITY CHECK — run.cpp exits 3
#                                                    unless E_GS = 207.18322156141 Ha
#                                                    is reproduced, proving the CSD3
#                                                    rebuild is the classical
#                                                    campaign's system.
#   2. gs04    ground state at dx = 0.40            PRODUCTION ground state (afterok
#                                                    gs05: no point building on an
#                                                    unverified system).
#   3. smoke   20 steps of wp/run.cpp               builds the production binary and
#                                                    runs the t=0 analytic gates
#                                                    (<p_z>, sigma_pz^2, T1, T1-T2,
#                                                    centroid, density std) which
#                                                    abort on a bad packet.
#   4. sweep   array 0-3: v = 2.0/2.5/3.0/3.5       the four production points, in
#                                                    PARALLEL, one A100 each.
#   5. vac     four CAP-only vacuum controls        the baseline that separates real
#                                                    stopping from CAP attrition
#                                                    (see plan section 6b item 5);
#                                                    runs concurrently with the sweep.
#   6. nb      build + EXECUTE the notebooks        one per velocity (density GIFs
#                                                    first, step-by-step S) plus the
#                                                    synthesis vs the classical curve.
#                                                    afterANY, so a single failed
#                                                    velocity still yields notebooks
#                                                    for the rest.
#
# WHY v = 4.0 AND 4.5 ARE NOT IN THE SWEEP: at dx = 0.40 they alias sigma_pz^2 by
# +17.9 % and +55.1 % (measured/modelled 2026-07-30). See run-wp-hd-wp.slurm.
#
# Plan: docs/plans/wavepacket-highdensity-sv-twin.md
set -euo pipefail

# repo root = two levels up from shared/bin/. SLURM_SUBMIT_DIR inside each job is
# taken from the CWD at sbatch time, so every job script resolves REPO_ROOT here.
cd "$(dirname "$0")/../.."

echo "Submitting the wavepacket high-density S(v) chain..."

GS05=$(sbatch --parsable shared/bin/run-wp-hd-gs.slurm 0.5)
echo "  1. gs05   (dx=0.50 fidelity check)      job $GS05"

GS04=$(sbatch --parsable --dependency=afterok:"$GS05" shared/bin/run-wp-hd-gs.slurm 0.4)
echo "  2. gs04   (dx=0.40 production GS)       job $GS04   [afterok $GS05]"

SMOKE=$(sbatch --parsable --dependency=afterok:"$GS04" shared/bin/run-wp-hd-wp.slurm smoke)
echo "  3. smoke  (build + t=0 gates)           job $SMOKE   [afterok $GS04]"

SWEEP=$(sbatch --parsable --dependency=afterok:"$SMOKE" --array=0-3 shared/bin/run-wp-hd-wp.slurm)
echo "  4. sweep  (v=2.0,2.5,3.0,3.5 parallel)  job $SWEEP   [afterok $SMOKE]"

VAC=$(sbatch --parsable --dependency=afterok:"$SMOKE" shared/bin/run-wp-hd-vac.slurm)
echo "  5. vac    (CAP-only baselines)          job $VAC   [afterok $SMOKE]"

NB=$(sbatch --parsable --dependency=afterany:"$SWEEP",afterany:"$VAC" \
      shared/bin/run-wp-hd-notebooks.slurm all)
echo "  6. nb     (notebooks + synthesis)         job $NB   [afterany $SWEEP,$VAC]"

echo ""
echo "Chain submitted. Monitor with:"
echo "    squeue -u \$USER"
echo "    tail -f wp-hd-gs-$GS05.out"
echo ""
echo "If the dx=0.50 fidelity check FAILS (exit 3), the whole chain stops and"
echo "nothing else runs — that is deliberate: it would mean the rebuilt ground"
echo "state is not the classical campaign's system."
echo ""
echo "To extend any completed point to more steps (never recompute from step 0):"
echo "    LJ_RESUME=1 sbatch shared/bin/run-wp-hd-wp.slurm <idx>   # with a larger N"
