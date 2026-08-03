#!/bin/bash
# submit-channeling-sic.sh — the FULL validated chain for the SIC campaign.
# Plan: docs/plans/wp-self-interaction-correction.md (reviewed 2026-08-02).
#
#   ./shared/bin/submit-channeling-sic.sh          # submit the whole chain
#
# Chain (each stage gates the next with afterok):
#
#   1. chan-tests   library gate incl. the NEW test_wp_sic_engine (kick
#                   semantics, exact Q re-orthogonalisation, D1 xc consistency)
#   2. wp-si        Tier V, vacuum: 5 configs (noninteracting/hartree/lda +
#                   lda/sic_h + lda/sic_pzrun). sic_pzrun carries the HARD
#                   closed-form free-dispersion gates — the binary exits 4 on
#                   failure, failing the job and blocking everything below
#                   (the plan's decision rule, mechanised).
#   3. chan-sic tierb   200 production steps with SIC: bath-integrity gates
#                   (the projection's first contact with a real occupied
#                   manifold — Tier V has none, plan §0/D6).
#   4. chan-sic prod    the 1500-step production run (third leg of the
#                   classical / wp / wp+SIC comparison).
#
# The classical and uncorrected-WP references are the COMPLETED channeling_twin
# runs; nothing re-runs them.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

echo "== channeling SIC chain =="

j_tests=$(sbatch --parsable shared/bin/run-chan-tests.slurm)
echo " 1. chan-tests            : $j_tests"

j_tierv=$(sbatch --parsable --dependency=afterok:"$j_tests" shared/bin/run-wp-si.slurm prod)
echo " 2. wp-si (Tier V vacuum) : $j_tierv   afterok($j_tests)"

j_tierb=$(sbatch --parsable --dependency=afterok:"$j_tierv" shared/bin/run-chan-sic.slurm tierb)
echo " 3. chan-sic tierb        : $j_tierb   afterok($j_tierv)"

j_prod=$(sbatch --parsable --dependency=afterok:"$j_tierb" shared/bin/run-chan-sic.slurm prod)
echo " 4. chan-sic prod         : $j_prod   afterok($j_tierb)"

echo
echo "watch:  squeue -u $USER -o '%.10i %.12j %.2t %.10M %.20E'"
echo "logs :  chan-tests-$j_tests.out  wp-si-$j_tierv.out  chan-sic-$j_tierb.out  chan-sic-$j_prod.out"
