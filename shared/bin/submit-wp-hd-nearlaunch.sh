#!/bin/bash
# submit-wp-hd-nearlaunch.sh — the near-launch sigma_WP = 0.5 campaign.
#
# THE HYPOTHESIS (docs/plans/effective-sigma-near-launch.md). It is not the
# LAUNCH sigma that sets the wavepacket-slab interaction, it is the packet's
# width when it ARRIVES. The far-launch campaign starts at z = -24 and the
# sigma = 0.5 packet disperses to 4.7-8.1 Bohr before it reaches the slab, so its
# "sigma = 0.5" label describes a packet that no longer exists by the time the
# physics happens. Launching at z = -14 (1.5 Bohr outside the face at -12.5,
# just beyond the 1 Bohr erfc softening) delivers it essentially undispersed at
# sigma/sqrt2 = 0.354 Bohr. Same sigma, same everything else, different arrival
# width: if S(v) moves, the launch sigma was never the controlling parameter.
#
# SCAN RESULT (job 32528019, 2026-08-01) — already run, gating already passed:
#   z = -14.0 removes 0.109 % of the packet to orthogonalisation (criterion 3 %),
#   the k_z marginal keeps R^2 = 0.99995 against the analytic N(k0, sigma_p^2),
#   and the packet CORE is 0.3559 Bohr vs the Gaussian 0.3536 (+0.65 %).
#   The z = -24 regression trial reproduced the campaign's recorded
#   max_overlap = 3.691564855e-4 to 12 significant figures.
# So the chain below starts at the smoke stage; re-run run-wp-hd-scan.slurm first
# only if the GS, dx or sigma change.
#
# CHAIN
#   1 smoke   builds the binary (inqkit gained the ortho-loss fields) + t=0 gates
#   2 sweep   array 0-3, v = 2.0/2.5/3.0/3.5, the four production points
#   3 vac     CAP-only baselines at the SAME launch z (subtractable step for step)
#   4 nb      notebooks + synthesis
#
# Every stage takes LJ_LAUNCH_Z=auto, which reads the scan's own
# inject_scan/results/scan/accepted_launch_z.txt — no number is relayed by hand.
# Run names are prefixed nl_ so every existing far-launch run, notebook and
# summary CSV keeps resolving unchanged.
#
# N_STEPS is deliberately IDENTICAL to the far-launch campaign (3623/2898/2415/
# 2070). The path is 10 Bohr shorter, so this buys extra post-exit plateau for
# the deposit estimator rather than costing anything, and keeps the time budget
# comparable point for point.
#
# Usage:  cd <repo> && bash shared/bin/submit-wp-hd-nearlaunch.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

SCAN_Z_FILE="ResearchProject/systems/localised_jellium/scripts/wp_highdensity_sv/inject_scan/results/scan/accepted_launch_z.txt"
if [ ! -f "$SCAN_Z_FILE" ]; then
  echo "FATAL: no accepted launch z. Run: sbatch shared/bin/run-wp-hd-scan.slurm" >&2
  exit 2
fi
Z="$(tr -d '[:space:]' < "$SCAN_Z_FILE")"
echo "near-launch campaign at z = $Z (standoff $(awk "BEGIN{printf \"%.1f\", -12.5-($Z)}") Bohr)"

# 1 — smoke: BUILDS the binary, runs 20 steps, applies the t=0 gates.
J1=$(sbatch --parsable --export=ALL,LJ_LAUNCH_Z=auto \
       shared/bin/run-wp-hd-wp.slurm smoke)
echo "  1 smoke        : $J1"

# 2 — the four production velocities, in parallel, after the build succeeds.
J2=$(sbatch --parsable --dependency=afterok:"$J1" --array=0-3 \
       --export=ALL,LJ_LAUNCH_Z=auto shared/bin/run-wp-hd-wp.slurm)
echo "  2 sweep 0-3    : $J2  (afterok $J1)"

# 3 — vacuum CAP baselines at the same launch z. afterany: a failed velocity
#     must not cost the baselines for the ones that worked.
J3=$(sbatch --parsable --dependency=afterany:"$J2" \
       --export=ALL,LJ_LAUNCH_Z=auto shared/bin/run-wp-hd-vac.slurm)
echo "  3 vac controls : $J3  (afterany $J2)"

# 4 — notebooks + synthesis.
J4=$(sbatch --parsable --dependency=afterany:"$J3" \
       shared/bin/run-wp-hd-notebooks.slurm)
echo "  4 notebooks    : $J4  (afterany $J3)"

echo ""
echo "submitted. monitor with: squeue -u \$USER"
echo "kill the chain with   : scancel $J1 $J2 $J3 $J4"
