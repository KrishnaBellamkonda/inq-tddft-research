#!/bin/bash
# submit-ladder.sh — submit the WHOLE cylindrical proximity-ladder campaign as one
# dependency-chained graph, so it runs to completion unattended.
# Plan: docs/plans/cylindrical-proximity-ladder.md
#
#   Usage:  cd <repo> && shared/bin/submit-ladder.sh
#           LADDER_GS_JOB=32667942 LADDER_BC_JOB=32668084 shared/bin/submit-ladder.sh
#           LADDER_RUNGS="r08 r06" shared/bin/submit-ladder.sh
#           LADDER_DRYRUN=1 shared/bin/submit-ladder.sh      # print, do not submit
#
# THE GRAPH
#
#   ladder-gs ─────┐                     (4 ground states, sequential: shared build dir)
#                  ├──> smoke(half,rung) ──> prod(half,rung) ──┐
#   buildcheck ────┘     20 steps, gates      1500 steps        │
#                                                               ├──> notebooks
#                                                               │    (per-rung +
#                                                               └───  cross-rung)
#
# WHY SMOKE IS IN THE CHAIN AND NOT OPTIONAL. A 20-step run costs ~2 minutes and
# executes every t=0 gate: injection norm, <p_z>=k0, var(p), the circular centroid,
# the Rayleigh f_bore/f_wall check, and the interaction-energy closure. Wiring it as
# an afterok parent means a rung with a broken ground state or a mis-set geometry
# never burns 1.5 GPU-hours to discover it.
#
# WHY prod DEPENDS ONLY ON ITS OWN SMOKE. Every (half, rung) is independent — the
# binary picks its geometry from CJ_RUNG at runtime — so the eight production runs
# fan out across whatever GPUs are free rather than serialising.
#
# WHY notebooks USE afterany, NOT afterok. If one rung fails its gates, the other
# three still have complete data and must still be written up. The notebook stage
# reports which rungs it found and which are missing, rather than refusing to run.
# That is the difference between a chain that degrades and one that dies.
set -uo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
[[ -f "$REPO_ROOT/shared/bin/csd3-env.sh" ]] || {
  echo "ERROR: run from the repo root." >&2; exit 1; }
cd "$REPO_ROOT" || exit 1

RUNGS="${LADDER_RUNGS:-r08 r06 r04 r00}"
DRY="${LADDER_DRYRUN:-0}"

sub () {   # sub <dependency-or-empty> <script> [args...]  -> prints the job id
  local dep="$1"; shift
  local args=(--parsable)
  [ -n "$dep" ] && args+=(--dependency="$dep" --kill-on-invalid-dep=yes)
  if [ "$DRY" = "1" ]; then
    echo "DRYRUN sbatch ${args[*]} $*" >&2
    echo "000000"
  else
    sbatch "${args[@]}" "$@"
  fi
}

# ---- stage 0: ground states + the one build, reused if already in flight -----
GS_JOB="${LADDER_GS_JOB:-}"
BC_JOB="${LADDER_BC_JOB:-}"
[ -z "$GS_JOB" ] && GS_JOB=$(sub "" shared/bin/run-ladder-gs.slurm)
[ -z "$BC_JOB" ] && BC_JOB=$(sub "" shared/bin/run-ladder-buildcheck.slurm)
echo "stage 0  ground states = $GS_JOB   build check = $BC_JOB"

READY="afterok:${GS_JOB}:${BC_JOB}"

# ---- stages 1-2: smoke then production, per (half, rung) --------------------
ALL_PROD=""
for RUNG in $RUNGS; do
  for HALF in wp classical; do
    # Time limits are per STAGE, and they matter more than they look. SLURM
    # charges the account's GPU-minutes budget (AssocGrpGRESMinutes) against
    # REQUESTED wall time, not used, and backfill favours short jobs — so asking
    # 12 h for a 20-step smoke run both burns budget and loses scheduling gaps.
    # Smoke measured ~123 s (wp) / 78 s (classical) at 104 states, so 30 min is a
    # ~7x margin; production is 0.9-1.8 h by state count, so 6 h is ~3.3x.
    S=$(sub "$READY"       --time=00:30:00 shared/bin/run-ladder-rt.slurm "$HALF" "$RUNG" smoke)
    P=$(sub "afterok:${S}" --time=06:00:00 shared/bin/run-ladder-rt.slurm "$HALF" "$RUNG" prod)
    ALL_PROD="${ALL_PROD}:${P}"
    printf 'stage 1-2  %-9s %-4s  smoke=%s  prod=%s\n' "$HALF" "$RUNG" "$S" "$P"
  done
done

# ---- stage 3: notebooks + cross-rung comparison -----------------------------
# afterany (see header): write up whatever completed.
NB=$(sub "afterany${ALL_PROD}" shared/bin/run-ladder-notebooks.slurm)
echo "stage 3  notebooks = $NB"

echo
echo "Submitted. Watch with:  squeue -u $USER -o '%.10i %.22j %.9T %.8M %.16R'"
echo "Rungs: $RUNGS"
