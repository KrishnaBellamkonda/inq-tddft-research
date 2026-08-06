#!/bin/bash
#SBATCH --job-name=bulk-t0
#SBATCH --account=mphil-nikiforakis-skcb2-sl2-gpu
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=00:30:00
#SBATCH --output=%x-%A_%a-%j.out
#SBATCH --error=%x-%A_%a-%j.out
#
# run-bulk-t0.slurm — regenerate the t = 0 DENSITY FIELDS for the bulk-jellium
# twin pairs at sigma_WP = 1, 2, 3. Report-2 draft-1 figures B-setup and B-pauli.
#
#   Usage:  sbatch --array=0-5 shared/bin/run-bulk-t0.slurm
#           sbatch shared/bin/run-bulk-t0.slurm 0        # single task
#
# WHY THIS EXISTS. The production runs
# (ResearchProject/systems/jellium/scripts/bulk_ks_stopping_rs4{,_sigma1,_sigma3})
# DID write results/raw/vti/density_{total,wp,delta} at Cfg::WRITE_EVERY = 8, but
# those fields were pruned from disk afterwards; only results/raw/observables/
# survives. The energy analysis is therefore intact and is NOT re-run here — only
# the t = 0 fields are missing, and they are what the two density figures need.
#
# WHAT IT DOES, AND WHY IT IS CHEAP AND SAFE.
#   * Re-uses the ALREADY-BUILT binaries (<family>/<half>/run). Nothing is
#     compiled, and no run.cpp is touched, so run provenance is unchanged.
#   * BKS_N_STEPS=2. The VTI write is gated on `abs_step % WRITE_EVERY == 0` and
#     the callback fires at step 0, so step 0 is written by the UNMODIFIED binary.
#     Two steps rather than one so the propagator is genuinely exercised and the
#     step-0 fields are not a special case of an empty loop.
#   * The binaries write to RELATIVE paths ("results/raw/...") while the ground
#     state is loaded from an ABSOLUTE GS_DIR. Running from a scratch CWD therefore
#     redirects ALL output there and CANNOT clobber the production observables.
#     This is the property that makes the whole approach safe — verify it holds
#     before adapting this script to another family.
#
# COST. Each task is a ground-state load (482 orbitals, 80x80x160 grid) plus two
# ETRS steps. The GS load dominates. 30 minutes requested; note that the GPU
# budget is charged on REQUESTED wall time, so do not inflate this.
#
# OUTPUT. ResearchProject/systems/jellium/scripts/bulk_t0_density/<tag>/results/
#   raw/vti/density_total/      both halves
#   raw/vti/density_wp/         wavepacket half only
#   raw/vti/density_gs_system/  the bath before the projectile exists
# The wavepacket half's BATH density is density_total - density_wp; the classical
# half's bath density IS density_total (its projectile is an external potential
# and carries no charge in the electron density). That asymmetry is the whole
# point of figure B-pauli, so do not "symmetrise" it.
#
#   idx  family                        half
#     0  bulk_ks_stopping_rs4          wp          sigma_WP = 2
#     1  bulk_ks_stopping_rs4          classical
#     2  bulk_ks_stopping_rs4_sigma1   wp          sigma_WP = 1
#     3  bulk_ks_stopping_rs4_sigma1   classical
#     4  bulk_ks_stopping_rs4_sigma3   wp          sigma_WP = 3
#     5  bulk_ks_stopping_rs4_sigma3   classical
#
set -euo pipefail

STAGE="${1:-${SLURM_ARRAY_TASK_ID:-}}"
if [ -z "$STAGE" ]; then
  echo "FATAL: pass a task index 0-5 (or submit as --array=0-5)" >&2
  exit 2
fi

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"

. "$REPO_ROOT/shared/bin/csd3-env.sh"

export INQ_SOURCE="$REPO_ROOT/inq-study"
export INQ_SHARE_PATH="$REPO_ROOT/inq/install/share"
export PSEUDOPOD_SHARE_PATH="$REPO_ROOT/inq/install/share/pseudopod"

. "$REPO_ROOT/shared/config.sh"

export PATH="$REPO_ROOT/shared/bin:$PATH"

SCRIPTS="$REPO_ROOT/ResearchProject/systems/jellium/scripts"
OUTROOT="$SCRIPTS/bulk_t0_density"

FAM_LIST=( bulk_ks_stopping_rs4 bulk_ks_stopping_rs4 \
           bulk_ks_stopping_rs4_sigma1 bulk_ks_stopping_rs4_sigma1 \
           bulk_ks_stopping_rs4_sigma3 bulk_ks_stopping_rs4_sigma3 )
HALF_LIST=( wp classical wp classical wp classical )
TAG_LIST=( s2_wp s2_classical s1_wp s1_classical s3_wp s3_classical )

IDX="$STAGE"
if [ "$IDX" -lt 0 ] || [ "$IDX" -gt 5 ]; then
  echo "FATAL: task index $IDX out of range 0-5" >&2
  exit 2
fi

FAM="${FAM_LIST[$IDX]}"
HALF="${HALF_LIST[$IDX]}"
TAG="${TAG_LIST[$IDX]}"
BIN="$SCRIPTS/$FAM/$HALF/run"
OUT="$OUTROOT/$TAG"

if [ ! -x "$BIN" ]; then
  echo "FATAL: binary missing or not executable: $BIN" >&2
  echo "       These runs deliberately do NOT build. Build it in place first." >&2
  exit 2
fi

N_STEPS="${BKS_N_STEPS:-2}"

echo "=== bulk t=0 density regeneration"
echo "    task    : $IDX  ($TAG)"
echo "    family  : $FAM"
echo "    half    : $HALF"
echo "    binary  : $BIN"
echo "    out dir : $OUT"
echo "    n_steps : $N_STEPS  (step 0 is what the figures need)"
echo

mkdir -p "$OUT"
cd "$OUT"

env INQ_SHARE_PATH="$INQ_SHARE_PATH" \
    PSEUDOPOD_SHARE_PATH="$PSEUDOPOD_SHARE_PATH" \
    BKS_N_STEPS="$N_STEPS" \
    "$BIN"

echo
echo "=== wrote:"
find "$OUT/results/raw/vti" -name "*.vti" 2>/dev/null | sed "s|$OUT/||" | sort
echo "=== done: $TAG"
