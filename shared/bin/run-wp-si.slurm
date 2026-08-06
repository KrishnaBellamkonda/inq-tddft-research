#!/bin/bash
#SBATCH --job-name=wp-si
#SBATCH --account=mphil-nikiforakis-skcb2-sl2-gpu
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=06:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.out
#
# run-wp-si.slurm — quantify the wavepacket self-interaction in VACUUM.
# Plan: docs/plans/wp-self-interaction-correction.md
#
#   sbatch shared/bin/run-wp-si.slurm smoke     # build + 20 steps of each theory
#   sbatch shared/bin/run-wp-si.slurm prod      # the three production runs
#
# ONE JOB RUNS ALL FIVE CONFIGURATIONS, sequentially on the same GPU and the
# same binary. That is deliberate: the comparison is only meaningful if nothing
# but `WP_THEORY`/`WP_SIC` differs, and running them from one job on one build
# removes any possibility of a stale binary or a divergent environment.
#
#   noninteracting          no self-interaction at all   <- the reference
#   hartree                 + Hartree self-interaction
#   lda                     + Hartree AND LDA xc self-interaction
#   sic_h     (lda+SIC-H)   lda MINUS the Hartree self-term      <- Tier V (b)
#   sic_pzrun (lda+SIC-PZ)  lda MINUS Hartree AND xc self-terms  <- Tier V (c),
#                           HARD-GATED against closed-form free dispersion:
#                           the binary exits 4 on gate failure, failing this
#                           job, and afterok then blocks Tier B + production
#                           (plan §4 decision rule).
#
# Completed outputs (run_summary.txt says run_completed = true) are SKIPPED so
# a resubmission never recomputes; set WP_FORCE=1 to override.
#
# Each is cheap (~15 min measured): ONE electron on a 144^3 grid, versus 104
# orbitals on 80x80x120 for the channeling twin.
set -euo pipefail

STAGE="${1:-prod}"
case "$STAGE" in
  smoke|prod) ;;
  *) echo "FATAL: argument must be 'smoke' or 'prod'; got '$STAGE'" >&2; exit 2 ;;
esac

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
if [[ ! -f "$REPO_ROOT/shared/bin/csd3-env.sh" ]]; then
  echo "FATAL: REPO_ROOT='$REPO_ROOT' is not the repo root." >&2; exit 2
fi
# shellcheck source=./csd3-env.sh
. "$REPO_ROOT/shared/bin/csd3-env.sh"

export INQ_SOURCE="$REPO_ROOT/inq-study"
export INQ_SHARE_PATH="$REPO_ROOT/inq/install/share"
export PSEUDOPOD_SHARE_PATH="$REPO_ROOT/inq/install/share/pseudopod"
# shellcheck source=../config.sh
. "$REPO_ROOT/shared/config.sh"
export PATH="$REPO_ROOT/shared/bin:$PATH"
export INQ_BUILD_JOBS="${SLURM_CPUS_PER_TASK:-8}"
export TMPDIR="${TMPDIR:-$REPO_ROOT/.build_tmp}"
mkdir -p "$TMPDIR"

RUN_DIR="$REPO_ROOT/ResearchProject/systems/vacuum/scripts/wp_selfinteraction"
cd "$RUN_DIR"

# Shared physics for all three runs. Only WP_THEORY varies below.
export WP_SIGMA="${WP_SIGMA:-4.0}"      # matched to channeling_twin
export WP_K0="${WP_K0:-0.0}"            # stationary: spreading is frame-independent
export WP_L="${WP_L:-72.0}"
export WP_H="${WP_H:-0.5}"
export WP_DT="${WP_DT:-0.02}"

if [ "$STAGE" = "smoke" ]; then
  export WP_NSTEPS=20 WP_SAVE_EVERY=10 WP_WF_EVERY=10
  SUFFIX="_smoke"
else
  export WP_NSTEPS="${WP_NSTEPS:-1500}"
  export WP_SAVE_EVERY="${WP_SAVE_EVERY:-150}"
  export WP_WF_EVERY="${WP_WF_EVERY:-150}"
  SUFFIX=""
fi

# PRE-FLIGHT DISK GUARD. Three runs x (11 density frames @ 24 MB + 11 complex
# wavefunctions @ 48 MB + a 48 MB checkpoint) ~ 2.5 GB, plus headroom.
NEED_GB="${WP_NEED_GB:-6}"
hd=$(quota -s 2>/dev/null | awk '/user\/skcb2/ {gsub(/[!*]/,"",$2); printf "%.0f", $3-$2; exit}')
if [ -n "${hd:-}" ] && [ "$STAGE" = "prod" ] && [ "$hd" -lt "$NEED_GB" ]; then
  echo "FATAL: only ${hd} GB free on /rds; the three runs need ~${NEED_GB} GB." >&2
  exit 3
fi

echo "=============================================="
echo " job     : ${SLURM_JOB_ID:-none}  wp-si [$STAGE]"
echo " node    : $(hostname)"
echo " gpu     : $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>&1 | head -1)"
echo " engine  : $INQ_SOURCE"
echo " sigma   : $WP_SIGMA Bohr   k0 = $WP_K0   box ${WP_L}^3 @ h=$WP_H"
echo " steps   : $WP_NSTEPS x $WP_DT"
echo " disk    : ${hd:-unknown} GB free (need ~${NEED_GB} GB)"
echo " started : $(date)"
echo "=============================================="

rc=0
FIRST=1
for CFG in noninteracting:none:noninteracting \
           hartree:none:hartree \
           lda:none:lda \
           lda:h:sic_h \
           lda:pzrun:sic_pzrun; do
  TH="${CFG%%:*}"; rest="${CFG#*:}"; SICM="${rest%%:*}"; OUTN="${rest#*:}"
  echo
  echo "########## THEORY = $TH   SIC = $SICM   OUT = ${OUTN}${SUFFIX} ##########"
  if [ -z "${WP_FORCE:-}" ] && \
     grep -q "run_completed = true" "results/${OUTN}${SUFFIX}/run_summary.txt" 2>/dev/null; then
    echo "  [skip] results/${OUTN}${SUFFIX} already complete"
    continue
  fi
  export WP_THEORY="$TH" WP_SIC="$SICM" WP_OUT="${OUTN}${SUFFIX}"
  t0=$(date +%s)
  if [ "$FIRST" = "1" ]; then
    # inq-run BUILDS then execs. Only the first pass needs to compile; the
    # remaining passes exec the binary directly so they cannot race the CMake tree.
    inq-run || rc=1
    FIRST=0
  else
    ./run || rc=1
  fi
  echo "  [$TH/$SICM] exit=$? wall=$(( $(date +%s) - t0 )) s"
done

echo
echo "=============================================="
echo " exit status : $rc"
echo " finished    : $(date)"
if [ "$rc" -eq 0 ]; then
  echo " outputs:"
  for OUTN in noninteracting hartree lda sic_h sic_pzrun; do
    d="$RUN_DIR/results/${OUTN}${SUFFIX}"
    [ -f "$d/run_summary.txt" ] && echo "   $d  ($(du -sh "$d" 2>/dev/null | cut -f1))"
  done
  echo
  echo " --- SIC gate summaries ---"
  for OUTN in sic_h sic_pzrun; do
    s="$RUN_DIR/results/${OUTN}${SUFFIX}/run_summary.txt"
    [ -f "$s" ] && { echo " [$OUTN]"; grep -A20 "sic_gates:" "$s" | head -12; }
  done
fi
echo "=============================================="
exit $rc
