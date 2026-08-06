#!/bin/bash
#SBATCH --job-name=wp-si-sweep
#SBATCH --account=mphil-nikiforakis-skcb2-sl2-gpu
#SBATCH --partition=ampere
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=04:00:00
#SBATCH --array=0-5%3
#SBATCH --output=%x-%A_%a.out
#SBATCH --error=%x-%A_%a.out
#
# run-wp-si-sweep.slurm — how does the wavepacket self-interaction scale with sigma?
#
#   sbatch shared/bin/run-wp-si-sweep.slurm
#
# One array task per sigma_WP; four theory levels each, sequentially on one GPU.
#
# THE SCALING PROTOCOL — why every run costs exactly the same
# -----------------------------------------------------------
# The free Schrodinger equation is scale-invariant under r -> lambda r,
# t -> lambda^2 t. So if the box, the spacing, the time step and the duration are
# all scaled off sigma, EVERY sigma is the same discrete problem:
#
#     L  = 18 sigma          box half-width = 6 sigma_dens(t_end)
#     h  = 0.125 sigma       grid is ALWAYS 144^3
#     dt = 0.00125 sigma^2   n_steps is ALWAYS 1500
#     t_end = 1.875 sigma^2  sigma_dens(t_end) = 1.5024 sigma  (2.12x the initial width)
#
# At sigma = 4 this reduces EXACTLY to the completed run (L=72, h=0.5, dt=0.02,
# 1500 steps, t_end=30), so the sigma=4 point of this sweep is a reproduction
# check on the harness as well as a data point.
#
# WHAT THE SWEEP ACTUALLY MEASURES. Under that scaling the only thing that
# changes is the DIMENSIONLESS coupling. Kinetic energy scales as 1/sigma^2
# while both the Hartree self-energy (E_PP ~ 1/sigma) and the LDA exchange
# potential (v_x ~ n^(1/3) ~ 1/sigma) scale as 1/sigma. Their ratio to the
# kinetic term is therefore ~ sigma: a WIDER packet is MORE strongly
# self-interacting per unit kinetic energy, at fixed dimensionless time.
#
# Because Hartree and LDA exchange carry the SAME 1/sigma scaling, an
# exchange-only cancellation would be sigma-independent. LDA CORRELATION has a
# different density dependence and does not follow the similarity. The sweep
# therefore tests whether the ~82 % cancellation measured at sigma = 4 is a
# constant of the functional or an accident of that density. (Inference from the
# scaling of the LDA functional forms; the run decides it.)
#
# FOUR CONFIGURATIONS PER SIGMA
#   noninteracting          the reference: exact free dispersion
#   hartree                 + Hartree self-interaction
#   lda                     + Hartree AND LDA xc self-interaction
#   lda + sic_pzrun         full Perdew-Zunger removal (must return to the reference)
#
# Completed outputs are SKIPPED; set WP_FORCE=1 to override.
set -euo pipefail

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

# The sweep reuses the binary built by run-wp-si.slurm. Building here would race
# across array tasks sharing one CMake tree.
if [ ! -x ./run ]; then
  echo "FATAL: ./run not built. Run 'sbatch shared/bin/run-wp-si.slurm prod' first" >&2
  echo "       (or build once with 'inq-run' in $RUN_DIR)." >&2
  exit 2
fi
if [ run.cpp -nt ./run ]; then
  echo "FATAL: run.cpp is newer than ./run — rebuild before sweeping." >&2
  exit 2
fi

SIGMAS=(1.0 2.0 3.0 4.0 6.0 8.0)
TAG=(   1p0 2p0 3p0 4p0 6p0 8p0)
IDX="${SLURM_ARRAY_TASK_ID:-3}"
SIG="${SIGMAS[$IDX]}"
STAG="${TAG[$IDX]}"

# The scale-exact protocol. Reduces to the completed sigma=4 run identically.
export WP_SIGMA="$SIG"
export WP_K0=0.0                                       # spreading is boost-invariant
export WP_L=$(python3   -c "print(18.0*$SIG)")
export WP_H=$(python3   -c "print(0.125*$SIG)")
export WP_DT=$(python3  -c "print(0.00125*$SIG*$SIG)")
export WP_NSTEPS=1500
export WP_SAVE_EVERY=750                               # 3 density frames
export WP_WF_EVERY=1500                                # first + last wavefunction

echo "=============================================="
echo " job     : ${SLURM_JOB_ID:-none}  wp-si-sweep  task $IDX"
echo " node    : $(hostname)"
echo " gpu     : $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>&1 | head -1)"
echo " sigma   : $WP_SIGMA Bohr"
echo " box     : ${WP_L}^3 Bohr @ h=$WP_H  -> $(python3 -c "print(round(float('$WP_L')/float('$WP_H')))")^3 grid"
echo " time    : $WP_NSTEPS x $WP_DT = $(python3 -c "print(float('$WP_DT')*1500)") a.u."
echo " started : $(date)"
echo "=============================================="

rc=0
for CFG in noninteracting:none hartree:none lda:none lda:pzrun; do
  TH="${CFG%%:*}"; SICM="${CFG##*:}"
  OUTN="sweep_s${STAG}_${TH}"
  [ "$SICM" != "none" ] && OUTN="sweep_s${STAG}_sic_pzrun"
  echo
  echo "########## sigma=$SIG  THEORY=$TH  SIC=$SICM  OUT=$OUTN ##########"
  if [ -z "${WP_FORCE:-}" ] && \
     grep -q "run_completed = true" "results/${OUTN}/run_summary.txt" 2>/dev/null; then
    echo "  [skip] results/${OUTN} already complete"; continue
  fi
  export WP_THEORY="$TH" WP_SIC="$SICM" WP_OUT="$OUTN"
  t0=$(date +%s)
  ./run || rc=1
  echo "  [$TH/$SICM] wall=$(( $(date +%s) - t0 )) s"
done

echo
echo "=============================================="
echo " sigma $SIG exit status : $rc"
echo " finished    : $(date)"
echo "=============================================="
exit $rc
