#!/usr/bin/env bash
# launch_2gpu.sh — run ONE slab_sE_sweep half across BOTH GPUs via MPI.
#
# Mirrors orchestrate_sE.py's slab_env() env contract EXACTLY, but instead of
# pinning the job to a single CUDA_VISIBLE_DEVICES it launches 2 MPI ranks (one
# per GPU, bound by gpu_rank_wrapper.sh). This gives ~2x the usable VRAM via
# INQ's state/domain decomposition, so large-L_z boxes that UVM-thrash on one
# 24 GB A30 fit in real memory.
#
# Usage:   ./launch_2gpu.sh <E_eV> <L_z> <classical|wp>
# Example: ./launch_2gpu.sh 276 60 classical
#
# Resumable: if results/<half>/checkpoint exists and last_step < target, resumes
# (LJ_RESUME=1) per the final-timestep-checkpoint rule. Detaches via setsid so it
# survives this shell / session; logs to runs/E<E>_Lz<L>/<half>.2gpu.log.
set -euo pipefail

ROOT=/local/data/public/skcb2/tddft
LJ=$ROOT/ResearchProject/systems/localised_jellium
BASE=$LJ/scripts/slab_sE_sweep

EV=${1:?need E_eV};  LZ=${2:?need L_z};  HALF=${3:?need classical|wp}
case "$HALF" in classical) BIN=$BASE/classical/run ;; wp) BIN=$BASE/wp/run ;;
  *) echo "half must be classical|wp"; exit 1 ;; esac

# --- energy(eV) -> v0(a.u.) map (locked design) ---
declare -A V0=( [31]=1.5 [54]=2.0 [100]=2.71 [167]=3.5 [276]=4.5 )
K0=${V0[$EV]:?unknown energy $EV}

# --- design constants (== orchestrate_sE.py) ---
DT=0.04; SIZE_FAC=1.6; WP_DRAIN_AU=60.0
SPACING=0.5; SIGMA=2.0; EDGE_W=1; PER=2; LX=50; LY=50; CAP_ETA=-1.0; CAP_L=10.0
# N(L_z) = round(n0 * Lx * Ly * L_z),  n0 = 0.001288
N=$(python3 -c "print(round(0.001288*$LX*$LY*$LZ))")
HALFW=$(python3 -c "print($LZ/2.0)")
LZBOX=$(python3 -c "print($LZ+60)")
LAUNCH_Z=$(python3 -c "print(-$LZ/2.0-10.0)")
if [ "$HALF" = classical ]; then
  NSTEPS=$(python3 -c "import math;print(int(math.ceil($SIZE_FAC*($LZ+45.0)/$K0/$DT)))")
else
  NSTEPS=$(python3 -c "import math;print(int(math.ceil(($SIZE_FAC*($LZ+45.0)/$K0+$WP_DRAIN_AU)/$DT)))")
fi

# IMPORTANT: a 2-GPU run must load a GS saved in the SAME 2-domain-partition
# layout. INQ's cross-layout checkpoint load (1-partition GS -> 2-partition run)
# corrupts the orbitals (verified 2026-09-18: garbage E_total ~2785 vs 923 Ha).
# So point at the 2-GPU-generated GS. Override with LJ_GS_DIR_OVERRIDE if needed.
GS_DIR=${LJ_GS_DIR_OVERRIDE:-$LJ/shared_gs/slab_sE_sweep/Lz${LZ}_2gpu/checkpoint}
[ -d "$GS_DIR" ] || { echo "FATAL: 2-GPU GS missing: $GS_DIR
  (generate it with a fresh 2-GPU GS SCF; the 1-GPU Lz$LZ/checkpoint will CORRUPT
   under 2 GPUs via INQ cross-layout load.)"; exit 2; }

RUNDIR=$BASE/runs/E${EV}_Lz${LZ}
OUTDIR=$RUNDIR/results/$HALF
mkdir -p "$RUNDIR"

# --- resume decision ---
RESUME=0
if [ -d "$OUTDIR/checkpoint" ] && [ -f "$OUTDIR/rt_state.txt" ]; then
  LAST=$(sed -n 's/^last_step=//p' "$OUTDIR/rt_state.txt" | head -1)
  if [ -n "$LAST" ] && [ "$LAST" -lt "$NSTEPS" ]; then RESUME=1;
    echo "resume: last_step=$LAST -> $NSTEPS"
  elif [ -n "$LAST" ]; then echo "already at target ($LAST>=$NSTEPS); nothing to do"; exit 0; fi
fi

export INQ_SHARE_PATH=$ROOT/inq/install/share
export PSEUDOPOD_SHARE_PATH=$ROOT/inq/install/share/pseudopod
export LJ_LX=$LX LJ_LY=$LY LJ_LZ=$LZBOX LJ_HALF=$HALFW LJ_N=$N \
  LJ_EDGE_W=$EDGE_W LJ_PERIODICITY=$PER LJ_SPACING=$SPACING LJ_SIGMA=$SIGMA \
  LJ_LAUNCH_Z=$LAUNCH_Z LJ_K0=$K0 LJ_N_STEPS=$NSTEPS LJ_DT=$DT \
  LJ_CAP_ETA=$CAP_ETA LJ_CAP_L=$CAP_L LJ_SAVE_EVERY=0 LJ_CKPT_EVERY=200 \
  LJ_GS_DIR=$GS_DIR LJ_OUT=$HALF LJ_RESUME=$RESUME
[ "$HALF" = classical ] && export LJ_MASS=1.0 LJ_DELTA=0.1

LOG=$RUNDIR/$HALF.2gpu.log
{
  echo "===== $(date '+%F %T') 2-GPU MPI  E${EV}_Lz${LZ}/$HALF  N=$N box=${LX}x${LY}x${LZBOX} steps->$NSTEPS resume=$RESUME ====="
} >> "$LOG"

cd "$RUNDIR"
# --bind-to none: let both ranks use all CPU cores (no oversubscription warning);
# wrapper masks each rank to one GPU via OMPI_COMM_WORLD_LOCAL_RANK.
setsid mpirun.openmpi -np 2 --bind-to none \
  -x INQ_SHARE_PATH -x PSEUDOPOD_SHARE_PATH \
  -x LJ_LX -x LJ_LY -x LJ_LZ -x LJ_HALF -x LJ_N -x LJ_EDGE_W -x LJ_PERIODICITY \
  -x LJ_SPACING -x LJ_SIGMA -x LJ_LAUNCH_Z -x LJ_K0 -x LJ_N_STEPS -x LJ_DT \
  -x LJ_CAP_ETA -x LJ_CAP_L -x LJ_SAVE_EVERY -x LJ_CKPT_EVERY -x LJ_GS_DIR \
  -x LJ_OUT -x LJ_RESUME ${LJ_MASS:+-x LJ_MASS} ${LJ_DELTA:+-x LJ_DELTA} \
  "$BASE/gpu_rank_wrapper.sh" "$BIN" >> "$LOG" 2>&1 &
echo "launched PID $! -> $LOG"
