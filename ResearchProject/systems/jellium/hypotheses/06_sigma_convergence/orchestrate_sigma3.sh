#!/usr/bin/env bash
# ============================================================================
# σ=3.0 large-width probe — 7-velocity S(v) ladder (plan §12).
#
# Self-contained: runs WITHOUT any session alive. Runs the σ=3 ladder
# {3.0,2.0,1.3,1.0,0.8,0.6,0.2} split across 2 GPUs, then emails BOTH the
# velocity- and energy-axis cumulative S(v) plots (vs the single point-charge
# Lindhard, σ=3 overlaid alongside σ={0.15,0.25,0.35,0.5}) + the Method A/B
# table via sigma_sweep_report.py --email.
#
# v0=1.0 is the new peak-refinement point between v0=0.8 and v0=1.3.
# One dedicated build (launch_z=-13, the standard 4σ rule for σ=3).
#
# Launch (detached, survives session exit):
#   setsid nohup bash orchestrate_sigma3.sh </dev/null >orchestrate_sigma3.log 2>&1 &
# ============================================================================
set -u

JELLIUM="/local/data/public/skcb2/tddft/ResearchProject/systems/jellium"
HERE="$JELLIUM/hypotheses/06_sigma_convergence"
VENV="/local/data/public/skcb2/tddft/venv/bin/python3"
RUNDIR="$JELLIUM/run_classical_n162_L50_sv_sigma3p0"
BIN="$RUNDIR/run"
LOG="$HERE/orchestrate_sigma3.log"

export PATH="/local/data/public/skcb2/tddft/shared/bin:$PATH"
export INQ_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share
export PSEUDOPOD_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share/pseudopod

log(){ echo "[$(date -u +%FT%TZ)] $*" >> "$LOG"; }

# velocity tag -> "v0:nsteps" (N_STEPS reused from the existing sweep; v1p0=700)
declare -A VEL=( [v3p0]="3.0:300" [v2p0]="2.0:450" [v1p3]="1.3:700"
                 [v1p0]="1.0:700" [v0p8]="0.8:700" [v0p6]="0.6:700"
                 [v0p2]="0.2:1000" )
# GPU split, balanced by total steps (sum 4550): GPU0 2400, GPU1 2150
GPU0_VT=(v0p2 v1p3 v0p6)          # 1000+700+700 = 2400
GPU1_VT=(v0p8 v1p0 v2p0 v3p0)     #  700+700+450+300 = 2150

run_one(){  # $1=gpu $2=vtag
  local gpu="$1" vt="$2"
  local v0="${VEL[$vt]%%:*}" ns="${VEL[$vt]##*:}"
  local vti=$(( ns/5 )); [ "$vti" -lt 1 ] && vti=1
  local obs=$(( ns/80 )); [ "$obs" -lt 1 ] && obs=1
  local out="$RUNDIR/results/$vt"
  mkdir -p "$out"
  log "START gpu$gpu sigma=3.0 v=$v0 nsteps=$ns vti_every=$vti -> $out"
  CUDA_VISIBLE_DEVICES="$gpu" PROJ_V0="$v0" SV_N_STEPS="$ns" \
    SV_VTI_EVERY="$vti" SV_OBS_EVERY="$obs" \
    SV_OUT_ROOT="$out" "$BIN" >> "$out/run.log" 2>&1
  log "DONE  gpu$gpu sigma=3.0 v=$v0 rc=$?"
}

run_queue(){  # $1=gpu $2...=vtags
  local gpu="$1"; shift
  for vt in "$@"; do run_one "$gpu" "$vt"; done
}

email_blocked(){  # $1=reason
  "$VENV" - "$1" <<'PY' 2>>"$LOG" || true
import sys
sys.path.insert(0,"/local/data/public/skcb2/tddft/inq-stack/python")
from inqview.email import send_run_email
send_run_email("[sigma-convergence] BLOCKED — sigma=3 ladder aborted",
               "The σ=3 ladder could not run.\n\nReason: %s\n" % sys.argv[1],
               to="chiddukanna@gmail.com")
PY
}

log "=== sigma=3 orchestrator start ==="
if [ ! -x "$BIN" ]; then
  log "FATAL: binary $BIN missing/not executable"
  email_blocked "run binary missing at $BIN (build/smoke did not complete)"
  exit 2
fi

run_queue 0 "${GPU0_VT[@]}" &
P0=$!
run_queue 1 "${GPU1_VT[@]}" &
P1=$!
wait "$P0" "$P1"
log "--- sigma=3 ladder complete; building + emailing S(v) (both figures) ---"
cd "$HERE"
"$VENV" sigma_sweep_report.py --email "3.0" >> "$LOG" 2>&1 \
  || log "WARN sigma_sweep_report (email) rc=$?"
log "=== sigma=3 orchestrator finished ==="
