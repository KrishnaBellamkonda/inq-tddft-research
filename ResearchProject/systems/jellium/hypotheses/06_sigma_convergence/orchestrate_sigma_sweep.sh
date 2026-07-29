#!/usr/bin/env bash
# ============================================================================
# Autonomous σ-convergence stopping-power sweep.
#
# Self-contained: runs WITHOUT any Claude session alive. For each σ ∈
# {0.15,0.25,0.35} it runs the 6-velocity ladder split across 2 GPUs, then emails
# the cumulative S(v) plot (vs the single point-charge Lindhard) via
# sigma_sweep_report.py --email. Base observables (ADR-0006 minimum set) are saved
# for every run; the full derived-observable pipeline (analyse.py) is intentionally
# deferred to post-sweep (not in the overnight critical path).
#
# Launch (detached, fires at 4am):
#   secs=$(( $(date -d 'tomorrow 04:00' +%s) - $(date +%s) ))   # or today 04:00
#   setsid nohup bash -c "sleep $secs; exec .../orchestrate_sigma_sweep.sh" \
#       </dev/null >.../orchestrate.log 2>&1 &
# ============================================================================
set -u

JELLIUM="/local/data/public/skcb2/tddft/ResearchProject/systems/jellium"
HERE="$JELLIUM/hypotheses/06_sigma_convergence"
VENV="/local/data/public/skcb2/tddft/venv/bin/python3"
PSP="$JELLIUM/shared/pseudopotentials"
BIN="$JELLIUM/run_classical_n162_L50_sv_sigma0p15/run"   # one shared binary
LOG="$HERE/orchestrate.log"

export PATH="/local/data/public/skcb2/tddft/shared/bin:$PATH"
export INQ_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share
export PSEUDOPOD_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share/pseudopod

log(){ echo "[$(date -u +%FT%TZ)] $*" >> "$LOG"; }

# σ tag -> UPF
declare -A SIG_UPF=(
  [0p15]="$PSP/electron_gaussian_sigma0p15.upf"
  [0p25]="$PSP/electron_gaussian_sigma0p25.upf"
  [0p35]="$PSP/electron_gaussian_sigma0p35.upf"
)
declare -A SIG_DIR=(
  [0p15]="$JELLIUM/run_classical_n162_L50_sv_sigma0p15"
  [0p25]="$JELLIUM/run_classical_n162_L50_sv_sigma0p25"
  [0p35]="$JELLIUM/run_classical_n162_L50_sv_sigma0p35"
)
declare -A SIG_LABEL=( [0p15]="0.15" [0p25]="0.25" [0p35]="0.35" )

# velocity tag -> "v0:nsteps"
VTAGS=(v3p0 v2p0 v1p3 v0p8 v0p6 v0p2)
declare -A VEL=( [v3p0]="3.0:300" [v2p0]="2.0:450" [v1p3]="1.3:700"
                 [v0p8]="0.8:700" [v0p6]="0.6:700" [v0p2]="0.2:1000" )
# GPU split (balance total steps): GPU0 heavy-low-v, GPU1 the rest
GPU0_VT=(v0p2 v0p8 v3p0)   # 1000+700+300 = 2000
GPU1_VT=(v0p6 v1p3 v2p0)   # 700+700+450 = 1850

run_one(){  # $1=gpu $2=sigtag $3=vtag
  local gpu="$1" sig="$2" vt="$3"
  local v0="${VEL[$vt]%%:*}" ns="${VEL[$vt]##*:}"
  local vti=$(( ns/5 )); [ "$vti" -lt 1 ] && vti=1
  local obs=$(( ns/80 )); [ "$obs" -lt 1 ] && obs=1
  local out="${SIG_DIR[$sig]}/results/$vt"
  mkdir -p "$out"
  log "START gpu$gpu sigma=$sig v=$v0 nsteps=$ns vti_every=$vti -> $out"
  CUDA_VISIBLE_DEVICES="$gpu" PROJ_V0="$v0" SV_N_STEPS="$ns" \
    SV_VTI_EVERY="$vti" SV_OBS_EVERY="$obs" SV_PSEUDO="${SIG_UPF[$sig]}" \
    SV_OUT_ROOT="$out" "$BIN" >> "$out/run.log" 2>&1
  log "DONE  gpu$gpu sigma=$sig v=$v0 rc=$?"
}

run_queue(){  # $1=gpu $2=sigtag $3...=vtags
  local gpu="$1" sig="$2"; shift 2
  for vt in "$@"; do run_one "$gpu" "$sig" "$vt"; done
}

email_blocked(){  # $1=reason
  "$VENV" - "$1" <<'PY' 2>>"$LOG" || true
import sys
sys.path.insert(0,"/local/data/public/skcb2/tddft/inq-stack/python")
from inqview.email import send_run_email
send_run_email("[sigma-convergence] BLOCKED — sweep aborted",
               "The autonomous σ-sweep could not run.\n\nReason: %s\n" % sys.argv[1],
               to="chiddukanna@gmail.com")
PY
}

log "=== orchestrator start ==="
if [ ! -x "$BIN" ]; then
  log "FATAL: binary $BIN missing/not executable"
  email_blocked "run binary missing at $BIN (build/smoke did not complete)"
  exit 2
fi

for sig in 0p15 0p25 0p35; do
  log "--- sigma $sig ladder begin ---"
  run_queue 0 "$sig" "${GPU0_VT[@]}" &
  P0=$!
  run_queue 1 "$sig" "${GPU1_VT[@]}" &
  P1=$!
  wait "$P0" "$P1"
  log "--- sigma $sig ladder complete; building + emailing S(v) ---"
  cd "$HERE"
  "$VENV" sigma_sweep_report.py --email "${SIG_LABEL[$sig]}" >> "$LOG" 2>&1 \
    || log "WARN sigma_sweep_report (email) rc=$?"
done

log "=== all sigma ladders complete ==="
"$VENV" "$HERE/sigma_sweep_report.py" --email "ALL (final)" >> "$LOG" 2>&1 || true
log "=== orchestrator finished ==="
