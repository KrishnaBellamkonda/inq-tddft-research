#!/bin/bash
# ============================================================================
# PERPENDICULAR arm of the impact-parameter comparison: classical + WP projectile
# hitting a finite coronene C24H12 flake HEAD-ON (flake in the x-y plane, beam +z
# perpendicular). SAME box (20x22x60) / grid (50 Ha) / sim-time (N=1319) / CAP
# (L=20, eta=-0.5) / projectile (E=100 eV, sigma=1.47) as the grazing arm — the
# ONLY difference is the flake orientation. Directly comparable.
#
# Uses FROZEN binaries scripts/perp/{cl,wp}/run (built+smoked first).
# Default: b=0 (through the flake centre). Extra impact parameters can be added.
# GPU 1 by default. All CAP results PROVISIONAL until inq-study Task #7.
# ============================================================================
set -u
PP=/local/data/public/skcb2/tddft/ResearchProject/systems/graphene/scripts/perp
BASE=/local/data/public/skcb2/tddft/ResearchProject/systems/graphene/perp
VENV=/local/data/public/skcb2/tddft/venv/bin/python3
GPU=${GPU:-1}
BLIST=${BLIST:-0}                 # space-separated impact parameters (lateral x-offset)
export INQ_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share
export PSEUDOPOD_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share/pseudopod
export PATH=/local/data/public/skcb2/tddft/shared/bin:$PATH
mkdir -p "$BASE"

run(){ local mode=$1 b=$2; local bin="$PP/$mode/run"; local name="${mode}_b${b}"
  local out="$BASE/run_$name"; mkdir -p "$out"
  if [ ! -x "$bin" ]; then echo "[skip] $name: binary $bin missing"; return; fi
  echo "[$(date +%H:%M:%S)] START $name (perp, b=$b) -> $out"
  GR_OUTDIR="$out/results" GR_CAP=1 GR_CX=$b GR_CY=0 GR_SEED=0 GR_TAG=$name GR_E_EV=100 \
    CUDA_VISIBLE_DEVICES=$GPU "$bin" > "$out/run.log" 2>&1
  echo "[$(date +%H:%M:%S)] END $name rc=$?"
}

for b in $BLIST; do run cl $b; done
for b in $BLIST; do run wp $b; done
echo "[$(date +%H:%M:%S)] PERPENDICULAR arm COMPLETE"

PYTHONPATH=/local/data/public/skcb2/tddft/inq-stack/python "$VENV" \
  "$PP/../../hypotheses/grazing/build_grazing_report.py" 2>/dev/null \
  || echo "comparison notebook builder not present yet (non-fatal)"
