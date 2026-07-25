#!/bin/bash
# ============================================================================
# Grazing / impact-parameter b-scan: classical + WP projectile grazing a finite
# coronene C24H12 flake (reoriented into the y-z plane). Projectile fires +z at
# x = b (impact parameter, perpendicular distance to the flake plane). CAP on the
# z-faces (same as the perpendicular runs). 6 runs: {cl,wp} x b={1,3,6} Bohr.
#
# Uses FROZEN binaries (no rebuild per run):
#   scripts/grazing/cl/run  (classical, He-symbol z_val=-1 projectile)
#   scripts/grazing/wp/run   (wave packet)
# Run on a free GPU AFTER both binaries are built+smoked. ~1-1.5 h/run -> overnight.
# All CAP results PROVISIONAL until inq-study Task #7.
# ============================================================================
set -u
GZ=/local/data/public/skcb2/tddft/ResearchProject/systems/graphene/scripts/grazing
BASE=/local/data/public/skcb2/tddft/ResearchProject/systems/graphene/grazing
VENV=/local/data/public/skcb2/tddft/venv/bin/python3
GPU=${GPU:-0}
export INQ_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share
export PSEUDOPOD_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share/pseudopod
export PATH=/local/data/public/skcb2/tddft/shared/bin:$PATH
mkdir -p "$BASE"

run(){ local mode=$1 b=$2; local bin="$GZ/$mode/run"; local name="${mode}_b${b}"
  local out="$BASE/run_$name"; mkdir -p "$out"
  if [ ! -x "$bin" ]; then echo "[skip] $name: binary $bin missing"; return; fi
  # RESUMABLE: skip runs already finished (run_completed=true); a partial/killed
  # run has no such marker -> it reruns cleanly (writers use overwrite=true).
  if grep -q 'run_completed = true' "$out/results/run_summary.txt" 2>/dev/null; then
    echo "[skip] $name: already completed"; return; fi
  echo "[$(date +%H:%M:%S)] START $name (b=$b Bohr) -> $out"
  GR_OUTDIR="$out/results" GR_CAP=1 GR_CX=$b GR_CY=0 GR_SEED=0 GR_TAG=$name GR_E_EV=100 \
    CUDA_VISIBLE_DEVICES=$GPU "$bin" > "$out/run.log" 2>&1
  echo "[$(date +%H:%M:%S)] END $name rc=$?"
}

for b in 1 3 6; do run cl $b; done     # classical: clean dE/dx vs b
for b in 1 3 6; do run wp $b; done     # wave packet: dispersive comparison
echo "[$(date +%H:%M:%S)] GRAZING b-SCAN COMPLETE"

# auto-build the grazing study notebook (skill convention)
PYTHONPATH=/local/data/public/skcb2/tddft/inq-stack/python "$VENV" \
  "$GZ/../../hypotheses/grazing/build_grazing_report.py" 2>/dev/null \
  || echo "grazing notebook builder not present yet (non-fatal)"
