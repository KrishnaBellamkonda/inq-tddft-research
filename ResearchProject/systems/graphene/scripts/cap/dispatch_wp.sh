#!/bin/bash
# Sequential WP campaign on GPU 1, invoking the PREBUILT ./run (no rebuild).
# 4 WP runs: 2 no-CAP baselines + centroid&channeling with CAP.
set -u
CAPDIR=/local/data/public/skcb2/tddft/ResearchProject/systems/graphene/scripts/cap
BASE=/local/data/public/skcb2/tddft/ResearchProject/systems/graphene/cap_scattering
VENV=/local/data/public/skcb2/tddft/venv/bin/python3
BIN="$CAPDIR/run_wp"
HX=4.6655; HY=-2.6840           # channeling hollow-site target (Bohr)
export INQ_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share
export PSEUDOPOD_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share/pseudopod
export PATH=/local/data/public/skcb2/tddft/shared/bin:$PATH
cd "$CAPDIR"

run() {  # name cap cx cy
  local name=$1 cap=$2 cx=$3 cy=$4
  local out="$BASE/run_$name"
  mkdir -p "$out"
  echo "[$(date +%H:%M:%S)] START $name (CAP=$cap cx=$cx cy=$cy) -> $out"
  GR_OUTDIR="$out/results" GR_CAP=$cap GR_CX=$cx GR_CY=$cy GR_TAG=$name GR_E_EV=100 \
    CUDA_VISIBLE_DEVICES=1 "$BIN" > "$out/run.log" 2>&1
  local rc=$?
  echo "[$(date +%H:%M:%S)] END $name rc=$rc"
  PYTHONPATH=/local/data/public/skcb2/tddft/inq-stack/python "$VENV" "$CAPDIR/post_and_email.py" "$out" "$name" >> "$out/post.log" 2>&1
}

run wp_centroid_nocap    0 0    0
run wp_channeling_nocap  0 $HX  $HY
run wp_centroid_cap      1 0    0
run wp_channeling_cap    1 $HX  $HY
echo "[$(date +%H:%M:%S)] WP CAMPAIGN COMPLETE"

# auto-build: regenerate the sweep study notebook from all completed runs
PYTHONPATH=/local/data/public/skcb2/tddft/inq-stack/python "$VENV" \
  "$BASE/../hypotheses/cap_scattering/build_report.py" || echo "notebook build failed (non-fatal)"
