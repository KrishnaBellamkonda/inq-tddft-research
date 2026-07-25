#!/bin/bash
# 6 classical-ensemble runs on GPU 1 (run AFTER WP campaign frees the GPU and
# AFTER run_cl is built+frozen). centroid x3 + channeling x3, seeds 1/2/3, CAP on.
set -u
CAPDIR=/local/data/public/skcb2/tddft/ResearchProject/systems/graphene/scripts/cap
BASE=/local/data/public/skcb2/tddft/ResearchProject/systems/graphene/cap_scattering
VENV=/local/data/public/skcb2/tddft/venv/bin/python3
BIN="$CAPDIR/run_cl"
HX=4.6655; HY=-2.6840
export INQ_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share
export PSEUDOPOD_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share/pseudopod
export PATH=/local/data/public/skcb2/tddft/shared/bin:$PATH
cd "$CAPDIR"
run(){ local name=$1 cx=$2 cy=$3 seed=$4; local out="$BASE/run_$name"; mkdir -p "$out"
  echo "[$(date +%H:%M:%S)] START $name seed=$seed -> $out"
  GR_OUTDIR="$out/results" GR_CAP=1 GR_CX=$cx GR_CY=$cy GR_SEED=$seed GR_TAG=$name GR_E_EV=100 \
    CUDA_VISIBLE_DEVICES=1 "$BIN" > "$out/run.log" 2>&1
  echo "[$(date +%H:%M:%S)] END $name rc=$?"
  PYTHONPATH=/local/data/public/skcb2/tddft/inq-stack/python "$VENV" "$CAPDIR/post_and_email.py" "$out" "$name" >> "$out/post.log" 2>&1; }
run cl_centroid_s1   0   0    1
run cl_centroid_s2   0   0    2
run cl_centroid_s3   0   0    3
run cl_channeling_s1 $HX $HY  1
run cl_channeling_s2 $HX $HY  2
run cl_channeling_s3 $HX $HY  3
echo "[$(date +%H:%M:%S)] CLASSICAL CAMPAIGN COMPLETE"

# auto-build: regenerate the sweep study notebook (now incl. classical KE stopping)
PYTHONPATH=/local/data/public/skcb2/tddft/inq-stack/python "$VENV" \
  "$BASE/../hypotheses/cap_scattering/build_report.py" || echo "notebook build failed (non-fatal)"
