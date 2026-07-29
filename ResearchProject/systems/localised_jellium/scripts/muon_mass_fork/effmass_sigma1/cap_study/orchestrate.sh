#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Autonomous CAP-parameter study for the σ=1 chirped WP run (2026-07-11).
# Fable-5 advisor design; Opus-approved. Three variants of the SAME simulation,
# differing ONLY in the CAP, to diagnose a suspected CAP problem. Decision metric
# is the TOTAL electron number N(t)=∫n dV (user: absorption must be judged on the
# whole-cell electron count, not the WP orbital norm alone).
#
#   R1 cap_gap19p5 : η=-1.0, centre |z|=36, width 8  → region [32,40], gap 19.5  (protect the wake)
#   R2 cap_eta0p4  : η=-0.4, centre |z|=32.5, width 15 → region [25,40], gap 12.5 (weak-η branch)
#   R3 cap_eta2p0  : η=-2.0, centre |z|=32.5, width 15 → region [25,40], gap 12.5 (strong-η branch)
#
# Reuses the parametrised binary  effmass_sigma1/wp/run  (already built) and the
# shared GS. Runs on 2 GPUs: R1(gpu0)+R2(gpu1) concurrently, R3 on the first freed
# GPU. Then builds a run-notebook per variant with the updated builder. Detached
# (setsid) so it survives logout; checkpointed every 300 steps (extension-ready).
# ---------------------------------------------------------------------------
set -u
ROOT=/local/data/public/skcb2/tddft
WPDIR=$ROOT/ResearchProject/systems/localised_jellium/scripts/muon_mass_fork/effmass_sigma1/wp
NBDIR=$ROOT/ResearchProject/systems/localised_jellium/hypotheses/muon_mass_fork
BUILDER=$ROOT/.claude/skills/run-notebook/run_notebook_builder.py
PY=$ROOT/venv/bin/python3
LOG=$ROOT/ResearchProject/systems/localised_jellium/scripts/muon_mass_fork/effmass_sigma1/cap_study/orchestrate.log
export INQ_SOURCE=$ROOT/inq-study
export PYTHONPATH=$ROOT/inq-stack/python
export MPLBACKEND=Agg
COMMON="EM_DT=0.04 EM_N_STEPS=900 EM_CKPT_EVERY=300"

exec >>"$LOG" 2>&1
echo "==== CAP study orchestrator start $(date) ===="

run_variant () {  # $1=gpu  $2=EM_OUT  $3=eta  $4=center  $5=width
  local gpu=$1 out=$2 eta=$3 cen=$4 wid=$5
  echo ">> launch $out on GPU$gpu : eta=$eta center=$cen width=$wid  $(date)"
  ( cd "$WPDIR" && env CUDA_VISIBLE_DEVICES=$gpu $COMMON \
      EM_OUT="$out" EM_CAP_ETA="$eta" EM_CAP_CENTER_BOHR="$cen" EM_CAP_WIDTH_BOHR="$wid" \
      ./run > "rt_${out}.log" 2>&1 )
  echo ">> DONE $out (exit $?)  $(date)"
}

# ---- launch R1 (gpu0) + R2 (gpu1) concurrently -----------------------------
run_variant 0 cap_gap19p5 -1.0 36.0 8.0  &  P1=$!
run_variant 1 cap_eta0p4  -0.4 32.5 15.0 &  P2=$!

# ---- wait for the first to free a GPU, then launch R3 on it -----------------
FREE_GPU=""
while :; do
  if ! kill -0 $P1 2>/dev/null; then FREE_GPU=0; break; fi
  if ! kill -0 $P2 2>/dev/null; then FREE_GPU=1; break; fi
  sleep 30
done
echo ">> a GPU freed (gpu$FREE_GPU); launching R3  $(date)"
run_variant "$FREE_GPU" cap_eta2p0 -2.0 32.5 15.0 &  P3=$!

wait $P1 $P2 $P3
echo "==== all three runs finished $(date) ===="

# ---- build a run-notebook per variant --------------------------------------
build_nb () {  # $1=EM_OUT  $2=cap_inner_bohr  $3=nb_stem
  local out=$1 capinner=$2 stem=$3
  echo ">> notebook $stem (cap_inner=$capinner)  $(date)"
  "$PY" "$BUILDER" "$WPDIR/results/$out" "$NBDIR/$stem.ipynb" \
      --run-cpp "$WPDIR/run.cpp" --cap-inner "$capinner" --rs 5.684 --proj-sigma 0.7071 \
      --launch-z -16.5 --v0 2.7111 --e-gs-ha -36.9404590471 --l-slab 25.0 --gif-seconds 17
  echo ">> notebook $stem done (exit $?)"
}
build_nb cap_gap19p5 32.0 effmass_sigma1_cap_gap19p5_wp_run
build_nb cap_eta0p4  25.0 effmass_sigma1_cap_eta0p4_wp_run
build_nb cap_eta2p0  25.0 effmass_sigma1_cap_eta2p0_wp_run
echo "==== CAP study orchestrator COMPLETE $(date) ===="
