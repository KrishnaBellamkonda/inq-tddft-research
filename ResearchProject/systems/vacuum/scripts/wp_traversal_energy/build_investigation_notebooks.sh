#!/usr/bin/env bash
# Build the per-run notebooks for the completed GPU-1 vacuum investigation suite
# (eta-sweep exp1a, partial-absorption exp2, and nocap_long) — the 9 runs that
# had no report/run_report.ipynb. Same builder + one-sided CAP annotation as the
# corrected 80-Bohr vacuum runs (regen_notebooks.sh).
set -uo pipefail
cd "$(dirname "$0")"
ROOT=/local/data/public/skcb2/tddft
PY() { PYTHONPATH=$ROOT/inq-stack/python $ROOT/venv/bin/python3 "$@"; }
CAMP=$ROOT/ResearchProject/systems/localised_jellium/scripts/wp_cap_energy_plateau

# CAP runs: one-sided CAP z in [7.5, 22.5], dt = 0.01, sigma=3, E=400 eV
declare -A ETA=( [exp1a_eta-0.3]=-0.3 [exp1a_eta-0.7]=-0.7 [exp1a_eta-1.0]=-1.0 \
                 [exp1a_eta-2.0]=-2.0 [exp1a_eta-3.5]=-3.5 )
declare -A NABS=( [exp2_N0p1]=0.1 [exp2_N0p3]=0.3 [exp2_N0p5]=0.5 )

fail=0
for run in "${!ETA[@]}"; do
  echo "===== $run (eta=${ETA[$run]}) ====="
  PY $CAMP/analyse.py results/$run --label wp --dt 0.01 \
     --slab-face 1000000 --cap-inner 7.5 --cap-lines 7.5,22.5 \
     --title "vacuum WP CAP eta=${ETA[$run]} (sigma=3, E=400eV, one-sided CAP)" \
     --per-frame-norm-wp || { echo "FAILED: $run"; fail=1; }
done

for run in "${!NABS[@]}"; do
  echo "===== $run (target N_abs=${NABS[$run]}) ====="
  PY $CAMP/analyse.py results/$run --label wp --dt 0.01 \
     --slab-face 1000000 --cap-inner 7.5 --cap-lines 7.5,22.5 \
     --title "vacuum WP partial-absorption N_abs~${NABS[$run]} (sigma=3, E=400eV)" \
     --per-frame-norm-wp || { echo "FAILED: $run"; fail=1; }
done

echo "===== nocap_long (no CAP, 1600 steps) ====="
PY $CAMP/analyse.py results/nocap_long --label wp --dt 0.01 \
   --slab-face 1000000 --cap-inner 7.5 --cap-lines 7.5,22.5 \
   --title "vacuum WP no-CAP long (sigma=3, E=400eV, 1600 steps)" \
   --per-frame-norm-wp || { echo "FAILED: nocap_long"; fail=1; }

echo "===== DONE (fail=$fail) ====="
ls -la results/exp1a_eta-*/report/run_report.ipynb \
       results/exp2_N*/report/run_report.ipynb \
       results/nocap_long/report/run_report.ipynb 2>/dev/null
