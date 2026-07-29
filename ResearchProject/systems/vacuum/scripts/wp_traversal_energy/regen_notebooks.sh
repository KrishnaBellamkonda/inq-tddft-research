#!/usr/bin/env bash
# Regenerate all figures/notebooks for the corrected 80-Bohr vacuum runs:
# setup figure (one-sided CAP), both per-run notebooks (per-frame-norm WP GIF,
# one-sided CAP lines), and the no-CAP-vs-CAP comparison notebook.
set -euo pipefail
cd "$(dirname "$0")"
ROOT=/local/data/public/skcb2/tddft
PY() { PYTHONPATH=$ROOT/inq-stack/python $ROOT/venv/bin/python3 "$@"; }
CAMP=$ROOT/ResearchProject/systems/localised_jellium/scripts/wp_cap_energy_plateau

for tag in nocap cap; do
  echo "== per-run notebook: $tag =="
  mkdir -p results/$tag/report
  PY $CAMP/analyse.py results/$tag --label wp --dt 0.01 \
     --slab-face 1000000 --cap-inner 7.5 --cap-lines 7.5,22.5 \
     --title "vacuum WP $tag (low-spread σ=3, E=400eV, 30×30×45)" --per-frame-norm-wp
done

echo "== setup figure (one-sided CAP, real t=0 density) =="
PY make_setup_figure.py results/cap results/cap/report/setup_vacuum_cap.png \
   --cap-inner 7.5 --cap-outer 22.5 --wp-launch -7.5 \
   --title "Vacuum WP-CAP setup (low-spread σ=3, E=400eV, t=0)"

echo "== comparison notebook =="
PY compare_notebook.py results/nocap results/cap results/comparison \
   --dt 0.01 --cap-inner 7.5 --cap-outer 22.5

echo "== DONE =="
ls -la results/comparison/nocap_vs_cap_comparison.ipynb results/*/report/run_report.ipynb
