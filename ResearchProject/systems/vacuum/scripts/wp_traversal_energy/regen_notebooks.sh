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
  PY $CAMP/analyse.py results/$tag --label wp --dt 0.02 \
     --slab-face 1000000 --cap-inner 30 --cap-lines 30,40 \
     --title "vacuum WP $tag (80-Bohr, 10σ clear)" --per-frame-norm-wp
done

echo "== setup figure (one-sided CAP, real t=0 density) =="
PY make_setup_figure.py results/cap results/cap/report/setup_vacuum_cap.png \
   --cap-inner 30 --cap-outer 40 --wp-launch -30 \
   --title "Vacuum WP-CAP run — setup (t = 0, 80-Bohr box)"

echo "== comparison notebook =="
PY compare_notebook.py results/nocap results/cap results/comparison \
   --dt 0.02 --cap-inner 30 --cap-outer 40

echo "== DONE =="
ls -la results/comparison/nocap_vs_cap_comparison.ipynb results/*/report/run_report.ipynb
