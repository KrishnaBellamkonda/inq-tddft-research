#!/usr/bin/env bash
# Assemble docs/reports/report1/code/ from the manifest.
# Idempotent: safe to re-run. Pure copy operation, no edits.
set -euo pipefail

ROOT=/local/data/public/skcb2/tddft
STAGE="$ROOT/docs/reports/report1/code"

JELLIUM_RUNS=(
  run_wp_n162_L50_E20_sigma1_v2
  run_wp_n162_L50_E25_sigma1_v2
  run_wp_n162_L50_E50_sigma1_v2
  run_wp_n162_L50_E100_sigma1_v2
  run_wp_n162_L50_E200_sigma1_v2
  run_wp_n162_L50_E300_sigma1_v2
  run_wp_n162_L50_E50_v2
  run_wp_n162_L50_E100_v2
  run_wp_n162_L50_E300_v2
  run_wp_n162_L50_E600_v2
  run_wp_n162_L50_E100_sigma0p5
  run_wp_n162_L50_E100_sigma3
  run_wp_n162_L50_E100_sigma8
  run_classical_n162_L50_E20
  run_classical_n162_L50_E25
  run_classical_n162_L50_E50_v2
  run_classical_n162_L50_E100_v2
  run_classical_n162_L50_E600_v2
  run_free_wp_L50_E25_sigma1_v2
  run_plasmon_n162_L50_E15
  run_wp_n162_L30_E50_highdens_sigma1_v2
  run_wp_n162_L30_E100_highdens_sigma1_v2
  run_wp_n162_L30_E200_highdens_sigma1_v2
  run_wp_n162_L30_E300_highdens_sigma1_v2
  run_wp_n162_L30_E100_highdens
  run_classical_n162_L30_E50_highdens
  run_classical_n162_L30_E100_highdens
  run_classical_n162_L30_E200_highdens
  run_classical_n162_L30_E300_highdens
)

CORONENE_RUNS=(
  run_save_gs_paper_replica
  run_propagate_paper_replica
  run_cc_bond
)

# Lean-tier file copier for a single run dir.
# Copies: run.cpp, analyse.py (if present), results/run_summary.txt (if present),
#         REPORT.md (if present).
copy_run() {
  local src="$1"
  local dst="$2"
  mkdir -p "$dst"
  [[ -f "$src/run.cpp" ]] && cp "$src/run.cpp" "$dst/run.cpp"
  [[ -f "$src/analyse.py" ]] && cp "$src/analyse.py" "$dst/analyse.py"
  [[ -f "$src/REPORT.md" ]] && cp "$src/REPORT.md" "$dst/REPORT.md"
  if [[ -f "$src/results/run_summary.txt" ]]; then
    mkdir -p "$dst/results"
    cp "$src/results/run_summary.txt" "$dst/results/run_summary.txt"
  fi
}

# ----- inq-stack library (full copy of headers + python sources) -----
mkdir -p "$STAGE/inq-stack"
rsync -a --delete \
  --include='*/' \
  --include='*.hpp' --include='*.h' \
  --exclude='*' \
  "$ROOT/inq-stack/include/" "$STAGE/inq-stack/include/"

rsync -a --delete \
  --include='*/' \
  --include='*.py' --include='*.toml' --include='*.cfg' --include='setup.py' \
  --exclude='__pycache__' --exclude='*.pyc' --exclude='*.egg-info' \
  --exclude='*' \
  "$ROOT/inq-stack/python/" "$STAGE/inq-stack/python/"

# Top-level inq-stack metadata (if any)
for f in pyproject.toml setup.py setup.cfg README.md; do
  [[ -f "$ROOT/inq-stack/$f" ]] && cp "$ROOT/inq-stack/$f" "$STAGE/inq-stack/$f"
done

# ----- Jellium runs -----
JEL_SRC="$ROOT/ResearchProject/systems/jellium"
JEL_DST="$STAGE/ResearchProject/systems/jellium"
mkdir -p "$JEL_DST/shared/cpp" "$JEL_DST/shared/configs"

# Shared infrastructure
cp "$JEL_SRC/shared/cpp/run_template.hpp" "$JEL_DST/shared/cpp/"
[[ -f "$JEL_SRC/shared/cpp/eigenvalues_writer.hpp" ]] && cp "$JEL_SRC/shared/cpp/eigenvalues_writer.hpp" "$JEL_DST/shared/cpp/"
[[ -f "$JEL_SRC/shared/cpp/results_paths.hpp" ]] && cp "$JEL_SRC/shared/cpp/results_paths.hpp" "$JEL_DST/shared/cpp/"

# All cited Cfg headers (copy the whole configs dir for simplicity — small files, gives examiner full view)
cp "$JEL_SRC/shared/configs/"*.hpp "$JEL_DST/shared/configs/" 2>/dev/null || true

for run in "${JELLIUM_RUNS[@]}"; do
  if [[ -d "$JEL_SRC/$run" ]]; then
    copy_run "$JEL_SRC/$run" "$JEL_DST/$run"
  else
    echo "WARNING: jellium run not found: $run" >&2
  fi
done

# ----- Coronene runs -----
COR_SRC="$ROOT/ResearchProject/systems/coronene"
COR_DST="$STAGE/ResearchProject/systems/coronene"
mkdir -p "$COR_DST/shared/cpp" "$COR_DST/shared/configs"

cp "$COR_SRC/shared/cpp/run_template.hpp" "$COR_DST/shared/cpp/"
[[ -f "$COR_SRC/shared/cpp/eigenvalues_writer.hpp" ]] && cp "$COR_SRC/shared/cpp/eigenvalues_writer.hpp" "$COR_DST/shared/cpp/"
[[ -f "$COR_SRC/shared/cpp/results_paths.hpp" ]] && cp "$COR_SRC/shared/cpp/results_paths.hpp" "$COR_DST/shared/cpp/"
[[ -f "$COR_SRC/shared/cpp/leed_screen_layout.hpp" ]] && cp "$COR_SRC/shared/cpp/leed_screen_layout.hpp" "$COR_DST/shared/cpp/"

cp "$COR_SRC/shared/configs/cc_bond_35x35x60.hpp" "$COR_DST/shared/configs/" 2>/dev/null || true

for run in "${CORONENE_RUNS[@]}"; do
  if [[ -d "$COR_SRC/$run" ]]; then
    copy_run "$COR_SRC/$run" "$COR_DST/$run"
  else
    echo "WARNING: coronene run not found: $run" >&2
  fi
done

# ----- Draft5 figure scripts -> staging inqview/report1/ -----
DRAFT5_SCRIPTS="$ROOT/docs/reports/report1/drafts/draft5/scripts"
STAGE_REPORT1="$STAGE/inq-stack/python/inqview/report1"
mkdir -p "$STAGE_REPORT1"
if [[ -d "$DRAFT5_SCRIPTS" ]]; then
  find "$DRAFT5_SCRIPTS" -maxdepth 1 -name 'make_fig*.py' -exec cp {} "$STAGE_REPORT1/" \;
  # Also copy the shared style file if it's in drafts and not already in report1/
  [[ -f "$DRAFT5_SCRIPTS/_shared_style.py" ]] && cp "$DRAFT5_SCRIPTS/_shared_style.py" "$STAGE_REPORT1/"
fi

echo "Done. Staging at: $STAGE"
echo "Summary:"
find "$STAGE" -type f | wc -l | xargs -I {} echo "  Files: {}"
du -sh "$STAGE" | awk '{print "  Size:  " $1}'
