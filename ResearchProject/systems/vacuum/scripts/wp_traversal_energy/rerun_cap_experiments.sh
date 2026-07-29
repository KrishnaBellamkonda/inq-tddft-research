#!/usr/bin/env bash
# ============================================================================
# rerun_cap_experiments.sh — two CAP-reflection experiments, one per GPU.
# Common box 30x30x60 [-30,30], sigma0=3, E=400 eV (k0=5.421), h=0.4, dt=0.01.
# Binary already built (env overrides only; no rebuild).
#   GPU0  cap_better  : W=25 (long/adiabatic), eta=-0.7 (GENTLE -> low reflection,
#                       per reflectivity facts); CAP z in [5,30], launch -10.
#   GPU1  cap_fulllen : W=30 = the ENTIRE +z half (full-length one-sided CAP),
#                       eta=-1.0; CAP z in [0,30], launch -15.
# Each -> run notebook (analyse.py) + energy diagnostics (decomposed + compounded,
# with wrap-time marker) appended. Verifies the END energy residual
#   R = E_total(tF) - norm(tF)*E0   (clean CAP -> ~0; eta=-3.5 gave ~416 eV).
# ============================================================================
set -uo pipefail
cd "$(dirname "$0")"
ROOT=/local/data/public/skcb2/tddft
JCAMP=$ROOT/ResearchProject/systems/localised_jellium/scripts/wp_cap_energy_plateau
PY=$ROOT/venv/bin/python3
export INQ_SOURCE=$ROOT/inq-study
export INQ_SHARE_PATH=${INQ_SHARE_PATH:-$ROOT/inq/install/share}
export PSEUDOPOD_SHARE_PATH=${PSEUDOPOD_SHARE_PATH:-$ROOT/inq/install/share/pseudopod}
export PATH="$ROOT/shared/bin:$PATH"
export PYTHONPATH=$ROOT/inq-stack/python
LOG=rerun_cap_experiments.log; : > "$LOG"
say(){ echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }
notify(){ $PY "$JCAMP/notify.py" "$1" "$2" >>"$LOG" 2>&1 || true; }
fail(){ say "FAILED: $1"; notify "CAP experiments FAILED: $1" "See $(pwd)/$LOG"; exit 1; }

say "=== CAP reflection experiments (better vs full-length one-sided) ==="
notify "CAP experiments started" "GPU0 cap_better W=25 eta=-0.7; GPU1 cap_fulllen W=30 eta=-1.0; box 30x30x60; log $(pwd)/$LOG"
rm -rf results/cap_better results/cap_fulllen

say "launching both concurrently..."
CUDA_VISIBLE_DEVICES=0 WP_OUT=cap_better  WP_ETA=-0.7 WP_LZ=60 WP_CAP_L=25 WP_LAUNCH_Z=-10 WP_NSTEPS=1000 ./run > cap_better.run.log  2>&1 &
P0=$!
CUDA_VISIBLE_DEVICES=1 WP_OUT=cap_fulllen WP_ETA=-1.0 WP_LZ=60 WP_CAP_L=30 WP_LAUNCH_Z=-15 WP_NSTEPS=1200 ./run > cap_fulllen.run.log 2>&1 &
P1=$!
wait $P0; R0=$?; wait $P1; R1=$?
[ $R0 -eq 0 ] || fail "cap_better run (exit $R0; cap_better.run.log)"
[ $R1 -eq 0 ] || fail "cap_fulllen run (exit $R1; cap_fulllen.run.log)"

# notebooks + diagnostics per run
for tag in cap_better cap_fulllen; do
  zc=$([ "$tag" = cap_better ] && echo 5 || echo 0)
  say "notebook + diagnostics: $tag (CAP inner z=$zc)..."
  $PY "$JCAMP/analyse.py" results/$tag --label wp --dt 0.01 --slab-face 1000000 \
      --cap-inner $zc --cap-lines $zc,30 --title "vacuum $tag (E=400eV)" --per-frame-norm-wp >>"$LOG" 2>&1 || say "  analyse $tag warn"
  $PY energy_diagnostics.py results/$tag --append >>"$LOG" 2>&1 || say "  diag $tag warn"
done

say "verifying END energy residual R = E_total(tF) - norm(tF)*E0 ..."
$PY - <<'PYV' 2>&1 | tee -a "$LOG"
import pandas as pd, numpy as np
HA=27.211386
for tag,eta,W in (("cap_better",-0.7,25),("cap_fulllen",-1.0,30)):
    en=pd.read_csv(f"results/{tag}/raw/observables/energies.csv")
    rs=pd.read_csv(f"results/{tag}/raw/observables/wp_real_space_stats.csv",comment="#")
    E0=en.total.iloc[0]*HA; Ef=en.total.iloc[-1]*HA
    Nf=rs.norm_check.iloc[-1]
    R=Ef-Nf*E0
    print(f"{tag}: eta={eta} W={W}  E0={E0:.1f}  E_final={Ef:.2f} eV  norm_final={Nf:.2e}  "
          f"residual R={R:.2f} eV  (eta=-3.5 baseline was ~416 eV)")
PYV

say "=== CAP experiments COMPLETE ==="
notify "CAP experiments COMPLETE" "residuals in log. Notebooks:
  results/cap_better/report/run_report.ipynb
  results/cap_fulllen/report/run_report.ipynb
See $(pwd)/$LOG"
say "done."
