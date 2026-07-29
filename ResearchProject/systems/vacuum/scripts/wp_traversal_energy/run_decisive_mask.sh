#!/usr/bin/env bash
# ============================================================================
# run_decisive_mask.sh — Phase 3 decisive experiment (energy-normalization).
# Same geometry as the production one-sided CAP (sigma0=3, E=400eV, h=0.4, dt=0.01,
# LZ=45, mask/CAP band z in [7.5,22.5], launch -7.5). Binary already built.
#   3a  mask + ETRS  -> norm LOST  (should reproduce the "energy shoots up" artifact)
#   3b  mask + CN    -> norm HELD ~1 (CN renormalises each step) -> should show NO rise
# If 3b's E_reported still rises with norm~1, the normalization hypothesis is refuted.
# Runs on GPU 1. Per-step norm (WP_MOM_EVERY=1). Diagnostics appended to each notebook.
# ============================================================================
set -uo pipefail
cd "$(dirname "$0")"
ROOT=/local/data/public/skcb2/tddft
JCAMP=$ROOT/ResearchProject/systems/localised_jellium/scripts/wp_cap_energy_plateau
PY=$ROOT/venv/bin/python3
export INQ_SOURCE=$ROOT/inq-study CUDA_VISIBLE_DEVICES=1
export INQ_SHARE_PATH=${INQ_SHARE_PATH:-$ROOT/inq/install/share}
export PSEUDOPOD_SHARE_PATH=${PSEUDOPOD_SHARE_PATH:-$ROOT/inq/install/share/pseudopod}
export PATH="$ROOT/shared/bin:$PATH" PYTHONPATH=$ROOT/inq-stack/python
LOG=decisive_mask.log; : > "$LOG"
say(){ echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }
notify(){ $PY "$JCAMP/notify.py" "$1" "$2" >>"$LOG" 2>&1 || true; }

say "=== Phase 3 decisive: mask+ETRS (3a) vs mask+CN (3b), GPU 1 ==="
notify "decisive mask experiment started" "3a mask+ETRS (norm lost), 3b mask+CN (norm held); GPU1; log $(pwd)/$LOG"
rm -rf results/exp3a_mask_etrs results/exp3b_mask_cn results/smoke_me results/smoke_mc

COMMON="WP_SIGMA=3 WP_K0=5.421 WP_LZ=45 WP_LPERP=30 WP_H=0.4 WP_DT=0.01 WP_CAP_L=15 WP_LAUNCH_Z=-7.5 WP_NSTEPS=800 WP_WF_EVERY=20 WP_MOM_EVERY=1"

say "3a: mask + ETRS (norm-losing)..."
env $COMMON WP_ABS=mask WP_PROP=etrs WP_OUT=exp3a_mask_etrs ./run > exp3a.run.log 2>&1 || { say "3a FAILED"; notify "decisive 3a FAILED" "see exp3a.run.log"; exit 1; }
say "3b: mask + CN (norm-preserving)..."
env $COMMON WP_ABS=mask WP_PROP=cn   WP_OUT=exp3b_mask_cn  ./run > exp3b.run.log 2>&1 || { say "3b FAILED"; notify "decisive 3b FAILED" "see exp3b.run.log"; exit 1; }

say "verifying norm behaviour + energy rise..."
$PY - <<'PYV' 2>&1 | tee -a "$LOG"
import pandas as pd, numpy as np
HA=27.211386
for tag,lab in (("exp3a_mask_etrs","mask+ETRS (expect norm LOST, energy rises)"),
                ("exp3b_mask_cn","mask+CN (expect norm HELD ~1, NO rise)")):
    en=pd.read_csv(f"results/{tag}/raw/observables/energies.csv")
    ms=pd.read_csv(f"results/{tag}/raw/observables/wp_momentum_stats.csv",comment='#')
    N=ms.norm_check/ms.norm_check.iloc[0]
    E=en.total.values*HA; E0=E[0]
    print(f"{tag}: {lab}")
    print(f"   norm: t0={N.iloc[0]:.4f} tF={N.iloc[-1]:.4f} min={N.min():.4f}")
    print(f"   E_reported: t0={E0:.1f} tF={E[-1]:.1f} eV  (rise={E[-1]-E0:+.1f})")
    print(f"   E_ext=E*norm: t0={E0*N.iloc[0]:.1f} tF={E[-1]*N.iloc[-1]:.2f} eV")
print("\nHYPOTHESIS: 3a rises + norm->0 + E_ext->0 ; 3b norm~1 + E_reported flat (NO rise).")
PYV

for tag in exp3a_mask_etrs exp3b_mask_cn; do
  $PY "$JCAMP/analyse.py" results/$tag --label wp --dt 0.01 --slab-face 1000000 --cap-inner 7.5 --cap-lines 7.5,22.5 --title "$tag" --per-frame-norm-wp >>"$LOG" 2>&1 || say "  analyse $tag warn"
  $PY energy_diagnostics.py results/$tag --append >>"$LOG" 2>&1 || say "  diag $tag warn"
done

say "=== decisive experiment COMPLETE ==="
notify "decisive mask experiment COMPLETE" "verdict in $(pwd)/$LOG. Notebooks: results/{exp3a_mask_etrs,exp3b_mask_cn}/report/run_report.ipynb"
say "done."
