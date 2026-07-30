#!/usr/bin/env bash
# ============================================================================
# run_cap_sweeps.sh — Phase 1a (eta sweep) + Phase 2 (partial-absorption ladder)
# on GPU 0, in parallel with the decisive mask experiment on GPU 1.
# Same geometry as production one-sided CAP (sigma0=3, E=400eV, h=0.4, dt=0.01,
# LZ=45, CAP band z in [7.5,22.5], launch -7.5). CAP = perturbations::absorbing (ETRS).
#   Phase 1a: eta in {-0.3,-0.7,-1.0,-2.0,-3.5} at fixed W=15
#             -> residual (E_reported - norm*E0) must NOT scale with eta (else reflection).
#   Phase 2 : eta in {-0.125,-0.217,-0.416} -> norm ends ~0.5,0.3,0.1
#             -> E_ext/E0 must equal norm at every endpoint (continuity of the fix).
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
LOG=cap_sweeps.log; : > "$LOG"
say(){ echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }
notify(){ $PY "$JCAMP/notify.py" "$1" "$2" >>"$LOG" 2>&1 || true; }
COMMON="WP_SIGMA=3 WP_K0=5.421 WP_LZ=45 WP_LPERP=30 WP_H=0.4 WP_DT=0.01 WP_CAP_L=15 WP_LAUNCH_Z=-7.5 WP_NSTEPS=800 WP_WF_EVERY=40 WP_MOM_EVERY=1 WP_ABS=cap"

say "=== Phase 1a eta-sweep + Phase 2 ladder (GPU 0) ==="
notify "CAP sweeps started (GPU0)" "eta sweep + partial-absorption ladder; log $(pwd)/$LOG"

for eta in -0.3 -0.7 -1.0 -2.0 -3.5; do
  tag="exp1a_eta${eta}"
  say "1a eta=$eta ..."
  env $COMMON WP_ETA=$eta WP_OUT=$tag ./run > $tag.run.log 2>&1 || say "  $tag FAILED"
done
for pair in "-0.125:N0p5" "-0.217:N0p3" "-0.416:N0p1"; do
  eta=${pair%%:*}; nm=${pair##*:}
  tag="exp2_${nm}"
  say "2 eta=$eta ($nm) ..."
  env $COMMON WP_ETA=$eta WP_OUT=$tag ./run > $tag.run.log 2>&1 || say "  $tag FAILED"
done

say "=== summary: residual vs eta, and E_ext/E0 vs norm ==="
$PY - <<'PYV' 2>&1 | tee -a "$LOG"
import pandas as pd, glob, re
HA=27.211386
print(f"{'run':16} {'eta':>6} {'norm_f':>8} {'E_rep_f':>9} {'E_ext_f':>9} {'resid=E_rep-normE0':>18} {'E_ext/E0':>9} {'~norm?':>7}")
for d in sorted(glob.glob("results/exp1a_eta*")) + sorted(glob.glob("results/exp2_*")):
    try:
        en=pd.read_csv(f"{d}/raw/observables/energies.csv"); ms=pd.read_csv(f"{d}/raw/observables/wp_momentum_stats.csv",comment='#')
        N=(ms.norm_check/ms.norm_check.iloc[0]); Nf=float(N.iloc[-1])
        E=en.total.values*HA; E0=E[0]; Ef=E[-1]; Eext=Ef*Nf
        resid=Ef-Nf*E0
        m=re.search(r'eta(-?\d+\.?\d*)',d); eta=m.group(1) if m else '?'
        print(f"{d.split('/')[-1]:16} {eta:>6} {Nf:8.4f} {Ef:9.1f} {Eext:9.2f} {resid:18.1f} {Eext/E0:9.4f} {'yes' if abs(Eext/E0-Nf)<0.05 else 'NO':>7}")
    except Exception as e: print(d, "err", e)
print("\nPredicted: E_ext/E0 == norm (yes) for ALL; residual ~ E_reported*(1-norm), NOT scaling with eta as reflection.")
PYV
say "=== CAP sweeps COMPLETE ==="
notify "CAP sweeps COMPLETE (GPU0)" "residual-vs-eta + E_ext/E0-vs-norm table in $(pwd)/$LOG"
say "done."
