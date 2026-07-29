#!/bin/bash
# Grid-spacing x r_cut sweep of the projectile pseudopotential representation error.
# Runs the already-built eval_projpot binary directly (no rebuild) with no GS load
# (ideal/impl/gap depend only on the background + pseudopotential, not the electrons).
set -u
cd "$(dirname "$0")"
BIN=$(find build -maxdepth 3 -type f -name run -perm -u+x 2>/dev/null | head -1)
[ -z "$BIN" ] && BIN=$(find . -maxdepth 3 -type f -name run -perm -u+x 2>/dev/null | grep -v run.cpp | head -1)
echo "binary: $BIN"
export INQ_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share
export PSEUDOPOD_SHARE_PATH=/local/data/public/skcb2/tddft/inq/install/share/pseudopod
UPF50=/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/pseudopotentials/electron_gaussian_wpsigma0p5.upf
UPF120=/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/pseudopotentials/electron_gaussian_wpsigma0p5_rc120.upf
export LJ_LAUNCH_Z=-24.5 LJ_GS_DIR=""
for dx in 0.5 0.4 0.3 0.25; do
  for rc in 50 120; do
    upf=$([ "$rc" = 50 ] && echo "$UPF50" || echo "$UPF120")
    echo "########## spacing=$dx  r_cut=$rc ##########"
    LJ_SPACING=$dx LJ_PROJ_UPF="$upf" "$BIN" 2>&1 | grep -E "ideal|impl|gap|n_proj_norm"
  done
done
echo "SWEEP_DONE"
