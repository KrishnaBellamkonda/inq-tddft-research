#!/usr/bin/env bash
# LOCKED CORRECTNESS GATE — cap_fix campaign.
# The primary metric (artifact_rise_eV -> 0) can be gamed trivially by a CAP
# that absorbs nothing (eta -> 0 gives rise = 0 AND no absorption). This gate
# enforces "the absorber still absorbs":
#   absorbed_e >= 0.5   (at least half an electron removed by t = 28)
#   fallback if charge.csv missing: drain_eV <= -10
# usage: ./checks.sh <run_name>
set -euo pipefail
RUN_NAME="${1:?usage: checks.sh <run_name>}"
cd "$(dirname "$0")"

OUT=$(/local/data/public/skcb2/tddft/venv/bin/python3 run_metrics.py "results/${RUN_NAME}")
ABSORBED=$(echo "$OUT" | awk -F= '/^METRIC absorbed_e=/{print $2}')
DRAIN=$(echo "$OUT" | awk -F= '/^METRIC drain_eV=/{print $2}')

/local/data/public/skcb2/tddft/venv/bin/python3 - "$ABSORBED" "$DRAIN" <<'EOF'
import math, sys
absorbed, drain = float(sys.argv[1]), float(sys.argv[2])
if not math.isnan(absorbed):
    ok = absorbed >= 0.5
    print(f"checks: absorbed_e={absorbed:.3f} (>=0.5 required) -> {'PASS' if ok else 'FAIL'}")
else:
    ok = drain <= -10.0
    print(f"checks: absorbed_e unavailable; drain_eV={drain:.2f} (<=-10 required) -> {'PASS' if ok else 'FAIL'}")
sys.exit(0 if ok else 1)
EOF
