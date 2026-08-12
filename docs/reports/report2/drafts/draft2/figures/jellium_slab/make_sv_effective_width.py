"""Draft-2 S(v) relabelled by effective width (fig 15).

Wrapper that calls build_sv_effective_width_s6.draw() with the output
redirected to this draft-2 figures directory.
  Output: slab_sv_effective_width_s56.png (600 DPI, bbox_inches=None)

The colour scheme is per-sigma (not WP/classical binary) — each sigma_WP
gets one colour; WP vs classical is encoded by marker fill + linestyle.
This is correct for sigma-series figures and differs from the WP/classical
binary convention intentionally.

Run:
  /local/data/public/skcb2/tddft/venv/bin/python3 make_sv_effective_width.py
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[7]
HERE = Path(__file__).parent
SRC  = REPO / "ResearchProject/systems/localised_jellium/hypotheses/sigma56_sv"

sys.path.insert(0, str(SRC))
sys.path.insert(0, str(REPO / "docs/reports/report2/drafts/draft1/figures"))
sys.path.insert(0, str(REPO / "inq-stack/python"))

import matplotlib
matplotlib.use("Agg")

import build_sv_effective_width_s6 as B

# Override the report output path to draft2
B.REPORT_FIG = HERE / "slab_sv_effective_width_s56.png"


def main() -> None:
    snap = B.snapshot()
    out  = HERE / "_sv_eff_width_tmp.png"
    B.draw(False, out, snap)
    out.unlink(missing_ok=True)
    if B.REPORT_FIG.exists():
        print(f"Draft-2 figure (600 DPI): {B.REPORT_FIG}")
    else:
        print(f"WARNING: {B.REPORT_FIG} not written — check draw() logic")


if __name__ == "__main__":
    main()
