"""Draft-2 bath norm lost vs sigma_WP (fig 24).

Wrapper calling make_norm_loss.draw() with output redirected to draft-2.
  Output: slab_norm_loss_vs_sigma.png (600 DPI, bbox_inches=None via draw())

Colour scheme: per-velocity (Okabe-Ito), WP solid / classical dashed.

Run:
  /local/data/public/skcb2/tddft/venv/bin/python3 make_norm_loss.py
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

import make_norm_loss as N

# Override REPORT_FIG to draft-2 path
N.REPORT_FIG = HERE / "slab_norm_loss_vs_sigma.png"


def main() -> None:
    df = N.collect()
    df = df.sort_values(["half", "sigma_WP", "v"]).reset_index(drop=True)
    df.to_csv(HERE / "norm_loss_table.csv", index=False)
    print(df.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    N.draw(df, HERE / "norm_loss_vs_sigma_local.png")
    if N.REPORT_FIG.exists():
        print(f"Draft-2 figure: {N.REPORT_FIG}")


if __name__ == "__main__":
    main()
