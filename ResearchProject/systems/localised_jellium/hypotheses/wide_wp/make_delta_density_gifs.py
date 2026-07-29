#!/usr/bin/env python3
"""Add delta-density GIFs to the wide-WP run notebook-figs folders.

For each wide-WP run, render two xz-slice density animations that complement the
existing total-density `lead_density.gif`:
  * density_delta.gif                — Δn(t) = n(t) − n(0)          (delta density)
  * density_delta_instantaneous.gif  — Δn(t) = n(t) − n(t−Δt)       (instantaneous)

Rules-compliant (report-figures skill): each GIF has a LINEAR (left) + LOG (right)
panel on ONE fixed shared colour scale (SymLogNorm for the signed field), rendered
through the canonical `_density_views._gif`, whose frames come from
`inqview.load_vti` (physical order — NO fftshift). Clim = symmetric 99.5th
percentile of |Δn| over all frames (same recipe as the Emilio-deck triptych, so
the two stay comparable).

Writes into each run's `hypotheses/wide_wp/<run>_run_notebook_figs/` folder.
Run: /local/data/public/skcb2/tddft/venv/bin/python3 make_delta_density_gifs.py
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np

ROOT = Path("/local/data/public/skcb2/tddft")
LJ = ROOT / "ResearchProject/systems/localised_jellium"
HYP = LJ / "hypotheses"
sys.path.insert(0, str(HYP))                        # _density_views
sys.path.insert(0, str(ROOT / "inq-stack/python"))  # inqview (load_vti)
from _density_views import _load_xz_series, _gif      # noqa: E402

SLAB_HALF = 12.5          # slab half-width (Bohr) -> cyan dotted faces
DT, WEVERY, FPS = 0.04, 6, 12

RUNS = [
    (LJ / "scripts/wide_wp/wp/results/wp_pbc_E300/raw/vti/density_total",
     HYP / "wide_wp/wp_pbc_E300_run_notebook_figs"),
    (LJ / "scripts/wide_wp/wp/results/results/wp_per2_E300_long/raw/vti/density_total",
     HYP / "wide_wp/wp_per2_E300_long_run_notebook_figs"),
]


def main() -> None:
    for total_dir, figs in RUNS:
        figs.mkdir(parents=True, exist_ok=True)
        steps, x, z, tot = _load_xz_series(str(total_dir))
        extent = [x[0], x[-1], z[0], z[-1]]
        views = {
            "density_delta": (
                r"$\Delta n = n(t)-n(0)$",
                [s - tot[0] for s in tot]),
            "density_delta_instantaneous": (
                r"$\Delta n = n(t)-n(t-\Delta t)$",
                [np.zeros_like(tot[0])] +
                [tot[i] - tot[i - 1] for i in range(1, len(tot))]),
        }
        for name, (title, frames) in views.items():
            a = float(np.percentile(np.abs(np.stack(frames)), 99.5)) or 1e-12
            out = figs / f"{name}.gif"
            _gif(frames, extent, "RdBu_r", -a, a, title, steps, DT, WEVERY,
                 str(out), SLAB_HALF, FPS)
            print(f"wrote {out}  ({len(frames)} frames, clim +/-{a:.2e})")


if __name__ == "__main__":
    main()
