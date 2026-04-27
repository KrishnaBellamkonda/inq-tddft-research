"""Cross-run comparison ("hypothesis" mode).

Used by ``coronene_postprocess.py hypothesis ...`` to collate the
``results/`` trees of several runs into a single ``hypotheses/<NN>_*/``
folder of comparison artefacts.

The function takes:

  hypothesis_dir : path to ``hypotheses/<NN>_*/``
  runs : list of (label, run_dir) tuples — ``run_dir`` is a path containing
         a populated ``results/`` subtree. ``label`` is what appears in
         legends.

Outputs (best-effort; missing inputs skipped silently):

  * ``leed_total_grid_<labels>.png`` — side-by-side total LEED at each
    screen; one panel row per run.
  * ``peak_intensity_vs_label.png`` — peak total LEED intensity at the
    central screen (label index halfway through the screen list) vs run
    label; useful for E and b scans.
  * ``energy_spectrum_overlay.png`` — FFT(total_energy) curves overlaid.
  * ``current_z_overlay.png`` — J_z(t) curves overlaid.
  * ``observable_drift_summary.csv`` — per-run total-energy drift,
    norm-conservation, etc.
"""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import numpy as np

from . import _common


def run_hypothesis(
    hypothesis_dir: Path,
    runs: Sequence[tuple[str, Path]],
    *,
    rebuild: bool = False,
) -> dict:
    out = _common.ensure_dir(hypothesis_dir)
    notes: dict = {"hypothesis_dir": str(out), "runs": []}

    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        return {"error": f"matplotlib missing: {exc}"}

    # ---- 1. Total LEED grid (one row per run) ----------------------------
    rows = []
    for label, run_dir in runs:
        total_dir = Path(run_dir) / "results" / "raw" / "screens" / "total"
        files = sorted(total_dir.glob("screen_*.dat")) if total_dir.exists() else []
        patterns = []
        if files:
            from .. import load_leed_pattern
            patterns = [load_leed_pattern(f) for f in files]
        rows.append((label, patterns))
        notes["runs"].append({"label": label, "n_screens": len(patterns)})

    if any(p for _, p in rows):
        max_screens = max(len(p) for _, p in rows)
        fig, axes = plt.subplots(
            len(rows), max_screens,
            figsize=(2.2 * max_screens, 2.2 * len(rows)),
            dpi=120,
            squeeze=False,
        )
        for i, (label, patterns) in enumerate(rows):
            for j in range(max_screens):
                ax = axes[i][j]
                if j < len(patterns):
                    pat = patterns[j]
                    ax.imshow(pat.data, origin="lower", cmap="viridis",
                              extent=pat.extent_bohr, aspect="equal", vmin=0)
                    if i == 0:
                        ax.set_title(pat.label, fontsize="x-small")
                else:
                    ax.axis("off")
                if j == 0:
                    ax.set_ylabel(label, fontsize="x-small")
                ax.set_xticks([]); ax.set_yticks([])
        fig.suptitle("Total LEED comparison", fontsize="medium")
        fig.tight_layout()
        out_png = out / "leed_total_grid.png"
        if _common.need_rebuild(out_png, rebuild):
            fig.savefig(out_png)
        plt.close(fig)
        notes["leed_grid"] = str(out_png)

    # ---- 2. Peak central-screen intensity vs label -----------------------
    peaks: list[tuple[str, float]] = []
    for label, patterns in rows:
        if not patterns:
            continue
        mid = patterns[len(patterns) // 2]
        peaks.append((label, float(mid.data.max())))
    if peaks:
        fig, ax = plt.subplots(figsize=(6, 4), dpi=120)
        labels, vals = zip(*peaks)
        ax.bar(labels, vals, color="steelblue")
        ax.set_ylabel("peak ∫ρ·dt at central screen (bohr⁻³·a.u.)")
        ax.set_title("Peak LEED intensity by run")
        fig.tight_layout()
        out_png = out / "peak_intensity_vs_label.png"
        if _common.need_rebuild(out_png, rebuild):
            fig.savefig(out_png)
        plt.close(fig)
        notes["peak_intensity"] = str(out_png)

    # ---- 3. Energy/current overlays --------------------------------------
    from .. import load_observables
    fig_e, ax_e = plt.subplots(figsize=(7, 4), dpi=120)
    fig_j, ax_j = plt.subplots(figsize=(7, 4), dpi=120)
    plotted = 0
    for label, run_dir in runs:
        csv = Path(run_dir) / "results" / "raw" / "observables" / "observables.csv"
        if not csv.exists():
            continue
        df = load_observables(csv)
        if "energy_total" in df.columns:
            ax_e.plot(df["time_au"], df["energy_total"] - df["energy_total"].iloc[0],
                      label=label, linewidth=1.0)
        if "current_z" in df.columns:
            ax_j.plot(df["time_au"], df["current_z"], label=label, linewidth=1.0)
        plotted += 1

    if plotted:
        ax_e.set_xlabel("Time (a.u.)")
        ax_e.set_ylabel("E_total - E_total(0)  (Ha)")
        ax_e.set_title("Total energy drift")
        ax_e.legend(fontsize="x-small")
        fig_e.tight_layout()
        out_png = out / "energy_drift_overlay.png"
        if _common.need_rebuild(out_png, rebuild):
            fig_e.savefig(out_png)

        ax_j.set_xlabel("Time (a.u.)")
        ax_j.set_ylabel("J_z (a.u.)")
        ax_j.set_title("z-current overlay")
        ax_j.legend(fontsize="x-small")
        fig_j.tight_layout()
        out_png = out / "current_z_overlay.png"
        if _common.need_rebuild(out_png, rebuild):
            fig_j.savefig(out_png)

    import matplotlib.pyplot as _plt
    _plt.close("all")

    return notes
