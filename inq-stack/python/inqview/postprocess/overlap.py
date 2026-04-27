"""Phase: ``overlap`` — animated WP-overlap-with-GS-orbitals bar chart.

Reads:
* ``results/raw/overlap/index.csv`` — step,time_au,file
* ``results/raw/overlap/overlap_NNNNNN.csv`` — single row of n_ref values
  written by ``OrbitalOverlapMatrix::snapshot_wp_only`` (see
  ``inq-stack/include/inqkit/observables/orbital_overlap.hpp:90``).

Produces:
* ``results/analysis/overlap/wp_overlap_with_gs_orbitals.gif`` — animated
  bar chart over GS orbital index, fixed y-axis (0 .. max), per
  ``docs/visualisation-instructions-v1.md`` §5.

This is the only overlap visualisation; the full O_ij matrix is no longer
emitted by the C++ runs (see plan §4.6 item 9).
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from . import _common
from . import pipeline as _pipeline


def _read_overlap_csv(path: Path) -> np.ndarray:
    """Read a wp_only-format overlap CSV: header line then one comma-separated row."""
    with path.open() as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            return np.array([float(x) for x in line.strip().split(",")])
    raise ValueError(f"overlap CSV {path} has no data row")


def run(results_dir: Path, *, run_name: str, rebuild: bool, **_) -> dict:
    raw = results_dir / "raw" / "overlap"
    index_csv = raw / "index.csv"
    if not index_csv.exists():
        _pipeline.skip(f"overlap index missing at {index_csv}")

    out_dir = _common.ensure_dir(results_dir / "analysis" / "overlap")
    out_gif = out_dir / "wp_overlap_with_gs_orbitals.gif"
    if not _common.need_rebuild(out_gif, rebuild):
        return {"gif": str(out_gif), "cached": True}

    try:
        import imageio.v2 as imageio
        import matplotlib.pyplot as plt
    except ImportError as exc:
        _pipeline.skip(f"missing imageio / matplotlib: {exc}")

    # Read the index file
    rows: list[tuple[int, float, Path]] = []
    with index_csv.open() as f:
        reader = csv.DictReader(f)
        for r in reader:
            step = int(r["step"])
            t_au = float(r["time_au"])
            f_csv = raw / r["file"]
            rows.append((step, t_au, f_csv))
    if not rows:
        _pipeline.skip(f"overlap index empty: {index_csv}")

    overlaps = [_read_overlap_csv(p) for _, _, p in rows]
    n_ref = overlaps[0].size
    arr = np.stack([o for o in overlaps], axis=0)  # (n_steps, n_ref)
    y_max = max(1.0, float(arr.max()) * 1.05)

    tmp = _common.ensure_dir(out_dir / ".__tmp_wp_overlap")
    pngs: list[Path] = []
    last_step = rows[-1][0]
    indices = np.arange(n_ref)
    for (step, t_au, _p), row in zip(rows, overlaps):
        fig, ax = plt.subplots(figsize=(7, 4), dpi=120)
        ax.bar(indices, row, color="steelblue")
        ax.set_xlim(-0.5, n_ref - 0.5)
        ax.set_ylim(0, y_max)
        ax.set_xlabel("Ground-state KS orbital index i")
        ax.set_ylabel(r"|⟨ψ$_i^{GS}$ | ψ$_{wp}(t)$⟩|²")
        ax.set_title(_common.title(
            run_name, "WP overlap with GS KS orbitals",
            step=step, total_steps=last_step, time_au=t_au))
        fig.tight_layout()
        p = tmp / f"f_{step:06d}.png"
        fig.savefig(p)
        plt.close(fig)
        pngs.append(p)

    with imageio.get_writer(out_gif, mode="I", fps=8, loop=0) as wr:
        for p in pngs:
            wr.append_data(imageio.imread(p))
    for p in pngs:
        p.unlink(missing_ok=True)
    tmp.rmdir()

    return {"gif": str(out_gif), "n_frames": len(rows), "n_ref": int(n_ref)}
