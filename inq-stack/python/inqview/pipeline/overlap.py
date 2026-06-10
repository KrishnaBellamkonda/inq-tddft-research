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
    out_stem_lin = out_dir / "wp_overlap_with_gs_orbitals"
    out_stem_log = out_dir / "wp_overlap_with_gs_orbitals_log"
    if (not _common.need_rebuild(out_stem_lin.with_suffix(".gif"), rebuild)
            and not _common.need_rebuild(out_stem_log.with_suffix(".gif"), rebuild)):
        return {"gif": str(out_stem_lin.with_suffix(".gif")), "cached": True}

    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        _pipeline.skip(f"missing matplotlib: {exc}")

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
    arr = np.stack([o for o in overlaps], axis=0)   # (n_steps, n_ref)
    # Data-driven y range. Drop the legacy max(1.0, ...) floor — overlap
    # values are typically ≪ 1 for a forward-scattered WP, and clamping
    # to [0,1] hides the structure (TODO 1f).
    y_max_lin = float(arr.max()) * 1.10
    y_min_lin = -0.05 * y_max_lin if y_max_lin > 0 else -0.01
    # Log axis floor: small positive value so symlog shows tiny overlaps.
    if y_max_lin > 0:
        nonzero = arr[arr > 0]
        y_log_floor = (float(nonzero.min()) if nonzero.size else 1e-12) * 0.5
    else:
        y_log_floor = 1e-12
    y_log_top = max(y_max_lin, 10 * y_log_floor)

    last_step = rows[-1][0]
    indices = np.arange(n_ref)

    def _render(stem: Path, log_scale: bool) -> dict:
        if not _common.need_rebuild(stem.with_suffix(".gif"), rebuild):
            return {"gif": str(stem.with_suffix(".gif")), "cached": True}
        tmp = _common.ensure_dir(out_dir / f".__tmp_{stem.name}")
        pngs: list[Path] = []
        for (step, t_au, _p), row in zip(rows, overlaps):
            fig, ax = plt.subplots(figsize=(7, 4), dpi=120)
            ax.bar(indices, row, color="steelblue")
            ax.set_xlim(-0.5, n_ref - 0.5)
            if log_scale:
                ax.set_yscale("log")
                ax.set_ylim(y_log_floor, y_log_top)
            else:
                ax.set_ylim(y_min_lin, y_max_lin)
            ax.set_xlabel("Ground-state KS orbital index i")
            ax.set_ylabel(r"|⟨ψ$_i^{GS}$ | ψ$_{wp}(t)$⟩|²"
                          + (" (log)" if log_scale else ""))
            ax.set_title(_common.title(
                run_name,
                "WP overlap with GS KS orbitals" + (" (log)" if log_scale else ""),
                step=step, total_steps=last_step, time_au=t_au))
            fig.tight_layout()
            p = tmp / f"f_{step:06d}.png"
            fig.savefig(p)
            plt.close(fig)
            pngs.append(p)
        outs = _common.write_animation(stem, pngs, fps=8)
        for p in pngs:
            p.unlink(missing_ok=True)
        tmp.rmdir()
        return {"gif": str(outs["gif"]),
                "mp4": str(outs["mp4"]) if outs["mp4"] else None}

    return {
        "linear":    _render(out_stem_lin, log_scale=False),
        "log":       _render(out_stem_log, log_scale=True),
        "n_frames":  len(rows),
        "n_ref":     int(n_ref),
        "y_max_lin": float(y_max_lin),
    }
