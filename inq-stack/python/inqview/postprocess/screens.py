"""Phase: ``screens`` — total / instantaneous / time-windowed LEED + checks.

Outputs under ``results/analysis/screens/``:

* ``total/`` — one PNG (linear and log) per screen, plus a 4×5 grid PNG.
* ``instantaneous/`` — one GIF per screen (fixed colour scale through time).
* ``time_windowed/`` — flat: ``screen_NN_tAAAAAA_to_tBBBBBB.png`` (and ``_log.png``).
* ``coordinate_checks/`` — ``screen_<NN>_raw_index_plot.png`` +
  ``..._coordinate_mapped_plot.png`` for the first screen, to detect the
  ``four-corner-split`` failure mode discussed in the spec §17.6.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from . import _common
from . import pipeline as _pipeline

_INSTANT_RE = re.compile(r"^screen_(\d{2})_t(\d{6})\.dat$")
_WINDOWED_RE = re.compile(
    r"^screen_(\d{2})_t(\d{6})_to_t(\d{6})_(forward|back|paper)\.dat$"
)


def _load_pattern(p: Path):
    from .. import load_leed_pattern
    return load_leed_pattern(p)


def _save_total_panel(pattern, out_png: Path, run_name: str, *,
                      log_scale: bool):
    from .. import plot_leed_pattern
    plot_leed_pattern(pattern, out_png, log_scale=log_scale)


def run(results_dir: Path, *, run_name: str, rebuild: bool, **_) -> dict:
    raw = results_dir / "raw" / "screens"
    out_dir = _common.ensure_dir(results_dir / "analysis" / "screens")
    notes: dict = {"out_dir": str(out_dir)}

    try:
        import imageio.v2 as imageio
        import matplotlib.pyplot as plt
    except ImportError as exc:
        _pipeline.skip(f"missing imageio / matplotlib: {exc}")

    # ---- total/ ----------------------------------------------------------
    total_dir = raw / "total"
    if total_dir.exists():
        out_total = _common.ensure_dir(out_dir / "total")
        files = sorted(total_dir.glob("screen_*.dat"))
        patterns = [_load_pattern(f) for f in files]
        for f, pat in zip(files, patterns):
            base = out_total / f.stem
            png = Path(str(base) + ".png")
            png_log = Path(str(base) + "_log.png")
            if _common.need_rebuild(png, rebuild):
                _save_total_panel(pat, png, run_name, log_scale=False)
            if _common.need_rebuild(png_log, rebuild):
                _save_total_panel(pat, png_log, run_name, log_scale=True)
        # Grid panel
        if patterns:
            cols, rows = 5, 4
            fig, axes = plt.subplots(rows, cols, figsize=(3 * cols, 3 * rows),
                                     dpi=120)
            axes = axes.ravel()
            for i, pat in enumerate(patterns):
                ax = axes[i]
                im = ax.imshow(pat.data, origin="lower", cmap="viridis",
                               extent=pat.extent_bohr, aspect="equal",
                               vmin=0)
                ax.set_title(f"{pat.label}\nz={_common.sigfigs(pat.z_bohr)}",
                             fontsize="x-small")
                ax.set_xticks([]); ax.set_yticks([])
            for j in range(len(patterns), rows * cols):
                axes[j].axis("off")
            fig.suptitle(_common.title(run_name, "total LEED grid"),
                         fontsize="medium")
            fig.tight_layout()
            out = out_total / "all_screens_grid.png"
            if _common.need_rebuild(out, rebuild):
                fig.savefig(out)
            plt.close(fig)
        notes["total"] = str(_common.ensure_dir(out_total))

    # ---- instantaneous/ — one GIF per screen ------------------------------
    inst_dir = raw / "instantaneous"
    if inst_dir.exists():
        out_inst = _common.ensure_dir(out_dir / "instantaneous")
        # Group files by screen index
        groups: dict[str, list[tuple[int, Path]]] = {}
        for f in inst_dir.glob("screen_*_t*.dat"):
            m = _INSTANT_RE.match(f.name)
            if not m:
                continue
            groups.setdefault(m.group(1), []).append((int(m.group(2)), f))
        for sid, items in groups.items():
            items.sort(key=lambda kv: kv[0])
            patterns = [_load_pattern(f) for _, f in items]
            if not patterns:
                continue
            data = np.stack([p.data for p in patterns], axis=0)
            vmin_lin = float(data.min())
            vmax_lin = float(np.percentile(data, 99))
            if vmax_lin <= vmin_lin:
                vmax_lin = vmin_lin + 1.0
            vmin_log = 0.0
            vmax_log = float(np.log1p(vmax_lin))

            for scale_label, transform, vmin, vmax, cbar in [
                ("",     lambda a: a,  vmin_lin, vmax_lin, "intensity (a.u.)"),
                ("_log", np.log1p,     vmin_log, vmax_log, r"log$_{10}$(1 + intensity)"),
            ]:
                stem = out_inst / f"screen_{sid}_time_evolution{scale_label}"
                gif_path = stem.with_suffix(".gif")
                if not _common.need_rebuild(gif_path, rebuild):
                    continue
                tmp = _common.ensure_dir(out_inst / f".__tmp_{sid}{scale_label}")
                pngs: list[Path] = []
                for (step, _path), pat in zip(items, patterns):
                    fig, ax = plt.subplots(figsize=(5, 5), dpi=120)
                    im = ax.imshow(transform(pat.data), origin="lower",
                                   cmap="viridis",
                                   extent=pat.extent_bohr, aspect="equal",
                                   vmin=vmin, vmax=vmax)
                    plt.colorbar(im, ax=ax, label=cbar)
                    ax.set_xlabel("x (bohr)"); ax.set_ylabel("y (bohr)")
                    ax.set_title(_common.title(
                        run_name, f"instantaneous screen_{sid}{scale_label}",
                        step=step, total_steps=items[-1][0],
                        time_au=pat.total_time_au,
                    ))
                    fig.tight_layout()
                    p = tmp / f"f_{step:06d}.png"
                    fig.savefig(p)
                    plt.close(fig)
                    pngs.append(p)
                _common.write_animation(stem, pngs, fps=6)
                for p in pngs:
                    p.unlink(missing_ok=True)
                tmp.rmdir()
        notes["instantaneous"] = str(out_inst)

    # ---- time_windowed/ — flat PNGs per window-screen --------------------
    win_dir = raw / "time_windowed"
    if win_dir.exists():
        out_win = _common.ensure_dir(out_dir / "time_windowed")
        for f in sorted(win_dir.glob("screen_*.dat")):
            m = _WINDOWED_RE.match(f.name)
            if not m:
                continue
            stem = f.stem
            pat = _load_pattern(f)
            for log, suf in [(False, ""), (True, "_log")]:
                out = out_win / f"{stem}{suf}.png"
                if _common.need_rebuild(out, rebuild):
                    _save_total_panel(pat, out, run_name, log_scale=log)
        notes["time_windowed"] = str(out_win)

    # ---- coordinate_checks/ — first total screen as raw vs mapped --------
    if total_dir.exists():
        out_cc = _common.ensure_dir(out_dir / "coordinate_checks")
        first = next(iter(sorted(total_dir.glob("screen_*.dat"))), None)
        if first is not None:
            pat = _load_pattern(first)  # already fftshifted by load_leed_pattern
            # Raw-index plot: undo the fftshift so this image shows the
            # FFT-natural layout the C++ writer produces (peak at corner with
            # four-fold split — the failure mode the spec §17.6 warns about).
            raw = np.fft.ifftshift(pat.data)
            fig, ax = plt.subplots(figsize=(5, 5), dpi=120)
            ax.imshow(raw, origin="lower", cmap="viridis", aspect="equal")
            ax.set_xlabel("array index ix (FFT-natural; index 0 = cell centre)")
            ax.set_ylabel("array index iy")
            ax.set_title(_common.title(
                run_name, f"{pat.label} raw-index plot (no fftshift)"))
            fig.tight_layout()
            out = out_cc / f"{pat.label}_raw_index_plot.png"
            if _common.need_rebuild(out, rebuild):
                fig.savefig(out)
            plt.close(fig)
            # Coordinate-mapped plot (uses fftshifted data + extent in bohr)
            _save_total_panel(pat, out_cc / f"{pat.label}_coordinate_mapped_plot.png",
                              run_name, log_scale=False)
        notes["coordinate_checks"] = str(out_cc)

    return notes
