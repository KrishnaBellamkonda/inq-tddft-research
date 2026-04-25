from __future__ import annotations

"""
Analysis for run_07_paper_replica.

Standard observable suite from docs/observables_reference.md, plus the
overlap-matrix GIFs filtered to the *physically meaningful* submatrix
(rows = occupied GS orbitals, columns = occupied + WP). The C++
OrbitalOverlapMatrix records a (wp_idx) x (wp_idx + 1) matrix where rows
include unoccupied extra-state ground orbitals; the GIFs here only show
rows i in [0, n_occupied) and the submatrix columns
{0, ..., n_occupied-1, wp_idx} (occupied block + WP slot).

Outputs under results/:
  vti/density_gs/                 (corrected via post-fix writer)
  vti/density_gs_orbitals/...
  vti/density_rt_target/
  vti/density_rt_wp/
  visualisation/observables.png
  visualisation/spectrum_*.png
  visualisation/leed/screen_NN.png            (full-run avg)
  visualisation/leed_window/screen_NN.png     (paper t1..t2 window)
  visualisation/leed_evolution.gif            (screen 10, snapshot series)
  visualisation/density_*.gif
  visualisation/overlap/overlap_orbital_jJJJJ.gif   (filtered to occupied+WP)
  diagnostics.txt
"""

import csv
import re
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

RUN_DIR = Path(__file__).resolve().parent
RESULTS = RUN_DIR / "results"
VTI_ROOT = RESULTS / "vti"
VIS = RESULTS / "visualisation"

REPO_ROOT = RUN_DIR.parents[3] if len(RUN_DIR.parents) >= 4 else RUN_DIR
sys.path.insert(0, str(REPO_ROOT / "inq-stack" / "python"))

import imageio.v3 as iio
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from inqview import convert_real_series_to_vti
from inqview.data import load_real_field
from inqview.overlap import (
    iter_overlap_series,
    plot_overlap_column_gif,
    pick_meaningful_columns,
)


GIF_FPS = 20
OVERLAP_GIF_FPS = 12
META_SUFFIX = ".meta.txt"


def _natural_key(p: Path):
    parts = re.split(r"(\d+)", p.name)
    return [int(x) if x.isdigit() else x for x in parts]


def _read_summary(path: Path) -> dict:
    out = {}
    if not path.exists():
        return out
    for line in path.read_text().splitlines():
        if "=" not in line:
            continue
        k, v = line.split("=", 1)
        out[k.strip()] = v.strip()
    return out


def _convert_all_vti():
    if not RESULTS.exists():
        return
    for src_name in [
        "density_gs",
        "density_rt_target",
        "density_rt_wp",
    ]:
        src = RESULTS / src_name
        if not src.exists():
            continue
        out_dir = VTI_ROOT / src_name
        out_dir.mkdir(parents=True, exist_ok=True)
        metas = sorted(
            (p for p in src.glob(f"*{META_SUFFIX}")
             if (p.parent / (p.name[: -len(META_SUFFIX)] + ".raw")).exists()),
            key=_natural_key,
        )
        if not metas:
            continue
        convert_real_series_to_vti(metas, out_dir, array_name="density")
        print(f"  VTI: {src_name} -> {len(metas)} file(s)")

    gs_orb = RESULTS / "density_gs_orbitals"
    if gs_orb.exists():
        for orb_dir in sorted(gs_orb.iterdir()):
            if not orb_dir.is_dir():
                continue
            metas = sorted(
                (p for p in orb_dir.glob(f"*{META_SUFFIX}")
                 if (p.parent / (p.name[: -len(META_SUFFIX)] + ".raw")).exists()),
                key=_natural_key,
            )
            if not metas:
                continue
            out_dir = VTI_ROOT / "density_gs_orbitals" / orb_dir.name
            out_dir.mkdir(parents=True, exist_ok=True)
            convert_real_series_to_vti(metas, out_dir, array_name="density")
        print(f"  VTI: density_gs_orbitals converted")


def _plot_observables():
    obs_csv = RESULTS / "observables.csv"
    if not obs_csv.exists():
        print("  Skip observables: observables.csv not found")
        return
    rows = []
    with obs_csv.open() as fh:
        r = csv.DictReader(fh)
        cols = r.fieldnames or []
        for row in r:
            rows.append(row)
    if not rows:
        return
    t = np.array([float(row["time_au"]) for row in rows])
    VIS.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(3, 1, figsize=(8, 8), dpi=120, sharex=True)
    e_total = np.array([float(row["energy_total"]) for row in rows]) if "energy_total" in cols else None
    if e_total is not None:
        e0 = e_total[0]
        ax[0].plot(t, (e_total - e0), color="black")
        ax[0].set_ylabel(r"$E - E_0$ (Ha)")
        ax[0].set_title("Energy drift")

    for c, key in zip(["C0", "C1", "C2"], ["current_x", "current_y", "current_z"]):
        if key in cols:
            ax[1].plot(t, [float(row[key]) for row in rows], color=c, label=key)
    ax[1].set_ylabel("Current (a.u.)")
    ax[1].legend(loc="best", fontsize=8)

    for c, key in zip(["C0", "C1", "C2"], ["dipole_x", "dipole_y", "dipole_z"]):
        if key in cols:
            ax[2].plot(t, [float(row[key]) for row in rows], color=c, label=key)
    ax[2].set_xlabel("time (a.u.)")
    ax[2].set_ylabel("Dipole (a.u.)")
    ax[2].legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(VIS / "observables.png")
    plt.close(fig)
    print(f"  Wrote {VIS.relative_to(RUN_DIR)}/observables.png")

    # Dipole power spectrum (paper Eq. 9 analogue)
    if "dipole_x" in cols:
        try:
            dx = np.array([float(row["dipole_x"]) for row in rows])
            dx -= dx.mean()
            n = len(dx)
            dt_au = float(t[1] - t[0])
            freqs_au = np.fft.rfftfreq(n, d=dt_au)  # angular freq
            mag = np.abs(np.fft.rfft(dx))
            mag2 = mag ** 2
            # convert frequency to eV
            HA_TO_EV = 27.21138625
            ev = freqs_au * 2 * np.pi * HA_TO_EV  # omega in eV
            fig, axs = plt.subplots(figsize=(8, 4), dpi=120)
            axs.semilogy(ev, mag2, color="black", lw=0.8)
            axs.set_xlim(0, 30)
            axs.set_xlabel(r"$\omega$ (eV)")
            axs.set_ylabel(r"$|p_x(\omega)|^2$")
            axs.set_title("Dipole excitation spectrum (x)")
            fig.tight_layout()
            fig.savefig(VIS / "spectrum_dipole_x.png")
            plt.close(fig)
            print(f"  Wrote spectrum_dipole_x.png")
        except Exception as e:
            print(f"  WARN: spectrum failed: {e}")


def _plot_leed_panels():
    for src, label in [
        ("screens",             "leed"),
        ("screens_leed_window", "leed_window"),
    ]:
        src_dir = RESULTS / src
        if not src_dir.exists():
            continue
        out_dir = VIS / label
        out_dir.mkdir(parents=True, exist_ok=True)
        for dat in sorted(src_dir.glob("screen_*.dat"), key=_natural_key):
            try:
                arr = np.loadtxt(dat)
            except Exception:
                continue
            fig, ax = plt.subplots(figsize=(5, 5), dpi=120)
            im = ax.imshow(arr, cmap="inferno", origin="lower")
            ax.set_title(f"{label} {dat.stem}")
            fig.colorbar(im, ax=ax, fraction=0.046)
            fig.tight_layout()
            fig.savefig(out_dir / f"{dat.stem}.png")
            plt.close(fig)
        print(f"  Wrote {label} panels: {out_dir.relative_to(RUN_DIR)}")


def _plot_density_gifs():
    VIS.mkdir(parents=True, exist_ok=True)
    for series_name in ["density_rt_target", "density_rt_wp"]:
        src = RESULTS / series_name
        if not src.exists():
            continue
        metas = sorted(
            (p for p in src.glob(f"*{META_SUFFIX}")
             if (p.parent / (p.name[: -len(META_SUFFIX)] + ".raw")).exists()),
            key=_natural_key,
        )
        if not metas:
            continue
        for axis_name, ax_idx in [("xz", 1), ("yz", 0), ("xy", 2)]:
            frames = []
            for m in metas:
                f = load_real_field(meta_path=m)
                arr = np.asarray(f.array)
                # pick mid-slice along ax_idx
                mid = arr.shape[ax_idx] // 2
                if ax_idx == 0:
                    sl = arr[mid, :, :]
                elif ax_idx == 1:
                    sl = arr[:, mid, :]
                else:
                    sl = arr[:, :, mid]
                fig, ax = plt.subplots(figsize=(4, 4), dpi=80)
                ax.imshow(sl.T, cmap="inferno", origin="lower")
                ax.set_title(f"{series_name} {axis_name} t={f.meta.time_au or 0:.2f}")
                fig.tight_layout()
                fig.canvas.draw()
                frames.append(np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy())
                plt.close(fig)
            iio.imwrite(VIS / f"{series_name}_{axis_name}.gif",
                        frames, fps=GIF_FPS, loop=0)
        print(f"  density GIFs for {series_name}")


def _plot_overlap_gifs():
    summary = _read_summary(RESULTS / "run_summary.txt")
    n_occupied = int(summary.get("n_occupied", 0))
    wp_idx = int(summary.get("wp_state_index", -1))
    if n_occupied <= 0 or wp_idx < 0:
        print("  Skip overlap: cannot read n_occupied/wp_state_index from run_summary.txt")
        return

    overlap_dir = RESULTS / "overlap"
    if not (overlap_dir / "index.csv").exists():
        print("  Skip overlap: results/overlap/index.csv not found")
        return

    series = list(iter_overlap_series(overlap_dir))
    if not series:
        print("  Skip overlap: no matrices loaded")
        return

    n_ref_raw, n_evolved_raw = series[0].matrix.shape
    print(f"  Overlap raw shape = {n_ref_raw} x {n_evolved_raw}, "
          f"n_occupied = {n_occupied}, wp_idx = {wp_idx}")

    # Filter to physically meaningful submatrix:
    # rows = occupied (0..n_occupied-1)
    # columns = occupied (0..n_occupied-1) plus WP (wp_idx, mapped to last column).
    cols_to_plot = pick_meaningful_columns(n_occupied, wp_idx, n_evolved_raw)

    out_dir = VIS / "overlap"
    out_dir.mkdir(parents=True, exist_ok=True)
    for j_global, j_label in cols_to_plot:
        col_data = [s.matrix[:n_occupied, j_global] for s in series]
        times    = [s.time_au for s in series]
        plot_overlap_column_gif(
            out_path=out_dir / f"overlap_orbital_{j_label:04d}.gif",
            n_ref_rows=n_occupied,
            col_data=col_data,
            times=times,
            title_prefix=f"Overlap (j={j_label}{' = WP' if j_global == wp_idx else ''})",
            fps=OVERLAP_GIF_FPS,
        )
    print(f"  Wrote {len(cols_to_plot)} overlap GIFs to {out_dir.relative_to(RUN_DIR)}")


def _diagnostics():
    summary_src = RESULTS / "run_summary.txt"
    out = RESULTS / "diagnostics.txt"
    parts = []
    summary = _read_summary(summary_src)
    if summary:
        parts.append("# run summary")
        for k, v in summary.items():
            parts.append(f"{k} = {v}")

    # N-electron conservation on density_rt_target + density_rt_wp
    target = sorted((RESULTS / "density_rt_target").glob(f"*{META_SUFFIX}"), key=_natural_key)
    wp     = sorted((RESULTS / "density_rt_wp"    ).glob(f"*{META_SUFFIX}"), key=_natural_key)
    if target and wp:
        n_check = min(20, len(target), len(wp))
        idxs = np.linspace(0, len(target) - 1, n_check, dtype=int)
        ne_vals = []
        for i in idxs:
            ft = load_real_field(meta_path=target[min(i, len(target) - 1)])
            fw = load_real_field(meta_path=wp[    min(i, len(wp)     - 1)])
            dx, dy, dz = ft.meta.spacing_bohr
            n_total = float((np.asarray(ft.array) + np.asarray(fw.array)).sum() * dx * dy * dz)
            ne_vals.append(n_total)
        ne = np.array(ne_vals)
        parts.append("\n# N-electron conservation (target+wp integrated)")
        parts.append(f"n_samples = {len(ne)}")
        parts.append(f"mean = {ne.mean():.6f}")
        parts.append(f"std  = {ne.std():.6e}")
        parts.append(f"min  = {ne.min():.6f}")
        parts.append(f"max  = {ne.max():.6f}")
        parts.append(f"drift_pct = {100*ne.std()/max(ne.mean(),1e-30):.4f}")

    out.write_text("\n".join(parts) + "\n")
    print(f"  Wrote {out.relative_to(RUN_DIR)}")


def main() -> int:
    if not RESULTS.exists():
        print(f"No results at {RESULTS}")
        return 1
    print("=== run_07 analysis ===")
    print("[1] VTI conversion"); _convert_all_vti()
    print("[2] Observables plot + spectrum"); _plot_observables()
    print("[3] LEED panels"); _plot_leed_panels()
    print("[4] Density GIFs"); _plot_density_gifs()
    print("[5] Overlap GIFs"); _plot_overlap_gifs()
    print("[6] Diagnostics"); _diagnostics()
    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
