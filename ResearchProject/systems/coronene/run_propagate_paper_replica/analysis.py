from __future__ import annotations

"""
Analysis for run_propagate_paper_replica.

Inputs (under results/):
  density_rt_target/density_t<step>.vti    (binary VTI, every 10 steps)
  density_rt_wp/density_t<step>.vti        (binary VTI, every 10 steps)
  observables.csv                          (every step)
  overlap/overlap_<step>.csv               (every step, wp-only single row)
  overlap/index.csv
  screens/screen_NN.dat                    (full-run accumulator)
  screens_leed_window/screen_NN.dat        (paper window accumulator)
  screens_snapshots/step_<step>/screen_NN.dat
  run_summary.txt

Outputs (written to results/visualisation/):
  observables.png             energy drift + currents + dipoles
  spectrum_dipole_x.png       FFT of dipole_x
  leed/screen_NN.png          full-run LEED panels
  leed_window/screen_NN.png   paper-window LEED panels
  leed_grid.png               20-panel grid for the paper window
  leed_screen_10_evolution.gif snapshot timeseries for screen 10
  density_rt_target_{xy,xz,yz}.gif
  density_rt_wp_{xy,xz,yz}.gif
  overlap_wp_heatmap.png      |<GS_i|WP(t)>|^2 vs time, all i
  overlap_wp_top.png          top-N dominant orbitals as line plot
  diagnostics.txt             N-electron conservation, energy drift, summary

This run uses the native VTI writer (no .raw/.meta.txt). The reader below
parses our binary VTI files directly via XML + base64; it does not depend
on a system VTK install.
"""

import base64
import csv
import re
import struct
import sys
from pathlib import Path
from xml.etree import ElementTree as ET

import imageio.v3 as iio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

RUN_DIR = Path(__file__).resolve().parent
RESULTS = RUN_DIR / "results"
VIS = RESULTS / "visualisation"

GIF_FPS = 20

HA_TO_EV = 27.21138625


# ─────────────────────────────────────────────────────────────────────────────
# VTI reader (matches inqkit::io::VTIImageDataWriter binary output)
# ─────────────────────────────────────────────────────────────────────────────

def _natural_key(p):
    name = p.name if hasattr(p, "name") else str(p)
    parts = re.split(r"(\d+)", name)
    return [int(x) if x.isdigit() else x for x in parts]


def read_vti(path: Path) -> dict:
    """Parse a .vti file (ASCII or binary) into nx,ny,nz / origin / spacing /
    arrays-by-name. Array shape = (nx, ny, nz) — VTK x-fastest stream is
    transposed back to our internal x-slowest convention."""
    tree = ET.parse(path)
    root = tree.getroot()
    img = root.find("ImageData")
    e = [int(x) for x in img.attrib["WholeExtent"].split()]
    nx = e[1] - e[0] + 1
    ny = e[3] - e[2] + 1
    nz = e[5] - e[4] + 1
    origin = tuple(float(x) for x in img.attrib["Origin"].split())
    spacing = tuple(float(x) for x in img.attrib["Spacing"].split())

    arrays = {}
    for da in img.find("Piece").find("PointData").findall("DataArray"):
        name = da.attrib["Name"]
        fmt = da.attrib["format"]
        text = da.text or ""
        if fmt == "ascii":
            values = np.fromstring(text, sep=" ", dtype=np.float64)
        elif fmt == "binary":
            blob = base64.b64decode("".join(text.split()))
            (n_bytes,) = struct.unpack("<Q", blob[:8])
            values = np.frombuffer(blob[8:8 + n_bytes], dtype="<f8").copy()
        else:
            raise ValueError(f"{path}: unknown format {fmt}")
        # VTK stream is x-fastest: index = ix + nx*(iy + ny*iz). Reshape
        # (nz, ny, nx) and transpose to our (nx, ny, nz) convention.
        arr = values.reshape((nz, ny, nx)).transpose(2, 1, 0).copy()
        arrays[name] = arr
    return {
        "nx": nx, "ny": ny, "nz": nz,
        "origin": origin, "spacing": spacing,
        "arrays": arrays,
    }


def _vti_meta_time(path: Path):
    """Best-effort extract t (au) from filename `density_tNNNNNN.vti`."""
    m = re.search(r"_t(\d+)", path.stem)
    return int(m.group(1)) if m else None


def _load_screen_centred(path: Path) -> np.ndarray:
    """
    Load a screen .dat file and FFT-shift it so the diffraction peak lands
    at the array centre.

    LeedPatternAccumulator (inq-stack/include/inqkit/screens/...) writes the
    screen in INQ's FFT-natural order (array index (0,0) corresponds to
    physical origin (x=0, y=0); positive coordinates first, then wrapped
    negative). Without np.fft.fftshift, matplotlib's imshow would put the
    diffraction peak at a corner with its 4-way symmetric tails distributed
    to the other three corners.

    The same convention bit run_06's density writer; the fix there was
    applied in C++. Here we keep the .dat layout untouched (so downstream
    tools see the canonical FFT-natural ordering) and shift only for
    visualisation.
    """
    arr = np.loadtxt(path)
    return np.fft.fftshift(arr)


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


# ─────────────────────────────────────────────────────────────────────────────
# 1. Observables + dipole spectrum
# ─────────────────────────────────────────────────────────────────────────────

def plot_observables() -> None:
    obs_csv = RESULTS / "observables.csv"
    if not obs_csv.exists():
        print("  Skip observables: observables.csv not found")
        return
    rows = []
    with obs_csv.open() as fh:
        r = csv.DictReader(fh)
        cols = r.fieldnames or []
        for row in r:
            if any(row.get(c) in (None, "") for c in cols):
                continue
            rows.append(row)
    if not rows:
        return
    t = np.array([float(row["time_au"]) for row in rows])
    VIS.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(3, 1, figsize=(8, 8), dpi=120, sharex=True)
    e_total = np.array([float(row["energy_total"]) for row in rows])
    e0 = e_total[0]
    ax[0].plot(t, (e_total - e0), color="black")
    ax[0].set_ylabel(r"$E - E_0$ (Ha)")
    ax[0].set_title("Energy drift")

    for c, key in zip(["C0", "C1", "C2"],
                      ["current_x", "current_y", "current_z"]):
        if key in cols:
            ax[1].plot(t, [float(row[key]) for row in rows],
                       color=c, label=key)
    ax[1].set_ylabel("Current (a.u.)")
    ax[1].legend(loc="best", fontsize=8)

    for c, key in zip(["C0", "C1", "C2"],
                      ["dipole_x", "dipole_y", "dipole_z"]):
        if key in cols:
            ax[2].plot(t, [float(row[key]) for row in rows],
                       color=c, label=key)
    ax[2].set_xlabel("time (a.u.)")
    ax[2].set_ylabel("Dipole (a.u.)")
    ax[2].legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(VIS / "observables.png")
    plt.close(fig)
    print(f"  Wrote observables.png")

    # Dipole power spectrum
    if "dipole_x" in cols:
        try:
            dx = np.array([float(row["dipole_x"]) for row in rows])
            dx -= dx.mean()
            n = len(dx)
            dt_au = float(t[1] - t[0])
            freqs_au = np.fft.rfftfreq(n, d=dt_au)
            mag2 = np.abs(np.fft.rfft(dx)) ** 2
            ev = freqs_au * 2 * np.pi * HA_TO_EV
            fig, axs = plt.subplots(figsize=(8, 4), dpi=120)
            axs.semilogy(ev, mag2, color="black", lw=0.8)
            axs.set_xlim(0, 30)
            axs.set_xlabel(r"$\omega$ (eV)")
            axs.set_ylabel(r"$|p_x(\omega)|^2$")
            axs.set_title("Dipole excitation spectrum (x)")
            fig.tight_layout()
            fig.savefig(VIS / "spectrum_dipole_x.png")
            plt.close(fig)
            print("  Wrote spectrum_dipole_x.png")
        except Exception as e:
            print(f"  WARN spectrum: {e}")


# ─────────────────────────────────────────────────────────────────────────────
# 2. LEED screens (full-run + paper window) + 20-panel grid + snapshot GIF
# ─────────────────────────────────────────────────────────────────────────────

def plot_leed_panels() -> None:
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
                arr = _load_screen_centred(dat)
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

    # 20-panel grid for the paper-window accumulator (most physically
    # interesting). Layout 4x5.
    win = RESULTS / "screens_leed_window"
    if win.exists():
        dats = sorted(win.glob("screen_*.dat"), key=_natural_key)
        if dats:
            arrs = [_load_screen_centred(p) for p in dats]
            n = len(arrs)
            cols = 5
            rows = (n + cols - 1) // cols
            fig, axs = plt.subplots(rows, cols, figsize=(3 * cols, 3 * rows),
                                    dpi=110)
            axs = np.atleast_2d(axs)
            vmax = max(a.max() for a in arrs)
            for k, (ax, arr, p) in enumerate(zip(axs.ravel(), arrs, dats)):
                ax.imshow(arr, cmap="inferno", origin="lower", vmin=0, vmax=vmax)
                ax.set_title(p.stem, fontsize=9)
                ax.set_xticks([]); ax.set_yticks([])
            for ax in axs.ravel()[n:]:
                ax.axis("off")
            fig.suptitle("LEED screens — paper window [T1,T2]")
            fig.tight_layout()
            fig.savefig(VIS / "leed_grid.png")
            plt.close(fig)
            print("  Wrote leed_grid.png (paper-window 20-panel view)")

    # Per-step snapshot GIF for one canonical screen index (10 = mid plane)
    snap_root = RESULTS / "screens_snapshots"
    if snap_root.exists():
        step_dirs = sorted(
            (d for d in snap_root.iterdir() if d.is_dir()),
            key=_natural_key,
        )
        if step_dirs:
            screen_name = "screen_10.dat"
            arrs, labels = [], []
            for d in step_dirs:
                p = d / screen_name
                if not p.exists():
                    continue
                try:
                    arrs.append(_load_screen_centred(p))
                    labels.append(d.name)
                except Exception:
                    continue
            if arrs:
                vmax = max(a.max() for a in arrs)
                frames = []
                for arr, lab in zip(arrs, labels):
                    fig, ax = plt.subplots(figsize=(4, 4), dpi=80)
                    ax.imshow(arr, cmap="inferno", origin="lower",
                              vmin=0, vmax=vmax)
                    ax.set_title(f"screen_10 {lab}")
                    fig.tight_layout()
                    fig.canvas.draw()
                    frames.append(np.asarray(fig.canvas.buffer_rgba())
                                  [:, :, :3].copy())
                    plt.close(fig)
                iio.imwrite(VIS / "leed_screen_10_evolution.gif", frames,
                            fps=GIF_FPS, loop=0)
                print("  Wrote leed_screen_10_evolution.gif")


# ─────────────────────────────────────────────────────────────────────────────
# 3. Density mid-slice GIFs (xy / xz / yz) for both target and WP series
# ─────────────────────────────────────────────────────────────────────────────

def plot_density_gifs() -> None:
    VIS.mkdir(parents=True, exist_ok=True)
    for series_name in ["density_rt_target", "density_rt_wp"]:
        src = RESULTS / series_name
        if not src.exists():
            continue
        vti_paths = sorted(src.glob("*.vti"), key=_natural_key)
        if not vti_paths:
            continue

        # Pre-load arrays once (memory: 61 frames × 23 MB = 1.4 GB — fine).
        print(f"  Loading {len(vti_paths)} VTI frames for {series_name}...")
        arrays, times = [], []
        for p in vti_paths:
            v = read_vti(p)
            arr = next(iter(v["arrays"].values()))  # only one PointData array
            arrays.append(arr)
            times.append(_vti_meta_time(p) or 0)
        global_max = max(a.max() for a in arrays)

        for axis_name, ax_idx in [("xz", 1), ("yz", 0), ("xy", 2)]:
            frames = []
            for arr, step in zip(arrays, times):
                mid = arr.shape[ax_idx] // 2
                if ax_idx == 0:
                    sl = arr[mid, :, :]
                elif ax_idx == 1:
                    sl = arr[:, mid, :]
                else:
                    sl = arr[:, :, mid]
                fig, ax = plt.subplots(figsize=(4, 4), dpi=80)
                ax.imshow(sl.T, cmap="inferno", origin="lower",
                          vmin=0, vmax=global_max)
                ax.set_title(f"{series_name} {axis_name} step={step}")
                fig.tight_layout()
                fig.canvas.draw()
                frames.append(np.asarray(fig.canvas.buffer_rgba())
                              [:, :, :3].copy())
                plt.close(fig)
            iio.imwrite(VIS / f"{series_name}_{axis_name}.gif",
                        frames, fps=GIF_FPS, loop=0)
        print(f"  Wrote {series_name} GIFs (xy, xz, yz)")


# ─────────────────────────────────────────────────────────────────────────────
# 4. WP-only overlap analysis: |<GS_i | WP(t)>|^2 vs time
# ─────────────────────────────────────────────────────────────────────────────

def plot_overlap_wp_only() -> None:
    overlap_dir = RESULTS / "overlap"
    index = overlap_dir / "index.csv"
    if not index.exists():
        print("  Skip overlap: index.csv not found")
        return

    rows = []
    with index.open() as fh:
        r = csv.DictReader(fh)
        for row in r:
            try:
                step = int(row["step"])
                t = float(row["time_au"])
            except (KeyError, ValueError):
                continue
            mat_path = overlap_dir / row["file"]
            if not mat_path.exists():
                continue
            arr = np.loadtxt(mat_path, delimiter=",", comments="#", ndmin=1)
            if arr.ndim == 0:
                arr = np.array([float(arr)])
            rows.append((step, t, arr))
    if not rows:
        print("  Skip overlap: no rows loaded")
        return

    rows.sort(key=lambda x: x[0])
    times = np.array([r[1] for r in rows])
    matrix = np.array([r[2] for r in rows])  # (n_steps, n_ref)
    n_steps, n_ref = matrix.shape
    print(f"  Overlap series: {n_steps} steps × {n_ref} GS reference orbitals")

    # Determine n_occupied for the line plot
    summary = _read_summary(RESULTS / "run_summary.txt")
    n_occupied = int(summary.get("n_occupied", 54))

    VIS.mkdir(parents=True, exist_ok=True)

    # Heatmap (orbital index vs time)
    fig, ax = plt.subplots(figsize=(10, 6), dpi=120)
    pcm = ax.pcolormesh(
        times, np.arange(n_ref), matrix.T,
        cmap="viridis", shading="auto",
        norm=matplotlib.colors.LogNorm(vmin=1e-12,
                                       vmax=max(matrix.max(), 1e-10)),
    )
    ax.axhline(n_occupied - 0.5, color="white", ls="--", lw=0.7,
               label=f"HOMO (i={n_occupied - 1})")
    ax.set_xlabel("time (a.u.)")
    ax.set_ylabel("GS reference orbital index")
    ax.set_title(r"$|\langle\psi_i^{GS}|\psi_{WP}(t)\rangle|^2$")
    fig.colorbar(pcm, ax=ax, label="overlap²")
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(VIS / "overlap_wp_heatmap.png")
    plt.close(fig)
    print("  Wrote overlap_wp_heatmap.png")

    # Top-N dominant orbitals
    peak = matrix.max(axis=0)
    top_n = 8
    top_idx = np.argsort(peak)[::-1][:top_n]
    fig, ax = plt.subplots(figsize=(10, 5), dpi=120)
    for i in top_idx:
        ax.plot(times, matrix[:, i], lw=1.0, label=f"i={i}")
    ax.set_yscale("log")
    ax.set_xlabel("time (a.u.)")
    ax.set_ylabel(r"$|\langle\psi_i^{GS}|\psi_{WP}(t)\rangle|^2$")
    ax.set_title(f"Top {top_n} GS orbitals overlapping with WP")
    ax.legend(ncol=2, fontsize=8, loc="best")
    fig.tight_layout()
    fig.savefig(VIS / "overlap_wp_top.png")
    plt.close(fig)
    print(f"  Wrote overlap_wp_top.png (top {top_n} orbitals)")


# ─────────────────────────────────────────────────────────────────────────────
# 5. Diagnostics text: N-electron conservation + summary
# ─────────────────────────────────────────────────────────────────────────────

def write_diagnostics() -> None:
    summary_src = RESULTS / "run_summary.txt"
    out = RESULTS / "diagnostics.txt"
    parts = []
    summary = _read_summary(summary_src)
    if summary:
        parts.append("# run summary")
        for k, v in summary.items():
            parts.append(f"{k} = {v}")

    target_dir = RESULTS / "density_rt_target"
    wp_dir = RESULTS / "density_rt_wp"
    if target_dir.exists() and wp_dir.exists():
        ts = sorted(target_dir.glob("*.vti"), key=_natural_key)
        ws = sorted(wp_dir.glob("*.vti"), key=_natural_key)
        if ts and ws:
            n_check = min(20, len(ts), len(ws))
            idxs = np.linspace(0, min(len(ts), len(ws)) - 1, n_check, dtype=int)
            ne_vals, t_vals = [], []
            for i in idxs:
                vt = read_vti(ts[i])
                vw = read_vti(ws[i])
                arr_t = next(iter(vt["arrays"].values()))
                arr_w = next(iter(vw["arrays"].values()))
                dx, dy, dz = vt["spacing"]
                n_total = float((arr_t + arr_w).sum() * dx * dy * dz)
                ne_vals.append(n_total)
                t_vals.append(_vti_meta_time(ts[i]))
            ne = np.array(ne_vals)
            parts.append("\n# N-electron conservation (target+wp integrated)")
            parts.append(f"n_samples = {len(ne)}")
            parts.append(f"mean = {ne.mean():.6f}")
            parts.append(f"std  = {ne.std():.6e}")
            parts.append(f"min  = {ne.min():.6f}")
            parts.append(f"max  = {ne.max():.6f}")
            parts.append(f"drift_pct = "
                         f"{100*ne.std()/max(ne.mean(),1e-30):.4f}")

    obs_csv = RESULTS / "observables.csv"
    if obs_csv.exists():
        with obs_csv.open() as fh:
            r = csv.DictReader(fh)
            es = [float(row["energy_total"])
                  for row in r if row.get("energy_total")]
        if es:
            parts.append("\n# Energy drift")
            parts.append(f"E_initial = {es[0]:.10f} Ha")
            parts.append(f"E_final   = {es[-1]:.10f} Ha")
            parts.append(f"drift     = {(es[-1] - es[0]):+.6e} Ha")
            parts.append(f"drift_per_step = "
                         f"{(es[-1] - es[0]) / max(len(es)-1, 1):+.3e} Ha")

    out.write_text("\n".join(parts) + "\n")
    print(f"  Wrote diagnostics.txt")


# ─────────────────────────────────────────────────────────────────────────────

def main() -> int:
    if not RESULTS.exists():
        print(f"No results at {RESULTS}")
        return 1
    print("=== run_propagate_paper_replica analysis ===")
    print("[1] Observables + spectrum");  plot_observables()
    print("[2] LEED panels + grid + GIF"); plot_leed_panels()
    print("[3] Density GIFs");             plot_density_gifs()
    print("[4] WP overlap plots");         plot_overlap_wp_only()
    print("[5] Diagnostics");              write_diagnostics()
    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
