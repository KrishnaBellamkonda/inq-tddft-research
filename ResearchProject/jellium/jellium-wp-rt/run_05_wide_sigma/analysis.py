"""analysis.py — jellium run_05_wide_sigma (N=38, 200 eV, σ=2.0 Å, +z)."""
from __future__ import annotations

import csv
import io
import sys
from pathlib import Path

import imageio.v3 as iio
import matplotlib.pyplot as plt
import numpy as np

RUN_DIR  = Path(__file__).parent
RESULTS  = RUN_DIR / "results"
VISDIR   = RESULTS / "visualisation"
VISDIR.mkdir(parents=True, exist_ok=True)

REPO_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO_ROOT / "inq-stack" / "python"))
from inqview import (
    AnimationSpec, FieldSeries, FourierTransform, ParaViewPipeline,
    VolumeRenderSpec, convert_real_series_to_vti, load_leed_pattern,
    load_observables, load_real_field, plot_density_slice,
    plot_leed_pattern, plot_observables_summary, plot_spectrum_summary,
)

PV_EXE = REPO_ROOT / "ParaView-6.1.0-MPI-Linux-Python3.12-x86_64" / "bin" / "pvbatch"

RUN_LABEL        = "05_wide_sigma"
N_SCREENS        = 20
N_VTI_FRAMES     = 50
GIF_FPS          = 6
OVERLAP_GIF_FPS  = 5

print(f"\n=== analysis: jellium run_{RUN_LABEL} ===")


def _subsample(files: list[Path], n: int) -> list[Path]:
    if len(files) <= n:
        return list(files)
    idx = np.round(np.linspace(0, len(files) - 1, n)).astype(int)
    return [files[i] for i in idx]


def _fig_to_rgb(fig) -> np.ndarray:
    buf = io.BytesIO()
    fig.savefig(buf, format="png")
    buf.seek(0)
    img = iio.imread(buf)
    return img[..., :3] if img.ndim == 3 and img.shape[2] == 4 else img


def _make_slice_gif(series_dir: Path, label: str, axis: int,
                    axis_name: str, output_path: Path) -> None:
    meta_files = sorted(series_dir.glob("*.meta.txt"))
    if not meta_files:
        print(f"  Skip GIF {label} ({axis_name}): no data")
        return
    selected = _subsample(meta_files, N_VTI_FRAMES)
    tmp_dir = output_path.parent / f"_tmp_{output_path.stem}"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    png_paths: list[Path] = []
    for i, mp in enumerate(selected):
        field = load_real_field(meta_path=mp)
        png = tmp_dir / f"frame_{i:04d}.png"
        plot_density_slice(field, png, axis=axis)
        png_paths.append(png)
    frames = [iio.imread(p) for p in png_paths]
    iio.imwrite(output_path, frames, fps=GIF_FPS, loop=0)
    for p in png_paths:
        p.unlink()
    tmp_dir.rmdir()
    print(f"  GIF {label} ({axis_name}): {output_path.name}")


def _load_overlap_matrix(csv_path: Path):
    if not csv_path.exists():
        return None
    rows = []
    with csv_path.open() as fh:
        for row in csv.reader(fh):
            try:
                rows.append([float(v) for v in row])
            except ValueError:
                pass
    return np.array(rows) if rows else None


print("\n[1] N-electron conservation")
total_dir = RESULTS / "density_rt_total"
if total_dir.exists():
    meta_files = sorted(total_dir.glob("*.meta.txt"))
    sampled = _subsample(meta_files, 20)
    nelec_vals = []
    for mp in sampled:
        f = load_real_field(meta_path=mp)
        dx, dy, dz = f.meta.spacing_bohr
        nelec_vals.append(float(f.array.sum() * dx * dy * dz))
    arr = np.array(nelec_vals)
    pct = 100.0 * arr.std() / arr.mean() if arr.mean() > 0 else float("inf")
    print(f"  N_elec mean={arr.mean():.4f}  std={arr.std():.4f}  "
          f"min={arr.min():.4f}  max={arr.max():.4f}  drift={pct:.4f}%")
    print("  PASS: conserved (<0.1%)" if pct < 0.1 else f"  WARN: N_elec drifts {pct:.3f}%")
else:
    print("  Skip: density_rt_total not found")

print("\n[2] Density consistency: total ≈ jellium + wp")
jell_dir = RESULTS / "density_rt_jellium"
wp_dir   = RESULTS / "density_rt_wp"
if total_dir.exists() and jell_dir.exists() and wp_dir.exists():
    t_files = sorted(total_dir.glob("*.meta.txt"))
    j_files = sorted(jell_dir.glob("*.meta.txt"))
    w_files = sorted(wp_dir.glob("*.meta.txt"))
    if t_files and j_files and w_files:
        mid = len(t_files) // 2
        ft = load_real_field(meta_path=t_files[mid])
        fj = load_real_field(meta_path=j_files[min(mid, len(j_files) - 1)])
        fw = load_real_field(meta_path=w_files[min(mid, len(w_files) - 1)])
        diff = np.abs(ft.array - fj.array - fw.array).max()
        ratio = diff / max(ft.array.max(), 1e-30)
        print(f"  max|total - jellium - wp| = {diff:.3e}  ratio = {ratio:.2e}")
        print("  PASS" if ratio < 1e-6 else "  WARN: inconsistency > 1e-6")
else:
    print("  Skip: density series not found")

print("\n[3] Observables summary")
obs_csv = RESULTS / "observables.csv"
if obs_csv.exists():
    plot_observables_summary(obs_csv, VISDIR / "observables_summary.png")
    print("  Saved: observables_summary.png")
else:
    print("  Skip: observables.csv not found")

print("\n[4] FFT spectra")
if obs_csv.exists():
    df = load_observables(obs_csv)
    ft_obj = FourierTransform()
    results_fft = []
    if "energy_total" in df.columns:
        results_fft.append(ft_obj.transform_energy(df))
    for col in ("current_x", "current_y", "current_z"):
        if col in df.columns:
            results_fft.append(ft_obj.transform_column(df, col))
    if results_fft:
        plot_spectrum_summary(results_fft, VISDIR / "spectra.png")
        print(f"  Saved: spectra.png ({len(results_fft)} panels)")
else:
    print("  Skip: no observables.csv")

print("\n[5] VTI conversion")
for series_name in ("density_rt_total", "density_rt_jellium", "density_rt_wp"):
    sdir = RESULTS / series_name
    if not sdir.exists():
        print(f"  Skip {series_name}: not found"); continue
    all_files = sorted(sdir.glob("*.meta.txt"))
    if not all_files:
        print(f"  Skip {series_name}: no frames"); continue
    subset = _subsample(all_files, N_VTI_FRAMES)
    tag = series_name.replace("density_rt_", "")
    series_obj = FieldSeries(root=sdir, files=subset, field_name=tag)
    convert_real_series_to_vti(series_obj, VISDIR / "vti" / series_name)
    print(f"  {series_name}: {len(subset)} frames → vti/{series_name}/")

print("\n[6] Density slice GIFs")
gif_dir = VISDIR / "gifs"
gif_dir.mkdir(parents=True, exist_ok=True)
for series_name in ("density_rt_total", "density_rt_jellium", "density_rt_wp"):
    sdir = RESULTS / series_name
    tag  = series_name.replace("density_rt_", "")
    _make_slice_gif(sdir, tag, axis=1, axis_name="xz",
                    output_path=gif_dir / f"{tag}_xz.gif")
    _make_slice_gif(sdir, tag, axis=0, axis_name="yz",
                    output_path=gif_dir / f"{tag}_yz.gif")

print("\n[7] LEED pattern grid")
screens_dir = RESULTS / "screens"
if screens_dir.exists():
    dat_files = sorted(screens_dir.glob("*.dat"))
    if len(dat_files) == N_SCREENS:
        patterns = [load_leed_pattern(d) for d in dat_files]
        fig, axes = plt.subplots(4, 5, figsize=(20, 16), dpi=100)
        vmin = min(p.data.min() for p in patterns)
        vmax = max(p.data.max() for p in patterns)
        last_im = None
        for ax, p in zip(axes.flat, patterns):
            last_im = ax.imshow(p.data, origin="lower", extent=p.extent_bohr,
                                vmin=vmin, vmax=vmax, cmap="inferno", aspect="equal")
            ax.set_title(f"z={p.z_bohr:.1f} bohr", fontsize=7)
            ax.set_xticks([]); ax.set_yticks([])
        if last_im is not None:
            fig.colorbar(last_im, ax=axes.flat[-1], label="ρ·dt (bohr⁻³·a.u.)")
        fig.suptitle(f"LEED patterns — {RUN_LABEL}", fontsize=12)
        fig.tight_layout()
        fig.savefig(VISDIR / "leed_grid.png", dpi=100)
        plt.close(fig)
        print("  Saved: leed_grid.png")
        ind_dir = VISDIR / "screens_individual"
        ind_dir.mkdir(exist_ok=True)
        for p, d in zip(patterns, dat_files):
            plot_leed_pattern(p, ind_dir / (d.stem + ".png"))
        print(f"  Saved: {len(dat_files)} individual PNGs")
    else:
        print(f"  Found {len(dat_files)} screens (expected {N_SCREENS})")
else:
    print("  Skip: results/screens not found")

print("\n[8] LEED time-evolution GIF")
snap_root = RESULTS / "screens_snapshots"
if snap_root.exists():
    snap_dirs = sorted(snap_root.glob("step_*"))
    frames_leed = []
    for snap_step in snap_dirs:
        dat = snap_step / "screen_10.dat"
        if not dat.exists():
            continue
        p = load_leed_pattern(dat)
        fig, ax = plt.subplots(figsize=(5, 5), dpi=80)
        ax.imshow(p.data, origin="lower", extent=p.extent_bohr,
                  cmap="inferno", aspect="equal")
        ax.set_title(f"screen_10  z={p.z_bohr:.1f} bohr\n"
                     f"t={p.total_time_au:.2f} a.u.  n={p.n_accum}", fontsize=8)
        ax.set_xticks([]); ax.set_yticks([])
        fig.tight_layout(pad=0.5)
        frames_leed.append(_fig_to_rgb(fig))
        plt.close(fig)
    if frames_leed:
        out = gif_dir / "leed_screen10_evolution.gif"
        iio.imwrite(out, frames_leed, fps=GIF_FPS, loop=0)
        print(f"  Saved: {out.name} ({len(frames_leed)} frames)")
    else:
        print("  No screen_10 snapshots found")
else:
    print("  Skip: screens_snapshots not found")

print("\n[9] Overlap matrix GIFs")
overlap_dir = RESULTS / "overlap"
index_csv   = overlap_dir / "index.csv"
if index_csv.exists():
    index_rows = []
    with index_csv.open() as fh:
        for row in csv.DictReader(fh):
            index_rows.append(row)
    if index_rows:
        first_mat = _load_overlap_matrix(overlap_dir / index_rows[0]["file"])
        if first_mat is not None:
            n_ref, n_evolved = first_mat.shape
            print(f"  Matrix shape: {n_ref}×{n_evolved}  ({len(index_rows)} snapshots)")
            matrices, times = [], []
            for row in index_rows:
                m = _load_overlap_matrix(overlap_dir / row["file"])
                if m is not None and m.shape == (n_ref, n_evolved):
                    matrices.append(m)
                    times.append(float(row["time_au"]))
            if matrices:
                ov_dir = gif_dir / "overlap"
                ov_dir.mkdir(parents=True, exist_ok=True)
                for j in range(n_evolved):
                    col_data = [m[:, j] for m in matrices]
                    ymax = max(c.max() for c in col_data) * 1.05 + 1e-9
                    frames_ov = []
                    for t, col in zip(times, col_data):
                        fig, ax = plt.subplots(figsize=(6, 3), dpi=80)
                        ax.bar(np.arange(n_ref), col, color="steelblue", edgecolor="none")
                        ax.set_xlim(-0.5, n_ref - 0.5)
                        ax.set_ylim(0, ymax)
                        ax.set_xlabel("GS orbital index i")
                        ax.set_ylabel(r"$|O_{ij}|^2$")
                        ax.set_title(f"Overlap: evolved orbital j={j}  t={t:.2f} a.u.")
                        fig.tight_layout()
                        frames_ov.append(_fig_to_rgb(fig))
                        plt.close(fig)
                    iio.imwrite(ov_dir / f"overlap_orbital_{j:04d}.gif",
                                frames_ov, fps=OVERLAP_GIF_FPS, loop=0)
                print(f"  Saved {n_evolved} overlap GIFs → overlap/")
else:
    print("  Skip: results/overlap/index.csv not found")

print("\n[10] ParaView 3D renders")
if PV_EXE.exists():
    sdir = RESULTS / "density_rt_total"
    if sdir.exists():
        all_pv = sorted(sdir.glob("*.meta.txt"))
        if all_pv:
            subset_pv = _subsample(all_pv, 20)
            series_pv = FieldSeries(root=sdir, files=subset_pv, field_name="total")
            render_spec = VolumeRenderSpec(array_name="total")
            anim_spec   = AnimationSpec(output_frames_dir=VISDIR / "paraview" / "frames")
            try:
                pipeline = ParaViewPipeline(pv_executable=PV_EXE)
                pipeline.render_density_from_meta_series(
                    series_pv, VISDIR / "paraview" / "vti", render_spec, anim_spec)
                print("  ParaView render complete")
            except Exception as exc:
                print(f"  ParaView render failed: {exc}")
        else:
            print("  Skip: no density frames")
    else:
        print("  Skip: density_rt_total not found")
else:
    print("  Skip: pvbatch not found")

print(f"\nDone. Outputs in {VISDIR.relative_to(RUN_DIR)}/")
