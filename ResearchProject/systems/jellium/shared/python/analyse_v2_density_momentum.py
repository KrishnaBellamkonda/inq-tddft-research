#!/usr/bin/env python3
"""Reusable analysis functions for jellium density evolution and
momentum/trajectory observables.

Each function takes a *results_dir* (the run directory, e.g.
``run_wp_n162_L50_E100_sigma1_v2``) and a *run_name* label for titles.
Outputs land in ``<results_dir>/results/analysis/observables/``.

VTI files are loaded with vtk (NOT pyvista).

Interpreter: /local/data/public/skcb2/tddft/venv/bin/python3
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
HA_TO_EV = 27.211386245988

# ---------------------------------------------------------------------------
# VTI loading
# ---------------------------------------------------------------------------

def _load_vti(path):
    """Load a VTI file and return a 3D numpy array shaped (nx, ny, nz)."""
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy

    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    img = reader.GetOutput()
    nx, ny, nz = img.GetDimensions()
    flat = vtk_to_numpy(img.GetPointData().GetArray(0)).astype(np.float64)
    return flat.reshape((nz, ny, nx)).transpose(2, 1, 0)


def _load_vti_with_meta(path):
    """Load VTI → (cube[nx,ny,nz], meta dict with origin/spacing)."""
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy

    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    img = reader.GetOutput()
    nx, ny, nz = img.GetDimensions()
    flat = vtk_to_numpy(img.GetPointData().GetArray(0)).astype(np.float64)
    cube = flat.reshape((nz, ny, nx)).transpose(2, 1, 0)
    meta = {
        "nx": nx, "ny": ny, "nz": nz,
        "origin": tuple(img.GetOrigin()),
        "spacing": tuple(img.GetSpacing()),
    }
    return cube, meta


# ---------------------------------------------------------------------------
# run_summary.txt parser
# ---------------------------------------------------------------------------
_KV_RE = re.compile(r"^(\w[\w.]+)\s*=\s*(.+)$")


def _parse_run_summary(results_dir: Path) -> dict:
    """Parse run_summary.txt into a flat dict of str values.

    Searches in ``results_dir/results/run_summary.txt`` first, then
    ``results_dir/run_summary.txt``.
    """
    for candidate in [
        results_dir / "results" / "run_summary.txt",
        results_dir / "run_summary.txt",
    ]:
        if candidate.exists():
            out = {}
            for line in candidate.read_text().splitlines():
                m = _KV_RE.match(line.strip())
                if m:
                    out[m.group(1)] = m.group(2).strip()
            return out
    return {}


def _get_float(summary: dict, key: str, default: float) -> float:
    """Extract a float from run_summary dict, falling back to *default*."""
    val = summary.get(key)
    if val is None:
        return default
    # Handle "50^3 (cubic, periodic)" → take first token
    tok = val.split()[0].replace("^3", "")
    try:
        return float(tok)
    except ValueError:
        return default


# ---------------------------------------------------------------------------
# wp_config.txt parser (fallback for sigma, k0 when run_summary missing)
# ---------------------------------------------------------------------------

def _parse_wp_config(results_dir: Path) -> dict:
    """Parse wp_config.txt for sigma_bohr, k0, energy_ev."""
    p = results_dir / "results" / "raw" / "observables" / "wp_config.txt"
    if not p.exists():
        return {}
    out = {}
    for line in p.read_text().splitlines():
        m = _KV_RE.match(line.strip())
        if m:
            out[m.group(1)] = m.group(2).strip()
    return out


# ---------------------------------------------------------------------------
# VTI frame listing
# ---------------------------------------------------------------------------
_STEP_RE = re.compile(r"_t(\d{6})\.vti$")


def _list_vti_frames(vti_dir: Path):
    """Return (sorted file list, step array, time_au array) for a VTI dir.

    *dt_au* is inferred from the CSV timestamps rather than needing
    run_summary; if unavailable, returns step indices as "times".
    """
    if not vti_dir.exists():
        return [], np.array([]), np.array([])
    files = sorted(vti_dir.glob("*.vti"))
    steps = []
    for f in files:
        m = _STEP_RE.search(f.name)
        steps.append(int(m.group(1)) if m else -1)
    return files, np.array(steps, dtype=int), steps


def _steps_to_times(steps: np.ndarray, dt_au: float) -> np.ndarray:
    return steps.astype(np.float64) * dt_au


# ---------------------------------------------------------------------------
# Helpers for extracting run parameters
# ---------------------------------------------------------------------------

def _run_params(results_dir: Path) -> dict:
    """Return a dict with commonly needed run parameters."""
    summary = _parse_run_summary(results_dir)
    wp_cfg = _parse_wp_config(results_dir)

    dt_au = _get_float(summary, "dt_au", 0.01)
    spacing = _get_float(summary, "spacing_bohr", 0.40)
    L = _get_float(summary, "cell_bohr", 50.0)

    # sigma from run_summary or wp_config
    sigma_str = summary.get("wp_sigma_bohr", wp_cfg.get("wp_sigma_bohr", "1"))
    sigma = float(sigma_str.split()[0])

    # k0 from run_summary or wp_config
    k0_str = summary.get("wp_k0_bohr_inv", wp_cfg.get("wp_k0_bohr_inv", "0 0 2.711"))
    k0_parts = k0_str.split()
    k0 = float(k0_parts[-1]) if len(k0_parts) > 1 else float(k0_parts[0])

    ekin_str = summary.get("wp_energy_ev", wp_cfg.get("wp_energy_ev", "100"))
    ekin_ev = float(ekin_str.split()[0])

    write_every = int(_get_float(summary, "write_every", 3))

    return {
        "dt_au": dt_au,
        "spacing_bohr": spacing,
        "L_bohr": L,
        "sigma_bohr": sigma,
        "k0": k0,
        "ekin_ev": ekin_ev,
        "write_every": write_every,
    }


def _ensure_output_dir(results_dir: Path) -> Path:
    out = results_dir / "results" / "analysis" / "observables"
    out.mkdir(parents=True, exist_ok=True)
    return out


def _pick_closest(arr: np.ndarray, target: float) -> int:
    return int(np.argmin(np.abs(arr - target)))


# ===========================================================================
# 1. plot_density_z_profile_evolution
# ===========================================================================

def plot_density_z_profile_evolution(
    results_dir: str | Path,
    run_name: str,
) -> Optional[Path]:
    """Heatmap of xy-averaged δn(z) vs time from density_delta_coarse VTIs.

    Returns the output path, or None if data is missing.
    """
    results_dir = Path(results_dir)
    params = _run_params(results_dir)
    dt_au = params["dt_au"]

    vti_dir = results_dir / "results" / "raw" / "vti" / "density_delta_coarse"
    files, steps, _ = _list_vti_frames(vti_dir)
    if len(files) == 0:
        print(f"[density_z_profile_evolution] No VTI files in {vti_dir}")
        return None

    times = _steps_to_times(steps, dt_au)
    print(f"[density_z_profile_evolution] Loading {len(files)} coarse frames ...")

    # Load first frame to get grid info
    cube0, meta = _load_vti_with_meta(files[0])
    nx, ny, nz = cube0.shape
    oz = meta["origin"][2]
    dz = meta["spacing"][2]
    z_axis = oz + np.arange(nz) * dz

    # Build 2D array: z_profile[frame_idx, z_idx]
    z_profiles = np.empty((len(files), nz), dtype=np.float64)
    for i, f in enumerate(files):
        cube = _load_vti(f)
        # Average over x and y → shape (nz,)
        z_profiles[i, :] = cube.mean(axis=(0, 1))
        if (i + 1) % 50 == 0:
            print(f"  ... loaded {i + 1}/{len(files)}")

    print(f"  ... all {len(files)} frames loaded.")

    # Plot heatmap
    out_dir = _ensure_output_dir(results_dir)
    fig, ax = plt.subplots(figsize=(10, 6))

    vmax = float(np.percentile(np.abs(z_profiles), 99.5))
    if vmax == 0:
        vmax = 1e-10
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    im = ax.pcolormesh(
        z_axis, times, z_profiles,
        cmap="RdBu_r", norm=norm, shading="nearest",
    )
    cb = fig.colorbar(im, ax=ax, label=r"$\langle\delta n\rangle_{xy}$ (a.u.$^{-3}$)")

    ax.set_xlabel("z (Bohr)")
    ax.set_ylabel("Time (a.u.)")
    ax.set_title(f"Density perturbation z-profile evolution — {run_name}")

    out_path = out_dir / "density_z_profile_evolution.png"
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[density_z_profile_evolution] Wrote {out_path}")
    return out_path


# ===========================================================================
# 2. plot_delta_density_xz_snapshots
# ===========================================================================

def plot_delta_density_xz_snapshots(
    results_dir: str | Path,
    run_name: str,
    n_snapshots: int = 4,
) -> Optional[Path]:
    """1xN grid of delta-n(x, y=0, z) at equally-spaced times.

    Returns the output path, or None if data is missing.
    """
    results_dir = Path(results_dir)
    params = _run_params(results_dir)
    dt_au = params["dt_au"]

    vti_dir = results_dir / "results" / "raw" / "vti" / "density_delta"
    files, steps, _ = _list_vti_frames(vti_dir)
    if len(files) == 0:
        print(f"[delta_density_xz_snapshots] No VTI files in {vti_dir}")
        return None

    times = _steps_to_times(steps, dt_au)

    # Select equally-spaced times within the full time range
    t_min, t_max = float(times[0]), float(times[-1])
    target_times = np.linspace(t_min, t_max, n_snapshots + 2)[1:-1]  # exclude endpoints
    selected_idx = [_pick_closest(times, t) for t in target_times]

    print(f"[delta_density_xz_snapshots] Loading {n_snapshots} snapshots ...")

    # Load selected frames
    slices = []
    for idx in selected_idx:
        cube, meta = _load_vti_with_meta(files[idx])
        nx, ny, nz = cube.shape
        y_mid = ny // 2
        sl = cube[:, y_mid, :]  # shape (nx, nz)
        slices.append({
            "slice": sl,
            "time": float(times[idx]),
            "meta": meta,
        })

    # Determine global colorscale
    all_abs = np.concatenate([np.abs(s["slice"]).ravel() for s in slices])
    vmax = float(np.percentile(all_abs, 99.5))
    if vmax == 0:
        vmax = 1e-10
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    # Plot
    out_dir = _ensure_output_dir(results_dir)
    fig, axes = plt.subplots(1, n_snapshots, figsize=(4.0 * n_snapshots + 1.0, 4.5))
    if n_snapshots == 1:
        axes = [axes]

    for ax, s in zip(axes, slices):
        sl = s["slice"]
        m = s["meta"]
        nx, nz = sl.shape
        ox, oz = m["origin"][0], m["origin"][2]
        dx, dz = m["spacing"][0], m["spacing"][2]
        extent = [oz, oz + (nz - 1) * dz, ox, ox + (nx - 1) * dx]
        im = ax.imshow(
            sl, origin="lower", cmap="RdBu_r", norm=norm,
            extent=extent, aspect="equal",
        )
        ax.set_title(f"t = {s['time']:.2f} a.u.", fontsize=10)
        ax.set_xlabel("z (Bohr)")

    axes[0].set_ylabel("x (Bohr)")
    for ax in axes[1:]:
        ax.set_yticklabels([])

    fig.colorbar(im, ax=axes, shrink=0.85, label=r"$\delta n$ (a.u.$^{-3}$)")
    fig.suptitle(
        rf"$\delta n(x, y{{=}}0, z)$ — {run_name}",
        y=1.02, fontsize=12,
    )

    out_path = out_dir / "delta_density_xz_snapshots.png"
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[delta_density_xz_snapshots] Wrote {out_path}")
    return out_path


# ===========================================================================
# 3. plot_momentum_band_free_vs_jellium
# ===========================================================================

def plot_momentum_band_free_vs_jellium(
    results_dir: str | Path,
    run_name: str,
    free_results_dir: Optional[str | Path] = None,
) -> Optional[Path]:
    """2-panel momentum band plot: free WP (left) vs jellium WP (right).

    If *free_results_dir* is not given, the left panel uses analytical
    free-particle values (constant p_z = k0, sigma_p = 1/(sigma*sqrt(2))).

    Returns the output path, or None if data is missing.
    """
    results_dir = Path(results_dir)
    params = _run_params(results_dir)
    k0 = params["k0"]
    sigma_r = params["sigma_bohr"]
    L = params["L_bohr"]

    # Load jellium momentum data
    mom_path = results_dir / "results" / "raw" / "observables" / "wp_momentum_stats.csv"
    pos_path = results_dir / "results" / "raw" / "observables" / "wp_real_space_stats.csv"
    if not mom_path.exists():
        print(f"[momentum_band] Missing {mom_path}")
        return None

    df_mom = pd.read_csv(mom_path, comment="#")
    df_mom["sigma_pz"] = np.sqrt(df_mom["sigma_pz2"].clip(lower=0))

    df_pos = pd.read_csv(pos_path, comment="#") if pos_path.exists() else None

    # Build free-WP data (analytical or from a run)
    if free_results_dir is not None:
        free_dir = Path(free_results_dir)
        free_mom_path = free_dir / "results" / "raw" / "observables" / "wp_momentum_stats.csv"
        free_pos_path = free_dir / "results" / "raw" / "observables" / "wp_real_space_stats.csv"
        if free_mom_path.exists():
            df_free_mom = pd.read_csv(free_mom_path, comment="#")
            df_free_mom["sigma_pz"] = np.sqrt(df_free_mom["sigma_pz2"].clip(lower=0))
            df_free_pos = (
                pd.read_csv(free_pos_path, comment="#") if free_pos_path.exists() else None
            )
            free_analytical = False
        else:
            print(f"[momentum_band] No free momentum CSV at {free_mom_path}, using analytical")
            free_analytical = True
    else:
        free_analytical = True

    # Prepare panels
    out_dir = _ensure_output_dir(results_dir)
    fig, (ax_free, ax_jel) = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

    # --- Left panel: free WP ---
    if free_analytical:
        # Analytical free-particle: p_z = k0 (constant), sigma_p = 1/(sigma*sqrt(2))
        sigma_p0 = 1.0 / (sigma_r * np.sqrt(2.0))
        if df_pos is not None:
            # Use jellium run's time axis to build matching z axis
            t = df_mom["time_au"].to_numpy()
            z0 = float(df_pos["z_mean"].iloc[0]) if "z_mean" in df_pos.columns else -21.0
            z_free = z0 + k0 * t
        else:
            t = df_mom["time_au"].to_numpy()
            z_free = -21.0 + k0 * t
        pz_free = np.full_like(t, k0)
        sigma_pz_free = np.full_like(t, sigma_p0)

        ax_free.fill_between(
            z_free, pz_free - sigma_pz_free, pz_free + sigma_pz_free,
            color="C0", alpha=0.25, label=r"$\langle p_z\rangle \pm \sigma_{p_z}$",
        )
        ax_free.plot(z_free, pz_free, color="C0", lw=1.6,
                     label=r"$\langle p_z\rangle$ (analytical)")
        ax_free.axhline(k0, color="gray", lw=0.8, ls="--", label=r"$p_0 = k_0$")
        ax_free.set_title("Free WP (analytical)")
    else:
        t = df_free_mom["time_au"].to_numpy()
        pz = df_free_mom["pz_mean"].to_numpy()
        sigma_pz = df_free_mom["sigma_pz"].to_numpy()
        if df_free_pos is not None and "z_mean" in df_free_pos.columns:
            z_c = np.interp(t, df_free_pos["time_au"].to_numpy(),
                            df_free_pos["z_mean"].to_numpy())
        else:
            z_c = -21.0 + pz[0] * t

        ax_free.fill_between(
            z_c, pz - sigma_pz, pz + sigma_pz,
            color="C0", alpha=0.25, label=r"$\langle p_z\rangle \pm \sigma_{p_z}$",
        )
        ax_free.plot(z_c, pz, color="C0", lw=1.6,
                     label=r"$\langle p_z\rangle(t)$")
        ax_free.axhline(float(pz[0]), color="gray", lw=0.8, ls="--",
                        label=r"$p_0$")
        ax_free.set_title("Free WP")

    # --- Right panel: jellium WP ---
    t_j = df_mom["time_au"].to_numpy()
    pz_j = df_mom["pz_mean"].to_numpy()
    sigma_pz_j = df_mom["sigma_pz"].to_numpy()

    if df_pos is not None and "z_mean" in df_pos.columns:
        z_j = np.interp(t_j, df_pos["time_au"].to_numpy(),
                        df_pos["z_mean"].to_numpy())
    else:
        z_j = -21.0 + k0 * t_j

    # Jellium region shading
    ax_jel.axvspan(-L / 2, +L / 2, alpha=0.06, color="gray", zorder=0)
    ax_jel.axvline(-L / 2, color="dimgray", lw=0.7, ls=":", zorder=1)
    ax_jel.axvline(+L / 2, color="dimgray", lw=0.7, ls=":", zorder=1)

    ax_jel.fill_between(
        z_j, pz_j - sigma_pz_j, pz_j + sigma_pz_j,
        color="C3", alpha=0.25, label=r"$\langle p_z\rangle \pm \sigma_{p_z}$",
    )
    ax_jel.plot(z_j, pz_j, color="C3", lw=1.6,
                label=r"$\langle p_z\rangle(t)$")
    ax_jel.axhline(float(pz_j[0]), color="gray", lw=0.8, ls="--",
                   label=r"$p_0$")
    ax_jel.set_title(f"Jellium WP — {run_name}")

    # Shared formatting
    for ax in (ax_free, ax_jel):
        ax.set_xlabel(r"Projectile centroid $\langle z\rangle$ (Bohr)")
        ax.legend(loc="best", fontsize=8.5)
        ax.grid(True, alpha=0.3)
    ax_free.set_ylabel(r"Momentum $p_z$ (a.u.)")

    fig.suptitle(
        f"Momentum band: free vs jellium — {run_name}",
        y=1.01, fontsize=12,
    )
    fig.tight_layout()
    out_path = out_dir / "momentum_band_free_vs_jellium.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[momentum_band] Wrote {out_path}")
    return out_path


# ===========================================================================
# 4. plot_sigma_xyz_vs_time
# ===========================================================================

def plot_sigma_xyz_vs_time(
    results_dir: str | Path,
    run_name: str,
) -> Optional[Path]:
    """Plot sigma_x, sigma_y, sigma_z vs time with analytical free-WP overlay.

    Returns the output path, or None if data is missing.
    """
    results_dir = Path(results_dir)
    params = _run_params(results_dir)
    sigma0 = params["sigma_bohr"]

    pos_path = results_dir / "results" / "raw" / "observables" / "wp_real_space_stats.csv"
    if not pos_path.exists():
        print(f"[sigma_xyz] Missing {pos_path}")
        return None

    df = pd.read_csv(pos_path, comment="#")
    t = df["time_au"].to_numpy()

    sigma_x = np.sqrt(df["sigma_x2"].clip(lower=0).to_numpy())
    sigma_y = np.sqrt(df["sigma_y2"].clip(lower=0).to_numpy())
    sigma_z = np.sqrt(df["sigma_z2"].clip(lower=0).to_numpy())

    # Analytical free-particle spreading: sigma_i(t) = (sigma0/sqrt(2)) * sqrt(1 + (t/sigma0^2)^2)
    # Note: sigma_x2 from the CSV is the variance <(x-<x>)^2>, so sqrt gives the std dev.
    # For a free Gaussian WP with initial half-width sigma0,
    # the wavefunction sigma is sigma0, so the density variance is sigma0^2/2.
    # The density std dev at t=0 is sigma0/sqrt(2).
    # Time evolution: sigma_density(t) = (sigma0/sqrt(2)) * sqrt(1 + (t/sigma0^2)^2)
    t_fine = np.linspace(0, float(t[-1]), 500)
    sigma_free = (sigma0 / np.sqrt(2.0)) * np.sqrt(1.0 + (t_fine / sigma0**2) ** 2)

    out_dir = _ensure_output_dir(results_dir)
    fig, ax = plt.subplots(figsize=(9, 5.5))

    ax.plot(t, sigma_x, color="C0", lw=1.6, label=r"$\sigma_x$ (TDDFT)")
    ax.plot(t, sigma_y, color="C1", lw=1.6, label=r"$\sigma_y$ (TDDFT)")
    ax.plot(t, sigma_z, color="C3", lw=1.6, label=r"$\sigma_z$ (TDDFT)")
    ax.plot(
        t_fine, sigma_free,
        color="gray", lw=1.2, ls="--",
        label=rf"Free-WP analytical ($\sigma_0$ = {sigma0:.2g} Bohr)",
    )

    ax.set_xlabel("Time (a.u.)")
    ax.set_ylabel(r"WP width $\sigma_i = \sqrt{\langle (x_i - \langle x_i\rangle)^2 \rangle}$ (Bohr)")
    ax.set_title(rf"WP spreading: $\sigma_x, \sigma_y, \sigma_z$ vs time — {run_name}")
    ax.legend(loc="best", fontsize=9)
    ax.grid(True, alpha=0.3)

    out_path = out_dir / "sigma_xyz_vs_time.png"
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[sigma_xyz] Wrote {out_path}")
    return out_path


# ===========================================================================
# 5. plot_density_diff_vs_free
# ===========================================================================

def plot_density_diff_vs_free(
    results_dir: str | Path,
    run_name: str,
    free_results_dir: Optional[str | Path] = None,
    n_snapshots: int = 4,
) -> Optional[Path]:
    """2xN grid: top row = free WP delta-n, bottom = jellium delta-n at matched times.

    If *free_results_dir* is not provided, compute an analytical free-WP
    density and show the jellium delta-n on both rows (top = analytical free,
    bottom = TDDFT jellium).

    Returns the output path, or None if data is missing.
    """
    results_dir = Path(results_dir)
    params = _run_params(results_dir)
    dt_au = params["dt_au"]
    k0 = params["k0"]
    sigma_r = params["sigma_bohr"]

    # Jellium density_delta VTIs
    jel_vti_dir = results_dir / "results" / "raw" / "vti" / "density_delta"
    jel_files, jel_steps, _ = _list_vti_frames(jel_vti_dir)
    if len(jel_files) == 0:
        print(f"[density_diff_vs_free] No jellium delta VTIs in {jel_vti_dir}")
        return None
    jel_times = _steps_to_times(jel_steps, dt_au)

    # Select times (equally spaced, excluding very start/end)
    t_min, t_max = float(jel_times[1]), float(jel_times[-2])
    target_times = np.linspace(t_min, t_max, n_snapshots + 2)[1:-1]
    jel_idx = [_pick_closest(jel_times, t) for t in target_times]

    # Free WP data
    has_free_vti = False
    if free_results_dir is not None:
        free_dir = Path(free_results_dir)
        free_vti_dir = free_dir / "results" / "raw" / "vti" / "density_delta"
        if not free_vti_dir.exists():
            # Try density_rt_delta
            free_vti_dir = free_dir / "results" / "raw" / "vti" / "density_rt_delta"
        if free_vti_dir.exists():
            free_files, free_steps, _ = _list_vti_frames(free_vti_dir)
            if len(free_files) > 0:
                free_params = _run_params(free_dir)
                free_dt = free_params["dt_au"]
                free_times = _steps_to_times(free_steps, free_dt)
                has_free_vti = True

    print(f"[density_diff_vs_free] Loading {n_snapshots} snapshot pairs ...")

    # Load jellium slices
    jel_slices = []
    for idx in jel_idx:
        cube, meta = _load_vti_with_meta(jel_files[idx])
        nx, ny, nz = cube.shape
        sl = cube[:, ny // 2, :]
        jel_slices.append({
            "slice": sl, "time": float(jel_times[idx]), "meta": meta,
        })

    # Load or compute free slices
    free_slices = []
    if has_free_vti:
        for s in jel_slices:
            t_target = s["time"]
            fi = _pick_closest(free_times, t_target)
            cube, meta = _load_vti_with_meta(free_files[fi])
            nx, ny, nz = cube.shape
            sl = cube[:, ny // 2, :]
            free_slices.append({
                "slice": sl, "time": float(free_times[fi]), "meta": meta,
            })
    else:
        # Analytical free-WP delta-n: Gaussian in (x,y,z) with spreading
        # dn(x,y,z,t) = |psi(t)|^2 - |psi(0)|^2 evaluated at y=0
        # |psi|^2 = (2*pi*sigma_t^2)^{-3/2} * exp(-r^2/(2*sigma_t^2))
        # where sigma_t^2 = sigma0^2/2 + t^2/(2*sigma0^2)  [density variance]
        # and r^2 = (x-x0)^2 + y^2 + (z-z0-k0*t)^2
        for s in jel_slices:
            m = s["meta"]
            nx, nz = m["nx"], m["nz"]
            ny = m["ny"]
            ox, oz = m["origin"][0], m["origin"][2]
            dx, dz = m["spacing"][0], m["spacing"][2]

            x = ox + np.arange(nx) * dx
            z = oz + np.arange(nz) * dz
            X, Z = np.meshgrid(x, z, indexing="ij")  # (nx, nz)

            t_now = s["time"]
            # Initial center
            z0 = -21.0  # default launch z
            x0 = 0.0

            # Density variance at time t: (sigma0^2 / 2) * (1 + (t/(sigma0^2))^2)
            var_0 = sigma_r**2 / 2.0
            var_t = var_0 * (1.0 + (t_now / sigma_r**2) ** 2)

            # Current Gaussian at t
            r2_now = (X - x0)**2 + (Z - z0 - k0 * t_now)**2
            dens_now = (2 * np.pi * var_t) ** (-1.0) * np.exp(-r2_now / (2.0 * var_t))
            # Note: this is the 2D marginal (y=0 slice); the y-integral gives
            # (2*pi*var_t)^{-1} for 2D, but for proper 3D normalization at y=0:
            dens_now = (2 * np.pi * var_t) ** (-1.5) * np.exp(-r2_now / (2.0 * var_t))
            # (The y=0 slice of a 3D Gaussian has the full 3D prefactor.)

            # Gaussian at t=0
            r2_0 = (X - x0)**2 + (Z - z0)**2
            dens_0 = (2 * np.pi * var_0) ** (-1.5) * np.exp(-r2_0 / (2.0 * var_0))

            delta = dens_now - dens_0

            free_slices.append({
                "slice": delta, "time": t_now, "meta": m,
            })

    # Determine shared colorscale
    all_vals = np.concatenate(
        [np.abs(s["slice"]).ravel() for s in jel_slices + free_slices]
    )
    vmax = float(np.percentile(all_vals, 99.5))
    if vmax == 0:
        vmax = 1e-10
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    # Plot: 2 rows x n_snapshots cols
    out_dir = _ensure_output_dir(results_dir)
    fig = plt.figure(figsize=(4.0 * n_snapshots + 1.5, 8.5))
    gs = gridspec.GridSpec(
        2, n_snapshots + 1,
        width_ratios=[1] * n_snapshots + [0.05],
        figure=fig, wspace=0.08, hspace=0.15,
        left=0.06, right=0.94, top=0.90, bottom=0.08,
    )

    row_labels = [
        "Free WP" if (has_free_vti or free_results_dir is None) else "Free (analytical)",
        f"Jellium — {run_name}",
    ]
    all_slices = [free_slices, jel_slices]

    for r, (row_data, row_label) in enumerate(zip(all_slices, row_labels)):
        for c, s in enumerate(row_data):
            ax = fig.add_subplot(gs[r, c])
            sl = s["slice"]
            m = s["meta"]
            nx_s, nz_s = sl.shape
            ox, oz = m["origin"][0], m["origin"][2]
            dx, dz = m["spacing"][0], m["spacing"][2]
            extent = [oz, oz + (nz_s - 1) * dz, ox, ox + (nx_s - 1) * dx]
            im = ax.imshow(
                sl, origin="lower", cmap="RdBu_r", norm=norm,
                extent=extent, aspect="equal",
            )
            if r == 0:
                ax.set_title(f"t = {s['time']:.2f} a.u.", fontsize=10)
            if c == 0:
                ax.set_ylabel(f"{row_label}\nx (Bohr)")
            else:
                ax.set_yticklabels([])
            if r == 1:
                ax.set_xlabel("z (Bohr)")
            else:
                ax.set_xticklabels([])

    cax = fig.add_subplot(gs[:, -1])
    fig.colorbar(im, cax=cax, label=r"$\delta n$ (a.u.$^{-3}$)")
    fig.suptitle(
        rf"$\delta n(x, y{{=}}0, z)$ — free vs jellium — {run_name}",
        y=0.97, fontsize=12,
    )
    out_path = out_dir / "density_diff_vs_free.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[density_diff_vs_free] Wrote {out_path}")
    return out_path


# ===========================================================================
# CLI entry point
# ===========================================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Density & momentum analysis for jellium WP runs",
    )
    parser.add_argument("results_dir", type=Path, help="Run directory")
    parser.add_argument("--name", default=None, help="Run name for titles")
    parser.add_argument("--free-dir", default=None, type=Path,
                        help="Free-WP run directory for comparison plots")
    parser.add_argument(
        "--plots", nargs="*",
        default=["z_profile", "xz_snapshots", "momentum_band", "sigma_xyz", "diff_vs_free"],
        help="Which plots to produce",
    )
    args = parser.parse_args()

    run_name = args.name or args.results_dir.name

    if "z_profile" in args.plots:
        plot_density_z_profile_evolution(args.results_dir, run_name)
    if "xz_snapshots" in args.plots:
        plot_delta_density_xz_snapshots(args.results_dir, run_name)
    if "momentum_band" in args.plots:
        plot_momentum_band_free_vs_jellium(args.results_dir, run_name,
                                           free_results_dir=args.free_dir)
    if "sigma_xyz" in args.plots:
        plot_sigma_xyz_vs_time(args.results_dir, run_name)
    if "diff_vs_free" in args.plots:
        plot_density_diff_vs_free(args.results_dir, run_name,
                                  free_results_dir=args.free_dir)
