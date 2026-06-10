"""fig_leed_validation — LEED backscatter validation against analytics.

Four-panel comparison for the center-target backscattering screen:
  (a) Simulation: |FFT[screen]|² (k-space LEED from screen 14, step 330)
  (b) Analytical structure factor |F(q)|² from atom positions
  (c) |FFT[n_GS(x,y)]|² from ground-state density
  (d) Azimuthal intensity I(θ) comparison at the first Bragg ring

Quantitative metrics: peak position match, D6h symmetry check,
relative intensity comparison.

Run:
    python -m inqview.report1.fig_leed_validation
"""

from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from pathlib import Path

from inqview.io.leed import load_leed_pattern
from inqview.report1._shared_style import (
    apply_style,
    column_widths_in,
    panel_label,
    TufteCritic,
    palette_sweep5,
)

ANG_TO_BOHR = 1.8897259886

SCREEN_IDX = 14
STEP = 330

CENTER_SNAP = Path(
    "ResearchProject/systems/coronene/run_propagate_paper_replica/"
    "results/screens_snapshots"
)
CENTER_REF = Path(
    "ResearchProject/systems/coronene/run_propagate_paper_replica/"
    f"results/screens/screen_{SCREEN_IDX:02d}.dat"
)
GS_VTI = Path(
    "ResearchProject/systems/coronene/run_save_gs_paper_replica/"
    "results/density_gs/density_t000000.vti"
)
XYZ_FILE = Path(
    "ResearchProject/systems/coronene/configurations/"
    "tsubonoya_2014_paper_replica/coronene_centred.xyz"
)
OUT = Path("docs/reports/report1/figures/fig_leed_validation.png")


def fix_snapshot_dx(pattern, ref_screen_path: Path):
    ref = load_leed_pattern(ref_screen_path)
    pattern.dx_bohr = ref.dx_bohr
    pattern.dy_bohr = ref.dy_bohr
    pattern.origin_x_bohr = -0.5 * pattern.nx * ref.dx_bohr
    pattern.origin_y_bohr = -0.5 * pattern.ny * ref.dy_bohr
    return pattern


def parse_xyz(path: Path):
    lines = path.read_text().strip().split("\n")
    n_atoms = int(lines[0])
    species, positions_bohr = [], []
    for line in lines[2:2 + n_atoms]:
        parts = line.split()
        species.append(parts[0])
        positions_bohr.append([float(parts[i]) * ANG_TO_BOHR for i in (1, 2, 3)])
    return species, np.array(positions_bohr)


def structure_factor_2d(species, pos_bohr, kx, ky):
    """Compute |F(q)|² on a 2D k-grid.

    F(q) = Σ_j f_j exp(i q·R_j) where f_C=6, f_H=1 (Born approx).
    """
    f_map = {"C": 6.0, "H": 1.0}
    KX, KY = np.meshgrid(kx, ky)
    F = np.zeros_like(KX, dtype=complex)
    for s, r in zip(species, pos_bohr):
        f_j = f_map.get(s, 1.0)
        F += f_j * np.exp(1j * (KX * r[0] + KY * r[1]))
    return np.abs(F) ** 2


def load_gs_density_xy(vti_path: Path):
    import vtk
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(str(vti_path))
    reader.Update()
    data = reader.GetOutput()
    dims = data.GetDimensions()
    spacing = data.GetSpacing()
    arr = data.GetPointData().GetArray(0)
    flat = np.array([arr.GetValue(i) for i in range(arr.GetNumberOfTuples())])
    rho = flat.reshape(dims[2], dims[1], dims[0]).transpose(2, 1, 0)
    n2d = rho.sum(axis=2) * spacing[2]
    return n2d, spacing[0], spacing[1]


def screen_to_kspace(pattern):
    """FFT screen real-space data to k-space LEED pattern."""
    nx = pattern.nx
    fft = np.fft.fftshift(np.fft.fft2(pattern.data))
    power = np.abs(fft) ** 2
    dk = 2 * np.pi / (nx * pattern.dx_bohr)
    k = (np.arange(nx) - nx // 2) * dk
    return power, k


def suppress_dc(power, margin=2):
    p = power.copy()
    cx = p.shape[0] // 2
    p[cx - margin:cx + margin + 1, cx - margin:cx + margin + 1] = 0
    return p


def find_peaks_2d(data, k, threshold_frac=0.1):
    """Find peaks above threshold in 2D k-space data."""
    from scipy.ndimage import maximum_filter, label
    thresh = threshold_frac * data.max()
    local_max = maximum_filter(data, size=5)
    peaks = (data == local_max) & (data > thresh)
    labeled, n_features = label(peaks)
    peak_list = []
    for i in range(1, n_features + 1):
        ys, xs = np.where(labeled == i)
        # Intensity-weighted centroid
        intensities = data[ys, xs]
        cx = np.average(xs, weights=intensities)
        cy = np.average(ys, weights=intensities)
        kx_peak = np.interp(cx, np.arange(len(k)), k)
        ky_peak = np.interp(cy, np.arange(len(k)), k)
        total_intensity = intensities.sum()
        peak_list.append((kx_peak, ky_peak, total_intensity))
    return peak_list


def check_d6h_symmetry(peaks, tol_deg=5.0):
    """Check if peak arrangement has 6-fold rotational symmetry."""
    if len(peaks) < 6:
        return False, "fewer than 6 peaks"
    angles = np.array([np.degrees(np.arctan2(ky, kx)) for kx, ky, _ in peaks])
    angles = np.sort(angles % 360)
    diffs = np.diff(angles)
    if len(diffs) < 5:
        return False, f"only {len(diffs)+1} peaks"
    # For D6h, expect angular spacings near 30° or 60°
    near_60 = np.sum(np.abs(diffs - 60) < tol_deg)
    near_30 = np.sum(np.abs(diffs - 30) < tol_deg)
    if near_60 >= 5 or near_30 >= 10:
        return True, f"{near_60} gaps near 60°, {near_30} near 30°"
    return False, f"gaps: {diffs.round(1)}"


def azimuthal_profile(data, k, k_ring, dk_width=0.3):
    """Extract azimuthal intensity I(θ) at a specific |k| ring."""
    KX, KY = np.meshgrid(k, k)
    K_mag = np.sqrt(KX**2 + KY**2)
    theta = np.arctan2(KY, KX)

    mask = np.abs(K_mag - k_ring) < dk_width
    if mask.sum() == 0:
        return np.array([]), np.array([])

    theta_bins = np.linspace(-np.pi, np.pi, 73)
    theta_centers = 0.5 * (theta_bins[:-1] + theta_bins[1:])
    I_theta = np.zeros(len(theta_centers))

    for i in range(len(theta_centers)):
        ring_and_wedge = mask & (theta >= theta_bins[i]) & (theta < theta_bins[i + 1])
        if ring_and_wedge.sum() > 0:
            I_theta[i] = data[ring_and_wedge].mean()

    return np.degrees(theta_centers), I_theta


def main() -> None:
    apply_style()

    # === Load simulation screen and FFT ===
    print("Loading backscatter screen (center target)...")
    center_path = CENTER_SNAP / f"step_{STEP:06d}" / f"screen_{SCREEN_IDX:02d}.dat"
    center_pat = fix_snapshot_dx(load_leed_pattern(center_path), CENTER_REF)
    sim_power, sim_k = screen_to_kspace(center_pat)
    sim_nodc = suppress_dc(sim_power)
    print(f"  Screen {SCREEN_IDX}, step {STEP}, z={center_pat.z_bohr:.2f}, "
          f"nx={center_pat.nx}, dx={center_pat.dx_bohr:.4f}")

    # === Analytical structure factor ===
    print("\nComputing analytical structure factor...")
    species, pos_bohr = parse_xyz(XYZ_FILE)
    n_C = sum(1 for s in species if s == "C")
    n_H = sum(1 for s in species if s == "H")
    print(f"  {n_C} C + {n_H} H atoms")
    sf = structure_factor_2d(species, pos_bohr, sim_k, sim_k)
    sf_nodc = suppress_dc(sf)

    # === GS density FFT ===
    print("\nComputing GS density FFT...")
    n2d, gs_dx, gs_dy = load_gs_density_xy(GS_VTI)
    nx_gs = n2d.shape[0]
    fft_gs = np.fft.fftshift(np.fft.fft2(n2d))
    gs_power = np.abs(fft_gs) ** 2
    dk_gs = 2 * np.pi / (nx_gs * gs_dx)
    k_gs = (np.arange(nx_gs) - nx_gs // 2) * dk_gs
    gs_nodc = suppress_dc(gs_power)

    # === Find peaks and compare ===
    print("\n" + "=" * 60)
    print("QUANTITATIVE COMPARISON")
    print("=" * 60)

    sim_peaks = find_peaks_2d(sim_nodc, sim_k, threshold_frac=0.05)
    sf_peaks = find_peaks_2d(sf_nodc, sim_k, threshold_frac=0.05)
    gs_peaks = find_peaks_2d(gs_nodc, k_gs, threshold_frac=0.05)

    print(f"\nPeak count: simulation={len(sim_peaks)}, "
          f"structure factor={len(sf_peaks)}, GS FFT={len(gs_peaks)}")

    # D6h symmetry check
    for label_str, peaks in [("Simulation", sim_peaks),
                              ("Structure factor", sf_peaks),
                              ("GS FFT", gs_peaks)]:
        ok, msg = check_d6h_symmetry(peaks)
        status = "✓ D6h" if ok else "✗ NOT D6h"
        print(f"  {label_str}: {status} — {msg}")

    # Bragg peak positions
    a_cc = 2.68  # Bohr (1.42 Å)
    k_inner = (2 * np.pi / a_cc) / np.sqrt(3)
    k_outer = 2 * (2 * np.pi / a_cc) / np.sqrt(3)
    print(f"\nExpected Bragg peaks: k_inner={k_inner:.3f}, k_outer={k_outer:.3f} Bohr⁻¹")

    for label_str, peaks in [("Simulation", sim_peaks),
                              ("Structure factor", sf_peaks)]:
        if not peaks:
            continue
        k_mags = [np.sqrt(kx**2 + ky**2) for kx, ky, _ in peaks]
        # Group into inner/outer rings
        inner = [k for k in k_mags if abs(k - k_inner) < 0.5]
        outer = [k for k in k_mags if abs(k - k_outer) < 0.5]
        if inner:
            print(f"  {label_str} inner ring: |k|={np.mean(inner):.3f} ± {np.std(inner):.3f} "
                  f"(expected {k_inner:.3f}, dev {abs(np.mean(inner)-k_inner)/k_inner*100:.1f}%)")
        if outer:
            print(f"  {label_str} outer ring: |k|={np.mean(outer):.3f} ± {np.std(outer):.3f} "
                  f"(expected {k_outer:.3f}, dev {abs(np.mean(outer)-k_outer)/k_outer*100:.1f}%)")

    # Azimuthal profile at first Bragg ring
    theta_sim, I_sim = azimuthal_profile(sim_nodc, sim_k, k_inner)
    theta_sf, I_sf = azimuthal_profile(sf_nodc, sim_k, k_inner)
    theta_gs, I_gs = azimuthal_profile(gs_nodc, k_gs, k_inner)

    # Normalize for comparison
    if I_sim.max() > 0:
        I_sim_norm = I_sim / I_sim.max()
    else:
        I_sim_norm = I_sim
    if I_sf.max() > 0:
        I_sf_norm = I_sf / I_sf.max()
    else:
        I_sf_norm = I_sf
    if I_gs.max() > 0:
        I_gs_norm = I_gs / I_gs.max()
    else:
        I_gs_norm = I_gs

    # Cross-correlation as similarity metric
    if len(I_sim_norm) > 0 and len(I_sf_norm) > 0:
        corr_sim_sf = np.corrcoef(I_sim_norm, I_sf_norm)[0, 1]
        print(f"\nAzimuthal cross-correlation (inner ring, |k|≈{k_inner:.2f}):")
        print(f"  Simulation vs Structure factor: r = {corr_sim_sf:.4f}")
    if len(I_sim_norm) > 0 and len(I_gs_norm) > 0 and len(I_gs_norm) == len(I_sim_norm):
        corr_sim_gs = np.corrcoef(I_sim_norm, I_gs_norm)[0, 1]
        print(f"  Simulation vs GS FFT:           r = {corr_sim_gs:.4f}")

    # === Figure: 2x2 ===
    W = column_widths_in["full"]
    fig, axes = plt.subplots(2, 2, figsize=(W, W * 0.95),
                              gridspec_kw={"wspace": 0.30, "hspace": 0.30})

    k_lim = 5.5

    # (a) Simulation FFT
    vmax_sim = sim_nodc.max()
    axes[0, 0].imshow(sim_nodc, origin="lower",
                       extent=[sim_k[0], sim_k[-1], sim_k[0], sim_k[-1]],
                       aspect="equal", cmap="inferno", vmin=0, vmax=vmax_sim,
                       interpolation="nearest")
    axes[0, 0].set_xlim(-k_lim, k_lim)
    axes[0, 0].set_ylim(-k_lim, k_lim)
    axes[0, 0].set_xlabel(r"$k_x$ (Bohr$^{-1}$)")
    axes[0, 0].set_ylabel(r"$k_y$ (Bohr$^{-1}$)")
    panel_label(axes[0, 0], "(a)", x=0.03, y=0.97)

    # (b) Analytical structure factor
    vmax_sf = sf_nodc.max()
    axes[0, 1].imshow(sf_nodc, origin="lower",
                       extent=[sim_k[0], sim_k[-1], sim_k[0], sim_k[-1]],
                       aspect="equal", cmap="inferno", vmin=0, vmax=vmax_sf,
                       interpolation="nearest")
    axes[0, 1].set_xlim(-k_lim, k_lim)
    axes[0, 1].set_ylim(-k_lim, k_lim)
    axes[0, 1].set_xlabel(r"$k_x$ (Bohr$^{-1}$)")
    axes[0, 1].set_ylabel(r"$k_y$ (Bohr$^{-1}$)")
    panel_label(axes[0, 1], "(b)", x=0.03, y=0.97)

    # (c) GS density FFT
    vmax_gs = gs_nodc.max()
    axes[1, 0].imshow(gs_nodc, origin="lower",
                       extent=[k_gs[0], k_gs[-1], k_gs[0], k_gs[-1]],
                       aspect="equal", cmap="inferno", vmin=0, vmax=vmax_gs,
                       interpolation="nearest")
    axes[1, 0].set_xlim(-k_lim, k_lim)
    axes[1, 0].set_ylim(-k_lim, k_lim)
    axes[1, 0].set_xlabel(r"$k_x$ (Bohr$^{-1}$)")
    axes[1, 0].set_ylabel(r"$k_y$ (Bohr$^{-1}$)")
    panel_label(axes[1, 0], "(c)", x=0.03, y=0.97)

    # (d) Azimuthal intensity I(θ) at inner Bragg ring
    if len(theta_sim) > 0:
        axes[1, 1].plot(theta_sim, I_sim_norm, "-", color=palette_sweep5[0],
                         linewidth=1.0, label="Simulation")
    if len(theta_sf) > 0:
        axes[1, 1].plot(theta_sf, I_sf_norm, "--", color=palette_sweep5[2],
                         linewidth=1.0, label=r"$|F(\mathbf{q})|^2$")
    if len(theta_gs) > 0:
        axes[1, 1].plot(theta_gs, I_gs_norm, ":", color=palette_sweep5[4],
                         linewidth=1.0, label=r"$|\tilde{n}_0|^2$")
    axes[1, 1].set_xlabel(r"$\theta$ (deg)")
    axes[1, 1].set_ylabel(r"$I(\theta) / I_{\max}$")
    axes[1, 1].set_xlim(-180, 180)
    axes[1, 1].set_ylim(bottom=0)
    axes[1, 1].legend(fontsize=6, loc="upper right", frameon=False)
    panel_label(axes[1, 1], "(d)", x=0.03, y=0.97)

    critic = TufteCritic()
    issues = critic.critique(fig)
    if issues:
        print(f"\nTufteCritic: {len(issues)} issue(s)")
        for iss in issues:
            print(f"  {iss}")

    fig.savefig(str(OUT), dpi=600, bbox_inches="tight", pad_inches=0.03)
    print(f"\nSaved -> {OUT}")
    plt.close(fig)


if __name__ == "__main__":
    main()
