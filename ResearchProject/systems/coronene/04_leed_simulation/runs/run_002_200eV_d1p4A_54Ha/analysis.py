#!/usr/bin/env python3
"""
analysis.py — run_002 Coronene TDDFT LEED post-processing
==========================================================

Reads all output files produced by run.cpp and generates publication-quality
figures in each results/ subfolder.

Outputs
-------
results/energy/total_energy.png         — total energy vs time
results/energy/kinetic_energy.png       — kinetic energy vs time
results/energy/all_energies.png         — all energy components vs time
results/energy/energy_fft.png           — FFT of E_total fluctuations
results/ks_overlaps/projected_occ_heatmap.png  — KS occupation heatmap
results/wp_orbital/wp_orbital_snapshots.png    — WP orbital density at snapshots
results/wp_trajectory/density_trajectory.png   — z-profile vs time (image)
results/density_snapshots/fig1_density_snapshots.png  — total density grid
results/leed_pattern/fig2_leed_pattern.png     — LEED I(x,y)

Usage
-----
  python3 analysis.py              (from run_002 directory)
  python3 analysis.py --no-show   (non-interactive: skip plt.show calls)

Reference: Tsubonoya, Hu, Watanabe, PRB 90, 035416 (2014)
"""

import sys
import os
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")  # non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import LogNorm

# ── Unit conversions ──────────────────────────────────────────────────────────
AU_TO_FS = 0.024188843   # 1 a.u. = 0.02419 fs
HA_TO_EV = 27.21138625

RESULTS = "results"

def savefig(path, **kw):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    plt.savefig(path, dpi=150, bbox_inches="tight", **kw)
    plt.close()
    print(f"  Saved: {path}")


# ════════════════════════════════════════════════════════════════════════════════
# 1.  Energy vs time
# ════════════════════════════════════════════════════════════════════════════════

def load_energy(path="results/energy/energy_vs_time.csv"):
    if not os.path.exists(path):
        print(f"  [SKIP] {path} not found")
        return None
    data = np.genfromtxt(path, delimiter=",", comments="#",
                         names=["step","t_au","E_total","E_kinetic",
                                "E_hartree","E_xc","E_external",
                                "E_nonlocal","E_ion"])
    return data


def plot_total_energy(data):
    if data is None:
        return
    t_fs = data["t_au"] * AU_TO_FS
    E    = data["E_total"]
    drift = (E[-1] - E[0]) * 1e6  # μHa

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(t_fs, E, lw=0.8, color="steelblue")
    ax.set_xlabel("Time (fs)")
    ax.set_ylabel("Total energy (Ha)")
    ax.set_title(f"Total energy vs time — run_002\n"
                 f"drift = {drift:+.1f} μHa over {len(t_fs)} steps")
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    fig.tight_layout()
    savefig("results/energy/total_energy.png")


def plot_kinetic_energy(data):
    if data is None:
        return
    t_fs = data["t_au"] * AU_TO_FS

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(t_fs, data["E_kinetic"], lw=0.8, color="darkorange")
    ax.set_xlabel("Time (fs)")
    ax.set_ylabel("Kinetic energy (Ha)")
    ax.set_title("Kinetic energy vs time — run_002")
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    fig.tight_layout()
    savefig("results/energy/kinetic_energy.png")


def plot_all_energies(data):
    if data is None:
        return
    t_fs = data["t_au"] * AU_TO_FS

    components = [
        ("E_kinetic",  "Kinetic",   "C0"),
        ("E_hartree",  "Hartree",   "C1"),
        ("E_xc",       "XC",        "C2"),
        ("E_external", "External",  "C3"),
        ("E_nonlocal", "Non-local", "C4"),
        ("E_ion",      "Ion-ion",   "C5"),
    ]

    fig, axes = plt.subplots(3, 2, figsize=(12, 9), sharex=True)
    axes = axes.flatten()

    for ax, (col, label, color) in zip(axes, components):
        ax.plot(t_fs, data[col], lw=0.8, color=color)
        ax.set_ylabel(f"{label} (Ha)", fontsize=9)
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.set_title(label, fontsize=10)

    for ax in axes[-2:]:
        ax.set_xlabel("Time (fs)")

    fig.suptitle("Energy components vs time — run_002", fontsize=12)
    fig.tight_layout()
    savefig("results/energy/all_energies.png")


def plot_energy_fft(data):
    if data is None:
        return
    t_au = data["t_au"]
    E    = data["E_total"]
    dt   = t_au[1] - t_au[0]

    # Detrend: subtract linear fit to isolate fluctuations
    coeffs = np.polyfit(t_au, E, 1)
    E_fluct = E - np.polyval(coeffs, t_au)

    N    = len(E_fluct)
    freqs = np.fft.rfftfreq(N, d=dt)   # a.u.^{-1}
    amps  = np.abs(np.fft.rfft(E_fluct)) / N
    # Convert frequency to eV (1 a.u. energy = HA_TO_EV eV; time in a.u.)
    # omega [a.u.] * HA_TO_EV = omega [eV] (since hbar=1 a.u.)
    omega_eV = freqs * HA_TO_EV * (2 * np.pi)  # angular -> linear frequency in eV

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.semilogy(omega_eV[1:], amps[1:], lw=0.7, color="purple")
    ax.set_xlabel("Frequency (eV)")
    ax.set_ylabel("|FFT| (Ha)")
    ax.set_title("FFT of total energy fluctuations — run_002")
    ax.set_xlim(0, min(30, omega_eV[-1]))
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    fig.tight_layout()
    savefig("results/energy/energy_fft.png")


# ════════════════════════════════════════════════════════════════════════════════
# 2.  Projected KS occupation heatmap
# ════════════════════════════════════════════════════════════════════════════════

def plot_ks_overlap_heatmap(path="results/ks_overlaps/projected_occ_vs_time.csv"):
    if not os.path.exists(path):
        print(f"  [SKIP] {path} not found")
        return

    raw = np.genfromtxt(path, delimiter=",", comments="#")
    if raw.ndim < 2 or raw.shape[1] < 4:
        print(f"  [SKIP] {path}: not enough columns")
        return

    steps = raw[:, 0].astype(int)
    t_au  = raw[:, 1]
    occ   = raw[:, 2:]   # shape (n_times, n_states)
    t_fs  = t_au * AU_TO_FS
    n_states = occ.shape[1]

    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.imshow(occ.T, origin="lower", aspect="auto",
                   extent=[t_fs[0], t_fs[-1], 0, n_states],
                   cmap="inferno", vmin=0.0, vmax=1.0)
    cbar = fig.colorbar(im, ax=ax, label="Projected occupation")
    ax.set_xlabel("Time (fs)")
    ax.set_ylabel("KS state index")
    ax.set_title("Projected KS occupation vs time — run_002\n"
                 r"$\sum_j f_j |\langle\phi_i^{GS}|\phi_j(t)\rangle|^2$")
    # Mark GS occupied states (0–53) and WP state (56)
    ax.axhline(53.5, color="cyan",  lw=0.8, ls="--", label="occupied / extra boundary")
    ax.axhline(55.5, color="yellow",lw=0.8, ls="--", label="WP state (last)")
    ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    savefig("results/ks_overlaps/projected_occ_heatmap.png")


# ════════════════════════════════════════════════════════════════════════════════
# 3.  WP orbital density snapshots
# ════════════════════════════════════════════════════════════════════════════════

def load_slice(path):
    """Load a 2D density slice text file. Returns (time_au, z_bohr, array2d)."""
    with open(path) as f:
        header = f.readline()  # # t=... z=...
    parts = header.strip("# \n").split()
    t_au   = float(parts[0].split("=")[1])
    z_bohr = float(parts[1].split("=")[1])
    arr = np.loadtxt(path, comments="#")
    return t_au, z_bohr, arr


def plot_wp_orbital_snapshots():
    files = sorted(glob.glob("results/wp_orbital/wp_slice_t*.txt"))
    if not files:
        print("  [SKIP] no WP orbital slice files found")
        return

    n = len(files)
    ncols = min(5, n)
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(3.5*ncols, 3.0*nrows))
    if n == 1:
        axes = np.array([[axes]])
    elif nrows == 1:
        axes = axes[np.newaxis, :]
    axes_flat = axes.flatten()

    vmax_global = 0.0
    slices = []
    for fp in files:
        t_au, z_bohr, arr = load_slice(fp)
        slices.append((t_au, z_bohr, arr))
        vmax_global = max(vmax_global, arr.max())

    for i, (t_au, z_bohr, arr) in enumerate(slices):
        ax = axes_flat[i]
        t_fs = t_au * AU_TO_FS
        im = ax.imshow(arr, origin="lower", aspect="equal",
                       cmap="hot",
                       vmin=0, vmax=vmax_global if vmax_global > 0 else 1)
        ax.set_title(f"t={t_fs:.3f} fs", fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])

    for j in range(i+1, len(axes_flat)):
        axes_flat[j].set_visible(False)

    fig.suptitle("|ψ_WP|² at z=z_flake — run_002  (WP orbital only)", fontsize=11)
    fig.colorbar(im, ax=axes.ravel().tolist(), label=r"$|\psi_{WP}|^2$ (bohr$^{-3}$)",
                 shrink=0.6)
    savefig("results/wp_orbital/wp_orbital_snapshots.png")


# ════════════════════════════════════════════════════════════════════════════════
# 4.  Z-trajectory density plot
# ════════════════════════════════════════════════════════════════════════════════

def plot_z_trajectory(path="results/wp_trajectory/density_z_profile_vs_time.csv"):
    if not os.path.exists(path):
        print(f"  [SKIP] {path} not found")
        return

    # File has 2 header lines then: step, t_au, n(z_0), n(z_1), ...
    raw = np.genfromtxt(path, delimiter=",", comments="#")
    if raw.ndim < 2 or raw.shape[1] < 4:
        print(f"  [SKIP] {path}: not enough columns")
        return

    t_au   = raw[:, 1]
    nz_arr = raw[:, 2:]   # shape (n_snaps, Nz)
    t_fs   = t_au * AU_TO_FS
    Nz     = nz_arr.shape[1]

    # Z axis in bohr — we don't know dz exactly; use index as proxy
    # For labelling: Lz_bohr ~ 59.9, so dz = 59.9 / Nz
    LZ_BOHR = 59.904
    z_bohr  = np.linspace(0, LZ_BOHR, Nz, endpoint=False)
    z_flake = LZ_BOHR / 2.0      # 29.95 bohr

    vmax = np.percentile(nz_arr, 99)

    fig, ax = plt.subplots(figsize=(10, 5))
    im = ax.imshow(nz_arr.T, origin="lower", aspect="auto",
                   extent=[t_fs[0], t_fs[-1], z_bohr[0], z_bohr[-1]],
                   cmap="inferno",
                   vmin=0, vmax=max(vmax, 1e-10))
    ax.axhline(z_flake, color="cyan", lw=0.8, ls="--", label=f"z_flake={z_flake:.1f} bohr")
    ax.set_xlabel("Time (fs)")
    ax.set_ylabel("z (bohr)")
    ax.set_title("Electron density at cell centre (x=Lx/2, y=Ly/2) along z — run_002")
    ax.legend(fontsize=8, loc="upper right")
    fig.colorbar(im, ax=ax, label=r"$n(x_c,y_c,z,t)$ (bohr$^{-3}$)")
    fig.tight_layout()
    savefig("results/wp_trajectory/density_trajectory.png")


# ════════════════════════════════════════════════════════════════════════════════
# 5.  Total density snapshots (Fig. 1 analogue)
# ════════════════════════════════════════════════════════════════════════════════

def plot_density_snapshots():
    files = sorted(glob.glob("results/density_snapshots/snapshot_t*.txt"))
    if not files:
        print("  [SKIP] no density snapshot files found")
        return

    n = len(files)
    ncols = min(5, n)
    nrows = (n + ncols - 1) // ncols

    slices = [load_slice(fp) for fp in files]
    vmax_global = max(arr.max() for _, _, arr in slices)

    fig, axes = plt.subplots(nrows, ncols, figsize=(3.5*ncols, 3.0*nrows))
    if n == 1:
        axes = np.array([[axes]])
    elif nrows == 1:
        axes = axes[np.newaxis, :]
    axes_flat = axes.flatten()

    for i, (t_au, z_bohr, arr) in enumerate(slices):
        ax = axes_flat[i]
        t_fs = t_au * AU_TO_FS
        im = ax.imshow(arr, origin="lower", aspect="equal",
                       cmap="inferno",
                       vmin=0, vmax=vmax_global if vmax_global > 0 else 1)
        ax.set_title(f"t={t_fs:.3f} fs", fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])

    for j in range(i+1, len(axes_flat)):
        axes_flat[j].set_visible(False)

    fig.suptitle("Total electron density n(x,y,z=z_flake,t) — run_002\n"
                 "(all KS orbitals, coronene plane)", fontsize=11)
    fig.colorbar(im, ax=axes.ravel().tolist(),
                 label=r"$n$ (bohr$^{-3}$)", shrink=0.6)
    savefig("results/density_snapshots/fig1_density_snapshots.png")


# ════════════════════════════════════════════════════════════════════════════════
# 6.  LEED pattern (Fig. 2 analogue)
# ════════════════════════════════════════════════════════════════════════════════

def plot_leed_pattern(path="results/leed_pattern/leed_pattern.txt"):
    if not os.path.exists(path):
        print(f"  [SKIP] {path} not found")
        return

    leed = np.loadtxt(path, comments="#")
    if leed.ndim != 2:
        print(f"  [SKIP] {path}: unexpected shape {leed.shape}")
        return

    Ny, Nx = leed.shape
    # Cell: Lx = Ly = 18.4 Å
    LX_ANG = 18.4
    LY_ANG = 18.4
    x = np.linspace(-LX_ANG/2, LX_ANG/2, Nx, endpoint=False)
    y = np.linspace(-LY_ANG/2, LY_ANG/2, Ny, endpoint=False)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Linear scale
    ax = axes[0]
    im0 = ax.pcolormesh(x, y, leed, cmap="hot", shading="auto")
    ax.set_aspect("equal")
    ax.set_xlabel("x (Å)")
    ax.set_ylabel("y (Å)")
    ax.set_title("LEED pattern I(x,y) — linear scale")
    fig.colorbar(im0, ax=ax, label=r"$\int n\,dt$ (bohr$^{-3}\cdot$a.u.)")

    # Log scale (clip zeros)
    ax = axes[1]
    leed_pos = np.clip(leed, leed[leed > 0].min() if (leed > 0).any() else 1e-30, None)
    im1 = ax.pcolormesh(x, y, leed_pos, cmap="hot", shading="auto",
                        norm=LogNorm(vmin=leed_pos.min(), vmax=leed_pos.max()))
    ax.set_aspect("equal")
    ax.set_xlabel("x (Å)")
    ax.set_ylabel("y (Å)")
    ax.set_title("LEED pattern I(x,y) — log scale")
    fig.colorbar(im1, ax=ax, label=r"$\int n\,dt$ (log scale)")

    fig.suptitle("Coronene LEED, 200 eV — run_002\n"
                 "Tsubonoya, Hu, Watanabe PRB 90, 035416 (2014) Fig. 2", fontsize=11)
    fig.tight_layout()
    savefig("results/leed_pattern/fig2_leed_pattern.png")

    # Also save FFT of LEED (reciprocal space)
    leed_fft = np.abs(np.fft.fftshift(np.fft.fft2(leed))) ** 2
    fig2, ax2 = plt.subplots(figsize=(6, 5))
    ax2.imshow(np.log1p(leed_fft), origin="lower", aspect="equal", cmap="inferno")
    ax2.set_title("FFT² of LEED pattern (reciprocal space)")
    ax2.set_xlabel("kx (arb.)")
    ax2.set_ylabel("ky (arb.)")
    savefig("results/leed_pattern/fig2_leed_reciprocal.png")


# ════════════════════════════════════════════════════════════════════════════════
# Main
# ════════════════════════════════════════════════════════════════════════════════

def main():
    print("=== run_002 analysis.py ===\n")
    print("Working directory:", os.getcwd())

    # --- Energy section ---
    print("\n[1] Energy vs time")
    energy = load_energy()
    plot_total_energy(energy)
    plot_kinetic_energy(energy)
    plot_all_energies(energy)
    plot_energy_fft(energy)

    # --- KS overlaps ---
    print("\n[2] Projected KS occupation")
    plot_ks_overlap_heatmap()

    # --- WP orbital snapshots ---
    print("\n[3] WP orbital density snapshots")
    plot_wp_orbital_snapshots()

    # --- Z-trajectory ---
    print("\n[4] Z-trajectory density")
    plot_z_trajectory()

    # --- Total density snapshots ---
    print("\n[5] Total density snapshots (Fig. 1)")
    plot_density_snapshots()

    # --- LEED pattern ---
    print("\n[6] LEED pattern (Fig. 2)")
    plot_leed_pattern()

    print("\n=== Done. All figures saved to results/ subfolders. ===")


if __name__ == "__main__":
    main()
