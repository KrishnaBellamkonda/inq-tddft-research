"""fig_loss_function — Loss function L(q,ω) from density Fourier modes.

Two-panel figure:
  (a) 2D loss function heatmap L(q_z, ω) constructed from density Fourier
      modes n_q_m(t). Overlaid: Bohm-Gross plasmon dispersion, e-h
      continuum boundaries, Lindhard cutoff.
  (b) 1D cuts at selected q modes (m=1,2,3) with plasmon peaks and
      e-h transition energies marked.

The loss function is proportional to |FFT[δn_q(t)]|² / q², which
captures both collective (plasmon) and single-particle (e-h) excitations.

Data: run_plasmon_n162_L50_E15 (N=162, L=50, r_s≈5.69, 100k steps / 2000 a.u.).

Run:
    python -m applications.report1.fig_loss_function
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from applications.report1._shared_style import (
    apply_style,
    palette_sweep5,
    column_widths_in,
    panel_label,
    TufteCritic,
)

RUN = "ResearchProject/systems/jellium/run_plasmon_n162_L50_E15"
NQ_VS_TIME = f"{RUN}/results/analysis/observables/n_q_vs_time.csv"
EIGEN_CSV = f"{RUN}/results/raw/observables/eigenvalues/eigenvalues.csv"
OCC_CSV = f"{RUN}/results/raw/observables/eigenvalues/occupations.csv"

OUT = "docs/reports/report1/figures/fig_loss_function.png"
OUT_1COL = "docs/reports/report1/one-column-figures/fig_loss_function.png"

HA_TO_EV = 27.2114
L_BOHR = 50.0
N_EL = 162
N_DENSITY = N_EL / L_BOHR**3
AU_TO_FS = 0.02418884


def bohm_gross(q: np.ndarray) -> np.ndarray:
    """Bohm-Gross plasmon dispersion ω²=ω_p²+(3/5)v_F²q² in eV."""
    omega_p_ha = np.sqrt(4 * np.pi * N_DENSITY)
    v_F = (3 * np.pi**2 * N_DENSITY) ** (1 / 3)
    omega_ha = np.sqrt(omega_p_ha**2 + 3 / 5 * v_F**2 * q**2)
    return omega_ha * HA_TO_EV


def eh_continuum_boundaries(q: np.ndarray):
    """Lindhard e-h continuum: ω_± = q²/2 ± q·v_F (in eV)."""
    v_F = (3 * np.pi**2 * N_DENSITY) ** (1 / 3)
    omega_plus = (q**2 / 2 + q * v_F) * HA_TO_EV
    omega_minus = np.abs(q**2 / 2 - q * v_F) * HA_TO_EV
    return omega_minus, omega_plus


def main() -> None:
    apply_style()

    # Load time-domain density Fourier modes and FFT ourselves
    df_t = pd.read_csv(NQ_VS_TIME)
    modes = sorted(df_t["m"].unique())
    q_unit = 2 * np.pi / L_BOHR
    omega_max_eV = 25.0

    mode_data = {}
    for m in modes:
        sub = df_t[df_t["m"] == m].sort_values("time_au")
        t = sub["time_au"].values
        nq_complex = sub["re_n_q"].values + 1j * sub["im_n_q"].values
        dt = t[1] - t[0]
        N = len(t)
        q = m * q_unit

        # Apply Hann window to reduce spectral leakage
        window = np.hanning(N)
        nq_windowed = nq_complex * window

        fft_vals = np.fft.fft(nq_windowed)
        freqs_au = np.fft.fftfreq(N, d=dt)
        # Take positive frequencies only
        pos_mask = freqs_au >= 0
        fft_vals = fft_vals[pos_mask]
        freqs_eV = freqs_au[pos_mask] * 2 * np.pi * HA_TO_EV

        power = np.abs(fft_vals) ** 2
        loss = power / q**2

        mask = freqs_eV <= omega_max_eV
        mode_data[m] = {
            "omega": freqs_eV[mask],
            "loss": loss[mask],
            "q": q,
        }

    n_omega = len(mode_data[modes[0]]["omega"])
    omega_axis = mode_data[modes[0]]["omega"]
    q_axis = np.array([m * q_unit for m in modes])

    loss_2d = np.zeros((len(modes), n_omega))
    for i, m in enumerate(modes):
        loss_2d[i, :] = mode_data[m]["loss"]

    # Load eigenvalues for e-h transitions
    eig = pd.read_csv(EIGEN_CSV)
    occ = pd.read_csv(OCC_CSV)
    merged = eig.merge(occ, on="state_index")
    occ_mask = merged["occupation"] > 0.5
    e_occ = merged.loc[occ_mask, "eigenvalue_ev"].values
    e_unocc = merged.loc[~occ_mask, "eigenvalue_ev"].values

    # Compute e-h transition histogram
    eh_transitions = []
    for eo in e_occ:
        for eu in e_unocc:
            dE = eu - eo
            if 0 < dE < omega_max_eV:
                eh_transitions.append(dE)
    eh_transitions = np.array(eh_transitions)

    # Physics parameters
    omega_p_eV = np.sqrt(4 * np.pi * N_DENSITY) * HA_TO_EV
    v_F = (3 * np.pi**2 * N_DENSITY) ** (1 / 3)
    q_c = np.sqrt(4 * np.pi * N_DENSITY) / v_F

    print(f"ω_p = {omega_p_eV:.3f} eV")
    print(f"v_F = {v_F:.4f} a.u.")
    print(f"q_c (Lindhard cutoff) = {q_c:.3f} Bohr⁻¹ (m_c ≈ {q_c * L_BOHR / (2*np.pi):.1f})")
    print(f"E-h gap: {eh_transitions.min():.3f} eV")
    print(f"E-h max: {eh_transitions.max():.3f} eV")

    # ── Figure ──
    W = column_widths_in["full"]
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(W, W * 0.45),
                                      gridspec_kw={"width_ratios": [1.3, 1],
                                                   "wspace": 0.35})

    # ══ Panel (a): 2D loss function heatmap ══
    q_edges = np.concatenate([
        [q_axis[0] - 0.5 * q_unit],
        0.5 * (q_axis[:-1] + q_axis[1:]),
        [q_axis[-1] + 0.5 * q_unit]
    ])
    omega_edges = np.concatenate([
        [omega_axis[0]],
        0.5 * (omega_axis[:-1] + omega_axis[1:]),
        [omega_axis[-1]]
    ])

    vmax = np.percentile(loss_2d[loss_2d > 0], 99)
    vmin = max(vmax * 1e-4, loss_2d[loss_2d > 0].min())

    im = ax_a.pcolormesh(
        omega_edges, q_edges, loss_2d,
        cmap="inferno",
        norm=LogNorm(vmin=vmin, vmax=vmax),
        shading="flat",
        rasterized=True,
    )

    # Overlay: Bohm-Gross dispersion
    q_overlay = np.linspace(0.01, q_axis[-1] * 1.3, 200)
    omega_BG = bohm_gross(q_overlay)
    ax_a.plot(omega_BG, q_overlay, "--", color="white", linewidth=1.0,
              label="Bohm--Gross", alpha=0.9)

    # Overlay: e-h continuum boundaries
    omega_minus, omega_plus = eh_continuum_boundaries(q_overlay)
    ax_a.plot(omega_minus, q_overlay, ":", color="#66CCFF", linewidth=0.8,
              label=r"$e$-$h$ continuum", alpha=0.8)
    ax_a.plot(omega_plus, q_overlay, ":", color="#66CCFF", linewidth=0.8,
              alpha=0.8)

    # Lindhard cutoff
    ax_a.axhline(q_c, color="#FF6666", linewidth=0.6, linestyle="-.",
                 alpha=0.7, label=f"$q_c$ (Landau)")

    ax_a.set_xlabel(r"$\omega$ (eV)")
    ax_a.set_ylabel(r"$q_z$ (Bohr$^{-1}$)")
    ax_a.set_xlim(0, omega_max_eV)
    ax_a.set_ylim(0, q_axis[-1] * 1.15)
    ax_a.legend(fontsize=5.5, loc="upper right", frameon=True,
                facecolor="black", edgecolor="#404040",
                labelcolor="white", framealpha=0.7)

    cbar = fig.colorbar(im, ax=ax_a, fraction=0.046, pad=0.04)
    cbar.set_label(r"$L(q,\omega)$ (arb.\ units)", fontsize=7)
    cbar.ax.tick_params(labelsize=6)
    panel_label(ax_a, "(a)", x=0.03, y=0.95)

    # ══ Panel (b): 1D cuts at m=1,2,3 ══
    colors_mode = [palette_sweep5[0], palette_sweep5[2], palette_sweep5[4]]

    for i, m in enumerate([1, 2, 3]):
        if m not in mode_data:
            continue
        omega = mode_data[m]["omega"]
        loss = mode_data[m]["loss"]
        q = mode_data[m]["q"]

        # Normalize for visual comparison
        loss_norm = loss / loss.max() if loss.max() > 0 else loss

        ax_b.plot(omega, loss_norm, "-", color=colors_mode[i], linewidth=0.9,
                  label=rf"$m={m}$ ($q={q:.2f}$)", alpha=0.85)

        # Mark Bohm-Gross prediction
        bg = bohm_gross(np.array([q]))[0]
        ax_b.axvline(bg, color=colors_mode[i], linewidth=0.4,
                     linestyle=":", alpha=0.5)

    # Mark e-h transition energy range
    ax_b.axvspan(eh_transitions.min(), eh_transitions.max(),
                 alpha=0.08, color="#66CCFF", zorder=0,
                 label=r"$e$-$h$ transitions")

    # Mark plasmon frequency
    ax_b.axvline(omega_p_eV, color="#808080", linewidth=0.7, linestyle="--",
                 alpha=0.6)
    ax_b.text(omega_p_eV + 0.2, 0.92, r"$\omega_p$", fontsize=6,
              color="#808080", ha="left", va="top")

    ax_b.set_xlabel(r"$\omega$ (eV)")
    ax_b.set_ylabel(r"$L(q,\omega) / L_{\max}$")
    ax_b.set_xlim(0, omega_max_eV)
    ax_b.set_ylim(0, 1.05)
    ax_b.legend(fontsize=5.5, loc="upper right", frameon=False)
    panel_label(ax_b, "(b)", x=0.03, y=0.95)

    # ── Save ──
    critic = TufteCritic()
    issues = critic.critique(fig)
    if issues:
        print(f"TufteCritic: {len(issues)} issue(s)")
        for iss in issues:
            print(f"  {iss}")

    fig.savefig(OUT, dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"Saved -> {OUT}")

    # One-column version (same figure, narrower)
    import os
    os.makedirs("docs/reports/report1/one-column-figures", exist_ok=True)
    fig.savefig(OUT_1COL, dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"Saved -> {OUT_1COL}")

    plt.close(fig)


if __name__ == "__main__":
    main()
