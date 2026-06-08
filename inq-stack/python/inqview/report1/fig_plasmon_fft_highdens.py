"""fig_plasmon_fft_highdens — Plasmon FFT for high-density jellium.

Fourier power spectrum of density Fourier modes n_q_m(t) for the
high-density run (N=162, L=30, r_s≈3.41).

Run:
    python -m inqview.report1.fig_plasmon_fft_highdens
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from inqview.report1._shared_style import (
    apply_style, palette_sweep5, column_widths_in, panel_label, TufteCritic,
)

RUN = "ResearchProject/systems/jellium/run_wp_n162_L30_E100_highdens_sigma1"
NQ_CSV = f"{RUN}/results/analysis/observables/n_q_vs_time.csv"
OUT = "docs/reports/report1/figures/fig_plasmon_fft_highdens.png"

HA_TO_EV = 27.2114
L_BOHR = 30.0
N_DENSITY = 162 / L_BOHR**3


def bohm_gross(m: int) -> float:
    q = 2 * np.pi * m / L_BOHR
    omega_p_ha = np.sqrt(4 * np.pi * N_DENSITY)
    v_F = (3 * np.pi**2 * N_DENSITY) ** (1 / 3)
    return np.sqrt(omega_p_ha**2 + 3 / 5 * v_F**2 * q**2) * HA_TO_EV


def main() -> None:
    apply_style()

    df = pd.read_csv(NQ_CSV)
    omega_p_eV = np.sqrt(4 * np.pi * N_DENSITY) * HA_TO_EV
    omega_max = 20.0

    colors = [palette_sweep5[0], palette_sweep5[2], palette_sweep5[4]]
    modes_to_plot = [1, 2, 3]

    W = column_widths_in["single"]
    fig, ax = plt.subplots(figsize=(W, W * 0.65))

    for i, m in enumerate(modes_to_plot):
        sub = df[df["m"] == m].sort_values("time_au")
        t = sub["time_au"].values
        nq = sub["re_n_q"].values + 1j * sub["im_n_q"].values
        dt = t[1] - t[0]
        N = len(t)

        window = np.hanning(N)
        fft_vals = np.fft.fft(nq * window)
        freqs_au = np.fft.fftfreq(N, d=dt)
        pos = freqs_au >= 0
        freqs_eV = freqs_au[pos] * 2 * np.pi * HA_TO_EV
        amp = np.abs(fft_vals[pos])

        mask = freqs_eV <= omega_max
        ax.semilogy(freqs_eV[mask], amp[mask], linewidth=0.8,
                     color=colors[i], label=f"$m={m}$", alpha=0.85)

        bg = bohm_gross(m)
        ax.axvline(bg, color=colors[i], ls=":", lw=0.5, alpha=0.5)

    ax.axvline(omega_p_eV, color="#808080", ls="--", lw=0.7, alpha=0.6)
    ax.text(omega_p_eV + 0.2, ax.get_ylim()[1] * 0.3,
            rf"$\omega_p = {omega_p_eV:.1f}$ eV", fontsize=5.5, color="#808080")

    bg1 = bohm_gross(1)
    ax.text(bg1 + 0.2, ax.get_ylim()[1] * 0.08,
            f"BG $m{{=}}1$: {bg1:.2f} eV", fontsize=5, color="grey")

    ax.set_xlabel(r"$\omega$ (eV)")
    ax.set_ylabel(r"$|\mathrm{FFT}[n_{q_m}(t)]|$")
    ax.set_xlim(0, omega_max)
    ax.legend(fontsize=6, loc="upper right", frameon=False)

    critic = TufteCritic()
    for iss in critic.critique(fig):
        print(f"  TufteCritic: {iss}")

    fig.savefig(OUT, dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"Saved -> {OUT}")
    plt.close(fig)


if __name__ == "__main__":
    main()
