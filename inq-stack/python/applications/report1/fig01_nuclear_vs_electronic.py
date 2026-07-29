"""fig01 — Nuclear vs electronic stopping crossover.

SRIM-calibrated S_n and S_e for proton-on-carbon using ZBL analytical
forms. Establishes why only electronic stopping matters for electron
projectiles.

Run:
    python -m applications.report1.fig01_nuclear_vs_electronic
"""

from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from applications.report1._shared_style import (
    apply_style,
    palette_sweep5,
    column_widths_in,
    TufteCritic,
)


def nuclear_stopping(E_keV_per_u: np.ndarray) -> np.ndarray:
    E = np.maximum(E_keV_per_u, 1e-4)
    eps = E * 1e3 / 502.0
    sn_red = 0.5 * np.log(1 + 1.1383 * eps) / (
        eps + 0.01321 * eps**0.21226 + 0.19593 * np.sqrt(eps)
    )
    return 24.0 * sn_red


def electronic_stopping(E_keV_per_u: np.ndarray) -> np.ndarray:
    E = np.maximum(E_keV_per_u, 1e-4)
    S_low = 1.4 * np.sqrt(E)
    I_keV = 0.078
    arg = np.maximum(2.0 * E / I_keV, 1.01)
    S_high = 350.0 * np.log(arg) / E
    S_high = np.maximum(S_high, 0.01)
    return 1.0 / (1.0 / S_low + 1.0 / S_high)


def find_crossover(E, Sn, Se):
    diff = Sn - Se
    sign_changes = np.where(np.diff(np.sign(diff)))[0]
    if len(sign_changes) == 0:
        return np.nan
    idx = sign_changes[0]
    x0, x1 = E[idx], E[idx + 1]
    d0, d1 = diff[idx], diff[idx + 1]
    return x0 - d0 * (x1 - x0) / (d1 - d0)


def main() -> None:
    apply_style()

    E = np.geomspace(0.1, 1e4, 2000)
    Sn = nuclear_stopping(E)
    Se = electronic_stopping(E)
    E_cross = find_crossover(E, Sn, Se)

    W = column_widths_in["single"]
    fig, ax = plt.subplots(figsize=(W, W * 0.78))

    c_nuc = palette_sweep5[0]
    c_elec = palette_sweep5[4]

    ax.loglog(E, Sn, color=c_nuc, linewidth=1.2, zorder=3)
    ax.loglog(E, Se, color=c_elec, linewidth=1.2, zorder=3)

    S_floor = 0.03
    ax.fill_between(E, S_floor, Sn, where=Sn >= Se,
                    color=c_nuc, alpha=0.10, linewidth=0, zorder=1)
    ax.fill_between(E, S_floor, Se, where=Se > Sn,
                    color=c_elec, alpha=0.10, linewidth=0, zorder=1)

    # crossover marker — short text near the dashed line, no leader arrow
    if not np.isnan(E_cross):
        ax.axvline(E_cross, color="#808080", linestyle="--", linewidth=0.6,
                   zorder=2)
        ax.text(E_cross * 1.15, 0.045,
                fr"$E_{{\mathrm{{cross}}}} \approx {E_cross:.0f}$ keV/u",
                fontsize=6, color="#606060", va="bottom", ha="left")

    # curve labels in clear whitespace, offset from the curves
    ax.text(0.18, 5.5, r"$S_n$", fontsize=8, color=c_nuc, va="bottom")
    ax.text(300, 25, r"$S_e$", fontsize=8, color=c_elec, va="bottom")

    ax.set_xlabel(r"$E$ (keV/u)")
    ax.set_ylabel(r"$S$ (eV / $10^{15}$ atoms cm$^{-2}$)")
    ax.set_xlim(0.1, 1e4)
    ax.set_ylim(0.03, 200)

    critic = TufteCritic()
    issues = critic.critique(fig)
    if issues:
        print(f"TufteCritic: {len(issues)} issue(s)")
        for iss in issues:
            print(f"  {iss}")

    out = "docs/reports/report1/figures/fig01_nuclear_vs_electronic.png"
    fig.savefig(out, dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"Saved → {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
