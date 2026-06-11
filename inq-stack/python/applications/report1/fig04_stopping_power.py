"""fig04 — Stopping power two-panel: E-sweep and σ-sweep.

Panel (a): |S| vs v for σ=1 Bohr E-sweep (WP + classical).
Panel (b): |S| vs σ for E=100 eV σ-sweep.

Stopping power extracted from observables.csv as
  S = |ΔE_total| / |Δz_centroid|
where Δz is estimated from dipole_z change.

Run:
    python -m applications.report1.fig04_stopping_power
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

from applications.report1._shared_style import (
    apply_style,
    palette_sweep5,
    column_widths_in,
    panel_label,
    TufteCritic,
    references,
)

BASE = Path("ResearchProject/systems/jellium")

# ── run specifications ──────────────────────────────────────────────
# E-sweep at σ≈1 Bohr (r_s = 5.69)
E_SWEEP_WP = {
    20:  "run_wp_n162_L50_E20_sigma1",
    25:  "run_wp_n162_L50_E25_sigma1",
    50:  "run_wp_n162_L50_E50_sigma1",
    100: "run_wp_n162_L50_E100_sigma1",
    200: "run_wp_n162_L50_E200_sigma1",
    300: "run_wp_n162_L50_E300_sigma1",
}

E_SWEEP_CL = {
    20:  "run_classical_n162_L50_E20",
    50:  "run_classical_n162_L50_E50",
    100: "run_classical_n162_L50_E100_v2",
    300: "run_classical_n162_L50_E300",
    600: "run_classical_n162_L50_E600_v2",
}

# σ-sweep at E=100 eV
SIGMA_SWEEP = {
    0.5: "run_wp_n162_L50_E100_sigma0p5",
    1.0: "run_wp_n162_L50_E100_sigma1",
    3.0: "run_wp_n162_L50_E100_sigma3",
    8.0: "run_wp_n162_L50_E100_sigma8",
}

M_E_AU = 1.0       # electron mass in a.u.
RS = 5.69           # Wigner-Seitz radius for main runs


def _v_from_E(E_eV: float) -> float:
    """Projectile velocity in a.u. from kinetic energy in eV."""
    E_au = E_eV / 27.211  # eV → Hartree
    return np.sqrt(2 * E_au / M_E_AU)


def extract_stopping(run_dir: str) -> tuple[float, float]:
    """Extract stopping power |S| and its uncertainty from observables.csv.

    Returns (S_eV_per_bohr, delta_z_bohr).
    S = |ΔE_total| / |Δz| where Δz is from dipole_z change.
    """
    csv = BASE / run_dir / "results" / "raw" / "observables" / "observables.csv"
    if not csv.exists():
        return np.nan, np.nan

    df = pd.read_csv(csv)

    # energy change from initial
    E0 = df["energy_total"].iloc[0]
    E_end = df["energy_total"].iloc[-1]
    dE = E_end - E0  # Hartree; positive if system gained energy

    # WP displacement from dipole_z (center of charge)
    z0 = df["dipole_z"].iloc[0]
    z_end = df["dipole_z"].iloc[-1]
    dz = z_end - z0   # Bohr

    if abs(dz) < 0.1:
        return np.nan, dz

    # stopping power: energy deposited into bath per unit distance
    # dE is gained by the system → lost by the WP → S = dE/|dz|
    S_ha_per_bohr = abs(dE) / abs(dz)
    S_eV_per_bohr = S_ha_per_bohr * 27.211

    return S_eV_per_bohr, dz


def main() -> None:
    apply_style()

    W = column_widths_in["full"]
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(W, W * 0.38),
                                      gridspec_kw={"wspace": 0.35})

    c_wp = palette_sweep5[4]    # deep navy
    c_cl = palette_sweep5[0]    # wine red

    # ── Panel (a): S vs v (E-sweep) ────────────────────────────────
    v_wp, S_wp = [], []
    for E, run in sorted(E_SWEEP_WP.items()):
        S, dz = extract_stopping(run)
        if not np.isnan(S):
            v_wp.append(_v_from_E(E))
            S_wp.append(S)

    v_cl, S_cl = [], []
    for E, run in sorted(E_SWEEP_CL.items()):
        S, dz = extract_stopping(run)
        if not np.isnan(S):
            v_cl.append(_v_from_E(E))
            S_cl.append(S)

    if v_wp:
        ax_a.plot(v_wp, S_wp, "o-", color=c_wp, markersize=3.5,
                  linewidth=0.9, label=r"WP ($\sigma = 1$ Bohr)")
    if v_cl:
        ax_a.plot(v_cl, S_cl, "s--", color=c_cl, markersize=3.5,
                  linewidth=0.9, label=r"Classical (point)")

    ax_a.set_xlabel(r"$v$ (a.u.)")
    ax_a.set_ylabel(r"$|S|$ (eV/Bohr)")
    ax_a.set_xscale("log")
    ax_a.set_yscale("log")
    ax_a.legend(fontsize=6, loc="upper left", frameon=True, framealpha=0.9,
                edgecolor="#b0b0b0")
    panel_label(ax_a, "(a)")

    # ── Panel (b): S vs σ (σ-sweep at E=100) ──────────────────────
    sigma_vals, S_sig = [], []
    for sig, run in sorted(SIGMA_SWEEP.items()):
        S, dz = extract_stopping(run)
        if not np.isnan(S):
            sigma_vals.append(sig)
            S_sig.append(S)

    # classical reference at E=100
    S_cl_100, _ = extract_stopping(E_SWEEP_CL.get(100, ""))
    if not np.isnan(S_cl_100):
        ax_b.axhline(S_cl_100, **references["theory"], label="Classical (point)")

    if sigma_vals:
        ax_b.plot(sigma_vals, S_sig, "o-", color=c_wp, markersize=3.5,
                  linewidth=0.9, label=r"WP, $E = 100$ eV")

    # regime shading
    rs = RS
    ax_b.axvspan(0.3, 0.2 * rs, color=palette_sweep5[0], alpha=0.06)
    ax_b.axvspan(0.2 * rs, 0.6 * rs, color="#C0C000", alpha=0.06)
    ax_b.axvspan(0.6 * rs, 12, color=palette_sweep5[4], alpha=0.06)
    ax_b.text(0.5, 0.92, r"$\sigma \ll r_s$", fontsize=5,
              transform=ax_b.transAxes, ha="center", va="top", color="#808080")

    ax_b.set_xlabel(r"$\sigma$ (Bohr)")
    ax_b.set_ylabel(r"$|S|$ (eV/Bohr)")
    ax_b.set_xscale("log")
    ax_b.set_yscale("log")
    ax_b.legend(fontsize=6, loc="upper right", frameon=True, framealpha=0.9,
                edgecolor="#b0b0b0")
    panel_label(ax_b, "(b)")

    critic = TufteCritic()
    issues = critic.critique(fig)
    if issues:
        print(f"TufteCritic: {len(issues)} issue(s)")
        for iss in issues:
            print(f"  {iss}")

    out = "docs/reports/report1/figures/fig04_stopping_power.png"
    fig.savefig(out, dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"Saved → {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
