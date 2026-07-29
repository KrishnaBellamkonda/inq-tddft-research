"""fig_matched_pair — WP–classical matched-pair stopping power comparison.

(a) S(E) log-log from the pipeline: classical Ehrenfest + WP σ=1 + WP σ=5
    + Bethe–Bloch asymptote.  x-axis: E (eV).  Uses S₂ (KS orbital)
    definition because S₁ (momentum) is negative for σ ≤ 1 due to
    free-particle spreading.
(b) Bath-only Hartree energy change ΔE_H^bath(t) at E=100 eV.
    The WP orbital's self-Hartree E_H_self = 1/(2σ_r√π) is subtracted
    using the actual σ_r(t) from wp_real_space_stats.csv, so only the
    bath electron response remains.  Classical line is unmodified
    (the point charge is external, not part of the KS system).

Run:
    python -m applications.report1.fig_matched_pair
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.interpolate import interp1d

from applications.report1._shared_style import (
    apply_style,
    palette_sweep5,
    column_widths_in,
    panel_label,
    references,
    TufteCritic,
)
from applications.report1.stopping_power_data import (
    collect_L50_data,
    bethe_curve,
    AU_TO_FS,
    HA_TO_EV,
)

N_EL_L50 = 162 / 50.0**3

CL_E100_DIR = Path("ResearchProject/systems/jellium/run_classical_n162_L50_E100_v2")
WP_S1_E100_DIR = Path("ResearchProject/systems/jellium/run_wp_n162_L50_E100_sigma1_v2")


def _read_hartree_series(run_dir: Path, t_end_au: float) -> tuple[np.ndarray, np.ndarray]:
    """Read ΔE_H(t) from observables.csv, return (t_fs, dEH_eV)."""
    obs = pd.read_csv(run_dir / "results" / "raw" / "observables" / "observables.csv")
    mask = obs["time_au"] <= t_end_au
    obs = obs[mask]
    t_fs = obs["time_au"].values * AU_TO_FS
    dEH = (obs["energy_hartree"].values - obs["energy_hartree"].iloc[0]) * HA_TO_EV
    return t_fs, dEH


def _wp_self_hartree_correction(
    run_dir: Path, t_end_au: float, obs_times_au: np.ndarray,
) -> np.ndarray:
    """Compute ΔE_H_self(t) = E_H_self(t) - E_H_self(0) using actual σ_r.

    E_H_self = 1/(2 σ_r √π) for a normalised 3D Gaussian charge,
    where σ_r = geometric mean of the three density-width components
    from wp_real_space_stats.csv.

    Returns ΔE_H_self in eV on the obs_times_au grid.
    """
    rs = pd.read_csv(
        run_dir / "results" / "raw" / "observables" / "wp_real_space_stats.csv",
        comment="#",
    )
    mask = rs["time_au"] <= t_end_au * 1.01
    rs = rs[mask]

    sigma_x = np.sqrt(rs["sigma_x2"].values)
    sigma_y = np.sqrt(rs["sigma_y2"].values)
    sigma_z = np.sqrt(rs["sigma_z2"].values)
    sigma_r = (sigma_x * sigma_y * sigma_z) ** (1.0 / 3.0)

    E_self_ha = 1.0 / (2.0 * sigma_r * np.sqrt(np.pi))
    t_au = rs["time_au"].values

    f = interp1d(t_au, E_self_ha, kind="linear", fill_value="extrapolate")
    E_self_interp = f(obs_times_au)

    dE_self_eV = (E_self_interp - E_self_interp[0]) * HA_TO_EV
    return dE_self_eV


def main() -> None:
    apply_style()

    print("Collecting L=50 data from pipeline...")
    data = collect_L50_data()

    # ── Layout with proper spacing ──
    W = column_widths_in["full"]
    fig, (ax_a, ax_b) = plt.subplots(
        1, 2, figsize=(W, W * 0.42),
        gridspec_kw={"wspace": 0.38},
    )

    c_cl = "black"
    c_s1 = palette_sweep5[0]
    c_s5 = palette_sweep5[3]

    # ═══════════════════════════════════════════════════════════════════════
    # Panel (a): S(E) — classical + WP σ=1 + WP σ=5 (S₂ definition)
    # ═══════════════════════════════════════════════════════════════════════

    cl_E_v2, cl_S_v2, cl_err_v2 = [], [], []
    cl_E_v1, cl_S_v1 = [], []
    for r in sorted(data.classical, key=lambda x: x.energy_eV):
        S = r.S_classical
        if np.isfinite(S) and S > 0:
            if r.is_v2:
                cl_E_v2.append(r.energy_eV)
                cl_S_v2.append(S)
                cl_err_v2.append(r.S_classical_err)
            else:
                cl_E_v1.append(r.energy_eV)
                cl_S_v1.append(S)

    all_cl_E = cl_E_v2 + cl_E_v1
    all_cl_S = cl_S_v2 + cl_S_v1
    if all_cl_E:
        order = np.argsort(all_cl_E)
        ax_a.plot(np.array(all_cl_E)[order], np.array(all_cl_S)[order],
                  "-", color=c_cl, linewidth=0.8, zorder=3)
    if cl_E_v2:
        ax_a.errorbar(cl_E_v2, cl_S_v2, yerr=cl_err_v2,
                      fmt="s", color=c_cl, markersize=4,
                      markeredgecolor="white", markeredgewidth=0.3,
                      capsize=1.5, capthick=0.5, elinewidth=0.5,
                      label="Classical", zorder=4)
    if cl_E_v1:
        ax_a.plot(cl_E_v1, cl_S_v1, "s", color=c_cl, markersize=4,
                  markerfacecolor="white", markeredgecolor=c_cl,
                  markeredgewidth=0.8,
                  label="Classical (v1)", zorder=4)

    # WP σ=1 (S₂)
    s1_E_v2, s1_S_v2, s1_E_v1, s1_S_v1 = [], [], [], []
    for r in sorted(data.wp_sigma1, key=lambda x: x.energy_eV):
        S = r.S2_eV_per_bohr
        if np.isfinite(S) and S > 0:
            if r.is_v2:
                s1_E_v2.append(r.energy_eV)
                s1_S_v2.append(S)
            else:
                s1_E_v1.append(r.energy_eV)
                s1_S_v1.append(S)

    all_s1_E = s1_E_v2 + s1_E_v1
    all_s1_S = s1_S_v2 + s1_S_v1
    if all_s1_E:
        order = np.argsort(all_s1_E)
        ax_a.plot(np.array(all_s1_E)[order], np.array(all_s1_S)[order],
                  "-", color=c_s1, linewidth=0.8, zorder=3)
    if s1_E_v2:
        ax_a.plot(s1_E_v2, s1_S_v2, "o", color=c_s1, markersize=4.5,
                  markeredgecolor="white", markeredgewidth=0.3,
                  label=r"WP $\sigma{=}1$", zorder=5)
    if s1_E_v1:
        ax_a.plot(s1_E_v1, s1_S_v1, "o", color=c_s1, markersize=4.5,
                  markerfacecolor="white", markeredgecolor=c_s1,
                  markeredgewidth=0.8, zorder=5)

    # WP σ=5 (S₂)
    s5_E_v2, s5_S_v2, s5_E_v1, s5_S_v1 = [], [], [], []
    for r in sorted(data.wp_sigma5, key=lambda x: x.energy_eV):
        S = r.S2_eV_per_bohr
        if np.isfinite(S) and S > 0:
            if r.is_v2:
                s5_E_v2.append(r.energy_eV)
                s5_S_v2.append(S)
            else:
                s5_E_v1.append(r.energy_eV)
                s5_S_v1.append(S)

    all_s5_E = s5_E_v2 + s5_E_v1
    all_s5_S = s5_S_v2 + s5_S_v1
    if all_s5_E:
        order = np.argsort(all_s5_E)
        ax_a.plot(np.array(all_s5_E)[order], np.array(all_s5_S)[order],
                  "-", color=c_s5, linewidth=0.8, zorder=3)
    if s5_E_v2:
        ax_a.plot(s5_E_v2, s5_S_v2, "^", color=c_s5, markersize=4.5,
                  markeredgecolor="white", markeredgewidth=0.3,
                  label=r"WP $\sigma{=}5$", zorder=5)
    if s5_E_v1:
        ax_a.plot(s5_E_v1, s5_S_v1, "^", color=c_s5, markersize=4.5,
                  markerfacecolor="white", markeredgecolor=c_s5,
                  markeredgewidth=0.8, zorder=5)

    # Bethe asymptote
    E_grid = np.logspace(np.log10(15), np.log10(3000), 200)
    S_bethe = bethe_curve(E_grid, N_EL_L50)
    mask = S_bethe > 0
    ax_a.plot(E_grid[mask], S_bethe[mask], **references["asymptote"],
              label=r"Bethe $v^{-2}$", zorder=2)

    ax_a.set_xscale("log")
    ax_a.set_yscale("log")
    ax_a.set_xlabel(r"Projectile energy $E$ (eV)")
    ax_a.set_ylabel(r"Stopping power $S_2$ (eV/Bohr)")
    ax_a.legend(fontsize=5.5, loc="upper right", frameon=True,
                edgecolor="#b0b0b0", fancybox=False)
    panel_label(ax_a, "(a)")

    # ═══════════════════════════════════════════════════════════════════════
    # Panel (b): Bath-only ΔE_H(t) at E=100 eV
    #   WP: subtract self-Hartree of the WP orbital (spreading Gaussian)
    #   Classical: unmodified (point charge is external to KS system)
    # ═══════════════════════════════════════════════════════════════════════

    wp_100_result = None
    for r in data.wp_sigma1:
        if r.energy_eV == 100:
            wp_100_result = r
            break
    t_end_au = wp_100_result.window.t_end if wp_100_result else 8.27

    print(f"\nPanel (b): Bath Hartree energy at E=100 eV")
    print(f"  Window: [0, {t_end_au:.2f}] a.u. = [0, {t_end_au * AU_TO_FS:.3f}] fs")

    # Classical: ΔE_H as-is
    t_cl, dEH_cl = _read_hartree_series(CL_E100_DIR, t_end_au)

    # WP: ΔE_H minus WP self-Hartree
    obs_wp = pd.read_csv(
        WP_S1_E100_DIR / "results" / "raw" / "observables" / "observables.csv"
    )
    obs_wp = obs_wp[obs_wp["time_au"] <= t_end_au]
    t_wp_fs = obs_wp["time_au"].values * AU_TO_FS
    dEH_wp_total = (obs_wp["energy_hartree"].values - obs_wp["energy_hartree"].iloc[0]) * HA_TO_EV

    dE_self = _wp_self_hartree_correction(
        WP_S1_E100_DIR, t_end_au, obs_wp["time_au"].values,
    )
    dEH_wp_bath = dEH_wp_total - dE_self

    print(f"  Classical ΔE_H:       {dEH_cl[0]:.3f} → {dEH_cl[-1]:.3f} eV")
    print(f"  WP total  ΔE_H:      {dEH_wp_total[0]:.3f} → {dEH_wp_total[-1]:.3f} eV")
    print(f"  WP self-H ΔE_H_self: {dE_self[0]:.3f} → {dE_self[-1]:.3f} eV")
    print(f"  WP bath   ΔE_H_bath: {dEH_wp_bath[0]:.3f} → {dEH_wp_bath[-1]:.3f} eV")

    ax_b.plot(t_cl, dEH_cl, "-", color=c_cl, linewidth=1.0,
              label="Classical", zorder=3)
    ax_b.plot(t_wp_fs, dEH_wp_bath, "-", color=c_s1, linewidth=1.0,
              label=r"WP $\sigma{=}1$", zorder=3)
    ax_b.axhline(0, color="#808080", linewidth=0.4, zorder=1)

    ax_b.set_xlabel(r"Time $t$ (fs)")
    ax_b.set_ylabel(r"$\Delta E_H^{\mathrm{bath}}(t)$ (eV)")
    ax_b.legend(fontsize=5.5, loc="upper left", frameon=True,
                edgecolor="#b0b0b0", fancybox=False)
    panel_label(ax_b, "(b)")

    # ── Save ──
    critic = TufteCritic()
    issues = critic.critique(fig)
    if issues:
        for iss in issues:
            print(f"  TufteCritic: {iss}")

    out = "docs/reports/report1/figures/fig_matched_pair.png"
    fig.savefig(out, dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"\nSaved -> {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
