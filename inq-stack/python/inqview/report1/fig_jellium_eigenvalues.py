"""fig_jellium_eigenvalues — KS eigenvalue level diagram for jellium GS (N=162, L=50).

Shows occupied KS eigenvalues grouped by |G|^2 shell, with degeneracy visible
as the number of horizontal lines per group. Overlays analytical free-electron
eigenvalues E_G = |G|^2 (2pi/L)^2 / 2 for comparison.

Run:
    python -m inqview.report1.fig_jellium_eigenvalues
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

from inqview.report1._shared_style import (
    apply_style,
    column_widths_in,
    panel_label,
    TufteCritic,
    palette_sweep5,
    references,
)

# --- Paths -------------------------------------------------------------------
EIGEN_CSV = (
    "ResearchProject/systems/jellium/save_gs/gs_L50_cubic_N162_dx0p40/"
    "results/raw/observables/eigenvalues/eigenvalues.csv"
)
OCC_CSV = (
    "ResearchProject/systems/jellium/save_gs/gs_L50_cubic_N162_dx0p40/"
    "results/raw/observables/eigenvalues/occupations.csv"
)
OUT = "docs/reports/report1/figures/fig_jellium_eigenvalues.png"

# --- Physical constants ------------------------------------------------------
L_BOHR = 50.0
HA_TO_EV = 27.21138625


def free_electron_shells(L: float, max_g2: int = 12):
    """Compute free-electron energy levels for cubic box of side L.

    Returns dict mapping |G|^2 -> (energy_eV, degeneracy).
    """
    shells = {}
    k_unit = 2.0 * np.pi / L
    for nx in range(-5, 6):
        for ny in range(-5, 6):
            for nz in range(-5, 6):
                g2 = nx**2 + ny**2 + nz**2
                if g2 > max_g2:
                    continue
                if g2 not in shells:
                    shells[g2] = {"energy_ha": 0.5 * g2 * k_unit**2, "count": 0}
                shells[g2]["count"] += 1

    result = {}
    for g2 in sorted(shells.keys()):
        result[g2] = (shells[g2]["energy_ha"] * HA_TO_EV, shells[g2]["count"])
    return result


def assign_shells(eigenvalues_ev: np.ndarray, shells: dict) -> list[int]:
    """Assign each KS eigenvalue to a free-electron shell by proximity."""
    shell_list = sorted(shells.keys())
    shell_energies = np.array([shells[g2][0] for g2 in shell_list])
    assignments = []
    for e in eigenvalues_ev:
        idx = np.argmin(np.abs(e - shell_energies))
        assignments.append(shell_list[idx])
    return assignments


def main() -> None:
    apply_style()

    # --- Load data -----------------------------------------------------------
    df_e = pd.read_csv(EIGEN_CSV)
    df_o = pd.read_csv(OCC_CSV)
    merged = df_e.merge(df_o, on="state_index")

    eigen_ev = merged["eigenvalue_ev"].values
    occ = merged["occupation"].values

    # Shift so lowest eigenvalue = 0
    e_shift = eigen_ev - eigen_ev.min()

    # --- Free-electron shells ------------------------------------------------
    shells = free_electron_shells(L_BOHR, max_g2=12)
    shell_assignments = assign_shells(e_shift, shells)

    # --- Group states by shell -----------------------------------------------
    shell_groups = {}
    for i, g2 in enumerate(shell_assignments):
        if g2 not in shell_groups:
            shell_groups[g2] = {"energies": [], "occs": []}
        shell_groups[g2]["energies"].append(e_shift[i])
        shell_groups[g2]["occs"].append(occ[i])

    # --- Figure: energy level diagram ----------------------------------------
    W = column_widths_in["full"]
    fig, ax = plt.subplots(figsize=(W, 3.0))

    colors_occ = palette_sweep5[0]     # wine red for occupied
    colors_unocc = "#AAAAAA"           # grey for unoccupied

    # Each shell gets a horizontal band. Within each band, draw horizontal lines
    # representing the individual KS states, spread vertically for visibility.
    x_positions = {}  # g2 -> x-center for each shell column
    sorted_shells = sorted(shell_groups.keys())
    bar_width = 0.6

    for col_idx, g2 in enumerate(sorted_shells):
        x_center = col_idx
        x_positions[g2] = x_center
        energies = shell_groups[g2]["energies"]
        occs = shell_groups[g2]["occs"]
        n_states_in_shell = len(energies)

        # Draw individual state lines
        for j, (e, o) in enumerate(zip(energies, occs)):
            color = colors_occ if o > 0.5 else colors_unocc
            lw = 1.2 if o > 0.5 else 0.8
            x_lo = x_center - bar_width / 2
            x_hi = x_center + bar_width / 2
            ax.plot([x_lo, x_hi], [e, e], color=color, linewidth=lw,
                    solid_capstyle="round")

        # Free-electron level as a dashed line (wider)
        fe_e = shells[g2][0]
        ax.plot(
            [x_center - bar_width * 0.7, x_center + bar_width * 0.7],
            [fe_e, fe_e],
            color="#000000", linewidth=1.0, linestyle="--", alpha=0.6,
        )

        # Shell label at bottom
        degen = shells[g2][1]
        n_found = n_states_in_shell
        ax.text(x_center, -0.12, rf"$|G|^2\!=\!{g2}$",
                fontsize=6.5, ha="center", va="top", transform=ax.get_xaxis_transform())
        ax.text(x_center, -0.20, f"({degen})",
                fontsize=5.5, ha="center", va="top", color="#606060",
                transform=ax.get_xaxis_transform())

    # Fermi level — inside the plot
    occ_mask = occ > 0.5
    e_fermi = e_shift[occ_mask].max()
    ax.axhline(e_fermi, color=palette_sweep5[2], linewidth=0.9, linestyle=":",
               zorder=0, alpha=0.7)
    ax.text(0.97, e_fermi, rf"$\varepsilon_F$",
            fontsize=8, color=palette_sweep5[2], ha="right", va="bottom",
            transform=ax.get_yaxis_transform())

    # HOMO annotation (state 80 — last occupied)
    homo_idx = int(occ_mask.sum()) - 1
    homo_e = e_shift[homo_idx]
    homo_g2 = shell_assignments[homo_idx]
    homo_x = x_positions[homo_g2]
    ax.annotate("HOMO", xy=(homo_x + bar_width / 2 + 0.05, homo_e),
                fontsize=6, color=colors_occ, ha="left", va="center",
                arrowprops=dict(arrowstyle="-", color=colors_occ, lw=0.6),
                xytext=(homo_x + bar_width / 2 + 0.4, homo_e + 0.08))

    # LUMO annotation (state 81 — first unoccupied)
    lumo_idx = homo_idx + 1
    lumo_e = e_shift[lumo_idx]
    lumo_g2 = shell_assignments[lumo_idx]
    lumo_x = x_positions[lumo_g2]
    ax.annotate("LUMO", xy=(lumo_x + bar_width / 2 + 0.05, lumo_e),
                fontsize=6, color=colors_unocc, ha="left", va="center",
                arrowprops=dict(arrowstyle="-", color=colors_unocc, lw=0.6),
                xytext=(lumo_x + bar_width / 2 + 0.4, lumo_e - 0.08))

    # WP orbital annotation (state 100 — last extra state)
    wp_idx = len(e_shift) - 1
    wp_e = e_shift[wp_idx]
    wp_g2 = shell_assignments[wp_idx]
    wp_x = x_positions[wp_g2]
    ax.annotate("WP orbital", xy=(wp_x - bar_width / 2 - 0.05, wp_e),
                fontsize=6, color="#188048", ha="right", va="center",
                arrowprops=dict(arrowstyle="-", color="#188048", lw=0.6),
                xytext=(wp_x - bar_width / 2 - 0.6, wp_e + 0.06))

    ax.set_ylabel(r"$\varepsilon_i - \varepsilon_0$ (eV)")
    ax.set_xlim(-0.7, len(sorted_shells) - 0.3)
    ax.set_ylim(-0.05, e_shift.max() * 1.05)

    # Remove x-axis ticks (labels are custom text)
    ax.set_xticks([])
    ax.set_xlabel(r"Free-electron shell $|G|^2$")

    # Legend
    legend_elements = [
        Line2D([0], [0], color=colors_occ, linewidth=1.5,
               label=f"Occupied ({int(occ_mask.sum())} states)"),
        Line2D([0], [0], color=colors_unocc, linewidth=1.0,
               label=f"Unoccupied ({int((~occ_mask).sum())} states)"),
        Line2D([0], [0], color="black", linewidth=1.0, linestyle="--",
               alpha=0.6, label=r"Free-electron $E_{|G|^2}$"),
    ]
    ax.legend(handles=legend_elements, loc="upper left", fontsize=7)

    # --- Save ----------------------------------------------------------------
    critic = TufteCritic()
    issues = critic.critique(fig)
    if issues:
        print(f"TufteCritic: {len(issues)} issue(s)")
        for iss in issues:
            print(f"  {iss}")

    fig.savefig(OUT, dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"Saved -> {OUT}")
    plt.close(fig)


if __name__ == "__main__":
    main()
