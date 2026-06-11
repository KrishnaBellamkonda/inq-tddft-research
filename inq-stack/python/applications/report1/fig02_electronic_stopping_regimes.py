"""fig02 — Electronic stopping regime diagram.

Full-width regime classification on the (v/v_F, kappa) plane with
colored regions for each historical theory's validity domain and the
electron-projectile constraint line. Theory-only view (no run data).

Regime boundaries (dimensional analysis for free-electron gas):
  x-axis: v/v_F  (host-response axis)
  y-axis: kappa = 2|Z_1|/v  (projectile-scattering axis)

For Z_1=1 electron in r_s=5.69 jellium:
  v_F = k_F = 0.337 a.u., omega_p = 3.47 eV
  Constraint: kappa * (v/v_F) = 2/v_F = 5.93

Run:
    python -m applications.report1.fig02_electronic_stopping_regimes
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from applications.report1._shared_style import (
    apply_style,
    column_widths_in,
    TufteCritic,
)

# ---------------------------------------------------------------------------
# Physical constants: r_s = 5.69 jellium (a.u.)
# ---------------------------------------------------------------------------
HA_TO_EV = 27.2114
N_ELECTRONS = 1.296e-3   # electron density (a.u.)
K_F = 0.337
V_F = K_F                # v_F = k_F for free-electron gas
OMEGA_P = np.sqrt(4 * np.pi * N_ELECTRONS)  # = 0.1276 Ha = 3.47 eV

# Electron-projectile constraint: kappa * (v/v_F) = 2/v_F
KAPPA_TIMES_VOVF = 2.0 / V_F   # = 5.93

# ---------------------------------------------------------------------------
# Regime colours (consistent with meeting-report palette)
# ---------------------------------------------------------------------------
REGIME_COLOURS = {
    "friction": "#e74c3c",
    "bragg":    "#f39c12",
    "bohr":     "#8b4513",
    "bloch":    "#9b59b6",
    "bethe":    "#3498db",
    "tddft":    "#2c3e50",
}

def draw_regime_regions(ax, *, alpha: float = 0.22, fs: float = 7.5) -> None:
    """Draw colored regime regions on a log-log (v/v_F, kappa) plot.

    Regions from dimensional analysis:
      - Friction (Fermi-Teller, ENRA, nonlinear DFT): v/v_F <= 0.5
      - Bragg peak (Lindhard RPA): 0.5 < v/v_F < 5, kappa <= 5
      - Bohr classical (impact parameter): kappa >= 5, v/v_F >= 0.5
      - Bloch interpolation (kappa ~ 1): 0.5 < kappa < 5, v/v_F >= 3
      - Bethe (perturbative, Born): v/v_F >= 5, kappa <= 0.5
    """
    xlo, xhi = ax.get_xlim()
    ylo, yhi = ax.get_ylim()

    # Friction: v/v_F <= 0.5
    ax.axvspan(xlo, 0.5, alpha=alpha, color=REGIME_COLOURS["friction"], zorder=0)
    ax.text(0.18, 1.8, r"Friction" + "\n" + r"$S \propto v$",
            ha="center", va="center", fontsize=fs, color="darkred",
            fontweight="bold")

    # Bragg / Lindhard: 0.5 < v/v_F < 5, kappa <= 5
    ax.fill_between([0.5, 5], ylo, 5, alpha=alpha,
                    color=REGIME_COLOURS["bragg"], zorder=0)
    ax.text(2.2, 3.5, r"Bragg peak" + "\n" + r"(Lindhard RPA)",
            ha="center", va="center", fontsize=fs - 0.5, color="#b07800",
            fontweight="bold", zorder=4,
            bbox=dict(boxstyle="round,pad=0.2", fc="white",
                      ec="none", alpha=0.85))

    # Bohr classical: kappa >= 5, v/v_F >= 0.5
    ax.fill_between([0.5, xhi], 5, yhi, alpha=alpha,
                    color=REGIME_COLOURS["bohr"], zorder=0)
    ax.text(4.0, 30, r"Bohr classical" + "\n" + r"$\kappa \gg 1$",
            ha="center", va="center", fontsize=fs, color="#5a2d0c",
            fontweight="bold")

    # Bloch interpolation: 0.5 < kappa < 5, v/v_F >= 3
    ax.fill_between([3, xhi], 0.5, 5, alpha=alpha,
                    color=REGIME_COLOURS["bloch"], zorder=0)
    ax.text(18, 2.0, r"Bloch" + "\n" + r"$\kappa \sim 1$",
            ha="center", va="center", fontsize=fs, color="purple",
            fontweight="bold")

    # Bethe: v/v_F >= 5, kappa <= 0.5
    ax.fill_between([5, xhi], ylo, 0.5, alpha=alpha,
                    color=REGIME_COLOURS["bethe"], zorder=0)
    ax.text(25, 0.08, r"Bethe" + "\n" + r"$S \sim \ln v / v^2$",
            ha="center", va="center", fontsize=fs, color="darkblue",
            fontweight="bold")


def draw_boundary_lines(ax) -> None:
    """Draw kappa=1 and v=v_F reference boundaries."""
    xlo, xhi = ax.get_xlim()
    ylo, yhi = ax.get_ylim()

    ax.axhline(1.0, color="black", lw=0.6, ls=":", zorder=1, alpha=0.7)
    ax.axvline(1.0, color="black", lw=0.6, ls=":", zorder=1, alpha=0.7)

    # Labels at margins
    ax.text(xhi * 0.55, 1.15, r"$\kappa = 1$",
            fontsize=7, color="black", alpha=0.7)
    ax.text(1.08, ylo * 2.0, r"$v = v_F$",
            fontsize=7, color="black", alpha=0.7, rotation=90)


def draw_constraint_line(ax) -> None:
    """Draw the electron-projectile constraint hyperbola."""
    vf_grid = np.logspace(-1, 2, 200)
    kappa_grid = KAPPA_TIMES_VOVF / vf_grid
    ax.plot(vf_grid, kappa_grid, "k--", lw=1.4, zorder=3)


def draw_tddft_annotation(ax) -> None:
    """Small italic text box for the rt-TDDFT annotation."""
    ax.text(
        0.55, 55,
        r"\textit{nonlinear rt-TDDFT (this work)}",
        ha="left", va="center", fontsize=7,
        color=REGIME_COLOURS["tddft"],
        bbox=dict(boxstyle="round,pad=0.25", fc="white",
                  ec=REGIME_COLOURS["tddft"], lw=0.5, alpha=0.9),
        zorder=5,
    )


def main() -> None:
    apply_style()

    W = column_widths_in["full"]
    fig, ax = plt.subplots(figsize=(W, W * 0.58))

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.1, 100)
    ax.set_ylim(0.03, 100)

    # Layer 1: regime regions
    draw_regime_regions(ax)

    # Layer 2: boundary reference lines
    draw_boundary_lines(ax)

    # Layer 3: electron-projectile constraint line
    draw_constraint_line(ax)

    # Layer 4: rt-TDDFT annotation
    draw_tddft_annotation(ax)

    # Axis labels
    ax.set_xlabel(r"$v / v_F$  (host-response axis)")
    ax.set_ylabel(r"$\kappa = 2|Z_1|/v$  (projectile-scattering axis)")

    # Minimal grid — very faint
    ax.grid(False)

    ax.tick_params(direction="in", which="both")

    # Tufte critique
    critic = TufteCritic()
    issues = critic.critique(fig)
    if issues:
        print(f"TufteCritic: {len(issues)} issue(s)")
        for iss in issues:
            print(f"  {iss}")

    out = Path("docs/reports/report1/figures/fig02_electronic_stopping_regimes.png")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"Saved -> {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
