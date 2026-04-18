#!/usr/bin/env python3
"""
03_ecut_convergence/plot_convergence.py
========================================
Post-processing for the E_cut convergence sweep.

Reads:   results/ecut_convergence.csv
Writes:  results/ecut_convergence.png  — total energy vs E_cut (Ha and eV)
         results/ecut_delta_E.png      — ΔE between consecutive E_cut values (meV)

The convergence criterion used in Tsubonoya et al. (PRB 90, 035416, 2014) is a
grid spacing of 0.16 Å ≈ 0.302 bohr, corresponding to E_cut = (π/h)²/2 ≈ 54 Ha.
We independently verify convergence by checking ΔE < 1 meV between consecutive
E_cut values, and identify the converged cutoff from the plateau.

Usage:
    cd 03_ecut_convergence
    python3 plot_convergence.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

RESULTS = "results"
CSV_FILE = os.path.join(RESULTS, "ecut_convergence.csv")

HA_TO_EV  = 27.211386245988   # CODATA 2018
HA_TO_MEV = HA_TO_EV * 1000.0

def load_csv(path):
    """Load ecut_convergence.csv; return structured dict."""
    data = np.loadtxt(path, delimiter=",", comments="#")
    return {
        "ecut_Ha"   : data[:, 0],
        "E_total_Ha": data[:, 1],
        "E_total_eV": data[:, 2],
        "grid_pts"  : data[:, 3].astype(int),
        "scf_steps" : data[:, 4].astype(int),
    }


def plot_energy_vs_ecut(d):
    """Plot total energy (Ha) vs E_cut, and grid size on secondary axis."""
    ecut = d["ecut_Ha"]
    etot = d["E_total_Ha"]

    fig, ax1 = plt.subplots(figsize=(7, 4.5))

    ax1.plot(ecut, etot, "o-", color="steelblue", linewidth=2, markersize=7,
             label="$E_{\\mathrm{total}}$ (Ha)")
    ax1.set_xlabel("$E_{\\mathrm{cut}}$ (Ha)", fontsize=13)
    ax1.set_ylabel("Total energy (Ha)", fontsize=13, color="steelblue")
    ax1.tick_params(axis="y", labelcolor="steelblue")

    ax2 = ax1.twinx()
    ax2.plot(ecut, d["grid_pts"] / 1e6, "s--", color="darkorange",
             linewidth=1.5, markersize=6, alpha=0.7, label="Grid points (M)")
    ax2.set_ylabel("Grid points (millions)", fontsize=12, color="darkorange")
    ax2.tick_params(axis="y", labelcolor="darkorange")

    # Mark the paper's target E_cut
    ax1.axvline(54, color="red", linestyle=":", linewidth=1.5, alpha=0.8,
                label="Paper: 54 Ha (0.16 Å grid)")

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, fontsize=10, loc="lower right")

    ax1.set_title("Coronene C₂₄H₁₂ — E$_{\\mathrm{cut}}$ convergence (LDA ground state)",
                  fontsize=12)
    ax1.grid(True, alpha=0.3)
    plt.tight_layout()
    out = os.path.join(RESULTS, "ecut_convergence.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out}")


def plot_delta_E(d):
    """Plot |ΔE| between consecutive E_cut values in meV."""
    ecut = d["ecut_Ha"]
    etot = d["E_total_Ha"]

    delta_meV = np.abs(np.diff(etot)) * HA_TO_MEV
    ecut_mid  = 0.5 * (ecut[:-1] + ecut[1:])   # midpoint x values for clarity

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.semilogy(ecut_mid, delta_meV, "o-", color="steelblue", linewidth=2, markersize=8)

    # Convergence threshold
    ax.axhline(1.0, color="red", linestyle="--", linewidth=1.5, label="1 meV threshold")

    # Annotate each point with its Δ value
    for x, y in zip(ecut_mid, delta_meV):
        ax.annotate(f"{y:.2f}", (x, y), textcoords="offset points",
                    xytext=(4, 6), fontsize=8.5, color="steelblue")

    ax.set_xlabel("$E_{\\mathrm{cut}}$ midpoint (Ha)", fontsize=13)
    ax.set_ylabel("|ΔE| (meV)", fontsize=13)
    ax.set_title("Energy change between consecutive $E_{\\mathrm{cut}}$ values\n"
                 "Convergence criterion: |ΔE| < 1 meV", fontsize=12)
    ax.legend(fontsize=11)
    ax.grid(True, which="both", alpha=0.3)
    plt.tight_layout()
    out = os.path.join(RESULTS, "ecut_delta_E.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out}")


def print_summary(d):
    """Print a human-readable table and identify the converged E_cut."""
    ecut = d["ecut_Ha"]
    etot = d["E_total_Ha"]
    gpts = d["grid_pts"]
    steps= d["scf_steps"]

    print("\n=== E_cut convergence summary ===")
    print(f"{'E_cut (Ha)':>11}  {'E_total (Ha)':>14}  {'E_total (eV)':>14}"
          f"  {'ΔE (meV)':>10}  {'Grid pts':>10}  {'SCF steps':>10}")
    print("-" * 80)
    for i in range(len(ecut)):
        de_str = "        —"
        if i > 0:
            de_mev = abs(etot[i] - etot[i-1]) * HA_TO_MEV
            de_str = f"{de_mev:10.3f}"
        print(f"  {ecut[i]:8.1f}  {etot[i]:14.6f}  {etot[i]*HA_TO_EV:14.4f}"
              f"  {de_str}  {gpts[i]:10d}  {steps[i]:10d}")

    # Find converged E_cut: first consecutive pair with ΔE < 1 meV
    delta_mev = np.abs(np.diff(etot)) * HA_TO_MEV
    converged_ecut = None
    for i, de in enumerate(delta_mev):
        if de < 1.0:
            converged_ecut = ecut[i+1]
            print(f"\n  CONVERGED at E_cut = {converged_ecut:.0f} Ha"
                  f"  (|ΔE| = {de:.3f} meV < 1 meV between"
                  f" {ecut[i]:.0f} and {ecut[i+1]:.0f} Ha)")
            break
    if converged_ecut is None:
        print("\n  WARNING: convergence not reached in sweep range — extend E_cut.")

    # Grid spacing from E_cut: h = π / sqrt(2 E_cut)
    print("\n  Grid spacing implied by each E_cut:")
    print(f"  {'E_cut (Ha)':>11}  {'h (bohr)':>10}  {'h (Å)':>10}")
    BOHR_TO_ANG = 0.529177210903
    for ec in ecut:
        h_bohr = np.pi / np.sqrt(2.0 * ec)
        print(f"    {ec:8.1f}  {h_bohr:10.4f}  {h_bohr*BOHR_TO_ANG:10.4f}")


if __name__ == "__main__":
    if not os.path.exists(CSV_FILE):
        print(f"ERROR: {CSV_FILE} not found. Run inq-run first.")
        raise SystemExit(1)

    d = load_csv(CSV_FILE)

    print("=== Energy vs E_cut plot ===")
    plot_energy_vs_ecut(d)

    print("\n=== ΔE convergence plot ===")
    plot_delta_E(d)

    print_summary(d)
    print("\nDone.")
