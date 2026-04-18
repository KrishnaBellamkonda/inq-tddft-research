#!/usr/bin/env python3
"""
plot_spreading.py
─────────────────
Static publication-quality figure: σ(t) numerical vs. analytical.

References
----------
[1] Angelo, "Quantum Wave-Packet Preparation and Electron Dynamics in
    Jellium", Cavendish Lab Report, Candidate 3221L.
    Eq. (6): σ(t) = σ₀ √(1 + t²/(4σ₀⁴))
    Eq. (7): T_rev = L²/π
[3] Robinett, Phys. Rep. 392, 1–119 (2004). Wave packet revivals.

Usage
-----
    cd ResearchProject/systems/jellium
    python3 analysis/plot_spreading.py

Output
------
    results/sigma_comparison.pdf   — σ(t) figure
    results/energy_conservation.pdf — total energy vs time
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ─── Paths ────────────────────────────────────────────────────────────────────
SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(SCRIPT_DIR, "..", "results")
SIGMA_FILE  = os.path.join(RESULTS_DIR, "sigma_t.txt")
GRID_FILE   = os.path.join(RESULTS_DIR, "grid_info.txt")

# ─── Read grid_info.txt ───────────────────────────────────────────────────────
def read_grid_info(path):
    info = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            key = parts[0]
            try:
                val = float(parts[1])
                info[key] = val
            except (IndexError, ValueError):
                pass
    return info

# ─── Main ─────────────────────────────────────────────────────────────────────
def main():
    if not os.path.exists(SIGMA_FILE):
        print(f"ERROR: {SIGMA_FILE} not found. Run the simulation first.")
        sys.exit(1)

    # Load σ(t) data
    # Columns: t  sigma_num_z  sigma_ana  cx  cy  cz  E_total  N_elec
    data = np.loadtxt(SIGMA_FILE, comments="#")
    if data.ndim == 1:
        data = data[np.newaxis, :]

    t         = data[:, 0]
    sigma_num = data[:, 1]
    sigma_ana = data[:, 2]
    cx        = data[:, 3]
    cy        = data[:, 4]
    cz        = data[:, 5]
    E_total   = data[:, 6]
    N_elec    = data[:, 7]

    # Read parameters
    info = {}
    if os.path.exists(GRID_FILE):
        info = read_grid_info(GRID_FILE)
    L        = info.get("L",       40.0)
    sigma_wp = info.get("sigma_wp", 1.0)
    k0       = info.get("k0",       1.5)
    T_rev    = L**2 / np.pi

    # Dense analytical curve [1] Eq.(6)
    t_dense  = np.linspace(0, t[-1] * 1.05, 500)
    sig_dense = sigma_wp * np.sqrt(1.0 + t_dense**2 / (4.0 * sigma_wp**4))

    # ── Figure 1: σ(t) ────────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle(
        f"Gaussian Wave Packet Spreading in Periodic Box\n"
        f"L = {L:.0f} bohr,  σ₀ = {sigma_wp:.1f} a₀,  "
        f"k₀ = {k0:.1f} a₀⁻¹,  T_rev = {T_rev:.1f} atu",
        fontsize=13, y=0.98
    )

    # ── Panel (a): σ(t) linear scale ─────────────────────────────────────────
    ax = axes[0, 0]
    ax.plot(t_dense, sig_dense, "r-",  lw=2,   label=r"Analytical $\sigma_0\sqrt{1+t^2/4\sigma_0^4}$  [1] Eq.(6)")
    ax.plot(t,       sigma_num, "b.",  ms=6,   label="INQ numerical")
    ax.axvline(T_rev, color="gray", ls="--", lw=1, alpha=0.7, label=f"$T_{{rev}}$ = {T_rev:.0f} atu")
    if T_rev < t[-1]:
        ax.axvline(T_rev/2, color="gray", ls=":", lw=1, alpha=0.5, label=f"$T_{{rev}}/2$ = {T_rev/2:.0f} atu")
    ax.set_xlabel("Time $t$ [atu]")
    ax.set_ylabel(r"Width $\sigma(t)$ [a$_0$]")
    ax.set_title("(a)  Packet width — linear scale")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # ── Panel (b): σ(t) log-log ───────────────────────────────────────────────
    ax = axes[0, 1]
    mask = (t > 0) & (sigma_num > 0)
    ax.loglog(t_dense[t_dense > 0], sig_dense[t_dense > 0], "r-",  lw=2,
              label="Analytical")
    ax.loglog(t[mask], sigma_num[mask],  "b.",  ms=6, label="INQ numerical")
    # Show t/2 asymptote (long-time limit)
    t_asym = t_dense[t_dense > 2]
    ax.loglog(t_asym, t_asym / 2.0, "k--", lw=1, alpha=0.6,
              label=r"$t/2$ (long-time limit)")
    ax.set_xlabel("Time $t$ [atu]")
    ax.set_ylabel(r"Width $\sigma(t)$ [a$_0$]")
    ax.set_title("(b)  Packet width — log–log scale")
    ax.legend(fontsize=9)
    ax.grid(True, which="both", alpha=0.3)

    # ── Panel (c): Relative error (σ_num − σ_ana) / σ_ana ────────────────────
    ax = axes[1, 0]
    rel_err = (sigma_num - sigma_ana) / sigma_ana * 100.0
    ax.plot(t, rel_err, "g-o", ms=4, lw=1.2)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("Time $t$ [atu]")
    ax.set_ylabel(r"$(σ_{num} - σ_{ana})/σ_{ana}$ [%]")
    ax.set_title("(c)  Relative error in σ(t)")
    ax.grid(True, alpha=0.3)

    # ── Panel (d): Centre of mass ⟨z⟩(t) ─────────────────────────────────────
    ax = axes[1, 1]
    # Expected: ⟨z⟩ = k₀·t wrapped to [-L/2, L/2) (INQ symmetric coords)
    cz_expected = ((k0 * t + L/2) % L) - L/2
    ax.plot(t, cz,          "b-",  lw=1.5, label=r"$\langle z \rangle$ (numerical)")
    ax.plot(t, cz_expected, "r--", lw=1.5, label=r"$k_0 t$ mod $L$ (expected)")
    ax.set_xlabel("Time $t$ [atu]")
    ax.set_ylabel(r"$\langle z \rangle$ [a$_0$]")
    ax.set_title("(d)  Centre of mass $⟨z⟩(t)$")
    ax.set_ylim(0, L)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    out_path = os.path.join(RESULTS_DIR, "sigma_comparison.pdf")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {out_path}")
    plt.close(fig)

    # ── Figure 2: Energy conservation ─────────────────────────────────────────
    # Skip t=0 row (E_total=0 stored as placeholder before propagation starts)
    mask_e = E_total != 0.0
    t_e, E_e = t[mask_e], E_total[mask_e]
    fig2, ax2 = plt.subplots(figsize=(8, 4))
    E0 = E_e[0]
    drift_e = (E_e - E0) / abs(E0) * 100.0
    ax2.plot(t_e, drift_e, "m-", lw=1.5)
    ax2.axhline(0, color="k", lw=0.8, ls="--")
    ax2.set_xlabel("Time $t$ [atu]")
    ax2.set_ylabel(r"$(E(t) - E_0)/|E_0|$ [%]")
    ax2.set_title(
        "Energy conservation — non-interacting ETRS\n"
        "(should be < 0.01% for non-interacting system)"
    )
    ax2.grid(True, alpha=0.3)
    out_path2 = os.path.join(RESULTS_DIR, "energy_conservation.pdf")
    fig2.savefig(out_path2, dpi=150, bbox_inches="tight")
    print(f"Saved: {out_path2}")
    plt.close(fig2)

    # ── Print summary ─────────────────────────────────────────────────────────
    print("\n── Summary ──────────────────────────────────────────────")
    print(f"  t range          : {t[0]:.1f} – {t[-1]:.1f} atu")
    print(f"  σ(0) numerical   : {sigma_num[0]:.4f} a₀  (expected {sigma_wp:.4f})")
    print(f"  Max |Δσ/σ|       : {np.max(np.abs(rel_err)):.4f} %")
    print(f"  ΔE max           : {np.max(np.abs(drift_e)):.4e} %")
    print(f"  N_elec (avg)     : {np.mean(N_elec):.6f}  (expected 1.000000)")
    print("─────────────────────────────────────────────────────────\n")

if __name__ == "__main__":
    main()
