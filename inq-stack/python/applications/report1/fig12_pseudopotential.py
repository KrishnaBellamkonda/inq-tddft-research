"""fig12 — Pseudopotential modification: ONCV-H original, sign-inverted,
and Nazarov-Gross Gaussian-smoothed reference.

Run:
    python -m applications.report1.fig12_pseudopotential
"""

from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
import gzip

from applications.report1._shared_style import (
    apply_style,
    palette_sweep5,
    column_widths_in,
    TufteCritic,
    panel_label,
)

UPF_PATH = (
    "inq/install/share/pseudopod/pseudopotentials/"
    "quantum-simulation.org/sg15/H_ONCV_PBE-1.2.upf.gz"
)


def load_vloc(upf_gz_path: str):
    """Load V_loc(r) from a gzipped UPF file. Returns (r, V_loc) in Bohr, Hartree."""
    with gzip.open(upf_gz_path, "rt") as f:
        tree = ET.parse(f)
    root = tree.getroot()
    r = np.array(root.find(".//PP_MESH/PP_R").text.split(), dtype=float)
    vloc_ry = np.array(root.find(".//PP_LOCAL").text.split(), dtype=float)
    return r, vloc_ry / 2.0  # Ry → Ha


def main() -> None:
    apply_style()

    r, vloc = load_vloc(UPF_PATH)

    # sign-inverted (electron projectile version)
    vloc_inv = -vloc

    # Nazarov-Gross Gaussian-smoothed reference
    sigma_smooth = 0.5  # Bohr
    A = np.abs(np.interp(2.0, r, vloc))  # match asymptotic -1/r at r=2
    vloc_gauss = -A * np.exp(-r**2 / (2 * sigma_smooth**2))

    # r_cutoff from UPF (rc for l=0 from the input file header)
    r_cutoff = 1.14  # from PP_INPUTFILE: rc = 1.13748

    W = column_widths_in["single"]
    fig, ax = plt.subplots(figsize=(W, W * 0.72))

    mask = r < 5.0

    ax.plot(r[mask], vloc[mask], color=palette_sweep5[4], linewidth=1.0,
            label=r"ONCV-H $V_{\mathrm{loc}}(r)$")
    ax.plot(r[mask], vloc_inv[mask], color=palette_sweep5[0], linewidth=1.0,
            label=r"Sign-inverted $-V_{\mathrm{loc}}(r)$")
    ax.plot(r[mask], vloc_gauss[mask], color="#808080", linewidth=0.9,
            linestyle="--", alpha=0.7,
            label=r"Gaussian $-A\,e^{-r^2/2\sigma^2}$")

    # r_cutoff marker
    ax.axvline(r_cutoff, color="#b0b0b0", linestyle=":", linewidth=0.5)
    ax.text(r_cutoff + 0.08, -3.8, r"$r_c$", fontsize=6.5, color="#808080")

    # zero line
    ax.axhline(0, color="#d0d0d0", linewidth=0.4)

    # mass annotation
    ax.text(2.2, 2.5,
            r"$m_{\mathrm{proj}} = m_e$",
            fontsize=7, color=palette_sweep5[0],
            bbox=dict(facecolor="white", edgecolor="#c0c0c0",
                      linewidth=0.3, pad=2, alpha=0.9))

    ax.set_xlabel(r"$r$ (Bohr)")
    ax.set_ylabel(r"$V(r)$ (Ha)")
    ax.set_xlim(0.01, 5.0)
    ax.set_ylim(-5.0, 5.0)

    ax.legend(fontsize=6, loc="lower right", frameon=True, framealpha=0.9,
              edgecolor="#b0b0b0", handlelength=1.5)

    ax.text(0.97, 0.02, r"\textit{SCHEMATIC}",
            transform=ax.transAxes, fontsize=5, color="#a0a0a0",
            ha="right", va="bottom")

    critic = TufteCritic()
    issues = critic.critique(fig)
    if issues:
        print(f"TufteCritic: {len(issues)} issue(s)")
        for iss in issues:
            print(f"  {iss}")

    out = "docs/reports/report1/figures/fig12_pseudopotential.png"
    fig.savefig(out, dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"Saved → {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
