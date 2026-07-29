#!/usr/bin/env python3
"""Phase 1 (vacuum exit) dashboard for the classical-highdensity-sv campaign.

Reads results/exit_scan.csv + nproj_z*.vti from the vac_exit run and produces:
  1. integral(n_proj) vs z_center   (should decay through +42.5, no secondary rise)
  2. wrap-witness (max near-face density) vs z_center  (must stay ~0)
  3. montage of n_proj x-z slices at z_center = 40, 42.5, 45, 48

Loads VTIs via inqview.load_vti (physical order; NEVER fftshift). PNG only.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from inqview import load_vti

RUN = Path(
    "/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
    "scripts/classical_highdensity_sv/vac_exit/results"
)
OUT = Path(
    "/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
    "hypotheses/classical_highdensity_sv/phase1_vac_exit"
)
OUT.mkdir(parents=True, exist_ok=True)

FAR_FACE = 42.5   # +Lz/2

try:
    from inqview.visualisation import style
    style.apply_theme()
except Exception:
    pass


def vti_name(zc: float) -> str:
    # matches the C++ formatting: nproj_z%+06.1f with '.'->'p'
    s = f"nproj_z{zc:+06.1f}".replace(".", "p")
    return s + ".vti"


def main() -> int:
    df = pd.read_csv(RUN / "exit_scan.csv").sort_values("z_center").reset_index(drop=True)
    print(df.to_string(index=False))

    # ---- Panel 1: integral vs z_center ----
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(df.z_center, df.integral, "o-", color="C0")
    ax.axvline(FAR_FACE, color="k", ls="--", lw=1)
    ax.axhline(1.0, color="0.6", ls=":", lw=1)
    ax.axhline(0.5, color="0.6", ls=":", lw=1)
    ax.annotate("far face +42.5", xy=(FAR_FACE, 0.5), xytext=(FAR_FACE - 11, 0.55),
                fontsize=9)
    ax.set_xlabel("projectile z-center (Bohr)")
    ax.set_ylabel(r"$\int n_{\rm proj}\,dV$")
    ax.set_title("Gaussian charge clipped at the +z open face (no wrap)")
    fig.tight_layout()
    fig.savefig(OUT / "integral_vs_zcenter.png", dpi=150)
    plt.close(fig)

    # ---- Panel 2: wrap witness ----
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(df.z_center, df.wrap_witness_max, "s-", color="C3")
    ax.axvline(FAR_FACE, color="k", ls="--", lw=1)
    ax.set_xlabel("projectile z-center (Bohr)")
    ax.set_ylabel("max density in near face (z < -38)")
    ax.set_title("Wrap witness: must stay ~0 (any bump = wraparound = FAIL)")
    peak_gauss = float(df.integral.max())  # order-of-magnitude ref only
    fig.tight_layout()
    fig.savefig(OUT / "wrap_witness_vs_zcenter.png", dpi=150)
    plt.close(fig)

    # ---- Panel 3: montage of x-z slices ----
    targets = [40.0, 42.5, 45.0, 48.0]
    fig, axes = plt.subplots(1, len(targets), figsize=(4 * len(targets), 5),
                             constrained_layout=True)
    # shared colour scale from the deepest-interior frame present
    vmax_ref = 0.0
    loaded = {}
    for zc in targets:
        p = RUN / vti_name(zc)
        if p.exists():
            f = load_vti(p)
            loaded[zc] = f
            vmax_ref = max(vmax_ref, float(f.data.max()))
    for ax, zc in zip(axes, targets):
        p = RUN / vti_name(zc)
        if zc not in loaded:
            ax.text(0.5, 0.5, f"missing\n{p.name}", ha="center", va="center",
                    transform=ax.transAxes)
            ax.set_title(f"z_center = {zc}")
            continue
        f = loaded[zc]
        sl = f.xz_slice(y=0.0)  # (nz, nx): rows=z, cols=x
        extent = [f.x[0], f.x[-1], f.z[0], f.z[-1]]
        im = ax.imshow(sl, origin="lower", extent=extent, aspect="auto",
                       cmap="viridis", vmin=0.0, vmax=vmax_ref)
        ax.axhline(FAR_FACE, color="w", ls="--", lw=1)
        ax.axhline(-FAR_FACE, color="w", ls="--", lw=1)
        ax.set_title(f"z_center = {zc}")
        ax.set_xlabel("x (Bohr)")
        if ax is axes[0]:
            ax.set_ylabel("z (Bohr)")
    fig.colorbar(im, ax=axes, shrink=0.8, label=r"$n_{\rm proj}$")
    fig.suptitle("n_proj x-z slice (mid-y): clipping at +z face, nothing at -z face")
    fig.savefig(OUT / "montage_xz_slices.png", dpi=150)
    plt.close(fig)

    # ---- verdict numbers ----
    wrap_peak = float(df.wrap_witness_max.max())
    # integral at deep interior (reference full norm)
    interior = df[df.z_center <= 0.0]
    interior_int = float(interior.integral.mean()) if len(interior) else float("nan")
    # integral at exactly the face
    face_row = df.iloc[(df.z_center - FAR_FACE).abs().argmin()]
    print(f"\ninterior integral (z<=0 mean) = {interior_int:.6g}")
    print(f"integral at z~{face_row.z_center} = {face_row.integral:.6g}")
    print(f"wrap-witness peak (over ALL z_center) = {wrap_peak:.6g}")
    print(f"dashboard written to {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
