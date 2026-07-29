#!/usr/bin/env python3
"""M5 — quantitative coronene diffraction analysis from the 3D ΔP scattering map.

Uses results/analysis/momentum/momentum_scatter_arrays.npz:
  dP3[kz,kx,ky] = |ψ(k)|²_after − |ψ(k)|²_before  (3D momentum-difference map),
  k0 = elastic-shell radius |k|.

Produces (batch2_figures/):
  m5_diffraction_shell_polar.png   ΔP sampled on the elastic shell |k|=k0,
                                    shown as a polar (θ=radius, φ=angle) map —
                                    forward hemisphere; discrete lobes visible.
  m5_diffraction_profiles.png      I(φ) (azimuthal) and I(θ)/I(k_perp) (radial)
                                    quantitative diffraction profiles, with the
                                    graphene/coronene reciprocal-lattice orders
                                    |G| overlaid as predicted diffraction angles.
  m5_diffraction_3d.png            3D render of the strongest positive-ΔP voxels
                                    (scattered intensity directions).

Coronene is graphene-like: a = 2.46 Å = 4.649 a0, C–C = 1.42 Å.
First reciprocal vector |G1| = 4π/(√3 a) = 1.561 a0^-1 → on the |k|=k0 shell a
diffracted lobe sits at θ = asin(|G|/k0).  These are PREDICTIONS overlaid for
you to judge whether the forward lobes are lattice diffraction.

Known-case (printed): shell radius matches k0; ΣΔP ≈ 0 (norm-conserving).
"""
from __future__ import annotations
import numpy as np
from pathlib import Path
from scipy.interpolate import RegularGridInterpolator
import matplotlib; matplotlib.use("Agg")
matplotlib.rcParams["savefig.bbox"] = "standard"
import matplotlib.pyplot as plt

RUN = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/coronene/run_broadening_35x35x80")
OUT = Path("/local/data/public/skcb2/tddft/docs/presentations/storyline/tasks/batch2_figures")
A_GRAPHENE = 2.46 / 0.529177210903       # a0
G1 = 4 * np.pi / (np.sqrt(3) * A_GRAPHENE)   # 1.561 a0^-1 (first reciprocal vector)
G_ORDERS = [("|G₁|", G1), ("√3|G₁|", np.sqrt(3) * G1), ("2|G₁|", 2 * G1)]

d = np.load(RUN / "results/analysis/momentum/momentum_scatter_arrays.npz")
kz, kx, ky, dP3 = d["kz"], d["kx"], d["ky"], d["dP3"].astype(float)
k0 = float(d["k0"])
print(f"M5: dP3{dP3.shape} k0={k0:.3f}  ΣΔP={dP3.sum():.2e} (≈0 norm-conserving)")
print(f"    graphene a={A_GRAPHENE:.3f} a0  |G1|={G1:.3f} a0^-1  "
      f"-> shell angle θ(|G1|)={np.degrees(np.arcsin(min(G1/k0,1))):.1f}°")

interp = RegularGridInterpolator((kz, kx, ky), dP3, bounds_error=False, fill_value=0.0)

# ---- sample on the elastic shell |k|=k0 over (theta, phi), forward hemisphere ----
nth, nph = 90, 180
th = np.linspace(0, np.pi / 2, nth)          # 0=forward(+kz) .. 90 deg = equator
ph = np.linspace(-np.pi, np.pi, nph)
TH, PH = np.meshgrid(th, ph, indexing="ij")
KZ = k0 * np.cos(TH); KX = k0 * np.sin(TH) * np.cos(PH); KY = k0 * np.sin(TH) * np.sin(PH)
shell = interp(np.stack([KZ.ravel(), KX.ravel(), KY.ravel()], axis=1)).reshape(TH.shape)

# ---- Fig 1: polar map (theta radius, phi angle), forward hemisphere ----
fig = plt.figure(figsize=(7.5, 6.4)); ax = fig.add_subplot(111, projection="polar")
vmax = np.percentile(np.abs(shell), 99.5)
pc = ax.pcolormesh(ph, np.degrees(th), shell, cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="auto")
ax.set_theta_zero_location("E"); ax.set_rlabel_position(135)
ax.set_title(f"M5 — ΔP on elastic shell |k|=k0={k0:.2f}  (θ=radius°, φ=azimuth)\n"
             "forward hemisphere; red=gain (scattered into), blue=loss", fontsize=10)
for lbl, G in G_ORDERS:
    if G < k0:
        ax.plot(ph, np.full_like(ph, np.degrees(np.arcsin(G / k0))), "k--", lw=0.8)
fig.colorbar(pc, ax=ax, shrink=0.7, label="ΔP (shell)")
fig.tight_layout(); fig.savefig(OUT / "m5_diffraction_shell_polar.png", dpi=150); plt.close(fig)
print(f"wrote {OUT/'m5_diffraction_shell_polar.png'}")

# ---- Fig 2: I(phi) and I(theta)/I(kperp) profiles ----
solid = np.sin(TH)                          # shell area element ∝ sinθ
I_phi = (np.clip(shell, 0, None) * solid).sum(axis=0)      # gain only, azimuthal
I_th = (np.clip(shell, 0, None) * solid).sum(axis=1)       # gain only, radial
kperp = k0 * np.sin(th)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 4.6))
a1.plot(np.degrees(ph), I_phi, color="C0", lw=1.6)
for a in (0, 60, 120, 180, -60, -120):
    a1.axvline(a, color="0.8", lw=0.8)
for a in (45, 135, -45, -135):
    a1.axvline(a, color="orange", ls=":", lw=1.0)
a1.set_xlabel("azimuth φ (deg)"); a1.set_ylabel("I(φ) gain (Σ over θ, sinθ-wtd)")
a1.set_title("Azimuthal diffraction profile\n(grey=6-fold 0/60/120; orange=observed ±45/±135)")
a1.grid(alpha=0.3)
a2.plot(kperp, I_th, color="C3", lw=1.6)
a2b = a2.twiny(); a2b.set_xlim(0, np.degrees(np.arcsin(min(kperp.max(), k0) / k0)))
for lbl, G in G_ORDERS:
    if G < k0:
        a2.axvline(G, color="k", ls="--", lw=1.0)
        a2.text(G, a2.get_ylim()[1] * 0.92, lbl, rotation=90, va="top", fontsize=8)
a2.set_xlabel("transverse momentum k$_\\perp$ = k0 sinθ (a0$^{-1}$)")
a2.set_ylabel("I(k$_\\perp$) gain (Σ over φ, sinθ-wtd)")
a2.set_title("Radial diffraction profile\n(dashed = graphene reciprocal orders |G|)")
a2.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(OUT / "m5_diffraction_profiles.png", dpi=150); plt.close(fig)
print(f"wrote {OUT/'m5_diffraction_profiles.png'}")

# ---- Fig 3: 3D render of strongest positive-ΔP voxels ----
KZc, KXc, KYc = np.meshgrid(kz, kx, ky, indexing="ij")
flat = dP3.ravel()
thr = np.percentile(flat, 99.7)
sel = flat > thr
fig = plt.figure(figsize=(7.5, 6.5)); ax = fig.add_subplot(111, projection="3d")
sc = ax.scatter(KXc.ravel()[sel], KYc.ravel()[sel], KZc.ravel()[sel],
                c=flat[sel], cmap="hot", s=12, alpha=0.7)
# elastic sphere wire
u, v = np.mgrid[0:2*np.pi:24j, 0:np.pi:12j]
ax.plot_wireframe(k0*np.cos(u)*np.sin(v), k0*np.sin(u)*np.sin(v), k0*np.cos(v),
                  color="0.7", lw=0.3, alpha=0.4)
ax.set_xlabel("kx"); ax.set_ylabel("ky"); ax.set_zlabel("kz (beam)")
ax.set_title("M5 — strongest scattered-into momentum directions (top 0.3% ΔP)\n"
             "grey wire = elastic shell |k|=k0")
fig.colorbar(sc, ax=ax, shrink=0.6, label="ΔP gain")
fig.tight_layout(); fig.savefig(OUT / "m5_diffraction_3d.png", dpi=150); plt.close(fig)
print(f"wrote {OUT/'m5_diffraction_3d.png'}")
