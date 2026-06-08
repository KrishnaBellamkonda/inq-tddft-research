#!/usr/bin/env python3
"""T-b: WP-vs-classical wake comparability at E=20 eV — uses BATH = total - wp.

Verified by integration: density_rt_system and density_rt_total both = 163 e
(162 bath + 1 WP); density_wp = 1 e; (total - wp) = 162 e = true bath.
So the bath-only response MUST subtract the WP orbital. We use total - wp.

Panels:
  A  bath-only induced z-profile (system: total-wp), WP vs classical, matched t.
  B  WP-in-z-profile contamination: standard delta (total, WP IN) vs bath-only.
  C  z_system profile = (total - wp)(z) at several t   [absolute bath line density]
  D  delta z_profile  = z_system(t) - z_system(0)      [induced bath line density]

Bath for classical run = density_rt_total directly (ion projectile, not in e-density).
Baseline = first RT frame. Known-case: bath induced at t0 == 0 everywhere.
Output: tb_wake_comparability_E20.png
"""
import glob, re, numpy as np, vtk
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib; matplotlib.use("Agg")
matplotlib.rcParams["savefig.bbox"] = "standard"
import matplotlib.pyplot as plt

JB = "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium"
WP = f"{JB}/run_wp_n162_L50_E20_sigma1_v2"; DT_WP = 0.01
CL = f"{JB}/run_classical_n162_L50_E20";    DT_CL = 0.02
OUT = f"{JB}/_compare_sigma1_E100/tb_wake_comparability_E20.png"
TARGETS = [2.0, 4.0, 6.0, 8.0, 10.0]

def zline(path):
    r = vtk.vtkXMLImageDataReader(); r.SetFileName(path); r.Update()
    img = r.GetOutput(); nx, ny, nz = img.GetDimensions()
    oz, sz = img.GetOrigin()[2], img.GetSpacing()[2]
    sx, sy = img.GetSpacing()[0], img.GetSpacing()[1]
    a = vtk_to_numpy(img.GetPointData().GetArray(0)).reshape(nz, ny, nx)
    return oz + sz * np.arange(nz), a.sum(axis=(1, 2)) * sx * sy

def frames(run, sub):
    out = []
    for f in glob.glob(f"{run}/results/raw/vti/{sub}/*.vti"):
        m = re.search(r"_t(\d+)\.vti$", f)
        if m: out.append((int(m.group(1)), f))
    return sorted(out)

def at_t(run, sub, dt, t):
    fr = frames(run, sub)
    if not fr: return None
    steps = np.array([s for s, _ in fr]); i = int(np.argmin(np.abs(steps * dt - t)))
    return fr[i][1], steps[i] * dt, steps[i]

def wp_bath_z(t):
    """WP run bath line density = total - wp at matched step."""
    ft, tt, st = at_t(WP, "density_total", DT_WP, t)
    # find wp frame at same step (wp saved less often) -> nearest
    z, ptot = zline(ft)
    fw = frames(WP, "density_wp"); ws = np.array([s for s, _ in fw])
    j = int(np.argmin(np.abs(ws - st))); _, pwp = zline(fw[j][1])
    return z, ptot - pwp, tt

# baselines
zb, bath0, _ = wp_bath_z(0.0)
fc0, _, _ = at_t(CL, "density_total", DT_CL, 0.0); zc, cl0 = zline(fc0)
ft0, _, _ = at_t(WP, "density_total", DT_WP, 0.0); _, tot0 = zline(ft0)

# known-case
_, bchk, _ = wp_bath_z(0.0)
print(f"[known-case] WP bath induced at t0: max|.|={np.abs(bchk-bath0).max():.2e} (==0)")

fig, axs = plt.subplots(2, 2, figsize=(12, 8))
axA, axB, axC, axD = axs[0, 0], axs[0, 1], axs[1, 0], axs[1, 1]
cols = plt.cm.viridis(np.linspace(0.1, 0.9, len(TARGETS)))
wp_pk, cl_pk = [], []
for c, t in zip(cols, TARGETS):
    z, bw, tw = wp_bath_z(t); ind_wp = bw - bath0
    axA.plot(z, ind_wp, color=c, lw=1.6, label=f"WP {tw:.0f}")
    wp_pk.append(np.abs(ind_wp).max())
    fc, tc, _ = at_t(CL, "density_total", DT_CL, t); zc, sc = zline(fc); ind_cl = sc - cl0
    axA.plot(zc, ind_cl, color=c, lw=1.1, ls="--")
    cl_pk.append(np.abs(ind_cl).max())
    axC.plot(z, bw, color=c, lw=1.4, label=f"t={tw:.0f}")           # absolute bath
    axD.plot(z, ind_wp, color=c, lw=1.4, label=f"t={tw:.0f}")        # induced bath
axA.axhline(0, color="0.6", lw=0.6); axA.set_title("A: bath induced Δn(z), WP (solid) vs classical (dashed)")
axA.set_xlabel("z (Bohr)"); axA.set_ylabel("Δn(z) (e/Bohr)"); axA.legend(fontsize=6, ncol=2); axA.grid(alpha=0.3)

# Panel B contamination at t=8
TB = 8.0
z, bw, tw = wp_bath_z(TB); ind_bath = bw - bath0
ftb, _, _ = at_t(WP, "density_total", DT_WP, TB); _, tot_tb = zline(ftb); ind_std = tot_tb - tot0
axB.plot(z, ind_bath, color="#0072B2", lw=1.8, label="bath only (total−wp)")
axB.plot(z, ind_std, color="#D55E00", lw=1.3, ls="--", label="standard Δ (total, WP IN)")
axB.axhline(0, color="0.6", lw=0.6); axB.set_title(f"B: WP-in-z-profile check (t={tw:.0f})")
axB.set_xlabel("z (Bohr)"); axB.set_ylabel("Δn(z) (e/Bohr)"); axB.legend(fontsize=7); axB.grid(alpha=0.3)
ratio = np.abs(ind_std).max() / max(np.abs(ind_bath).max(), 1e-12)
axB.text(0.02, 0.02, f"std/bath peak ~ {ratio:.0f}x", transform=axB.transAxes, fontsize=8,
         bbox=dict(boxstyle="round", fc="white", ec="0.6"))

axC.axhline(0, color="0.6", lw=0.6); axC.set_title("C: z_system profile = (total−wp)(z), WP run")
axC.set_xlabel("z (Bohr)"); axC.set_ylabel("bath n(z) (e/Bohr)"); axC.legend(fontsize=7); axC.grid(alpha=0.3)
axD.axhline(0, color="0.6", lw=0.6); axD.set_title("D: delta z_profile = z_system(t) − z_system(0), WP run")
axD.set_xlabel("z (Bohr)"); axD.set_ylabel("Δn(z) (e/Bohr)"); axD.legend(fontsize=7); axD.grid(alpha=0.3)

fig.tight_layout(); fig.savefig(OUT, dpi=150); plt.close(fig)
print("WP bath wake peak per t:", [f"{v:.2e}" for v in wp_pk])
print("classical wake peak per t:", [f"{v:.2e}" for v in cl_pk])
print(f"contamination ratio at t={TB}: {ratio:.1f}x")
print(f"wrote {OUT}")
