#!/usr/bin/env python3
"""Q1(vii) / todo-10: final-time xz total-density map (log), with region demarcations.

Last RT frame, mid-y xz slice of density_total, LOG colour scale. Dashed lines mark the
jellium-slab faces (|z|=12.5), the slab<->CAP boundaries / CAP inner edges (|z|=17.5), and
the box edges. WP and classical side by side. VTIs loaded via inqview.load_vti (physical
order, NO fftshift — per the VTI coordinate rule).
Output: qa_vii_xz_logdensity_final.png.
"""
import glob
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from inqview import load_vti
from inqview.visualisation import style

style.apply_theme()
DT = 0.02
ROOT = "/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium"
WP = f"{ROOT}/scripts/fullsuite_wp/results/p5_wp/raw"
CL = f"{ROOT}/scripts/fullsuite_classical/results/p5_classical/raw"
OUT = f"{ROOT}/hypotheses/03_cap_stopping"


def ftime(p):
    return int(re.search(r"_t(\d+)\.vti", p).group(1)) * DT


def final_vf(raw):
    f = sorted(glob.glob(f"{raw}/vti/density_total/density_t*.vti"), key=ftime)[-1]
    return load_vti(f), ftime(f)


vf_wp, t_wp = final_vf(WP)
vf_cl, t_cl = final_vf(CL)

# common log scale across both panels (shared-colorbar rule)
both = np.concatenate([vf_wp.data.ravel(), vf_cl.data.ravel()])
vmax = float(both.max())
vmin = vmax / 1e4
norm = LogNorm(vmin=vmin, vmax=vmax)

fig, axes = plt.subplots(1, 2, figsize=(8.8, 4.6), sharey=True)
for ax, vf, tag, t in ((axes[0], vf_wp, "WP run", t_wp), (axes[1], vf_cl, "classical", t_cl)):
    img = vf.xz_slice(0.0)                                  # (nz, nx), rows = z
    ext = [vf.x.min(), vf.x.max(), vf.z.min(), vf.z.max()]
    im = ax.imshow(np.clip(img, vmin, None), origin="lower", extent=ext,
                   aspect="auto", cmap=style.cmap_for("sequential"), norm=norm)
    for z, c, ls in ((12.5, "cyan", "--"), (-12.5, "cyan", "--"),
                     (17.5, "lime", "--"), (-17.5, "lime", "--"),
                     (25.0, "white", ":"), (-25.0, "white", ":")):
        ax.axhline(z, color=c, ls=ls, lw=1.0)
    ax.set_title(f"{tag} — total density, t={t:.0f} a.u. (log)", fontsize=9)
    ax.set_xlabel("x (Bohr)")
axes[0].set_ylabel("z (Bohr)")
# region labels on the left panel
for zc, lab in ((21.25, "+z CAP"), (15.0, "slab↔CAP"), (0.0, "jellium slab"),
                (-15.0, "slab↔CAP"), (-21.25, "−z CAP")):
    axes[0].annotate(lab, xy=(vf_wp.x.min() + 1, zc), fontsize=6, color="white", va="center")
cb = fig.colorbar(im, ax=axes, fraction=0.046, pad=0.02)
cb.set_label("total electron density (a₀⁻³, log)")
fig.suptitle("Final-time total density (xz, mid-y) — cyan=slab faces ±12.5, lime=CAP edges ±17.5",
             fontsize=8.5)
fig.savefig(f"{OUT}/qa_vii_xz_logdensity_final.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"WP final t={t_wp:.0f}, classical final t={t_cl:.0f}; vmax={vmax:.3e}")
print("wrote qa_vii_xz_logdensity_final.png")
