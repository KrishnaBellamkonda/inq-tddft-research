#!/usr/bin/env python3
"""Phase 5 (screening) GIF builder — campaign localised-jellium-dynamics-analysis.

From the WP and classical RT frames (at rest, r=12), build per-frame density
differences and render two GIFs:
  total_diff(t)   = n_WP_total(t) − n_CL_total(t)          (includes the WP orbital)
  induced_diff(t) = bath_WP(t) − n_CL(t),  bath_WP = n_WP_total − n_WP_orbital
                    (bath-only screening response; n_GS cancels in the difference)
Renders the y=0 (x–z) slice, shared symmetric colorbar across frames. Uses
inqview.load_vti (physical order — NEVER fftshift a VTI). Robust to missing runs.
"""
import sys, glob, re
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

sys.path.insert(0, "/local/data/public/skcb2/tddft/inq-stack/python")
from inqview import load_vti
try:
    from inqview.visualisation import style as _st; _st.apply()
except Exception: pass

LJ = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium")
DYN = LJ/"scripts/localised_jellium_dynamics"
P5 = DYN/"runs/p5"
OUT = LJ/"hypotheses/localised_jellium_dynamics"; OUT.mkdir(parents=True, exist_ok=True)

def frames(sub):
    fs = sorted(glob.glob(str(P5/sub/"*.vti")))
    return {int(re.search(r"_t(\d+)\.vti", f).group(1)): f for f in fs}

wp_tot = frames("wp/results/wp/frames/total")
wp_orb = frames("wp/results/wp/frames/wp")
cl_tot = frames("cl/results/cl/frames/total")
steps = sorted(set(wp_tot) & set(cl_tot))
if not steps:
    print("Phase 5 GIF: no matching WP/CL frames yet — skipping."); sys.exit(0)

def slice_xz(path):
    v = load_vti(path, expect_centered_axis="z")   # physical order; no fftshift
    ny = v.data.shape[1]
    return v.data[:, ny//2, :], v.x, v.z           # (x,z) at y=0

def build(kind):
    imgs, vmax = [], 0.0
    data = {}
    for s in steps:
        tot_wp, x, z = slice_xz(wp_tot[s])
        tot_cl, _, _ = slice_xz(cl_tot[s])
        if kind == "total":
            diff = tot_wp - tot_cl
        else:  # induced: subtract WP orbital from WP total → bath
            orb, _, _ = slice_xz(wp_orb[s]) if s in wp_orb else (0*tot_wp, x, z)
            diff = (tot_wp - orb) - tot_cl
        data[s] = (diff, x, z)
        vmax = max(vmax, np.percentile(np.abs(diff), 99.5))
    try:
        import imageio.v2 as imageio
    except Exception:
        print("imageio unavailable — saving first/last PNG panels only")
        imageio = None
    for s in steps:
        diff, x, z = data[s]
        fig, ax = plt.subplots(figsize=(6.4, 4.2))
        im = ax.imshow(diff.T, origin="lower", aspect="auto",
                       extent=[x[0], x[-1], z[0], z[-1]], cmap="RdBu_r",
                       vmin=-vmax, vmax=vmax)
        ax.axhspan(-12.5, 12.5, color="0.5", alpha=0.12)
        ax.set_xlabel("x (Bohr)"); ax.set_ylabel("z (Bohr)")
        ax.set_title(f"{kind} density diff  (step {s})")
        fig.colorbar(im, ax=ax, label="Δn (e/Bohr³)")
        fig.tight_layout()
        png = OUT/f"_frame_{kind}_{s:06d}.png"; fig.savefig(png, dpi=110); plt.close(fig)
        imgs.append(str(png))
    if imageio is not None and imgs:
        gif = OUT/f"screening_{kind}.gif"
        with imageio.get_writer(str(gif), mode="I", duration=0.15) as w:
            for p in imgs: w.append_data(imageio.imread(p))
        print("wrote", gif)
    # keep a t=0 vs final static panel
    if imgs:
        print(f"{kind}: {len(imgs)} frames, |Δn|max≈{vmax:.2e}")

for kind in ("total", "induced"):
    try: build(kind)
    except Exception as e: print(f"{kind} GIF failed: {e}")
