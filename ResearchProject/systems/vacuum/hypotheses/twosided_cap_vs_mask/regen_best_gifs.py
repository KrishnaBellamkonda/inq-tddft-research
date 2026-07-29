#!/usr/bin/env python3
"""FB-004 — regenerate the best-CAP / best-mask density GIFs with the corrected
y-axis label (n_WP, was an unrendered '⟂' glyph). Minimal standalone replica of
the gif cell in build_twosided_report.py (which only ASSEMBLES the notebook, so
the cell never executes on a plain run). Source of truth for the labels.
"""
from __future__ import annotations
import re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import vtk
from vtk.util.numpy_support import vtk_to_numpy

HERE = Path(__file__).resolve().parent          # hypotheses/ dir (gifs land here)
SWEEP = HERE.parent.parent / "twosided_cap_vs_mask"   # actual run dirs (systems/vacuum/...)
ANCHOR_E, ETA_STAR = 10, -0.50


def parse(p):
    out = {}
    for ln in Path(p).read_text().splitlines():
        k, _, v = ln.partition(" ")
        try:
            out[k] = float(v)
        except ValueError:
            out[k] = v
    return out


def zprof(p):
    rd = vtk.vtkXMLImageDataReader(); rd.SetFileName(str(p)); rd.Update()
    img = rd.GetOutput(); dm = img.GetDimensions()
    a = vtk_to_numpy(img.GetPointData().GetArray(0)).reshape(dm[2], dm[1], dm[0])
    sp = img.GetSpacing(); org = img.GetOrigin()
    z = org[2] + sp[2] * np.arange(dm[2])
    return z, a.sum(axis=(1, 2))


recs = []
for mode in ("cap", "mask"):
    for d in sorted(SWEEP.glob(f"run_{mode}_*")):
        f = d / "results/epsilon.txt"
        if f.exists():
            r = parse(f); r["name"] = d.name; recs.append(r)
cap = [r for r in recs if r["mode"] == "cap"]
mask = [r for r in recs if r["mode"] == "mask"]
capL = [r for r in cap if abs(r["eta_Ha"] - ETA_STAR) < 1e-6]


def best(data):
    pool = [r for r in data if abs(r["E_eV"] - ANCHOR_E) < 0.6] or data
    return min(pool, key=lambda r: r["epsilon"]) if pool else None


def sci(x):
    """x10^n mathtext (FB-003) for the gif title."""
    if x == 0:
        return "0"
    import math
    e = int(math.floor(math.log10(abs(x))))
    return rf"{x/10**e:.1f}\times10^{{{e}}}"


for data, out, lab in [(capL, "fig_best_cap_density.gif", "CAP"),
                       (mask, "fig_best_mask_density.gif", "mask")]:
    sc = best(data)
    if not sc:
        print(f"{lab}: no run"); continue
    d = SWEEP / sc["name"]
    vti = d / "results/raw/vti/density_wp"
    frames = sorted(vti.glob("density_wp_t*.vti"),
                    key=lambda p: int(re.search(r"_t(\d+)", p.name).group(1)))
    if not frames:
        print(f"{lab}: no frames"); continue
    z0, _ = zprof(frames[0])
    profs = [zprof(p)[1] for p in frames]
    ymax = max(p.max() for p in profs) * 1.05
    r = parse(d / "results/epsilon.txt")
    zin, Lh = r["z_in"], r["Lhalf"]
    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    (line,) = ax.plot(z0, profs[0])
    for sgn in (+1, -1):
        ax.axvspan(sgn * zin, sgn * (zin + Lh), color="C3", alpha=0.15)
    ax.set_xlabel("z (Bohr)")
    ax.set_ylabel(r"$n_{\rm WP}$")                      # FB-004
    ax.set_ylim(0, ymax)
    ax.set_title(rf'best {lab}: L={int(r["L_total"])}  $\varepsilon={sci(r["epsilon"])}$')
    ttl = ax.text(0.02, 0.92, "", transform=ax.transAxes, fontsize=8)

    def upd(k, line=line, profs=profs, ttl=ttl):
        line.set_ydata(profs[k]); ttl.set_text(f"{k+1}/{len(profs)}")
        return line, ttl

    animation.FuncAnimation(fig, upd, frames=len(profs), interval=120,
                            blit=False).save(HERE / out, writer="pillow", dpi=90)
    plt.close(fig)
    print(f"wrote {out}  ({(HERE/out).stat().st_size} bytes, {len(profs)} frames)")
