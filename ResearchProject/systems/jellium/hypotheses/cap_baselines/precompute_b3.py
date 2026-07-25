#!/usr/bin/env python3
"""Precompute the heavy B3 (wavepacket) artefacts for the cap_baselines notebook.

One strided read pass over the aligned B3 and B1 density_system VTI frames →
  fig_b3_density_xz.gif   total electronic density, mid-y xz slice, animated
  fig_b3_wake_xz.gif      WP-induced density  (B3 − B1),  mid-y xz slice
  fig_b3_efield_xz.gif    E_z from the B3 total density (FFT-Poisson), mid-y slice
  cap_b3_wp_centroid.csv  WP centroid z(t) from the positive part of the wake

Done in the builder process (not the notebook kernel) so the executed notebook
stays fast. PROVISIONAL until Task #7.
"""
import re
import sys
import csv
from pathlib import Path

import numpy as np

sys.path.insert(0, "/local/data/public/skcb2/tddft/inq-stack/python")
import vtk  # noqa: E402
from vtk.util.numpy_support import vtk_to_numpy  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.animation as animation  # noqa: E402
from inqview.analysis.efield import electric_field  # noqa: E402

HERE = Path(__file__).resolve().parent
RES = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/"
           "scripts/cap_baselines/results")
DT = 0.02
MAXF = 100


def vol(path):
    r = vtk.vtkXMLImageDataReader(); r.SetFileName(str(path)); r.Update()
    img = r.GetOutput(); nx, ny, nz = img.GetDimensions()
    o, s = img.GetOrigin(), img.GetSpacing()
    a = vtk_to_numpy(img.GetPointData().GetArray(0)).reshape(nz, ny, nx)  # [nz,ny,nx]
    return o, s, a


def frame_map(sub):
    d = RES / sub / "raw" / "vti" / "density_system"
    out = {}
    for f in d.glob("density_t*.vti"):
        out[int(re.search(r"_t(\d+)", f.name).group(1))] = f
    return out


def main():
    import os
    T_MAX = os.environ.get("T_MAX")
    T_MAX = float(T_MAX) if T_MAX else None
    SUF = os.environ.get("SUFFIX", "")
    b3, b1 = frame_map("b3_wp_E100"), frame_map("b1_eta0p50")
    steps = sorted(set(b3) & set(b1))
    if T_MAX is not None:
        steps = [s for s in steps if s*DT <= T_MAX]   # restrict to the window first
    if len(steps) > MAXF:
        steps = steps[::len(steps)//MAXF + 1]
    print(f"{len(steps)} aligned frames"
          + (f" in window t<=[{T_MAX}] suffix='{SUF}'" if T_MAX else " (full)"))

    dens, wake, ez = [], [], []
    cod_rows = []
    x = z = None
    midy = None
    for st in steps:
        o, s, a3 = vol(b3[st])
        _, _, a1 = vol(b1[st])
        nz, ny, nx = a3.shape
        if midy is None:
            midy = ny // 2
            x = o[0] + s[0]*np.arange(nx)
            z = o[2] + s[2]*np.arange(nz)
        dens.append(a3[:, midy, :])             # [nz,nx]
        w = a3 - a1                              # WP-induced density (3D)
        wake.append(w[:, midy, :])
        # E_z from the B3 total density (efield wants [nx,ny,nz])
        E = electric_field(np.transpose(a3, (2, 1, 0)), (s[0], s[1], s[2]))
        ez.append(np.transpose(E.ez, (2, 1, 0))[:, midy, :])
        # WP centroid from the positive part of the wake (robust)
        wp = np.maximum(w, 0.0)
        wsum = wp.sum()
        zc = float((wp.sum(axis=(1, 2)) * z).sum() / wsum) if wsum > 0 else float("nan")
        cod_rows.append((st*DT, zc, float(wsum)))
    times = [st*DT for st in steps]

    with open(HERE/f"cap_b3_wp_centroid{SUF}.csv", "w", newline="") as fh:
        wri = csv.writer(fh); wri.writerow(["t_au", "wp_centroid_z", "wp_weight"])
        wri.writerows(cod_rows)
    print(f"wrote cap_b3_wp_centroid{SUF}.csv")

    def gif(slabs, out, title, cmap, clim, label="n (e/Bohr$^3$)"):
        # SHARED, FIXED colour scale: clim=(vmin,vmax) is passed in (computed once
        # from the full window) so full vs transit (and B1 vs B3) GIFs of the same
        # quantity use an identical colorbar gradation. See the shared-colorbar rule.
        vmin, vmax = clim
        fig, ax = plt.subplots(figsize=(5.4, 4.6))
        im = ax.imshow(slabs[0], origin="lower", aspect="auto",
                       extent=[x[0], x[-1], z[0], z[-1]], vmin=vmin, vmax=vmax, cmap=cmap)
        for ze in (-15, 15):
            ax.axhline(ze, color="0.2", ls="--", lw=1)
        ax.set_xlabel("x (Bohr)"); ax.set_ylabel("z (Bohr)")
        fig.colorbar(im, ax=ax, label=label)
        ttl = ax.set_title("")
        def upd(i):
            im.set_data(slabs[i]); ttl.set_text(f"{title}   t = {times[i]:6.1f} a.u.")
            return im, ttl
        animation.FuncAnimation(fig, upd, frames=len(slabs), interval=120, blit=False
                                ).save(HERE/out, writer="pillow", dpi=90)
        plt.close(fig); print(f"wrote {out} ({len(slabs)} frames)")

    # Shared clim: compute ONCE from the full window, persist, and reuse for the
    # transit window (and the B1 density GIF) so every comparable GIF matches.
    import json
    clim_file = HERE/"cap_b3_clim.json"
    if T_MAX is None:
        d, w, e = np.array(dens), np.abs(np.array(wake)), np.abs(np.array(ez))
        clim = {"density": [0.0, float(np.percentile(d, 99.5))],
                "wake":   [-float(np.percentile(w, 99.0)), float(np.percentile(w, 99.0))],
                "efield": [-float(np.percentile(e, 99.0)), float(np.percentile(e, 99.0))]}
        json.dump(clim, open(clim_file, "w"), indent=2)
        print(f"wrote {clim_file}: {clim}")
    else:
        clim = json.load(open(clim_file))   # reuse the full-window scale
        print(f"reusing clim from {clim_file}")

    gif(dens, f"fig_b3_density_xz{SUF}.gif", "B3 total density", "inferno", clim["density"])
    gif(wake, f"fig_b3_wake_xz{SUF}.gif", "B3 WP-induced density (B3-B1)", "RdBu_r",
        clim["wake"], label=r"$\delta n$ (e/Bohr$^3$)")
    gif(ez, f"fig_b3_efield_xz{SUF}.gif", "B3 $E_z$ (FFT-Poisson)", "RdBu_r",
        clim["efield"], label=r"$E_z$ (a.u.)")


if __name__ == "__main__":
    main()
