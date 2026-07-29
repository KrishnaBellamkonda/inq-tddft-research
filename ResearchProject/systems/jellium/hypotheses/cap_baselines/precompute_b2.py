#!/usr/bin/env python3
"""Precompute the heavy B2 (classical projectile) artefacts for cap_baselines.

B2 = CAP + a *classical* sigma=0.5 electron projectile (in `ions`, Ehrenfest),
100 eV, launched z0=-13 moving +z. Because the projectile is a classical ion it
is NOT part of the electron density, so:
  * B2 density is BATH-ONLY (no wavepacket spike), and
  * B2 - B1 is the pure bath response to the projectile's bare Coulomb field.

One strided read pass over the aligned B2 and B1 density_system VTI frames ->
  fig_b2_density_xz.gif   total electronic (bath) density, mid-y xz slice
  fig_b2_wake_xz.gif      projectile-induced bath density (B2 - B1), mid-y slice
  fig_b2_efield_xz.gif    E_z from the B2 total density (FFT-Poisson), mid-y slice

The exact classical projectile z(t) (from electron_track.csv), folded into the
periodic box [-25,25), is overlaid as a marker on every frame.

Shared colour scale: reuses cap_b3_clim.json (written by precompute_b3.py on its
full run) so B1/B2/B3 GIFs of the same quantity share one identical colorbar
gradation (report-figures shared-colorbar rule, production rule 7). Window-aware
via T_MAX/SUFFIX env, exactly like precompute_b3.py.

Done in the builder process (not the notebook kernel). PROVISIONAL until Task #7.
"""
import re
import os
import sys
import csv
import json
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
LBOX = 50.0   # periodic box side (Bohr) for folding the projectile track


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


def load_track():
    """Classical projectile (t_au, z_abs) from electron_track.csv (dedup step 0)."""
    t, z = [], []
    seen = set()
    with open(RES/"b2_classical_E100"/"raw/observables/electron_track.csv") as fh:
        for row in csv.DictReader(fh):
            st = int(row["step"])
            if st in seen:
                continue
            seen.add(st)
            t.append(float(row["time_au"])); z.append(float(row["z"]))
    return np.array(t), np.array(z)


def fold(z):
    """Fold an absolute z into the centred periodic box [-LBOX/2, LBOX/2)."""
    return ((z + LBOX/2) % LBOX) - LBOX/2


def main():
    T_MAX = os.environ.get("T_MAX")
    T_MAX = float(T_MAX) if T_MAX else None
    SUF = os.environ.get("SUFFIX", "")
    b2, b1 = frame_map("b2_classical_E100"), frame_map("b1_eta0p50")
    steps = sorted(set(b2) & set(b1))
    if T_MAX is not None:
        steps = [s for s in steps if s*DT <= T_MAX]
    if len(steps) > MAXF:
        steps = steps[::len(steps)//MAXF + 1]
    print(f"{len(steps)} aligned frames"
          + (f" in window t<=[{T_MAX}] suffix='{SUF}'" if T_MAX else " (full)"))

    # exact projectile z(t), folded, sampled at the chosen frame times
    tt, zt = load_track()
    times = [st*DT for st in steps]
    z_proj = [fold(float(np.interp(t, tt, zt))) for t in times]

    dens, wake, ez = [], [], []
    x = z = None
    midy = None
    for st in steps:
        o, s, a2 = vol(b2[st])
        _, _, a1 = vol(b1[st])
        nz, ny, nx = a2.shape
        if midy is None:
            midy = ny // 2
            x = o[0] + s[0]*np.arange(nx)
            z = o[2] + s[2]*np.arange(nz)
        dens.append(a2[:, midy, :])
        wake.append((a2 - a1)[:, midy, :])
        E = electric_field(np.transpose(a2, (2, 1, 0)), (s[0], s[1], s[2]))
        ez.append(np.transpose(E.ez, (2, 1, 0))[:, midy, :])

    def gif(slabs, out, title, cmap, clim, label="n (e/Bohr$^3$)"):
        # SHARED, FIXED colour scale (clim passed in, from cap_b3_clim.json) so the
        # B1/B2/B3 GIFs of the same quantity share one colorbar gradation, and the
        # full vs transit windows match. See report-figures production rule 7.
        vmin, vmax = clim
        fig, ax = plt.subplots(figsize=(5.4, 4.6))
        im = ax.imshow(slabs[0], origin="lower", aspect="auto",
                       extent=[x[0], x[-1], z[0], z[-1]], vmin=vmin, vmax=vmax, cmap=cmap)
        for ze in (-15, 15):
            ax.axhline(ze, color="0.2", ls="--", lw=1)
        # classical projectile marker (exact, folded into the box)
        pm, = ax.plot([0.0], [z_proj[0]], marker="_", ms=22, mew=2.5,
                      color="lime", ls="none")
        ax.set_xlabel("x (Bohr)"); ax.set_ylabel("z (Bohr)")
        fig.colorbar(im, ax=ax, label=label)
        ttl = ax.set_title("")
        def upd(i):
            im.set_data(slabs[i]); pm.set_data([0.0], [z_proj[i]])
            ttl.set_text(f"{title}   t = {times[i]:6.1f} a.u.")
            return im, pm, ttl
        animation.FuncAnimation(fig, upd, frames=len(slabs), interval=120, blit=False
                                ).save(HERE/out, writer="pillow", dpi=90)
        plt.close(fig); print(f"wrote {out} ({len(slabs)} frames)")

    clim = json.load(open(HERE/"cap_b3_clim.json"))   # reuse the shared B1/B3 scale
    print(f"reusing shared clim: {clim}")

    gif(dens, f"fig_b2_density_xz{SUF}.gif", "B2 bath density (+ classical e$^-$)",
        "inferno", clim["density"])
    gif(wake, f"fig_b2_wake_xz{SUF}.gif", "B2 induced bath density (B2-B1)", "RdBu_r",
        clim["wake"], label=r"$\delta n$ (e/Bohr$^3$)")
    gif(ez, f"fig_b2_efield_xz{SUF}.gif", "B2 $E_z$ (FFT-Poisson)", "RdBu_r",
        clim["efield"], label=r"$E_z$ (a.u.)")


if __name__ == "__main__":
    main()
