#!/usr/bin/env python3
"""Region-resolved bath-drainage analysis for the CAP-in-jellium Baseline 1 runs.

Reads the density_system VTI series for the two B1 runs (eta=-0.5, eta=-0.10),
integrates the bath electron number in the FREE region |z|<15 vs the CAP SLABS
|z|>=15, and the linear profile n(z,t). The decisive question: how long does the
free region survive (the wake window for B2/B3 is ~14-30 a.u.)?

Outputs (this folder):
  cap_b1_region_drainage.csv   t, N_free, N_slab, N_total per run
  fig_b1_region_drainage.png   N_free/N_slab/N_total vs t, both eta
  fig_b1_density_carpet.png    n(z,t) carpet, eta=-0.5

Pure post-processing (VTK + numpy + matplotlib). PROVISIONAL until Task #7.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/local/data/public/skcb2/tddft/inq-stack/python")
import vtk  # noqa: E402
from vtk.util.numpy_support import vtk_to_numpy  # noqa: E402

BASE = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/"
            "scripts/cap_baselines/results")
DT = 0.020
FREE_HALF = 15.0  # |z| < 15 Bohr is the free region; >= is the CAP slab
RUNS = [("b1_eta0p50", "eta=-0.5"), ("b1_eta0p10", "eta=-0.10")]


def read_vti(path):
    """Return (density[nz,ny,nx], origin(x,y,z), spacing(x,y,z))."""
    r = vtk.vtkXMLImageDataReader()
    r.SetFileName(str(path))
    r.Update()
    img = r.GetOutput()
    nx, ny, nz = img.GetDimensions()
    origin = img.GetOrigin()
    spacing = img.GetSpacing()
    a = vtk_to_numpy(img.GetPointData().GetArray(0)).reshape(nz, ny, nx)
    return a, origin, spacing


def frame_list(sub):
    import glob
    import re
    base = BASE / sub / "raw" / "vti" / "density_system"
    out = []
    for f in glob.glob(str(base / "*.vti")):
        m = re.search(r"_t(\d+)\.vti$", f)
        if m:
            out.append((int(m.group(1)), f))
    return sorted(out)


def region_numbers(sub):
    """Per frame: (t, N_free, N_slab, N_total) and the (z, n(z,t)) carpet."""
    frames = frame_list(sub)
    ts, nf, ns, nt = [], [], [], []
    zprof = []
    zc = None
    for step, path in frames:
        d, origin, spacing = read_vti(path)  # d[nz,ny,nx]
        dx, dy, dz = spacing
        nz = d.shape[0]
        if zc is None:
            zc = origin[2] + np.arange(nz) * dz  # physical z per slab index
        # linear density lambda(z) = (sum over x,y) * dx*dy  [electrons/Bohr]
        lam = d.sum(axis=(1, 2)) * dx * dy
        free_mask = np.abs(zc) < FREE_HALF
        N_free = float((lam[free_mask]).sum() * dz)
        N_slab = float((lam[~free_mask]).sum() * dz)
        ts.append(step * DT)
        nf.append(N_free)
        ns.append(N_slab)
        nt.append(N_free + N_slab)
        zprof.append(lam)
    return (np.array(ts), np.array(nf), np.array(ns), np.array(nt),
            zc, np.array(zprof))


def main():
    import csv
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    try:
        from inqview.visualisation import style
        style.apply()
    except Exception:
        pass

    here = Path(__file__).resolve().parent
    data = {}
    for sub, label in RUNS:
        t, nf, ns, ntot, zc, zprof = region_numbers(sub)
        data[sub] = dict(t=t, nf=nf, ns=ns, nt=ntot, zc=zc, zprof=zprof, label=label)
        # survival markers
        N0free = nf[0]
        def t_at(frac):
            below = np.where(nf <= frac * N0free)[0]
            return float(t[below[0]]) if len(below) else float("nan")
        print(f"\n{label}: N_free(0)={N0free:.1f}  "
              f"t(90%)={t_at(0.90):.1f}  t(50%)={t_at(0.50):.1f}  "
              f"t(10%)={t_at(0.10):.1f} a.u.")
        for tt in (5, 10, 15, 20, 30):
            i = int(np.argmin(np.abs(t - tt)))
            print(f"   t={tt:>3} a.u.: N_free={nf[i]:6.1f} ({100*nf[i]/N0free:5.1f}%)  "
                  f"N_slab={ns[i]:6.1f}")

    # combined CSV
    with open(here / "cap_b1_region_drainage.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["run", "t_au", "N_free", "N_slab", "N_total"])
        for sub, _ in RUNS:
            d = data[sub]
            for i in range(len(d["t"])):
                w.writerow([sub, f"{d['t'][i]:.4f}", f"{d['nf'][i]:.4f}",
                            f"{d['ns'][i]:.4f}", f"{d['nt'][i]:.4f}"])

    # figure 1: region drainage
    fig, ax = plt.subplots(1, 1, figsize=(7.0, 4.2))
    colors = {"b1_eta0p50": "C3", "b1_eta0p10": "C0"}
    for sub, _ in RUNS:
        d = data[sub]
        ax.plot(d["t"], d["nf"], color=colors[sub], lw=2,
                label=f"{d['label']}: free |z|<15")
        ax.plot(d["t"], d["ns"], color=colors[sub], lw=1, ls="--",
                label=f"{d['label']}: slab |z|>=15")
    ax.axvspan(0, 30, color="gray", alpha=0.12)
    ax.text(15, ax.get_ylim()[1]*0.92, "wake window\n(~14-30 a.u.)",
            ha="center", va="top", fontsize=8, color="0.3")
    ax.set_xlabel("time (a.u.)")
    ax.set_ylabel("bath electrons in region")
    ax.set_title("Baseline 1: CAP drains the FREE region too, within the wake window")
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(here / "fig_b1_region_drainage.png", dpi=150)
    print(f"\nwrote {here/'fig_b1_region_drainage.png'}")

    # figure 2: n(z,t) carpet for eta=-0.5
    d = data["b1_eta0p50"]
    fig2, ax2 = plt.subplots(1, 1, figsize=(7.0, 4.2))
    im = ax2.pcolormesh(d["zc"], d["t"], d["zprof"], shading="auto", cmap="inferno")
    for zedge in (-15, 15):
        ax2.axvline(zedge, color="cyan", ls="--", lw=1)
    fig2.colorbar(im, ax=ax2, label="linear density n(z) [e/Bohr]")
    ax2.set_xlabel("z (Bohr)")
    ax2.set_ylabel("time (a.u.)")
    ax2.set_title("Baseline 1 eta=-0.5: n(z,t) — bath collapses inward (CAP edges dashed)")
    fig2.tight_layout()
    fig2.savefig(here / "fig_b1_density_carpet.png", dpi=150)
    print(f"wrote {here/'fig_b1_density_carpet.png'}")


if __name__ == "__main__":
    main()
