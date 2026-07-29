#!/usr/bin/env python3
"""build_grazing_report.py — grazing / impact-parameter study notebook.

A FINITE coronene C24H12 graphene flake reoriented into the y-z plane (flake
normal = x); a +z projectile grazes it at impact parameter b = x-offset. This
builder is PARTIAL-TOLERANT: it renders whatever runs exist under
``graphene/grazing/run_{cl,wp}_b*`` and is re-run by the dispatcher at batch end.

Reuses the validated kernels inqview.analysis.{planar_density,diffraction} and the
canonical visualisation theme. All CAP results PROVISIONAL until inq-study Task #7.

Run (venv): PYTHONPATH=.../inq-stack/python python3 build_grazing_report.py
"""
from __future__ import annotations

import glob
import os

import nbformat as nbf
from nbconvert.preprocessors import ExecutePreprocessor
from nbformat.v4 import new_code_cell, new_markdown_cell, new_notebook

HERE = os.path.dirname(os.path.abspath(__file__))
SYS = os.path.normpath(os.path.join(HERE, "..", ".."))   # systems/graphene
FIG = os.path.join(HERE, "figs"); os.makedirs(FIG, exist_ok=True)

C: list = []
def md(s): C.append(new_markdown_cell(s))
def code(s): C.append(new_code_cell(s))

md(r"""# Grazing / impact-parameter scattering off a finite graphene flake

**System.** A finite, H-passivated graphene flake (**coronene C$_{24}$H$_{12}$**)
reoriented into the **y–z plane** (flake normal = $x$). A projectile travels
along **+z**, *parallel* to the sheet, and grazes it at **impact parameter
$b$ = perpendicular $x$-offset**. Two-sided sin² CAP on the $z$-faces ($L=20$)
absorbs the projectile + scattered flux. Ground state: closed-shell,
HOMO–LUMO gap $\approx2.76$ eV (LDA), $E=-150.77$ Ha, 108 e.

**Why a finite flake (not periodic bulk).** A beam *parallel* to an infinite
periodic sheet can never enter/exit it, and graphene's lattice is incommensurate
with the 60-Bohr vacuum needed for the CAP. A finite flake is the natural
"fly-past" geometry and is what the user specified.

**Projectiles.** *Classical* (Gaussian $-1$ ion, the now-fixed `z_valence=-1`
He-symbol UPF — clean trajectory, no dispersion) and *wave packet* (disperses
over the long parallel path — included for comparison). $E=100$ eV, $\sigma=1.47$
Bohr. Scan $b\in\{1,3,6\}$ Bohr.

> **PROVISIONAL** — every CAP number awaits the inq-study engine regression
> (Task #7). Feasibility replica (reduced flake/ensemble), not converged values.""")

code(f"""import os, glob, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from IPython.display import Image, display
SYS = {SYS!r}; FIG = {FIG!r}
HA = 27.211386245988
try:
    from inqview.visualisation import style as _S
    _S.apply_theme() if hasattr(_S, "apply_theme") else None
except Exception as e:
    print("theme not applied:", e)

def read_summary(p):
    d = {{}}
    if not os.path.exists(p): return d
    for ln in open(p):
        if "=" in ln:
            k, v = ln.split("=", 1); d[k.strip()] = v.strip()
    return d

def read_cols(p):
    if not os.path.exists(p): return {{}}
    import csv
    rows = list(csv.DictReader(open(p)))
    if not rows: return {{}}
    out = {{k: [] for k in rows[0]}}
    for r in rows:
        for k, v in r.items():
            try: out[k].append(float(v))
            except (ValueError, TypeError): out[k].append(np.nan)
    return {{k: np.array(v) for k, v in out.items()}}

def b_of(name):
    import re
    m = re.search(r"_b([0-9p.]+)", name);
    return float(m.group(1).replace("p", ".")) if m else np.nan""")

# --- classical: KE stopping vs impact parameter b ---
md(r"""## 1. Classical projectile — energy loss vs impact parameter $b$

Projectile kinetic energy from `electron_track.csv`; **stopping = $KE(t_0)-KE(t_f)$**
(the run.cpp now records `KE_initial_eV` from the actual launch $v_0$). As $b$
decreases the projectile grazes denser regions of the flake $\Rightarrow$ more
energy loss expected.""")

code(r"""cls = sorted(glob.glob(os.path.join(SYS, "grazing", "run_cl_b*")))
pts = []
fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
for d in cls:
    tk = read_cols(os.path.join(d, "results/raw/observables/electron_track.csv"))
    if not tk or "vz" not in tk or len(tk["vz"]) < 3: continue
    ke = 0.5 * (tk["vx"]**2 + tk["vy"]**2 + tk["vz"]**2) * HA
    b = b_of(os.path.basename(d))
    ax[0].plot(tk["time_au"], ke, lw=1.4, label=f"b={b:g} Bohr")
    pts.append((b, float(ke[0] - ke[-1])))
ax[0].set(xlabel="time (a.u.)", ylabel="projectile KE (eV)", title="Classical KE(t) by b")
ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)
if pts:
    pts.sort()
    bb = [p[0] for p in pts]; dd = [p[1] for p in pts]
    ax[1].plot(bb, dd, "o-", color="#b03030")
    ax[1].set(xlabel="impact parameter b (Bohr)", ylabel="ΔKE (eV)",
              title="Stopping vs impact parameter")
    ax[1].grid(alpha=0.3)
    print("b (Bohr)   ΔKE (eV)")
    for b, dl in pts: print(f"  {b:6.2f}   {dl:8.3f}")
else:
    ax[1].text(0.5, 0.5, "no classical grazing runs yet", ha="center")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_grazing_stopping_vs_b.png"), dpi=130)
plt.show()""")

# --- WP: survival / absorption vs b ---
md(r"""## 2. Wave-packet projectile — survival $\varepsilon$ vs $b$

The WP grazes the flake; the inner-region survival $\varepsilon(b)$ and absorbed
fraction quantify how strongly each impact parameter scatters/absorbs the packet.
(The WP disperses over the long parallel path — read qualitatively.)""")

code(r"""wps = sorted(glob.glob(os.path.join(SYS, "grazing", "run_wp_b*")))
rows = []
for d in wps:
    s = read_summary(os.path.join(d, "results/run_summary.txt"))
    if s.get("run_completed") != "true": continue
    b = b_of(os.path.basename(d))
    eps = s.get("epsilon") or s.get("wp_norm_tau") or s.get("N0")
    rows.append((b, s))
    print(f"b={b:g}: " + " ".join(f"{k}={s[k]}" for k in ("N0","wp_norm_tau","epsilon","absorbed_fraction") if k in s))
if not rows:
    print("no completed WP grazing runs yet — re-run after the b-scan finishes")""")

# --- field diagnostics per b (planar dn + LEED), partial-tolerant ---
md(r"""## 3. Field diagnostics — planar $\Delta n(z,t)$ and LEED per $b$

Reuses the validated `inqview.analysis` kernels. Planar $\Delta n(z,t)$ shows the
density disturbance the grazing projectile drags along the sheet; the LEED screens
show the scattered/diffracted flux. Rendered for whichever runs exist.""")

code(r"""from inqview.analysis.planar_density import planar_delta_map
from inqview.analysis.diffraction import diffraction_pattern
from inqview.io.leed import load_leed_pattern

def _load_vti(path):
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy
    r = vtk.vtkXMLImageDataReader(); r.SetFileName(str(path)); r.Update()
    img = r.GetOutput(); nx, ny, nz = img.GetDimensions()
    flat = vtk_to_numpy(img.GetPointData().GetArray(0)).astype(np.float64)
    return flat.reshape((nz, ny, nx)).transpose(2, 1, 0), img.GetOrigin(), img.GetSpacing(), (nx, ny, nz)

def planar_dn(run_dir, cat="density_rt_system"):
    cdir = os.path.join(run_dir, "results/raw/vti", cat)
    fs = sorted(glob.glob(cdir + "/*.vti"),
                key=lambda f: int(os.path.basename(f).split("_t")[1].split(".")[0]))
    if not fs: return None
    cubes, times, meta = [], [], None
    for f in fs:
        c, org, sp, shp = _load_vti(f); cubes.append(c); meta = (org, sp, shp)
        times.append(int(os.path.basename(f).split("_t")[1].split(".")[0]) * 0.02)
    org, sp, shp = meta; z = org[2] + np.arange(shp[2]) * sp[2]
    return planar_delta_map(cubes, np.array(times), z, axis=2, cell_area=sp[0] * sp[1])

runs = sorted(glob.glob(os.path.join(SYS, "grazing", "run_cl_b*"))) or \
       sorted(glob.glob(os.path.join(SYS, "grazing", "run_wp_b*")))
shown = 0
for d in runs:
    dm = planar_dn(d)
    if dm is None: continue
    fig, ax = plt.subplots(figsize=(6.8, 4.0))
    vmax = np.percentile(np.abs(dm.dn), 99.5) or 1e-12
    im = ax.imshow(dm.dn, origin="lower", aspect="auto", cmap="RdBu_r",
                   vmin=-vmax, vmax=vmax, extent=dm.extent)
    ax.axhline(0, color="k", lw=0.6, ls=":")
    for zc in (20, -20): ax.axhline(zc, color="0.4", lw=0.6, ls="--")
    ax.set(xlabel="time (a.u.)", ylabel="z (bohr)",
           title=f"planar Δn(z,t) — {os.path.basename(d)}")
    fig.colorbar(im, ax=ax, label="Δn (e/bohr)")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, f"fig_dn_{os.path.basename(d)}.png"), dpi=120)
    plt.show(); shown += 1
    if shown >= 3: break
if shown == 0:
    print("no VTI frames yet for planar Δn")""")

md(r"""## 4. Takeaway (PROVISIONAL)

- Finite-flake grazing geometry established: coronene C$_{24}$H$_{12}$ ∥ beam,
  impact parameter $b=x$-offset, classical (clean trajectory) + WP projectiles.
- Classical stopping vs $b$ (§1) is the headline observable; expect monotonic
  rise as $b\to0$. WP survival vs $b$ (§2) corroborates qualitatively.
- All CAP numbers **PROVISIONAL** until the inq-study engine regression (Task #7);
  feasibility replica, not converged paper values.""")

nb = new_notebook(cells=C, metadata={"kernelspec": {"name": "python3", "display_name": "Python 3"}})
ep = ExecutePreprocessor(timeout=2400, kernel_name="python3")
ep.preprocess(nb, {"metadata": {"path": HERE}})
out = os.path.join(HERE, "grazing_study.ipynb")
nbf.write(nb, out)
print("wrote + executed", out)
