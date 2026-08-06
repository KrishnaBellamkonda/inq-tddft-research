#!/usr/bin/env python3
"""Build (and execute) the PHASE-SPACE comparison notebook for one twin pair.

    venv/bin/python build_phase_notebook.py bulk_ks_stopping
    venv/bin/python build_phase_notebook.py --all

Plan: docs/plans/bulk-jellium-ks-stopping.md

WHAT THIS NOTEBOOK IS FOR. The per-pair `energy_component_comparison_*.ipynb`
compares the two halves in the ENERGY ledger. This one compares them in PHASE
SPACE: the (z, v) portrait of the projectile, how the two trajectories part
company, and how the local stopping S(z) = -dT/dz differs along the path. Same
pair, different question -- the energy notebook asks "where did the energy go",
this one asks "how did the projectile actually move".

WHY (z, v) IS THE RIGHT PLANE. Both halves are light free-Ehrenfest projectiles,
so they DECELERATE strongly rather than traversing at fixed velocity
(.claude/rules/light-projectile-stopping.md). A single S per run therefore hides
the physics: one run sweeps a whole velocity range. The phase portrait shows that
sweep directly -- each run traces a curve through (z, v), and the classical and
wavepacket curves peeling apart IS the quantum effect, plotted.

CONFIG PROVENANCE. The geometry/window constants are IMPORTED from each family's
own build_run_notebook.py rather than retyped here, so the two notebooks for a
family can never disagree about the box, the launch point, or the fit window.

Composes: ks_stopping.py (all arithmetic), inqview (density GIFs, per
.claude/rules/notebook-density-gif.md).
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import nbformat as nbf
from nbclient import NotebookClient

HERE = Path(__file__).resolve().parent
HYP = HERE.parent
REPO = HERE.parents[4]
SCRIPTS = REPO / "ResearchProject/systems/jellium/scripts"

FAMILIES = ["bulk_ks_stopping", "bulk_ks_stopping_sigma3",
            "bulk_ks_stopping_rs4", "bulk_ks_stopping_rs4_sigma3"]


def load_cfg(family: str) -> dict:
    """Import CFG from that family's build_run_notebook.py.

    Loaded under a unique module name so several families can be built in one
    process without the second import silently returning the first one's module.
    """
    path = HYP / family / "build_run_notebook.py"
    if not path.exists():
        raise FileNotFoundError(f"no build_run_notebook.py for {family} at {path}")
    spec = importlib.util.spec_from_file_location(f"_brn_{family}", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)           # module level is CFG + defs; main() is guarded
    return dict(mod.CFG)


def md(text: str) -> nbf.NotebookNode:
    return nbf.v4.new_markdown_cell(text)


def code(src: str) -> nbf.NotebookNode:
    return nbf.v4.new_code_cell(src)


def build(family: str) -> nbf.NotebookNode:
    cfg = load_cfg(family)
    out_dir = HYP / family
    cells: list[nbf.NotebookNode] = []

    cells.append(md(f"""# Phase-space comparison — **{family}**

**Classical projectile vs. wavepacket, in the $(z, v)$ plane.**

$\\sigma_{{\\rm WP}} = {cfg['SIGMA_WP']}$ Bohr &nbsp;·&nbsp;
$r_s = {cfg['R_S']}$ &nbsp;·&nbsp;
$E_{{\\rm kin}}(0) = {cfg['EKIN_EV']}$ eV &nbsp;·&nbsp;
box ${cfg['LX']:.0f}\\times{cfg['LY']:.0f}\\times{cfg['LZ']:.0f}$ Bohr &nbsp;·&nbsp;
$N_e = {cfg['N_E']}$

---

Both projectiles are **light** ($m = m_e$) and propagate under **free Ehrenfest
dynamics**, so neither travels at constant velocity — they decelerate strongly as
they deposit kinetic energy into the bath
(`.claude/rules/light-projectile-stopping.md`). A single stopping power per run
therefore *averages over a whole velocity range* and hides what actually happened.

The phase portrait puts that back: each run traces one curve through $(z, v)$, and
**the gap that opens between the two curves is the quantum effect, drawn directly.**

### What is comparable to what

The classical projectile has one kinetic energy, $T_{{\\rm cl}} = \\tfrac12 m v^2$.
The wavepacket has two, and only one of them is its counterpart:

| | | |
|---|---|---|
| $T_2 = \\langle p\\rangle^2/2m$ | **drift** KE | the classical analogue |
| $T_1 = \\langle p^2\\rangle/2m$ | total orbital KE | drift + internal |
| $T_1 - T_2$ | **internal** (momentum-width) energy | **no classical counterpart** |

Comparing $T_{{\\rm cl}}$ against $T_1$ would charge the wavepacket for its own
zero-point spread — about 2.6 eV at $\\sigma_{{\\rm WP}}=2$ Bohr *before it has moved
at all*. Everything below therefore compares $(z,v)$ and $(z,T_2)$, and reports
$T_1-T_2$ separately as the quantum-only channel.

Velocities need no conversion: the projectile is an electron, $m=1$ in atomic
units, so $v = \\langle p_z\\rangle$ exactly."""))

    cells.append(code(f"""import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import Image, display

HERE    = Path("{HERE}")
HYP     = Path("{HYP}")
SCRIPTS = Path("{SCRIPTS}")
FAMILY  = "{family}"
OUT     = Path("{out_dir}")
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(Path("{REPO}") / "inq-stack/python"))

import ks_stopping as K
from ks_stopping import HA_TO_EV

CFG = {cfg!r}

plt.rcParams.update({{"figure.dpi": 120, "figure.figsize": (9, 4),
                     "axes.grid": True, "grid.alpha": 0.25,
                     "axes.spines.top": False, "axes.spines.right": False}})

pair = K.load_pair(SCRIPTS, FAMILY, box_length_z=CFG["LZ"], z0=CFG["Z0"])
CL, WP = pair.cl, pair.wp
print(f"classical : {{len(CL.t):4d}} steps,  t = {{CL.t[0]:.2f}} .. {{CL.t[-1]:.2f}} a.u.")
print(f"wavepacket: {{len(WP.t):4d}} steps,  t = {{WP.t[0]:.2f}} .. {{WP.t[-1]:.2f}} a.u.")
print(f"v0        : {{pair.v0:.6f}} Bohr/a.u.   (both halves launched at the same k0)")
print(f"interactions: wp={{pair.ix_wp is not None}}  classical={{pair.ix_cl is not None}}")
if pair.ix_cl is not None and np.isfinite(pair.ix_cl.clip_time):
    print(f"clipping onset: t = {{pair.ix_cl.clip_time:.2f}} a.u. "
          f"(hard upper bound on any fit window)")"""))

    # ---------------------------------------------- density GIF (always-on rule)
    cells.append(md("""---
## 1. Visual intuition — the density matrix

`.claude/rules/notebook-density-gif.md`: the animated $n(x,z,t)$ on the mid-$y$
propagation slice, rows {classical, wavepacket, WP−classical} × columns
{density, induced, instantaneous}. This is the most direct picture of what the two
projectiles do differently; everything quantitative below is a projection of it.

*Bulk note:* there is no slab, so the dashed lines mark the **cell faces** — the
boundary the classical cloud eventually runs into."""))

    cells.append(code("""from inqview.visualisation import make_twin_density_matrix

gifs = make_twin_density_matrix(
    classical_dir=str(SCRIPTS / FAMILY / "classical" / "results"),
    wp_dir=str(SCRIPTS / FAMILY / "wp" / "results"),
    out_dir=str(OUT / f"phase_gifs_{FAMILY}"),
    dt=CFG["DT"], slab_face=CFG["LZ"] / 2.0,
    frames_max=30, fps=8,
    total_subpath="raw/vti/density_total",
)
print(f"{len(gifs)} density GIFs" if gifs else "no density frames found")"""))

    cells.append(code("""for row, col, path, ttl in gifs:
    print(f"--- {row} / {col} ---")
    display(Image(filename=path))"""))

    # ---------------------------------------------------------- phase portrait
    cells.append(md(r"""---
## 2. The phase portrait

The headline figure: both projectiles in the $(z, v_z)$ plane. Time runs left to
right (both launch at the same $z_0$ with the same $v_0$). Markers are placed at
equal time intervals, so **where the markers bunch up, the projectile is slow**.

The shaded band is the stopping fit window. The red dotted line, where present, is
where the classical projectile's Gaussian charge cloud starts leaving the grid —
a hard bound on how far right anything may be read."""))

    # NOTE: raw string. These cells contain LaTeX like \rangle and \rm; in a
    # non-raw builder string Python eats the \r as a carriage return and the
    # emitted cell is a syntax error before it is ever executed.
    cells.append(code(r"""fig, ax = plt.subplots(figsize=(9.5, 5.5))
ax.axvspan(np.interp(CFG["FIT_T0"], CL.t, CL.z),
           np.interp(CFG["FIT_T1"], CL.t, CL.z), color="0.88", zorder=0,
           label="fit window")

ax.plot(CL.z, CL.vz, lw=2.0, color="C1", label="classical  $v_z$", zorder=3)
ax.plot(WP.s3, WP.pz, lw=2.0, color="C0", label=r"wavepacket  $\langle p_z\rangle$",
        zorder=3)

# Equal-TIME markers: spacing encodes speed, which a plain line cannot show.
for t_, z_, v_, c in ((CL.t, CL.z, CL.vz, "C1"), (WP.t, WP.s3, WP.pz, "C0")):
    tm = np.arange(t_[0], t_[-1], 2.0)
    ax.plot(np.interp(tm, t_, z_), np.interp(tm, t_, v_), "o", ms=3.2,
            color=c, zorder=4)

d = pair.divergence(frac=0.05)
if np.isfinite(d["t"]):
    ax.plot([d["z_cl"], d["z_wp"]], [d["v_cl"], d["v_wp"]], "-", color="0.35",
            lw=1.0, zorder=5)
    ax.annotate(f"histories part by 5% of $v_0$\nt = {d['t']:.1f} a.u.",
                xy=(d["z_cl"], d["v_cl"]), xytext=(0.05, 0.18),
                textcoords="axes fraction", fontsize=8.5, color="0.25",
                arrowprops=dict(arrowstyle="->", color="0.35", lw=1.0))

if pair.ix_cl is not None and np.isfinite(pair.ix_cl.clip_time):
    ax.axvline(np.interp(pair.ix_cl.clip_time, CL.t, CL.z), color="C3", ls=":",
               lw=1.4, label="classical cloud clips the face")

ax.axhline(pair.v0, color="0.6", lw=0.8, ls="--")
ax.text(CL.z[0], pair.v0, "  $v_0$", va="bottom", ha="left", fontsize=8, color="0.4")
ax.set_xlabel("projectile position  $z$  [Bohr]")
ax.set_ylabel(r"velocity  $v_z$  [Bohr / a.u.]")
ax.set_title(f"Phase portrait — {FAMILY}   (markers every 2 a.u. of time)")
ax.legend(fontsize=8.5, frameon=False, loc="upper right")
ax.margins(x=0.02)
fig.tight_layout()
fig.savefig(OUT / f"phase_portrait_{FAMILY}.png", dpi=150, bbox_inches="tight")
plt.show()"""))

    cells.append(md(r"""**How to read it.** A projectile losing energy moves
*down* as it moves *right*. The steeper the descent, the harder it is being
stopped. Two features carry the physics:

1. **Which curve sits higher** at a given $z$ — that projectile has been stopped
   less by the time it got there.
2. **Where the curves separate** — before that point the two representations are
   dynamically equivalent, and the wavepacket is behaving like a classical
   particle. After it, they are not."""))

    # --------------------------------------------------------- velocity vs time
    cells.append(md(r"""---
## 3. Velocity histories, and the gap between them

The same data against time rather than position, plus the difference
$\Delta v = v_{\rm WP} - v_{\rm cl}$ that the phase portrait shows only implicitly."""))

    cells.append(code(r"""fig, ax = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
for a in ax:
    a.axvspan(CFG["FIT_T0"], CFG["FIT_T1"], color="0.88", zorder=0)
    if pair.ix_cl is not None and np.isfinite(pair.ix_cl.clip_time):
        a.axvline(pair.ix_cl.clip_time, color="C3", ls=":", lw=1.3)

ax[0].plot(CL.t, CL.vz, lw=1.8, color="C1", label="classical")
ax[0].plot(WP.t, WP.pz, lw=1.8, color="C0", label="wavepacket")
ax[0].axhline(pair.v0, color="0.6", lw=0.8, ls="--")
ax[0].set_ylabel(r"$v_z$  [Bohr / a.u.]")
ax[0].set_title("Velocity history")
ax[0].legend(fontsize=8.5, frameon=False)

v_wp_on_cl = np.interp(CL.t, WP.t, WP.pz)
dv = v_wp_on_cl - CL.vz
ax[1].plot(CL.t, dv, lw=1.8, color="C4")
ax[1].axhline(0.0, color="0.6", lw=0.8, ls="--")
ax[1].fill_between(CL.t, 0.0, dv, color="C4", alpha=0.15)
ax[1].set_xlabel("time  [a.u.]")
ax[1].set_ylabel(r"$v_{\rm WP} - v_{\rm cl}$  [Bohr / a.u.]")
ax[1].set_title("Velocity gap — positive means the wavepacket is being stopped LESS")
for a in ax: a.margins(x=0)
fig.tight_layout()
fig.savefig(OUT / f"phase_velocity_{FAMILY}.png", dpi=150, bbox_inches="tight")
plt.show()

i1 = np.searchsorted(CL.t, CFG["FIT_T1"])
print(f"  v/v0 at end of fit window : classical {CL.vz[i1]/pair.v0:.4f}   "
      f"wavepacket {v_wp_on_cl[i1]/pair.v0:.4f}")
print(f"  velocity gap there        : {dv[i1]:+.5f} Bohr/a.u. "
      f"({100*dv[i1]/pair.v0:+.2f}% of v0)")"""))

    # ------------------------------------------------------------- KE and S(z)
    cells.append(md(r"""---
## 4. Kinetic energy along the path, and the internal channel

Left: the comparable pair — classical $T_{\rm cl}=\tfrac12mv^2$ against the
wavepacket **drift** energy $T_2=\langle p\rangle^2/2m$.

Right: the wavepacket's **internal** energy $T_1-T_2$, which has no classical
counterpart at all. For a free Gaussian it starts at $3/(8\sigma^2)$ and grows as
the packet is deformed by the bath. Any apparent "stopping" that actually lives in
this channel is momentum-space broadening, not drift-momentum loss."""))

    cells.append(code(r"""fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.6))

ax[0].axvspan(np.interp(CFG["FIT_T0"], CL.t, CL.z),
              np.interp(CFG["FIT_T1"], CL.t, CL.z), color="0.88", zorder=0)
ax[0].plot(CL.z, CL.T * HA_TO_EV, lw=1.8, color="C1",
           label=r"classical $\frac{1}{2}mv^2$")
ax[0].plot(WP.s3, WP.T2 * HA_TO_EV, lw=1.8, color="C0",
           label=r"wavepacket drift $\langle p\rangle^2/2m$")
ax[0].set_xlabel("$z$  [Bohr]"); ax[0].set_ylabel("kinetic energy  [eV]")
ax[0].set_title("Drift kinetic energy — the comparable pair")
ax[0].legend(fontsize=8.5, frameon=False)

ax[1].axvspan(CFG["FIT_T0"], CFG["FIT_T1"], color="0.88", zorder=0)
ax[1].plot(WP.t, WP.localisation_energy, lw=1.8, color="C2")
sigma0 = CFG["SIGMA_WP"]
ax[1].axhline(3.0 / (8.0 * sigma0**2) * HA_TO_EV, color="0.5", lw=1.0, ls="--")
ax[1].text(WP.t[-1], 3.0 / (8.0 * sigma0**2) * HA_TO_EV,
           r"  $3/(8\sigma^2)$", va="center", ha="right", fontsize=8, color="0.4")
ax[1].set_xlabel("time  [a.u.]")
ax[1].set_ylabel(r"$T_1 - T_2$  [eV]")
ax[1].set_title("Internal (momentum-width) energy — quantum only")
for a in ax: a.margins(x=0)
fig.tight_layout()
fig.savefig(OUT / f"phase_kinetic_{FAMILY}.png", dpi=150, bbox_inches="tight")
plt.show()"""))

    cells.append(md(r"""---
## 5. Local stopping power $S(z) = -\mathrm{d}T/\mathrm{d}z$

A rolling least-squares slope of $T(z)$ — the same estimator the windowed fit
uses, evaluated locally instead of over one wide window, so the local curve and
the reported single number cannot disagree.

**Read the interior only.** The first and last few points repeat the nearest
interior value (the estimator runs out of samples), so the flat shoulders at both
ends are an artefact of the method, not impact physics."""))

    cells.append(code(r"""HW = 12
S_cl = K.local_stopping(CL.z,  CL.T,  half_width=HW)
S_wp = K.local_stopping(WP.s3, WP.T2, half_width=HW)

fig, ax = plt.subplots(figsize=(9.5, 5.0))
ax.axvspan(np.interp(CFG["FIT_T0"], CL.t, CL.z),
           np.interp(CFG["FIT_T1"], CL.t, CL.z), color="0.88", zorder=0,
           label="fit window")
ax.plot(CL.z[HW:-HW],  S_cl[HW:-HW], lw=1.8, color="C1", label="classical")
ax.plot(WP.s3[HW:-HW], S_wp[HW:-HW], lw=1.8, color="C0",
        label="wavepacket (drift KE)")
if pair.ix_cl is not None and np.isfinite(pair.ix_cl.clip_time):
    ax.axvline(np.interp(pair.ix_cl.clip_time, CL.t, CL.z), color="C3", ls=":",
               lw=1.4, label="cloud clips the face")
ax.axhline(0.0, color="0.6", lw=0.8, ls="--")
ax.set_xlabel("$z$  [Bohr]")
ax.set_ylabel(r"$S(z) = -\mathrm{d}T/\mathrm{d}z$  [eV / Bohr]")
ax.set_title(f"Local stopping power along the path — {FAMILY}  "
             f"(rolling OLS, {2*HW+1} samples)")
ax.legend(fontsize=8.5, frameon=False)
ax.margins(x=0.02)
fig.tight_layout()
fig.savefig(OUT / f"phase_local_S_{FAMILY}.png", dpi=150, bbox_inches="tight")
plt.show()

w_cl = (CL.t >= CFG["FIT_T0"]) & (CL.t <= CFG["FIT_T1"])
w_wp = (WP.t >= CFG["FIT_T0"]) & (WP.t <= CFG["FIT_T1"])
print(f"  mean S in fit window : classical {np.nanmean(S_cl[w_cl]):.4f}   "
      f"wavepacket {np.nanmean(S_wp[w_wp]):.4f}  eV/Bohr")
r = np.nanmean(S_cl[w_cl]) / np.nanmean(S_wp[w_wp]) if np.nanmean(S_wp[w_wp]) else np.nan
print(f"  classical / wavepacket ratio : {r:.3f}")

# CROSS-CHECK. The mean of local slopes and a single global OLS slope are
# DIFFERENT estimators and need not agree when T(s) is curved. Comparing them is
# therefore a real test: if they agree, the local curve can be trusted to
# decompose the same number the windowed fit reports. S_24 is the WP fit built
# on the DRIFT kinetic energy T2, so it is the like-for-like partner of S_wp.
gl_cl = K.fit_classical(CL, CFG["FIT_T0"], CFG["FIT_T1"])
gl_wp = K.fit_all_wp(WP, CFG["FIT_T0"], CFG["FIT_T1"])["S_24"]
print()
print("  cross-check, local mean vs global OLS over the same window:")
for nm, loc, glo in (("classical", np.nanmean(S_cl[w_cl]), gl_cl.S_ev_per_bohr),
                     ("wavepacket", np.nanmean(S_wp[w_wp]), gl_wp.S_ev_per_bohr)):
    dev = 100.0 * abs(loc - glo) / abs(glo) if glo else float("nan")
    print(f"    {nm:<11}: local {loc:.4f}   global {glo:.4f}   "
          f"differ by {dev:.1f}%  [{'OK' if dev < 15 else 'CHECK CURVATURE'}]")
print(f"    global ratio: {gl_cl.S_ev_per_bohr / gl_wp.S_ev_per_bohr:.3f}")"""))

    # ----------------------------------------------------------- E_PP on z axis
    cells.append(md(r"""---
## 6. Self-Hartree $E_{PP}$ on the same position axis

$E_{PP}$ is the projectile self-Hartree — the uncancelled LDA self-interaction a
wavepacket has *because it is an occupied KS orbital*, and a classical external
potential simply does not have. Putting it on the **same $z$ axis** as the phase
portrait answers the question the energy notebook cannot: *where along the path*
does the packet shed it, and does that coincide with where the trajectories part?

For the classical half this must be a **horizontal line** until the cloud clips
the face — a rigid translating charge cannot change its own self-energy. That
flatness is a free end-to-end validation of the whole measurement.

**Gauge:** absolute $E_{PP}$ carries the charged-cell $G{=}0$ gauge. Only the
*difference* between the halves, in this same cell, is gauge-clean."""))

    cells.append(code(r"""if pair.ix_wp is None or pair.ix_cl is None:
    print("interactions.csv missing on one half — section skipped")
else:
    z_wp, epp_wp = pair.epp_on_z("wp")
    z_cl, epp_cl = pair.epp_on_z("classical")

    fig, ax = plt.subplots(2, 1, figsize=(9.5, 7.5), sharex=True)
    for a in ax:
        a.axvspan(np.interp(CFG["FIT_T0"], CL.t, CL.z),
                  np.interp(CFG["FIT_T1"], CL.t, CL.z), color="0.88", zorder=0)
        if np.isfinite(pair.ix_cl.clip_time):
            a.axvline(np.interp(pair.ix_cl.clip_time, CL.t, CL.z), color="C3",
                      ls=":", lw=1.4)

    ax[0].plot(z_cl, epp_cl, lw=1.8, color="C1", label="classical (rigid cloud)")
    ax[0].plot(z_wp, epp_wp, lw=1.8, color="C0", label="wavepacket (disperses)")
    ax[0].set_ylabel(r"$E_{PP}$  [eV]")
    ax[0].set_title(r"Projectile self-Hartree along the path")
    ax[0].legend(fontsize=8.5, frameon=False)

    # Gauge-clean difference, on the classical z axis (both halves, same cell).
    d_epp = np.interp(CL.z, z_wp, epp_wp) - np.interp(CL.z, z_cl, epp_cl)
    ax[1].plot(CL.z, d_epp, lw=1.8, color="C4")
    ax[1].axhline(0.0, color="0.6", lw=0.8, ls="--")
    ax[1].fill_between(CL.z, 0.0, d_epp, color="C4", alpha=0.15)
    if np.isfinite(d["t"]):
        ax[1].axvline(d["z_cl"], color="0.35", lw=1.0, ls="-.")
        ax[1].text(d["z_cl"], ax[1].get_ylim()[0], " trajectories part here",
                   fontsize=8, color="0.25", va="bottom", rotation=90)
    ax[1].set_xlabel("$z$  [Bohr]")
    ax[1].set_ylabel(r"$\Delta E_{PP}$  [eV]")
    ax[1].set_title(r"Gauge-clean $\Delta E_{PP} = E_{PP}^{\rm WP} - E_{PP}^{\rm cl}$")
    for a in ax: a.margins(x=0)
    fig.tight_layout()
    fig.savefig(OUT / f"phase_epp_{FAMILY}.png", dpi=150, bbox_inches="tight")
    plt.show()

    # Windowed to [FIT_T0, min(FIT_T1, clipping onset)] -- the honest bound.
    m = pair.ix_cl.in_window(CFG["FIT_T0"], CFG["FIT_T1"])
    t_end = pair.ix_cl.t[m][-1]
    e_w = np.interp(t_end, pair.ix_wp.t, pair.ix_wp.e_pp) * HA_TO_EV
    e_c = np.interp(t_end, pair.ix_cl.t, pair.ix_cl.e_pp) * HA_TO_EV
    e_w0 = pair.ix_wp.e_pp[0] * HA_TO_EV
    e_c0 = pair.ix_cl.e_pp[0] * HA_TO_EV
    print(f"  window            : t = {CFG['FIT_T0']:.2f} .. {t_end:.2f} a.u.")
    print(f"  E_PP(0)           : wp {e_w0:.4f}   classical {e_c0:.4f}  eV")
    print(f"  dE_PP(0)          : {e_w0-e_c0:+.2e} eV   "
          f"[{'PASS' if abs(e_w0-e_c0) < 1e-3 else 'FAIL'} sigma-matching gate]")
    print(f"  dE_PP(end of win) : {e_w-e_c:+.4f} eV")
    print(f"  classical drift   : {e_c-e_c0:+.2e} eV   "
          f"[{'PASS' if abs(e_c-e_c0) < 1e-4 else 'FAIL'} rigid-cloud gate]")"""))

    # ------------------------------------------------------------------ summary
    cells.append(md("""---
## 7. Takeaway

The numbers this notebook establishes, and what they do and do not mean."""))

    cells.append(code("""print(f"=== {FAMILY} ===")
print(f"  sigma_WP = {CFG['SIGMA_WP']} Bohr,  r_s = {CFG['R_S']},  "
      f"E_kin(0) = {CFG['EKIN_EV']} eV")
print()
print(f"  v0                          : {pair.v0:.5f} Bohr/a.u.")
print(f"  v/v0 at end of fit window   : classical {CL.vz[i1]/pair.v0:.4f}, "
      f"wavepacket {v_wp_on_cl[i1]/pair.v0:.4f}")
if np.isfinite(d["t"]):
    print(f"  trajectories part (5% v0)   : t = {d['t']:.2f} a.u., "
          f"z = {d['z_cl']:+.2f} Bohr")
else:
    print("  trajectories part (5% v0)   : never within this run")
print(f"  mean S in window            : classical {np.nanmean(S_cl[w_cl]):.4f}, "
      f"wavepacket {np.nanmean(S_wp[w_wp]):.4f} eV/Bohr  (ratio {r:.3f})")
print(f"  internal T1-T2              : {WP.localisation_energy[0]:.3f} -> "
      f"{WP.localisation_energy[-1]:.3f} eV")
if pair.ix_cl is not None and np.isfinite(pair.ix_cl.clip_time):
    print(f"  cloud clipping onset        : t = {pair.ix_cl.clip_time:.2f} a.u. "
          f"(fit ends {CFG['FIT_T1']:.2f} -> "
          f"{'CLEAR' if CFG['FIT_T1'] < pair.ix_cl.clip_time else 'CONTAMINATED'})")"""))

    cells.append(md(f"""**Limitations, in order of importance.**

1. **Not a steady-state stopping power.** The plasma period
   ({cfg['PLASMA_PERIOD']} a.u.) exceeds the run
   ({cfg['DT']*cfg['N_STEPS']:.1f} a.u.), so no steady wake forms. Every $S$ here
   is an *initial-drag* number, and $S(z)$ is a local slope within that regime.
2. **$S(z)$ is a rolling estimator**, not a measurement at a point. Its
   resolution is the window width; structure narrower than that is smoothed away,
   and the outer $\\pm$12 samples are filled, not fitted.
3. **The divergence marker is a reading aid.** "5% of $v_0$" is a threshold chosen
   for legibility, not a physically derived scale. Its *value* is arbitrary; only
   the *ordering* it reveals is meaningful.
4. **Classical cloud clipping** bounds the usable range on the right — it is a
   property of the trajectory meeting the box, not of the physics, and it is not
   recorded in any config header.
5. **One velocity, one width, one density per notebook.** The width-vs-density
   contrast needs the 2x2 across all four pairs
   (`summarise_epp_across_pairs.py`), not this notebook.

**Provenance.** All arithmetic by
`hypotheses/bulk_ks_stopping/ks_stopping.py` (23 known-case tests in
`tests/test_ks_stopping.py`). Geometry and window constants imported directly from
`hypotheses/{family}/build_run_notebook.py`, so this notebook and the run notebooks
cannot disagree. Density GIFs by `inqview.visualisation.make_twin_density_matrix`
(physical-order VTIs, no fftshift — `.claude/rules/vti-coordinate-mapping.md`).
"""))

    nb = nbf.v4.new_notebook(cells=cells)
    nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python",
                                 "name": "python3"}
    return nb


def check_cells_compile(nb: nbf.NotebookNode) -> list[str]:
    """Compile every code cell before the notebook is written.

    Guards a bug class that cost a full build cycle on 2026-08-01: the cell
    sources are emitted from triple-quoted strings in THIS file, so a non-raw
    string containing LaTeX like ``\\rangle`` or ``\\rm`` has its ``\\r`` eaten as a
    carriage return and the emitted cell is a SyntaxError -- which only surfaces
    after the (slow) GIF stage has already run. Compiling here fails in
    milliseconds instead.
    """
    bad = []
    for i, c in enumerate(nb.cells):
        if c.cell_type != "code":
            continue
        try:
            compile(c.source, f"<cell {i}>", "exec")
        except SyntaxError as e:
            bad.append(f"cell {i}: {e.msg} (line {e.lineno}): {(e.text or '').strip()[:70]}")
    return bad


def build_one(family: str, execute: bool = True) -> int:
    out = HYP / family / f"run_pair_phase_{family}.ipynb"
    nb = build(family)
    bad = check_cells_compile(nb)
    if bad:
        print(f"FATAL {family}: {len(bad)} cell(s) do not compile — not written:")
        for b in bad:
            print("   ", b)
        return 1
    out.write_text(nbf.writes(nb))
    print(f"wrote {out} ({len(nb.cells)} cells)" + ("; executing..." if execute else ""))
    if not execute:
        return 0
    client = NotebookClient(nb, timeout=3600, kernel_name="python3",
                            resources={"metadata": {"path": str(HYP / family)}},
                            allow_errors=True)
    client.execute()
    errs = sum(1 for c in nb.cells
               if any(o.get("output_type") == "error" for o in c.get("outputs", [])))
    out.write_text(nbf.writes(nb))
    print(f"  executed with {errs} error(s) -> {out}")
    return 1 if errs else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("family", nargs="?", help="one of: " + ", ".join(FAMILIES))
    ap.add_argument("--all", action="store_true", help="build every family")
    ap.add_argument("--no-execute", action="store_true", help="write without running")
    a = ap.parse_args()

    if a.all:
        targets = FAMILIES
    elif a.family in FAMILIES:
        targets = [a.family]
    else:
        ap.error(f"family must be one of {FAMILIES} (got {a.family!r}), or use --all")

    rc = 0
    for f in targets:
        rc |= build_one(f, execute=not a.no_execute)
    return rc


if __name__ == "__main__":
    sys.exit(main())
