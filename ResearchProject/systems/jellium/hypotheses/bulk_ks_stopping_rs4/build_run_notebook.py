#!/usr/bin/env python3
"""Build (and execute) the run notebook for one half of the DENSITY-REPLICA
(r_s = 3.99) bulk-jellium KS-stopping twin pair.

    venv/bin/python build_run_notebook.py wp
    venv/bin/python build_run_notebook.py classical

Plan: docs/plans/bulk-jellium-ks-stopping.md

House narrative (notebook-making skill): context -> formulas (every term defined)
-> full reconstructable setup -> results -> takeaway. Plus the always-on rule
.claude/rules/notebook-density-gif.md: the density-matrix GIFs go at the TOP and
must be DISPLAYED inline, not merely written to disk.

The notebook narrates; ks_stopping.py does the arithmetic. Nothing is recomputed
here by hand, so the notebook and any later cross-run comparison cannot disagree.
"""
from __future__ import annotations

import sys
from pathlib import Path

import nbformat as nbf
from nbclient import NotebookClient

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]
SYSTEM = REPO / "ResearchProject/systems/jellium"

# Geometry / physics constants, mirrored from
# shared/configs/bulk_ks_stopping_L40x40x80_rs4.hpp. Kept as literals so the
# notebook is readable standalone; the run_summary.txt printed in the setup cell
# is the authority if they ever disagree.
CFG = dict(
    LX=40.0, LY=40.0, LZ=80.0, N_E=482, DX=0.50,
    SIGMA_WP=2.0, EKIN_EV=100.0, K0=2.7111, Z0=-32.0,
    DT=0.04, N_STEPS=646,
    # NOTE the INVERSION vs the r_s=5.702 pair: at L_xy=40 the TRANSVERSE limit
    # (18.43) binds before the longitudinal one (18.97), so FIT_T1 comes from
    # the periodic images in x/y, not from the +z face.
    FIT_T0=4.0, FIT_T1=18.4270,
    T_IFW=18.9719, T_TRANSVERSE=18.4270,
    R_S=3.9874, HW_P_EV=5.9194, V_F=0.4813, PLASMA_PERIOD=28.88,
)


def md(text: str) -> nbf.NotebookNode:
    return nbf.v4.new_markdown_cell(text)


def code(src: str) -> nbf.NotebookNode:
    return nbf.v4.new_code_cell(src)


def build(half: str) -> nbf.NotebookNode:
    is_wp = half == "wp"
    run_dir = SYSTEM / "scripts/bulk_ks_stopping_rs4" / half
    title = "Wavepacket" if is_wp else "Classical"
    cells: list[nbf.NotebookNode] = []

    # ---------------------------------------------------------------- header
    cells.append(md(f"""# Bulk-jellium KS stopping power — **{title} projectile**

A {CFG['EKIN_EV']:.0f} eV electron crossing a fully periodic bulk jellium bath,
represented { "as a Gaussian wavepacket occupying its own Kohn-Sham orbital"
              if is_wp else
              "as a classical Gaussian charge moving under free Ehrenfest dynamics" }.

This is one half of a **twin pair**: both runs share the same ground state, the
same box, the same projectile energy and the same charge distribution, and differ
*only* in how the projectile is represented. Everything below is therefore
directly comparable with the other notebook.

| | |
|---|---|
| Cell | {CFG['LX']:.0f} x {CFG['LY']:.0f} x {CFG['LZ']:.0f} Bohr, periodic in **x, y and z** |
| Bath | N = {CFG['N_E']} electrons, r_s = {CFG['R_S']}, hw_p = {CFG['HW_P_EV']} eV, v_F = {CFG['V_F']} |
| Projectile | {CFG['EKIN_EV']:.0f} eV, k0 = v = {CFG['K0']} a.u., sigma = {CFG['SIGMA_WP']} Bohr, launched at z = {CFG['Z0']} |
| Grid / time | dx = {CFG['DX']} Bohr, dt = {CFG['DT']} a.u., {CFG['N_STEPS']} steps (t = {CFG['DT']*CFG['N_STEPS']:.2f} a.u.) |
| Fit window | t in [{CFG['FIT_T0']}, {CFG['FIT_T1']:.2f}] a.u. |

**Plan:** `docs/plans/bulk-jellium-ks-stopping.md`

---

### Read this before quoting any number

**The bath cannot complete one plasma oscillation during this run.**
2*pi/w_p = {CFG['PLASMA_PERIOD']} a.u., against a {CFG['DT']*CFG['N_STEPS']:.1f} a.u. run and a
{CFG['FIT_T1']-CFG['FIT_T0']:.1f} a.u. fit window. This is forced, not an oversight: a light
{CFG['EKIN_EV']:.0f} eV electron crosses this box in ~26 a.u., so the steady-state wake
criterion (5 plasma periods ~ 247 a.u.) is geometrically unreachable. Per
`.claude/rules/light-projectile-stopping.md` we therefore extract S as an
**initial-drag** stopping power. It is *not* a converged steady-state S(v), and
should not be compared with Lindhard or Bethe values as though it were.
"""))

    # ------------------------------------------------------- setup / imports
    cells.append(code(f"""import sys, json, math
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from IPython.display import Image, display

HERE = Path("{HERE}")
RUN  = Path("{run_dir}")
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(Path("{REPO}") / "inq-stack/python"))
# The stopping-power-extraction skill ships the reference kernels (Correa 2018
# Eq.10 free-intercept fit, the fixed-20%-time transient cut, the guards). We
# import them rather than re-implement so the cross-check uses the same code
# every other run in this project is measured with.
sys.path.insert(0, str(Path("{REPO}") / ".claude/skills/stopping-power-extraction"))

import ks_stopping as K
from ks_stopping import HA_TO_EV

CFG = {CFG!r}
HALF = "{half}"

plt.rcParams.update({{"figure.dpi": 120, "figure.figsize": (9, 4),
                     "axes.grid": True, "grid.alpha": 0.25,
                     "axes.spines.top": False, "axes.spines.right": False}})

print("run directory:", RUN)
summary = (RUN / "results/run_summary.txt")
print(summary.read_text() if summary.exists() else "run_summary.txt not found")"""))

    # ------------------------------------------------- density GIFs (RULE)
    cells.append(md(f"""---
## 1. Visual intuition — the density evolution

Per `.claude/rules/notebook-density-gif.md` this comes **first**: the animated
n(r,t) in the x-z plane (mid-y slice) is the most direct picture of what the
projectile does, and every quantitative result below is an attempt to compress
these animations into one number.

Three *kinds* are shown for each *category*:

- **density** — n(x,z,t) itself
- **delta0** — Dn = n(t) - n(0), the **induced** density: the wake the projectile
  leaves behind, with the uniform bath subtracted off
- **dstep** — Dn = n(t+dt) - n(t), the instantaneous change, which makes the
  propagating disturbance visible even where the accumulated wake is flat

{"and three categories: the **total** electron system, the **wavepacket** orbital "
 "|psi_wp|^2 alone, and the **bath** (total minus wavepacket)."
 if is_wp else
 "The classical run has no projectile orbital, so only the **total** system is "
 "available — the projectile itself is a pseudo-ion, visible only through the "
 "hole it carves in the electron density."}

VTIs are written in **physical order** by inqkit (`.claude/rules/vti-coordinate-mapping.md`),
so they are loaded through `inqview.load_vti` and never `fftshift`ed."""))

    cells.append(code(f"""from inqview.visualisation import make_density_gif_battery

OUT = RUN / "results" / "report"
OUT.mkdir(parents=True, exist_ok=True)

# Bulk jellium: no slab faces and no absorbing potential, so no guide lines.
gifs, vmax = make_density_gif_battery(
    str(RUN / "results"), str(OUT),
    run_label=HALF, dt=CFG["DT"] * CFG.get("WRITE_EVERY", 2),
    slab_face=0.0, cap_inner=0.0, cap_lines=(),
    frames_max=40, fps=10,
    run_title="{title} projectile — bulk jellium",
)
print(f"built {{len(gifs)}} gifs, shared density scale vmax = {{vmax:.3e}}")
for cat, kind, path, ttl in gifs:
    print(f"  {{cat:6s}} {{kind:8s}} {{Path(path).name}}")"""))

    cells.append(md("""### The animations

Displayed inline (base64-embedded into the stored outputs), so they animate when
this notebook is reopened without needing the sidecar `.gif` files."""))

    cells.append(code("""for cat, kind, path, ttl in gifs:
    print(f"=== {ttl} ===")
    display(Image(filename=path))"""))

    # ----------------------------------------------------------- formulas
    if is_wp:
        cells.append(md(r"""---
## 2. What we are computing, and why there is more than one answer

The stopping power is the energy the projectile loses per unit distance
travelled:

$$ S \;=\; -\frac{dT}{ds} $$

For a *classical* projectile both $T$ and $s$ are unambiguous. For a **quantum**
projectile neither is: the packet has a spread in momentum and a spread in
position, so "the kinetic energy" and "the position" each admit more than one
reasonable definition. This run measures all of them.

### Two kinetic energies

$$ T_1 \;=\; \frac{\langle \hat p^2\rangle}{2m}
   \;=\; \tfrac12\!\int\! k^2\,|\tilde\psi_{\rm wp}(\mathbf k,t)|^2\,d^3k $$

$$ T_2 \;=\; \frac{\langle \hat{\mathbf p}\rangle^2}{2m}
   \;=\; \tfrac12\Big(\langle p_x\rangle^2+\langle p_y\rangle^2+\langle p_z\rangle^2\Big) $$

$T_1$ is the full kinetic energy of the orbital — the quantity INQ itself tracks.
$T_2$ keeps only the **drift** momentum and throws away everything in the spread.
Their difference is exactly the momentum-width energy:

$$ T_1 - T_2 \;=\; \tfrac12\sum_d \sigma_{p_d}^2
   \;\xrightarrow[\ t=0\ ]{}\; \frac{3}{4\sigma_{\rm WP}^2} \;=\; 5.102\ {\rm eV}
   \quad (\sigma_{\rm WP}=2\ {\rm Bohr}) $$

where $\sigma_{p_d}^2 = \langle p_d^2\rangle - \langle p_d\rangle^2$. At $t=0$
this is pure **localisation** (zero-point) energy, fixed by the packet width. As
the run proceeds it also collects **elastic/angular scattering**: a collision
that deflects momentum without removing it increases the spread while reducing
$\langle p_z\rangle$. So:

> $S_1$ counts scattered momentum as *still kinetic*; $S_2$ counts it as *lost*.
> **Their difference measures how much of the apparent stopping is angular
> scattering rather than genuine drift-momentum loss.** That is the central
> question this twin pair was built to answer.

### Two positions

$$ s_3 \;=\; \langle z\rangle_{\rm circ}
   \;=\; \frac{L_z}{2\pi}\arg\big\langle \psi\big|e^{\,i 2\pi z/L_z}\big|\psi\big\rangle
   \qquad\text{(density centroid)} $$

$$ s_4 \;=\; z_0 + \int_0^t \frac{\langle p_z\rangle(t')}{m}\,dt'
   \qquad\text{(integrated drift momentum)} $$

$s_3$ is the **circular** centroid, not the naive $\int z|\psi|^2 d^3r$. In a
periodic cell the naive integral is discontinuous across a face and — worse —
slides *smoothly to a wrong value* while the packet straddles it, so the failure
is invisible to the eye. The phase estimator above is exact for any cell. Both
are plotted below so you can see where they part company.

### The identity you should check first

This run has **no ions and no absorbing potential**, so the Kohn-Sham
Hamiltonian is purely local (kinetic + Hartree + ALDA). Ehrenfest's theorem then
gives, exactly,

$$ \frac{d\langle z\rangle}{dt} \;=\; \frac{\langle p_z\rangle}{m}
   \qquad\Longrightarrow\qquad s_3 \equiv s_4 . $$

**So $s_3$ and $s_4$ are not two independent measurements — their agreement is a
validation of the propagation.** Any drift between them localises a problem:
either the packet has wrapped and the naive centroid crept in, or the wavepacket
orbital is leaking norm into the bath. (Contrast the earlier `qsp5` runs, where an
absorbing potential made the propagation non-unitary and broke this identity at
$t\approx5$ a.u.)

The four reported numbers are then $S_{ij} = -dT_i/ds_j$ for $i\in\{1,2\}$,
$j\in\{3,4\}$ — of which $S_{13}\!\approx\!S_{14}$ and $S_{23}\!\approx\!S_{24}$
by the identity, leaving $S_{1\ast}$ vs $S_{2\ast}$ as the real physical
contrast."""))
    else:
        cells.append(md(r"""---
## 2. What we are computing

$$ S_{\rm cl} \;=\; -\frac{dT}{dz},\qquad T=\tfrac12 m v_z^2 $$

For the classical projectile there is no ambiguity: it is a point-like Gaussian
charge with a definite position $z(t)$ and a definite velocity $v_z(t)$, both
integrated by INQ's Ehrenfest ion propagator. This single number is the
**reference** against which the wavepacket run's four KS-orbital definitions are
judged.

Two independent routes to $T$ are recorded and must agree:

* `ke_ion_ha` in `electron_track.csv` — computed from $v$ and the ion mass;
* `energy_ion_kinetic` in `observables.csv` — INQ's own bookkeeping.

Their agreement is a free consistency check on the ion integrator.

### This projectile decelerates, by design

It has the mass of an electron and only 100 eV of kinetic energy, so the
electronic drag slows it substantially over the traversal. Per
`.claude/rules/light-projectile-stopping.md` a regression over the *whole*
trajectory would average $S$ over every velocity between $v_0$ and whatever it
ends at — which is not $S$ at $v_0$. Two windows are therefore reported:

* the **shared** window $t\in[4,19]$ a.u., for comparability with the wavepacket twin;
* the **initial-drag** window $v_z \ge 0.85\,v_0$, which is the honest $S(v_0)$.

### Matching the wavepacket

The projectile is a Gaussian pseudopotential whose *charge* standard deviation is
$\sigma_{\rm WP}/\sqrt2 = 1.4142$ Bohr. That is exactly the density width of the
wavepacket at $t=0$, so both projectiles present the identical charge cloud
$\exp(-r^2/\sigma_{\rm WP}^2)$ to the bath and both are labelled $\sigma=2$
(`.claude/rules/sigma-wp-convention.md`)."""))

    # ------------------------------------------------------ load kinematics
    if is_wp:
        cells.append(code("""run = K.load_wp_run(RUN, box_length_z=CFG["LZ"], z0=CFG["Z0"])
print(f"{len(run.t)} steps, t = {run.t[0]:.2f} .. {run.t[-1]:.2f} a.u.")

t0v = dict(
    norm=run.norm[0], pz=run.pz[0],
    T1_eV=run.T1[0]*HA_TO_EV, T2_eV=run.T2[0]*HA_TO_EV,
    localisation_eV=run.localisation_energy[0],
    s3=run.s3[0], s4=run.s4[0], sigma_z=run.sigma_z[0],
)
print("\\n--- t = 0 against analytic expectations ---")
print(f"  norm (real)   = {t0v['norm']:.6f}   (expect 1)")
print(f"  norm drift    = {abs(run.norm[-1]-run.norm[0]):.3e}   (orbital stays normalised)")
print(f"  Parseval sum  = {run.parseval[0]:.4e}   (momentum-space FFT constant, NOT 1;"
      f" varies by {100*np.ptp(run.parseval)/run.parseval[0]:.2e} % over the run)")
print(f"  <p_z>         = {t0v['pz']:.6f}   (expect {CFG['K0']})")
print(f"  T1            = {t0v['T1_eV']:.3f} eV (expect 105.10)")
print(f"  T2            = {t0v['T2_eV']:.3f} eV (expect {CFG['EKIN_EV']:.2f})")
print(f"  T1-T2         = {t0v['localisation_eV']:.3f} eV (expect 5.102 = 3/(4 sigma^2))")
print(f"  sigma_z(circ) = {t0v['sigma_z']:.4f} Bohr (expect {CFG['SIGMA_WP']/np.sqrt(2):.4f})")
print(f"  s3, s4        = {t0v['s3']:.3f}, {t0v['s4']:.3f} Bohr (expect {CFG['Z0']})")"""))
    else:
        cells.append(code("""run = K.load_classical_run(RUN, box_length_z=CFG["LZ"])
print(f"{len(run.t)} steps, t = {run.t[0]:.2f} .. {run.t[-1]:.2f} a.u.")
print(f"  z  : {run.z[0]:.3f} -> {run.z[-1]:.3f} Bohr")
print(f"  v_z: {run.vz[0]:.4f} -> {run.vz[-1]:.4f} a.u.  ({run.v_fraction[-1]*100:.1f}% of v0)")
print(f"  T  : {run.T[0]*HA_TO_EV:.2f} -> {run.T[-1]*HA_TO_EV:.2f} eV"
      f"   (lost {(run.T[0]-run.T[-1])*HA_TO_EV:.2f} eV)")"""))

    # ---------------------------------------------------- KE terms vs time
    cells.append(md("""---
## 3. The kinetic-energy terms, individually

Each term on its own axes first, so nothing is hidden by a shared scale."""))

    if is_wp:
        cells.append(code("""fig, ax = plt.subplots(3, 1, figsize=(9, 9), sharex=True)
W = dict(color="0.85", zorder=0)

ax[0].axvspan(CFG["FIT_T0"], CFG["FIT_T1"], **W)
ax[0].plot(run.t, run.T1*HA_TO_EV, lw=1.4, color="C0")
ax[0].set_ylabel("$T_1=\\\\langle p^2\\\\rangle/2m$  [eV]")
ax[0].set_title("Kinetic-energy definitions, separately")

ax[1].axvspan(CFG["FIT_T0"], CFG["FIT_T1"], **W)
ax[1].plot(run.t, run.T2*HA_TO_EV, lw=1.4, color="C1")
ax[1].set_ylabel("$T_2=\\\\langle p\\\\rangle^2/2m$  [eV]")

ax[2].axvspan(CFG["FIT_T0"], CFG["FIT_T1"], **W)
ax[2].plot(run.t, run.localisation_energy, lw=1.4, color="C3")
ax[2].axhline(5.102, ls=":", color="k", lw=1)
ax[2].annotate("$3/4\\\\sigma^2$ = 5.102 eV (t=0 localisation energy)",
               xy=(run.t[-1], 5.102), ha="right", va="bottom", fontsize=8)
ax[2].set_ylabel("$T_1-T_2$  [eV]")
ax[2].set_xlabel("time  [a.u.]")
for a in ax: a.margins(x=0)
fig.tight_layout(); plt.show()

print(f"T1: {run.T1[0]*HA_TO_EV:8.3f} -> {run.T1[-1]*HA_TO_EV:8.3f} eV"
      f"   (change {(run.T1[-1]-run.T1[0])*HA_TO_EV:+.3f})")
print(f"T2: {run.T2[0]*HA_TO_EV:8.3f} -> {run.T2[-1]*HA_TO_EV:8.3f} eV"
      f"   (change {(run.T2[-1]-run.T2[0])*HA_TO_EV:+.3f})")
print(f"T1-T2: {run.localisation_energy[0]:.3f} -> {run.localisation_energy[-1]:.3f} eV")"""))

        cells.append(md("""**How to read the third panel.** If $T_1-T_2$ stays flat, the packet's momentum
spread is unchanged and any energy loss is pure drift deceleration — the quantum
and classical pictures should then agree. If it *grows*, the bath is broadening
the momentum distribution (elastic/angular scattering), and $S_1$ and $S_2$ will
separate by exactly that growth rate per unit path."""))
    else:
        cells.append(code("""fig, ax = plt.subplots(3, 1, figsize=(9, 9), sharex=True)
W = dict(color="0.85", zorder=0)

ax[0].axvspan(CFG["FIT_T0"], CFG["FIT_T1"], **W)
ax[0].plot(run.t, run.T*HA_TO_EV, lw=1.4, color="C0")
ax[0].set_ylabel("$T=\\\\frac{1}{2} m v_z^2$  [eV]")
ax[0].set_title("Classical projectile kinematics")

ax[1].axvspan(CFG["FIT_T0"], CFG["FIT_T1"], **W)
ax[1].plot(run.t, run.vz, lw=1.4, color="C1")
ax[1].axhline(0.85*run.vz[0], ls=":", color="k", lw=1)
ax[1].annotate("0.85 $v_0$ (initial-drag window edge)",
               xy=(run.t[-1], 0.85*run.vz[0]), ha="right", va="bottom", fontsize=8)
ax[1].set_ylabel("$v_z$  [a.u.]")

ax[2].axvspan(CFG["FIT_T0"], CFG["FIT_T1"], **W)
ax[2].plot(run.t, run.z, lw=1.4, color="C2")
ax[2].axhline(CFG["LZ"]/2, ls="--", color="C3", lw=1)
ax[2].annotate("+z face", xy=(run.t[0], CFG["LZ"]/2), va="bottom", fontsize=8)
ax[2].set_ylabel("$z$  [Bohr]"); ax[2].set_xlabel("time  [a.u.]")
for a in ax: a.margins(x=0)
fig.tight_layout(); plt.show()"""))

        cells.append(code("""# Cross-check T against INQ's own energy_ion_kinetic bookkeeping.
en = K.load_energies(RUN)
if "energy_ion_kinetic" in en.columns:
    m = np.interp(run.t, en["time_au"], en["energy_ion_kinetic"])
    dev = np.max(np.abs(m - run.T)) * HA_TO_EV
    print(f"max |ke_ion_ha - energy_ion_kinetic| = {dev:.3e} eV"
          f"  -> {'AGREE' if dev < 1e-3 else 'DISAGREE — investigate'}")
else:
    print("energy_ion_kinetic not in observables.csv")"""))

    # ------------------------------------------------ position terms vs time
    if is_wp:
        cells.append(md("""---
## 4. The position terms, individually

`s3` is the circular centroid, `s3_naive` the raw $\\int z|\\psi|^2$, `s4` the
integrated drift momentum."""))

        cells.append(code("""fig, ax = plt.subplots(3, 1, figsize=(9, 9), sharex=True)
W = dict(color="0.85", zorder=0)

ax[0].axvspan(CFG["FIT_T0"], CFG["FIT_T1"], **W)
ax[0].plot(run.t, run.s3, lw=1.4, color="C0", label="$s_3$ circular centroid")
ax[0].plot(run.t, run.s3_naive, lw=1.0, ls="--", color="C3",
           label=r"naive $\\int z|\\psi|^2$")
ax[0].axhline(CFG["LZ"]/2, ls=":", color="k", lw=1)
ax[0].annotate("+z face", xy=(run.t[0], CFG["LZ"]/2), va="bottom", fontsize=8)
ax[0].set_ylabel("position  [Bohr]"); ax[0].legend(fontsize=8, frameon=False)
ax[0].set_title("Position definitions, separately")

ax[1].axvspan(CFG["FIT_T0"], CFG["FIT_T1"], **W)
ax[1].plot(run.t, run.s4, lw=1.4, color="C1")
ax[1].set_ylabel(r"$s_4=z_0+\\int\\langle p_z\\rangle dt$  [Bohr]")

ax[2].axvspan(CFG["FIT_T0"], CFG["FIT_T1"], **W)
ax[2].plot(run.t, run.ehrenfest_residual, lw=1.4, color="C3")
ax[2].axhline(0, color="k", lw=0.8)
ax[2].set_ylabel("$s_3-s_4$  [Bohr]"); ax[2].set_xlabel("time  [a.u.]")
for a in ax: a.margins(x=0)
fig.tight_layout(); plt.show()

fitm = (run.t >= CFG["FIT_T0"]) & (run.t <= CFG["FIT_T1"])
res = np.max(np.abs(run.ehrenfest_residual[fitm]))
print(f"max |s3 - s4| inside the fit window = {res:.4f} Bohr")
print("  Ehrenfest identity", "HOLDS" if res < 0.1 else "VIOLATED — investigate")
print(f"max |s3 - s3_naive| over whole run  = "
      f"{np.max(np.abs(run.s3 - run.s3_naive)):.3f} Bohr")"""))

        cells.append(md("""**The naive-vs-circular gap is the whole reason the circular estimator exists.**
While the packet is well inside the cell the two curves are indistinguishable.
Once the leading tail reaches the $+z$ face the naive integral starts averaging
the wrapped part of the packet against the rest and bends *away* from the true
trajectory — smoothly, with no discontinuity to warn you."""))

        cells.append(md("""### Packet width — is the interference-free window where we claimed?

The fit window's upper edge was derived from free-Gaussian dispersion,
$\\sigma_d(t)=\\sqrt{\\sigma^2/2+t^2/2\\sigma^2}$, as the time when the leading
$3\\sigma_d$ tail reaches the $+z$ face. The run reports its own periodic-safe
width (`sigma_z_circ`), so that prediction is directly testable."""))

        cells.append(code("""sd_free = np.sqrt(CFG["SIGMA_WP"]**2/2 + run.t**2/(2*CFG["SIGMA_WP"]**2))

fig, ax = plt.subplots(figsize=(9, 4))
ax.axvspan(CFG["FIT_T0"], CFG["FIT_T1"], color="0.85", zorder=0)
ax.plot(run.t, run.sigma_z, lw=1.5, color="C0", label=r"measured $\\sigma_z$ (circular)")
ax.plot(run.t, sd_free, lw=1.2, ls="--", color="k", label=r"free dispersion $\\sigma_d(t)$")
ax.axhline(CFG["LX"]/6, ls=":", color="C3")
ax.annotate("$L_{xy}/6$: transverse images touch",
            xy=(run.t[0], CFG["LX"]/6), va="bottom", fontsize=8, color="C3")
ax.axvline(CFG["T_IFW"], ls="-.", color="C1")
ax.annotate("longitudinal IFW", xy=(CFG["T_IFW"], 1), rotation=90,
            fontsize=8, color="C1", ha="right")
ax.set_xlabel("time  [a.u.]"); ax.set_ylabel(r"$\\sigma_z$  [Bohr]")
ax.legend(fontsize=8, frameon=False); ax.margins(x=0)
fig.tight_layout(); plt.show()

i = np.argmin(np.abs(run.t - CFG["FIT_T1"]))
print(f"at t = {run.t[i]:.2f} a.u. (fit-window edge):")
print(f"   measured sigma_z = {run.sigma_z[i]:.2f} Bohr")
print(f"   free dispersion  = {sd_free[i]:.2f} Bohr")
print(f"   ratio measured/free = {run.sigma_z[i]/sd_free[i]:.3f}"
      "   (<1 => the bath CONFINES the packet relative to vacuum)")"""))

    # ----------------------------------------------------- stopping power
    cells.append(md("""---
## 5. Extracting the stopping power, step by step

The recipe, in order:

1. **Restrict to the analysis window.** Drop $t<4$ a.u. (the launch transient:
   the packet is orthogonalised against the occupied manifold at $t=0$ and the
   bath has not yet responded) and $t>t_{\\rm IFW}$ (where the projectile begins
   to interfere with its own periodic images).
2. **Plot $T$ against $s$** — not against $t$. The stopping power is energy per
   unit *path*, and the projectile does not move at constant speed.
3. **Fit a straight line** by ordinary least squares. The slope is $dT/ds$ in
   Ha/Bohr; $S=-\\,dT/ds$, converted to eV/Bohr.
4. **Check linearity** ($r^2$). If $T(s)$ is visibly curved, a single $S$ does not
   describe the run and the curvature should be reported instead.
5. **Price the window choice.** Both window edges are judgement calls, so they
   are moved independently by $\\pm3$ a.u. and half the resulting spread in $S$ is
   quoted as a systematic. This is normally much larger than the fit's own
   standard error, and it is the number that should be believed."""))

    if is_wp:
        cells.append(code("""fits = K.fit_all_wp(run, CFG["FIT_T0"], CFG["FIT_T1"])
for k, f in fits.items():
    print(f.summary())"""))

        cells.append(code("""fig, axes = plt.subplots(2, 2, figsize=(11, 8))
for ax, (key, f) in zip(axes.flat, fits.items()):
    ax.plot(f.s_fit, f.T_fit*HA_TO_EV, ".", ms=3, color="C0", label="data")
    ax.plot(f.s_fit, f.T_model*HA_TO_EV, "-", lw=1.6, color="C3", label="OLS fit")
    ax.set_title(f.label, fontsize=9)
    ax.set_xlabel("path $s$  [Bohr]"); ax.set_ylabel("$T$  [eV]")
    ax.annotate(f"S = {f.S_ev_per_bohr:.2f} $\\\\pm$ {f.uncertainty:.2f} eV/Bohr\\n"
                f"$r^2$ = {f.r2:.4f}",
                xy=(0.04, 0.06), xycoords="axes fraction", fontsize=8,
                bbox=dict(fc="w", ec="0.7", alpha=0.9))
    ax.legend(fontsize=7, frameon=False)
fig.suptitle("Stopping power: four KS-orbital definitions", y=1.00)
fig.tight_layout(); plt.show()"""))

        cells.append(code("""rows = []
for k, f in fits.items():
    rows.append(dict(combination=k, definition=f.label.split("(")[1].rstrip(")"),
                     S_eV_per_Bohr=round(f.S_ev_per_bohr, 2),
                     stat=round(f.stderr, 2), syst=round(f.window_syst, 2),
                     total_unc=round(f.uncertainty, 2), r2=round(f.r2, 4),
                     n=f.n_points))
tab = pd.DataFrame(rows)
display(tab)

print("\\n--- the two comparisons that matter ---")
print(f"position definitions (should agree, Ehrenfest):")
print(f"   S_13 vs S_14 : {fits['S_13'].S_ev_per_bohr:.3f} vs "
      f"{fits['S_14'].S_ev_per_bohr:.3f}  "
      f"(diff {abs(fits['S_13'].S_ev_per_bohr-fits['S_14'].S_ev_per_bohr):.3f})")
print(f"   S_23 vs S_24 : {fits['S_23'].S_ev_per_bohr:.3f} vs "
      f"{fits['S_24'].S_ev_per_bohr:.3f}  "
      f"(diff {abs(fits['S_23'].S_ev_per_bohr-fits['S_24'].S_ev_per_bohr):.3f})")
print(f"\\nKE definitions (the PHYSICS contrast):")
d = fits['S_13'].S_ev_per_bohr - fits['S_23'].S_ev_per_bohr
print(f"   S_1* - S_2*  : {d:+.3f} eV/Bohr")
print("   => " + ("momentum spread GROWS along the path: part of the apparent\\n"
                  "      stopping is angular scattering, not drift loss."
                  if d < 0 else
                  "momentum spread SHRINKS/flat: the loss is drift deceleration."))"""))

        # ------------------------- 5b. the four S values on one plot --------
        cells.append(md(r"""---
### 5b. The four stopping powers on one plot — and what they are measured against

The `stopping-power-extraction` skill fixes the headline method by run geometry:
for a **continuous traversal** through a homogeneous medium the answer is
*Method A* — the slope of the deposited energy $\Delta E_{\rm total}(x)$ after a
fixed 20%-of-time transient cut — and $-dT/dx$ is only the conservation
cross-check.

**That hierarchy inverts here, and the reason is worth stating precisely.** In the
wavepacket run the projectile *is* an occupied Kohn–Sham orbital, so its energy is
already inside `energy_total`. The system is closed, nothing absorbs, and
$E_{\rm total}$ is constant to $2.6\times10^{-4}$ eV over the whole run — there is
no deposit curve to fit. Method A is not merely inaccurate here, it is
**undefined**: the quantity it fits is identically zero by construction.

So for this half the four $-dT_i/ds_j$ slopes *are* the measurement, which is
exactly why the KS-orbital-dependent definitions were requested. The external
reference comes from the **classical twin**, where the projectile sits outside
`energy_total` and Method A applies normally. Below we compute that reference with
the skill's own kernels, on the *same* fit window, so the comparison is like for
like."""))

        cells.append(code(r"""import stopping_power as SP   # the skill's reference kernels

CLASSICAL = RUN.parent / "classical"
ref = {}
try:
    cl  = K.load_classical_run(CLASSICAL, box_length_z=CFG["LZ"])
    enc = K.load_energies(CLASSICAL)

    # Method A signal: deposited electronic energy vs projectile path.
    t_e = enc["time_au"].to_numpy()
    dE  = (enc["energy_total"] - enc["energy_total"].iloc[0]).to_numpy() * HA_TO_EV
    x_e = np.interp(t_e, cl.t, cl.z) - cl.z[0]

    # (i) the skill's locked default: fixed 20%-of-time transient cut, free intercept
    fA = SP.fixed_time_fraction(t_e, x_e, dE, frac=0.20)

    # (ii) the same deposit fit restricted to the WP fit window, for like-for-like
    mE = (t_e >= CFG["FIT_T0"]) & (t_e <= CFG["FIT_T1"])
    fW = SP.free_fit(x_e[mE], dE[mE], x_e[mE].min(), x_e[mE].max())

    # (iii) sanity channel: -dKE_ion/dx over the same window (independent CSV)
    mk   = (cl.t >= CFG["FIT_T0"]) & (cl.t <= CFG["FIT_T1"])
    xk   = cl.z[mk] - cl.z[0]
    dKE  = (cl.T[mk][0] - cl.T[mk]) * HA_TO_EV
    fK   = SP.free_fit(xk, dKE, xk.min(), xk.max())

    ref = dict(methodA=fA["S"], methodA_se=fA["se"], methodA_r2=fA["r2"],
               window=fW["S"], window_se=fW["se"], window_r2=fW["r2"],
               kinetic=fK["S"], kinetic_se=fK["se"])

    print("CLASSICAL TWIN — the external reference (skill kernels)")
    print(f"  Method A, fixed 20% time cut  S = {fA['S']:.3f} +/- {fA['se']:.3f} eV/Bohr"
          f"   (r2 = {fA['r2']:.4f}, x0 = {fA['x0']:.1f} Bohr, status {fA['status']})")
    print(f"  Method A on the WP window     S = {fW['S']:.3f} +/- {fW['se']:.3f} eV/Bohr"
          f"   (r2 = {fW['r2']:.4f})")
    print(f"  sanity: -dKE_ion/dx           S = {fK['S']:.3f} +/- {fK['se']:.3f} eV/Bohr")
    agree = 100.0 * abs(fW["S"] - fK["S"]) / fW["S"]
    print(f"\n  deposit vs kinetic channel agree to {agree:.2f}%  ->  "
          + ("energy conservation confirmed; the classical reference is trustworthy."
             if agree < 5 else "DISAGREEMENT — investigate before using this reference."))
except FileNotFoundError as exc:
    print("classical twin not available, reference band omitted:", exc)

# The deposit method on THIS run, to show it is empty rather than merely small.
enw = K.load_energies(RUN)
dEw = (enw["energy_total"] - enw["energy_total"].iloc[0]).to_numpy() * HA_TO_EV
print(f"\nWAVEPACKET RUN, deposit channel: max |dE_total| = {np.max(np.abs(dEw)):.2e} eV")
print(f"  compare the projectile's own KE loss over the window: "
      f"{(run.T2[(run.t>=CFG['FIT_T0'])&(run.t<=CFG['FIT_T1'])][0] - run.T2[(run.t>=CFG['FIT_T0'])&(run.t<=CFG['FIT_T1'])][-1])*HA_TO_EV:.3f} eV")
print("  => Method A is UNDEFINED for this half (closed system). Use -dT/ds.")"""))

        cells.append(code(r"""fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8),
                         gridspec_kw=dict(width_ratios=[1.15, 1.0]))

# ---- left: the slopes you are actually fitting, all four on one axes -------
ax = axes[0]
STYLE = {"S_13": ("C0", "-"), "S_14": ("C0", "--"),
         "S_23": ("C3", "-"), "S_24": ("C3", "--")}
for key, f in fits.items():
    c, ls = STYLE[key]
    s_rel = f.s_fit - f.s_fit[0]
    dT    = (f.T_fit - f.T_fit[0]) * HA_TO_EV
    ax.plot(s_rel, dT, ls, color=c, lw=1.0, alpha=0.45)
    ax.plot(s_rel, (f.T_model - f.T_model[0]) * HA_TO_EV, ls, color=c, lw=2.0,
            label=f"{key}:  S = {f.S_ev_per_bohr:.3f}")
if ref:
    s_ax = np.linspace(0, max(f.s_fit[-1]-f.s_fit[0] for f in fits.values()), 50)
    ax.plot(s_ax, -ref["window"]*s_ax, ":", color="k", lw=2.0,
            label=f"classical twin: {ref['window']:.3f}")
ax.set_xlabel("path travelled inside the window,  $s-s_0$  [Bohr]")
ax.set_ylabel(r"$\Delta T$  [eV]")
ax.set_title("The fitted slopes ($S=-dT/ds$), all definitions")
ax.legend(fontsize=8, frameon=False, loc="lower left")

# ---- right: the four S values with their honest error bars ----------------
ax = axes[1]
keys = ["S_13", "S_14", "S_23", "S_24"]
lbl  = [r"$S_{13}$" "\n" r"$\langle p^2\rangle/2m$" "\n" "centroid",
        r"$S_{14}$" "\n" r"$\langle p^2\rangle/2m$" "\n" r"$\int\langle p\rangle dt$",
        r"$S_{23}$" "\n" r"$\langle p\rangle^2/2m$" "\n" "centroid",
        r"$S_{24}$" "\n" r"$\langle p\rangle^2/2m$" "\n" r"$\int\langle p\rangle dt$"]
y    = [fits[k].S_ev_per_bohr for k in keys]
yerr = [fits[k].uncertainty   for k in keys]
cols = ["C0", "C0", "C3", "C3"]
xpos = np.arange(4)
ax.bar(xpos, y, color=cols, alpha=0.35, width=0.62)
ax.errorbar(xpos, y, yerr=yerr, fmt="o", color="k", ms=4, lw=1.4, capsize=4)
for xi, (yi, ei) in enumerate(zip(y, yerr)):
    ax.annotate(f"{yi:.3f}\n$\\pm${ei:.3f}", xy=(xi, yi + ei), ha="center",
                va="bottom", fontsize=8)
if ref:
    ax.axhline(ref["window"], ls=":", color="k", lw=1.6)
    ax.annotate(f"classical twin, same window: {ref['window']:.3f} eV/Bohr",
                xy=(3.42, ref["window"]), ha="right", va="bottom", fontsize=8)
ax.set_xticks(xpos); ax.set_xticklabels(lbl, fontsize=7.5)
ax.set_ylabel("$S$  [eV/Bohr]")
ax.set_title("Stopping power — four KS-orbital definitions")
ax.set_ylim(0, max(max(np.array(y)+np.array(yerr)),
                   ref.get("window", 0)) * 1.28)
fig.tight_layout()
(HERE / "figures").mkdir(exist_ok=True)
fig.savefig(HERE / "figures" / "four_stopping_powers.png",
            dpi=150, bbox_inches="tight")
plt.show()

print("saved -> figures/four_stopping_powers.png")
if ref:
    print(f"\nratio classical / S_2* = {ref['window']/fits['S_23'].S_ev_per_bohr:.1f}x")
    print(f"ratio classical / S_1* = {ref['window']/fits['S_13'].S_ev_per_bohr:.1f}x")"""))

        cells.append(md(r"""**Reading the two panels.**

*Left* — the four $\Delta T(s)$ series, shifted to a common origin so the slopes
can be compared directly. The two blue curves ($T_1=\langle p^2\rangle/2m$) lie
almost on top of each other, as do the two red ($T_2=\langle p\rangle^2/2m$): the
**choice of position definition changes nothing**, which is the Ehrenfest identity
$d\langle z\rangle/dt=\langle p_z\rangle/m$ showing up in the final answer. The
blue–red separation is the real result, and it is a factor of ~4.

*Right* — the same four numbers with the uncertainty that should actually be
quoted (fit standard error and window-sensitivity systematic in quadrature).

**Why $S_1 \ne S_2$.** $T_1-T_2=\tfrac32\sum_d\sigma_{p_d}^2$ is the momentum-width
energy. $T_1$ counts it, $T_2$ does not. As the packet propagates its momentum
distribution broadens, so $T_1$ falls *more slowly* than $T_2$ — part of the drift
energy has been converted into spread rather than lost to the bath. $S_2$ is
therefore the **drift deceleration**, and $S_1$ is the drift deceleration *minus*
the growth of the localisation term. Neither is wrong; they answer different
questions, and quoting one without saying which is the error to avoid.

**Caveat on $r^2$.** The $S_{1j}$ fits have $r^2\approx0.66$ against
$r^2\approx0.999$ for $S_{2j}$. $T_1(s)$ is not straight over this window, so its
single-slope value is a window-averaged number rather than a well-defined
gradient — read it with that in mind, and prefer the $T_2$ channel for the
headline drift stopping."""))

    else:
        cells.append(code("""fit_shared = K.fit_classical(run, CFG["FIT_T0"], CFG["FIT_T1"])
fit_early  = K.fit_classical_early(run)
print(fit_shared.summary())
print(fit_early.summary())

fig, axes = plt.subplots(1, 2, figsize=(11, 4))
for ax, f in zip(axes, (fit_shared, fit_early)):
    ax.plot(f.s_fit, f.T_fit*HA_TO_EV, ".", ms=3, color="C0", label="data")
    ax.plot(f.s_fit, f.T_model*HA_TO_EV, "-", lw=1.6, color="C3", label="OLS fit")
    ax.set_title(f.label, fontsize=9)
    ax.set_xlabel("path $z$  [Bohr]"); ax.set_ylabel("$T$  [eV]")
    ax.annotate(f"S = {f.S_ev_per_bohr:.2f} $\\\\pm$ {f.uncertainty:.2f} eV/Bohr\\n"
                f"$r^2$ = {f.r2:.4f},  mean v = {f.mean_v:.3f}",
                xy=(0.04, 0.06), xycoords="axes fraction", fontsize=8,
                bbox=dict(fc="w", ec="0.7", alpha=0.9))
    ax.legend(fontsize=7, frameon=False)
fig.tight_layout(); plt.show()"""))

    # -------------------------------------------------- energy decomposition
    cells.append(md("""---
## 6. Energy decomposition and the conservation check

Every component INQ tracks. The first eight sum to the total; `nvxc` and
`eigenvalues` are diagnostics that sit outside it.

There is no absorbing potential anywhere in this study, so nothing removes energy
and the bookkeeping must close. What "close" means differs between the halves:

* **Wavepacket run** — the projectile *is* an electron orbital, so it sits inside
  `energy_total`. The system is closed and the total must be **constant**; its
  drift is pure numerical error, and every stopping number inherits it.
* **Classical run** — INQ leaves `energy_ion_kinetic` at zero (verified from the
  run output), so `energy_total` is the **electronic** energy alone. It is
  *supposed* to rise, by exactly the kinetic energy the projectile gives up.
  Reading that rise as "drift" would be a mistake; the real test is the
  **closure** dE_electronic = -dT_projectile."""))

    cells.append(code("""en = K.load_energies(RUN)
if HALF == "classical":
    # energy_ion_kinetic is left at ZERO by INQ (verified), so energy_total here
    # is the ELECTRONIC energy only. It is SUPPOSED to rise, by exactly the
    # kinetic energy the projectile gives up: the test is closure, not flatness.
    cons = K.conservation_check(en, projectile_ke_loss_ev=(run.T[0]-run.T[-1])*HA_TO_EV)
    for k, v in cons.items():
        print(f"  {k:26s} {v: .6f}")
    mm = cons["closure_mismatch_pct"]
    verdict = ("EXCELLENT" if mm < 1 else "acceptable" if mm < 5 else
               "POOR — energy is not accounted for")
    print(f"\\n  electron/ion energy closure: {verdict} ({mm:.2f}% of the transfer)")
else:
    # The wavepacket IS inside energy_total, so the system is closed and the
    # total must be CONSTANT. Its drift is pure numerical error.
    cons = K.conservation_check(en)
    for k, v in cons.items():
        print(f"  {k:26s} {v: .6f}")
    verdict = ("EXCELLENT" if abs(cons["drift_ev"]) < 0.01 else
               "acceptable" if abs(cons["drift_ev"]) < 0.1 else
               "POOR — dt or dx too coarse; treat S with caution")
    print(f"\\n  energy conservation: {verdict}")

comp = [c for c in ("energy_kinetic","energy_hartree","energy_xc","energy_external",
                    "energy_nonlocal","energy_ion","energy_ion_kinetic",
                    "energy_exact_exchange") if c in en.columns]

fig, ax = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
ax[0].axvspan(CFG["FIT_T0"], CFG["FIT_T1"], color="0.85", zorder=0)
ax[0].plot(en["time_au"], en["delta_e_total_ev"], lw=1.4, color="k")
ax[0].set_ylabel(r"$E_{\\rm total}(t)-E_{\\rm total}(0)$  [eV]")
ax[0].set_title("Total electronic energy relative to t=0")

ax[1].axvspan(CFG["FIT_T0"], CFG["FIT_T1"], color="0.85", zorder=0)
for c in comp:
    y = (en[c] - en[c].iloc[0]) * HA_TO_EV
    if np.max(np.abs(y)) > 1e-6:
        ax[1].plot(en["time_au"], y, lw=1.2, label=c.replace("energy_", ""))
ax[1].set_xlabel("time  [a.u.]"); ax[1].set_ylabel("component $-$ its $t{=}0$ value  [eV]")
ax[1].legend(fontsize=8, ncol=2, frameon=False)
ax[1].set_title("Energy decomposition, each relative to its own start")
for a in ax: a.margins(x=0)
fig.tight_layout(); plt.show()"""))

    # ------------------------------- interaction-energy section (decomposed-interaction-energies rule)
    cells.append(md(r"""---
## 7. Pairwise interaction energies — the comparison INQ's scalars cannot make

Section 6 showed every energy component INQ tracks. None of them can be compared
against the twin, because the two representations put the projectile in
**different ledger terms**:

| | classical | wavepacket |
|---|---|---|
| projectile enters as | external potential | occupied KS orbital |
| `energy_external` | non-zero | identically 0 |
| `energy_hartree` | $E_{SS}$ | $E_{SS}+E_{PS}+E_{PP}$ |

So a raw `energy_hartree` comparison puts a **net** quantity against a **gross**
one. The pairwise decomposition below is representation-independent and *is*
comparable. Splitting the charge into P (projectile), S (bath electrons) and
B (neutralising background):

$$E_{SS}=\tfrac12\!\int n_S\phi_S \qquad
  E_{PP}=\tfrac12\!\int n_P\phi_P \qquad
  E_{PS}=\int n_S\phi_P$$

* $E_{PS}$ is the projectile–bath interaction — **the term that stops it**.
* $E_{PP}$ is the projectile **self-Hartree**. A wavepacket is an occupied KS
  orbital, so in LDA it feels its own Hartree field with no exact-exchange
  cancellation. A classical projectile is an external potential and has no such
  term — for it $E_{PP}$ is a *rigid constant of the motion*, and its constancy
  is a free validation of the whole pipeline.

**This is bulk.** The background is uniform, so $\phi_+$ is pure $G{=}0$, which
INQ drops: $\phi_+\equiv0$ and $E_{SB}=E_{PB}=E_{BB}=0$ identically. They are
written as columns anyway so the schema matches the slab systems. All the physics
here is in $E_{SS}$, $E_{PP}$, $E_{PS}$.

**Gauge caveat.** Absolute $E_{PP}$ carries the charged-cell $G{=}0$ gauge. Only
closure sums and WP-minus-classical differences *within the same cell* are
gauge-clean — never quote an absolute $E_{PP}$ across geometries."""))

    cells.append(code("""ix = K.load_interactions(RUN, HALF)

# CLOSURE GATE (.claude/rules/decomposed-interaction-energies.md). The terms must
# sum back to INQ's own Hartree energy; this is what makes them trustworthy.
resid = np.nanmax(np.abs(ix.closure))
print(f"  rows                : {len(ix.t)}  (steps {ix.step.min()}..{ix.step.max()})")
print(f"  closure vs INQ E_H  : max|resid| = {resid:.2e} Ha  "
      f"[{'PASS' if resid < 1e-9 else 'FAIL'}]")
for nm, v in (("E_SB", ix.e_sb), ("E_PB", ix.e_pb), ("E_BB", ix.e_bb)):
    print(f"  {nm} == 0 (bulk)     : max|{nm}| = {np.nanmax(np.abs(v)):.2e} Ha")

if HALF == "classical":
    # E_PP of a RIGID Gaussian cloud is a constant of the motion -- but only
    # while the cloud is fully inside the box. Near the +z face its tail is
    # clipped off the grid, charge is lost and E_PP decays. That is not physics,
    # and it is a hard upper bound on any fit window using this run.
    clean = ix.norm >= 1.0 - K.CLIP_TOL
    spread = np.nanmax(ix.e_pp[clean]) - np.nanmin(ix.e_pp[clean])
    print(f"  E_PP constant       : spread = {spread:.2e} Ha over {clean.sum()}"
          f"/{len(ix.e_pp)} unclipped rows  [{'PASS' if spread < 1e-9 else 'FAIL'}]")
    if np.isfinite(ix.clip_time):
        print(f"  cloud clipping onset: t = {ix.clip_time:.2f} a.u. "
              f"(norm_proj -> {np.nanmin(ix.norm):.6f})")
        print(f"                        fit window ends {CFG['FIT_T1']:.2f} -> "
              f"{'CLEAR' if CFG['FIT_T1'] < ix.clip_time else 'CONTAMINATED'}")
    else:
        print("  cloud clipping onset: never (projectile stays inside the box)")
else:
    print(f"  norm_wp             : {ix.norm.min():.6f} .. {ix.norm.max():.6f} "
          f"(no CAP in this study, so it must stay ~1)")"""))

    cells.append(md(r"""### The three active terms

Each shown relative to its own $t{=}0$ value, so the gauge offset cancels and
only the *change* — the physics — is on the axis."""))

    # NOTE: raw string. This cell has an f-string containing \n; in a NON-raw
    # builder string Python turns that into a REAL newline and the emitted cell
    # is an unterminated string literal. Cost one full array job (2026-08-01).
    cells.append(code(r"""fig, ax = plt.subplots(3, 1, figsize=(9, 9), sharex=True)
for a in ax:
    a.axvspan(CFG["FIT_T0"], CFG["FIT_T1"], color="0.85", zorder=0)
    if HALF == "classical" and np.isfinite(ix.clip_time):
        a.axvline(ix.clip_time, color="C3", ls=":", lw=1.4)

for c, lab, col in (("e_ss", r"$E_{SS}$ bath-bath", "C0"),
                    ("e_ps", r"$E_{PS}$ projectile-bath", "C1"),
                    ("e_pp", r"$E_{PP}$ projectile self", "C2")):
    y = getattr(ix, c)
    ax[0].plot(ix.t, (y - y[0]) * HA_TO_EV, lw=1.4, color=col, label=lab)
ax[0].set_ylabel(r"$E - E(0)$  [eV]")
ax[0].set_title("Pairwise terms, each relative to its own start")
ax[0].legend(fontsize=8, frameon=False)

ax[1].plot(ix.t, ix.e_pp * HA_TO_EV, lw=1.5, color="C2")
ax[1].set_ylabel(r"$E_{PP}$  [eV]")
ax[1].set_title(r"$E_{PP}$ absolute — flat for a rigid classical cloud, "
                "collapsing for a dispersing packet")

ax[2].plot(ix.t, ix.norm, lw=1.4, color="k")
ax[2].axhline(1.0, color="0.6", lw=0.8, ls="--")
ax[2].set_ylabel("norm$_{wp}$" if HALF == "wp" else "norm$_{proj}$")
ax[2].set_xlabel("time  [a.u.]")
ax[2].set_title("Projectile charge on the grid — the clipping diagnostic")

if HALF == "classical" and np.isfinite(ix.clip_time):
    ax[2].annotate(f"cloud clips the +z face\nt = {ix.clip_time:.1f} a.u.",
                   xy=(ix.clip_time, 1.0), xytext=(0.62, 0.35),
                   textcoords="axes fraction", fontsize=8, color="C3",
                   arrowprops=dict(arrowstyle="->", color="C3", lw=1.0))
for a in ax: a.margins(x=0)
fig.tight_layout(); plt.show()"""))

    cells.append(md(r"""**How to read this.** The shaded band is the stopping fit
window; the dotted red line (classical only) is where the projectile's Gaussian
cloud starts leaving the grid.

* **Classical.** $E_{PP}$ must be a horizontal line until the clipping onset —
  the cloud is rigid and translation-invariant, so its self-energy cannot change.
  Any drift *inside* the window would mean a grid artefact (egg-box error)
  contaminating the run.
* **Wavepacket.** $E_{PP}$ *falls*, because the packet disperses: a free Gaussian
  spreads as $\sigma_d(t)=\sqrt{\sigma^2/2+t^2/2\sigma^2}$, and a more spread-out
  charge has less self-energy. That decay is the quantum channel with no
  classical counterpart, and it is what the per-pair comparison notebook
  quantifies as the gauge-clean $\Delta E_{PP}$."""))

    # ---------------------------------------------------------- takeaway
    cells.append(code(f"""summary_out = dict(half=HALF, cfg=CFG, conservation=cons)
{"summary_out['fits'] = {k: dict(S_eV_per_Bohr=f.S_ev_per_bohr, stat=f.stderr, syst=f.window_syst, r2=f.r2, n=f.n_points) for k, f in fits.items()}"
 if is_wp else
 "summary_out['fits'] = {'S_cl_shared': dict(S_eV_per_Bohr=fit_shared.S_ev_per_bohr, stat=fit_shared.stderr, syst=fit_shared.window_syst, r2=fit_shared.r2), 'S_cl_initial': dict(S_eV_per_Bohr=fit_early.S_ev_per_bohr, stat=fit_early.stderr, syst=fit_early.window_syst, r2=fit_early.r2, mean_v=fit_early.mean_v)}"}

{"summary_out['classical_reference'] = ref  # skill Method A + kinetic channel" if is_wp else ""}
out = HERE / f"{{HALF}}_stopping_summary.json"
out.write_text(json.dumps(summary_out, indent=2, default=float))
print("wrote", out)
print(json.dumps(summary_out["fits"], indent=2, default=float))"""))

    cells.append(md(f"""---
## 8. Takeaway and limitations

**What this notebook establishes.** {
"Four stopping powers from the two KS-orbital kinetic-energy definitions crossed "
"with the two position definitions, plus the Ehrenfest cross-check that collapses "
"the position axis to one degree of freedom. The physically meaningful spread is "
"between $S_1$ and $S_2$: it measures how much apparent stopping is momentum-space "
"broadening rather than drift-momentum loss."
if is_wp else
"The classical reference stopping power, in two windows: the shared window for "
"comparability with the wavepacket twin, and the initial-drag window that is the "
"honest $S(v_0)$ for a decelerating light projectile."}

**Limitations, in order of importance.**

1. **Not a steady-state stopping power.** The plasma period ({CFG['PLASMA_PERIOD']} a.u.)
   exceeds the run ({CFG['DT']*CFG['N_STEPS']:.1f} a.u.). The bath never completes one
   oscillation, so no steady wake forms. This is an *initial-drag* number.
2. **Finite box.** $L_{{xy}}={CFG['LX']:.0f}$ Bohr was chosen so transverse periodic
   images stay separated until $t={CFG['T_TRANSVERSE']:.1f}$ a.u., past the longitudinal
   limit of {CFG['T_IFW']:.1f} a.u. — but the projectile is still one member of an
   infinite transverse array, as it is in any supercell calculation.
3. **The window is a judgement call**, which is why the systematic from moving its
   edges is reported alongside every $S$ and is usually the dominant uncertainty.
4. **ALDA.** Adiabatic LDA has no memory; dynamical xc effects on stopping are
   outside this model.
5. **One velocity, one width.** $S$ at 100 eV and $\\sigma=2$ Bohr only. Nothing
   here constrains the $v$- or $\\sigma$-dependence.

**Provenance.** Arithmetic by `hypotheses/bulk_ks_stopping/ks_stopping.py`
(11 known-case tests, `tests/test_ks_stopping.py`). Geometry constants pinned by
`static_assert` in `shared/configs/bulk_ks_stopping_L40x40x80_rs4.hpp`. Ground state
validated against the analytic plane-wave result before either run started.
"""))

    nb = nbf.v4.new_notebook(cells=cells)
    nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python",
                                 "name": "python3"}
    return nb


def main() -> int:
    if len(sys.argv) != 2 or sys.argv[1] not in ("wp", "classical"):
        print(__doc__)
        return 2
    half = sys.argv[1]
    nb = build(half)
    out = HERE / f"bulk_ks_stopping_rs4_{half}.ipynb"
    nbf.write(nb, str(out))
    print(f"wrote {out} ({len(nb.cells)} cells)")

    print("executing...")
    client = NotebookClient(nb, timeout=3600, kernel_name="python3",
                            resources={"metadata": {"path": str(HERE)}},
                            allow_errors=True)
    client.execute()
    nbf.write(nb, str(out))

    n_err = sum(1 for c in nb.cells if c.cell_type == "code"
                for o in c.get("outputs", []) if o.get("output_type") == "error")
    print(f"executed with {n_err} error(s) -> {out}")
    return 1 if n_err else 0


if __name__ == "__main__":
    raise SystemExit(main())
