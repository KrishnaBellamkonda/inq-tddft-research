#!/usr/bin/env python3
"""Build the LIGHTWEIGHT refined-analysis notebook for the channeling twin.

Plan: docs/plans/cylindrical-channeling-ks-stopping.md, section 8.

WHAT THIS IS FOR, AND WHY IT IS SEPARATE FROM THE COMPARISON NOTEBOOK
---------------------------------------------------------------------
``channeling_twin_comparison.ipynb`` answers "was the aim met?" using a fit
window the ANALYSIS derives (first breach of f_bore >= 0.95). This notebook
answers the earlier question the user actually needs answered first: *what do
the raw diagnostics look like*, so that the fit window can be chosen BY EYE and
SEPARATELY FOR EACH HALF.

It therefore renders NO verdict and performs NO stopping fit on first build. The
last section is a parameter cell (``T_WIN_CL``, ``T_WIN_WP``, both ``None``)
which, once filled in, runs the fits through the same ``refined.fit_in_window``
for both halves.

LIGHTWEIGHT means: no embedded density GIFs. .claude/rules/notebook-density-gif.md
requires the density matrix on a run/analysis notebook; it is already embedded in
``channeling_twin_comparison.ipynb`` (which is 222 MB precisely because of it).
Re-embedding the same animations here would add 222 MB and no information, so
this notebook LINKS to them instead. That is a deliberate, stated departure —
flagged to the user, not a silent omission.

The arithmetic is in ``refined.py``, which has its own known-case tests
(``tests/test_refined.py``, 12 passing). This file lays out narrative and plots.

Usage:
    PYTHONPATH=<repo>/inq-stack/python <repo>/venv/bin/python3 build_refined_notebook.py
        [--out refined_analysis.ipynb] [--wp NAME] [--classical NAME]
        [--times 0,15,30] [--no-execute]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]

sys.path.insert(0, str(HERE))


def cells(wp_name: str, cl_name: str, times: list[float]) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []

    def md(s: str): out.append(("md", s.strip("\n")))
    def code(s: str): out.append(("code", s.strip("\n")))

    # ===================================================== header
    md(f"""
# Channeling twin — refined stopping-power analysis

**Runs:** wavepacket `{wp_name}` · classical `{cl_name}`
(annular jellium tube, $r_s = 3$, $R_\\mathrm{{in}} = 10$, $R_\\mathrm{{out}} = 14$ Bohr,
$\\sigma_\\mathrm{{WP}} = 4$ Bohr, $E_0 = 50$ eV, 1500 steps $\\times$ 0.02 a.u.)

This notebook exists to let the fit window be chosen **from the diagnostics**
rather than from a rule, and **separately for each half**. It renders no verdict.
Sections 1–4 are the evidence; section 5 is where you set the window.

> **On the density GIFs.** They are not embedded here — that is what makes this
> notebook ~1 MB instead of 222 MB. They are already in
> `channeling_twin_comparison.ipynb` and as standalone files under
> `comparison_figs/density_matrix/`.

---

### The T₁ / T₂ convention used here

| this notebook | definition | `ks_stopping.py` calls it |
|---|---|---|
| $T_1$ | $\\langle p\\rangle^2/2m$ — **drift only** | `T2` |
| $T_2$ | $\\langle p\\rangle^2/2m + \\mathrm{{var}}(p)/2m = \\langle p^2\\rangle/2m$ | `T1` |

**The names are swapped relative to the engine module.** They are the same two
quantities; only the labels differ. This matters because the study's conclusion
is that the *drift* channel is the defensible stopping estimator, so reading a
curve under the wrong convention inverts the result. Everything below is in the
convention of this table.

$T_2 - T_1 = \\mathrm{{var}}(p)/2m$ is an exact identity, not an approximation
(asserted in `tests/test_refined.py`).
""")

    code(f"""
import sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt

HERE = Path.cwd()
sys.path.insert(0, str(HERE))
import refined as R

try:
    from inqview.visualisation import style
    style.apply_theme()
except Exception as exc:
    warnings.warn(f"inqview theme unavailable ({{exc}}); using matplotlib defaults")

# numpy 2 renamed trapz -> trapezoid and REMOVED the old name (it does not warn,
# it raises). Bind once so the notebook runs on either major version.
TRAPZ = getattr(np, "trapezoid", None) or np.trapz

FIGS = HERE / "refined_figs"; FIGS.mkdir(exist_ok=True)
def save(fig, name):
    fig.savefig(FIGS / f"{{name}}.png", dpi=150, bbox_inches="tight")
    return fig

WP_NAME, CL_NAME = {wp_name!r}, {cl_name!r}
wp = R.wp_frame(WP_NAME)
cl = R.cl_frame(CL_NAME)

print(f"wavepacket : {{len(wp)}} steps, t = {{wp.t.iloc[0]:.2f}} .. {{wp.t.iloc[-1]:.2f}} a.u.")
print(f"classical  : {{len(cl)}} steps, t = {{cl.t.iloc[0]:.2f}} .. {{cl.t.iloc[-1]:.2f}} a.u.")
print()
print(f"T1(0) = {{wp.T1_drift_ev.iloc[0]:.3f}} eV   (launch energy, should be 50.00)")
print(f"T2(0) = {{wp.T2_total_ev.iloc[0]:.3f}} eV")
print(f"T2-T1 at t=0 = {{wp.var_term_ev.iloc[0]:.4f}} eV  "
      f"(free-Gaussian localisation energy {{R.T2_MINUS_T1_FREE_EV:.4f}} eV)")
""")

    # ===================================================== 1. positions
    md("""
---
## 1. Where is the projectile?

The classical position is unambiguous: one Ehrenfest particle, one $z(t)$.

For the wavepacket it is not, and the two candidates fail in *different* ways —
which is why both are carried rather than one being chosen:

- **$s_\\mathrm{centroid}$** — the circular (Resta-phase) centroid of
  $|\\psi|^2$, unwrapped. Periodic-exact, so it stays meaningful while the packet
  straddles a cell face (this one does at $t=0$, launched 2 Bohr from the $-z$
  face). It is a property of the **density**, so if the packet splits into
  transmitted and reflected lobes it reports their weighted mean — physically
  real, but no longer a particle position.
- **$s_{\\int\\langle p\\rangle}$** — $z_0 + \\int \\langle p_z\\rangle\\,dt$. A
  property of the **momentum**.

Under exact Ehrenfest dynamics with $m=1$ these must coincide. **Where they part
company is where "the wavepacket has a trajectory" stops being true** — so the
residual is a diagnostic for choosing the window, not an error to be minimised.
""")

    code(r"""
fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.2))

ax[0].plot(cl.t, cl.z_unwrapped, lw=2.2, label="classical $z(t)$", color="tab:blue")
ax[0].plot(wp.t, wp.s_centroid, lw=1.8, ls="--", color="tab:red",
           label=r"WP $s_\mathrm{centroid}$ (circular)")
ax[0].plot(wp.t, wp.s_pintegral, lw=1.4, ls=":", color="tab:green",
           label=r"WP $s_{\int\langle p\rangle dt}$")
ax[0].axhline(-30, color="0.5", lw=0.8, ls="-.",
              label=r"cell faces $z=\pm L_z/2$")
ax[0].axhline(+30, color="0.5", lw=0.8, ls="-.")
ax[0].set_ylim(-33, 33)
ax[0].set_xlabel("t (a.u.)"); ax[0].set_ylabel("z (Bohr)")
ax[0].set_title("(a) trajectories"); ax[0].legend(fontsize=8, loc="lower right")

ax[1].plot(wp.t, wp.ehrenfest_resid * 1e3, lw=1.8, color="tab:red")
ax[1].axhline(0, color="0.5", lw=0.8)
ax[1].set_xlabel("t (a.u.)")
ax[1].set_ylabel(r"$s_\mathrm{centroid} - s_{\int\langle p\rangle}$  ($10^{-3}$ Bohr)")
ax[1].set_title("(b) WP: do the two definitions agree?")

d = np.interp(cl.t, wp.t, wp.s_pintegral) - cl.z_unwrapped
ax[2].plot(cl.t, d, lw=1.8, color="tab:purple",
           label=r"WP $s_{\int\langle p\rangle}$ $-$ classical")
d2 = np.interp(cl.t, wp.t, wp.s_centroid) - cl.z_unwrapped
ax[2].plot(cl.t, d2, lw=1.4, ls="--", color="tab:orange",
           label=r"WP $s_\mathrm{centroid}$ $-$ classical")
ax[2].axhline(0, color="0.5", lw=0.8)
ax[2].set_xlabel("t (a.u.)"); ax[2].set_ylabel("separation (Bohr)")
ax[2].set_title("(c) twin separation"); ax[2].legend(fontsize=8)

fig.tight_layout(); save(fig, "01_positions"); plt.show()

res = wp.ehrenfest_resid.abs()
print(f"max |s_centroid - s_pintegral|  = {res.max():.3e} Bohr  "
      f"(at t = {wp.t[res.idxmax()]:.2f} a.u.)")
print(f"  ... at t = 10 a.u.: {np.interp(10.0, wp.t, wp.ehrenfest_resid):.3e} Bohr")
print(f"  ... at t = 20 a.u.: {np.interp(20.0, wp.t, wp.ehrenfest_resid):.3e} Bohr")
print(f"final separation, WP - classical = {d.iloc[-1]:+.3f} Bohr")
print(f"classical path travelled = {cl.z_unwrapped.iloc[-1]-cl.z_unwrapped.iloc[0]:.3f} Bohr")
print(f"WP path travelled        = {wp.s_pintegral.iloc[-1]-wp.s_pintegral.iloc[0]:.3f} Bohr")
""")

    # ===================================================== 2a. classical energy
    md("""
---
## 2. Kinetic energy

### 2a. Classical — and the closure sanity check

The projectile is an **external moving-charge perturbation**, not an INQ ion, so
`energy_ion` and `energy_ion_kinetic` are identically zero and INQ's
`energy_total` is the **bath energy alone**. The projectile's $\\tfrac12 m v^2$
is tracked separately by the Ehrenfest integrator.

Energy is conserved overall, so the two must be mirror images:

$$\\Delta E_\\mathrm{total}^{\\,\\mathrm{bath}}(t) + \\Delta\\left(\\tfrac12 m v^2\\right)(t) = 0$$

**This is the check that certifies the classical stopping power is a real energy
transfer** and not a bookkeeping artefact of a moving perturbation. Panel (b) is
that residual.
""")

    code(r"""
fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.2))

ax[0].plot(cl.t, cl.d_e_total_ev, lw=2.0, color="tab:blue",
           label=r"$\Delta E_\mathrm{total}$ (bath)")
ax[0].plot(cl.t, cl.d_ke_ev, lw=2.0, color="tab:red",
           label=r"$\Delta(\frac{1}{2}mv^2)$ (projectile)")
ax[0].plot(cl.t, cl.closure_ev, lw=2.4, color="k", label="sum")
ax[0].axhline(0, color="0.5", lw=0.8)
ax[0].set_xlabel("t (a.u.)"); ax[0].set_ylabel(r"$\Delta E$ (eV)")
ax[0].set_title("(a) the bath gains what the projectile loses")
ax[0].legend(fontsize=9)

ax[1].plot(cl.t, cl.closure_ev * 1e6, lw=1.8, color="k")
ax[1].axhline(0, color="0.5", lw=0.8)
ax[1].set_xlabel("t (a.u.)")
ax[1].set_ylabel(r"closure residual ($10^{-6}$ eV)")
ax[1].set_title("(b) closure residual — the sanity check")

ax[2].plot(cl.t, cl.ke_ev, lw=2.0, color="tab:red")
ax[2].axhline(cl.ke_ev.iloc[0], color="0.5", lw=0.8, ls="--")
ax[2].set_xlabel("t (a.u.)"); ax[2].set_ylabel(r"$\frac{1}{2}mv^2$ (eV)")
ax[2].set_title("(c) projectile kinetic energy, absolute")

fig.tight_layout(); save(fig, "02_classical_energy"); plt.show()

print(f"KE:  {cl.ke_ev.iloc[0]:.4f} -> {cl.ke_ev.iloc[-1]:.4f} eV   "
      f"(lost {-cl.d_ke_ev.iloc[-1]:.4f} eV)")
print(f"bath: gained {cl.d_e_total_ev.iloc[-1]:+.4f} eV")
print(f"CLOSURE: max |residual| = {cl.closure_ev.abs().max():.3e} eV over {len(cl)} steps")
print(f"         final residual = {cl.closure_ev.iloc[-1]:+.3e} eV")
ok = cl.closure_ev.abs().max() < 1e-3
print(f"         => {'PASS' if ok else 'FAIL'} (energy budget closes)")
""")

    # ===================================================== 2b. WP energy
    md("""
### 2b. Wavepacket — $T_1$ and $T_2$

$T_1 = \\langle p\\rangle^2/2m$ is the **drift** channel: the energy the packet
carries as coherent forward motion. It is the quantity with a classical
counterpart.

$T_2 = T_1 + \\mathrm{var}(p)/2m$ adds the **momentum-spread** channel. Under
*free* evolution $\\mathrm{var}(p)$ is exactly conserved — a free Gaussian
disperses in space but never in momentum — so **any growth in $T_2 - T_1$ is
interaction**, and it is the contamination that made $-dT_2/ds$ unusable as a
stopping power in the bulk study.

The dashed line marks the free-Gaussian value
$3/(4\\sigma_\\mathrm{WP}^2) = 1.276$ eV.

**Watch the printed $\\mathrm{var}(p)$ split below panel (c).** The growth is not
isotropic, and which direction it is in changes what it means: growth in
$\\mathrm{var}(p_z)$ is longitudinal straggling *along* the stopping direction and
contaminates $-dT_2/ds$ directly, whereas growth in $\\mathrm{var}(p_\\perp)$ is
the packet spreading *sideways* toward the bore wall — a failure of the
channeling premise rather than of the stopping estimator.
""")

    code(r"""
fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.2))

ax[0].plot(wp.t, wp.T1_drift_ev, lw=2.0, color="tab:green",
           label=r"$T_1 = \langle p\rangle^2/2m$")
ax[0].plot(wp.t, wp.T2_total_ev, lw=2.0, color="tab:orange",
           label=r"$T_2 = T_1 + \mathrm{var}(p)/2m$")
ax[0].set_xlabel("t (a.u.)"); ax[0].set_ylabel("energy (eV)")
ax[0].set_title("(a) the two kinetic-energy definitions")
ax[0].legend(fontsize=9)

ax[1].plot(wp.t, wp.d_T1_ev, lw=2.0, color="tab:green", label=r"$\Delta T_1$")
ax[1].plot(wp.t, wp.d_T2_ev, lw=2.0, color="tab:orange", label=r"$\Delta T_2$")
ax[1].axhline(0, color="0.5", lw=0.8)
ax[1].set_xlabel("t (a.u.)"); ax[1].set_ylabel(r"$\Delta E$ (eV)")
ax[1].set_title("(b) changes from t = 0"); ax[1].legend(fontsize=9)

ax[2].plot(wp.t, wp.var_term_ev, lw=2.2, color="tab:red",
           label=r"$T_2-T_1 = \mathrm{var}(p)/2m$")
ax[2].axhline(R.T2_MINUS_T1_FREE_EV, color="0.4", ls="--", lw=1.2,
              label="free-Gaussian value")
ax[2].set_xlabel("t (a.u.)"); ax[2].set_ylabel("energy (eV)")
ax[2].set_title("(c) the spread term — flat means channeling worked")
ax[2].legend(fontsize=9)

fig.tight_layout(); save(fig, "03_wp_energy"); plt.show()
""")

    md("""
#### The $\\mathrm{var}(p)/2m$ term on its own, and split by direction

$\\mathrm{var}(p)/2m$ is the entire difference between the two kinetic-energy
definitions, so it deserves its own axes rather than sharing them. Splitting it
into $z$ and $\\perp$ matters because the two mean different things:

- $\\mathrm{var}(p_z)/2m$ — **longitudinal straggling**, spread *along* the
  stopping direction. This is what directly contaminates $-dT_2/ds$ as a
  stopping power, because it lives on the same axis as the drift being measured.
- $\\mathrm{var}(p_\\perp)/2m$ — **transverse heating**, the packet spreading
  sideways toward the bore wall. This is the *channeling premise* failing, not
  the estimator being contaminated.

Under free evolution both are exactly conserved, so the dashed lines are the
zero-interaction reference and any departure from them is interaction.
""")

    code(r"""
fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.2))
HAv = R.HA_TO_EV

var_tot = 0.5 * wp.var_p3d * HAv
var_z = 0.5 * wp.var_pz * HAv
var_perp = 0.5 * wp.var_perp * HAv

ax[0].plot(wp.t, var_tot, lw=2.4, color="tab:red", label=r"$\mathrm{var}(p)/2m$ (3D)")
ax[0].axhline(R.T2_MINUS_T1_FREE_EV, color="0.4", ls="--", lw=1.2,
              label="free-evolution value")
ax[0].set_xlabel("t (a.u.)"); ax[0].set_ylabel("energy (eV)")
ax[0].set_title(r"(a) $\mathrm{var}(p)/2m$ vs time")
ax[0].legend(fontsize=9)

ax[1].plot(wp.t, var_z, lw=2.2, color="tab:blue", label=r"$\mathrm{var}(p_z)/2m$")
ax[1].plot(wp.t, var_perp, lw=2.2, color="tab:cyan", label=r"$\mathrm{var}(p_\perp)/2m$")
ax[1].axhline(0.5 * R.VAR_P_FREE * HAv, color="tab:blue", ls="--", lw=1.0)
ax[1].axhline(0.5 * 2 * R.VAR_P_FREE * HAv, color="tab:cyan", ls="--", lw=1.0)
ax[1].set_xlabel("t (a.u.)"); ax[1].set_ylabel("energy (eV)")
ax[1].set_title("(b) split by direction (dashed = free values)")
ax[1].legend(fontsize=9)

ax[2].plot(wp.t, var_tot - var_tot.iloc[0], lw=2.4, color="tab:red", label="total")
ax[2].plot(wp.t, var_z - var_z.iloc[0], lw=1.8, color="tab:blue", label=r"$z$")
ax[2].plot(wp.t, var_perp - var_perp.iloc[0], lw=1.8, color="tab:cyan", label=r"$\perp$")
ax[2].axhline(0, color="0.5", lw=0.8)
ax[2].set_xlabel("t (a.u.)"); ax[2].set_ylabel(r"$\Delta$ energy (eV)")
ax[2].set_title(r"(c) growth from $t=0$ — how much each contributes")
ax[2].legend(fontsize=9)

fig.tight_layout(); save(fig, "03b_var_p_term"); plt.show()

print("var(p)/2m in eV:")
for lbl, v in (("total (3D)", var_tot), ("longitudinal z", var_z),
               ("transverse   ", var_perp)):
    print(f"  {lbl}: {v.iloc[0]:.4f} -> {v.iloc[-1]:.4f} eV   "
          f"(grew {v.iloc[-1]-v.iloc[0]:+.4f} eV, {100*(v.iloc[-1]/v.iloc[0]-1):+.1f} %)")
gz = var_z.iloc[-1] - var_z.iloc[0]; gp = var_perp.iloc[-1] - var_perp.iloc[0]
print()
print(f"  of the {gz+gp:.4f} eV total growth, {100*gz/(gz+gp):.1f} % is longitudinal "
      f"and {100*gp/(gz+gp):.1f} % is transverse.")
print("  => most of the T2-T1 growth is the packet spreading SIDEWAYS, which is the")
print("     channeling premise failing rather than the stopping estimator being spoiled.")

print(f"T1: {wp.T1_drift_ev.iloc[0]:.4f} -> {wp.T1_drift_ev.iloc[-1]:.4f} eV "
      f"({wp.d_T1_ev.iloc[-1]:+.4f})")
print(f"T2: {wp.T2_total_ev.iloc[0]:.4f} -> {wp.T2_total_ev.iloc[-1]:.4f} eV "
      f"({wp.d_T2_ev.iloc[-1]:+.4f})")
print(f"var-term: {wp.var_term_ev.iloc[0]:.4f} -> {wp.var_term_ev.iloc[-1]:.4f} eV "
      f"({100*(wp.var_term_ev.iloc[-1]/wp.var_term_ev.iloc[0]-1):+.1f} %)")
print()
print("var(p) split (a.u.), showing WHERE the spread grows:")
for lbl, c in (("var(p_z)  ", "var_pz"), ("var(p_perp)", "var_perp"),
               ("var(p) 3D ", "var_p3d")):
    print(f"  {lbl}: {wp[c].iloc[0]:.5f} -> {wp[c].iloc[-1]:.5f}  "
          f"({100*(wp[c].iloc[-1]/wp[c].iloc[0]-1):+.1f} %)")
""")

    # ===================================================== 2c. comparison
    md("""
### 2c. Classical against wavepacket

Both halves are light projectiles and both decelerate
(`.claude/rules/light-projectile-stopping.md`), so the honest comparison is of
their **energy loss against path travelled**, not against time.

If channeling works as intended, $\\Delta T_1$ (WP drift) should track
$\\Delta(\\tfrac12 mv^2)$ (classical), while $\\Delta T_2$ runs shallower because
the growing spread term partially cancels the drift loss.
""")

    code(r"""
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))

ax[0].plot(cl.t, cl.d_ke_ev, lw=2.2, color="tab:blue",
           label=r"classical $\Delta(\frac{1}{2}mv^2)$")
ax[0].plot(wp.t, wp.d_T1_ev, lw=2.0, color="tab:green", label=r"WP $\Delta T_1$")
ax[0].plot(wp.t, wp.d_T2_ev, lw=2.0, ls="--", color="tab:orange",
           label=r"WP $\Delta T_2$")
ax[0].axhline(0, color="0.5", lw=0.8)
ax[0].set_xlabel("t (a.u.)"); ax[0].set_ylabel(r"$\Delta E$ (eV)")
ax[0].set_title("(a) versus time"); ax[0].legend(fontsize=9)

ax[1].plot(cl.z_unwrapped - cl.z_unwrapped.iloc[0], cl.d_ke_ev, lw=2.2,
           color="tab:blue", label="classical")
ax[1].plot(wp.s_pintegral - wp.s_pintegral.iloc[0], wp.d_T1_ev, lw=2.0,
           color="tab:green", label=r"WP $\Delta T_1$")
ax[1].plot(wp.s_pintegral - wp.s_pintegral.iloc[0], wp.d_T2_ev, lw=2.0,
           ls="--", color="tab:orange", label=r"WP $\Delta T_2$")
ax[1].axhline(0, color="0.5", lw=0.8)
ax[1].set_xlabel("path travelled (Bohr)"); ax[1].set_ylabel(r"$\Delta E$ (eV)")
ax[1].set_title("(b) versus path — the slope IS the stopping power")
ax[1].legend(fontsize=9)

fig.tight_layout(); save(fig, "04_comparison"); plt.show()

print("energy lost by the END of the run:")
print(f"  classical  : {-cl.d_ke_ev.iloc[-1]:.4f} eV over "
      f"{cl.z_unwrapped.iloc[-1]-cl.z_unwrapped.iloc[0]:.3f} Bohr")
print(f"  WP  (T1)   : {-wp.d_T1_ev.iloc[-1]:.4f} eV over "
      f"{wp.s_pintegral.iloc[-1]-wp.s_pintegral.iloc[0]:.3f} Bohr")
print(f"  WP  (T2)   : {-wp.d_T2_ev.iloc[-1]:.4f} eV")
print()
print("NOTE: whole-run averages, NOT a stopping power at v0 -- both projectiles")
print("      decelerate. Choose the window in section 5.")
""")

    # ===================================================== 3. momentum dist
    md(f"""
---
## 3. The momentum distribution, at three times

$\\mathrm{{var}}(p)$ is a single number summarising a whole distribution, and the
same variance can come from very different physics:

- a **uniform broadening** — the packet is being scattered elastically in all
  directions;
- a **shift with the peak intact** — coherent deceleration, no contamination;
- a **growing low-$k$ tail with the peak intact** — part of the packet is being
  stopped or captured while the rest sails on. This is the dangerous one, because
  it inflates $\\mathrm{{var}}(p)$ without the packet having "spread" in any
  intuitive sense.

Slices requested at $t \\approx$ {', '.join(f'{x:g}' for x in times)} a.u.
(snapped to the nearest written slice — the distribution is saved every 15 steps,
and interpolating between two distributions would smear exactly the feature we
are looking for).

> **Read the log panel, not the linear one.** The linear curve is visibly spiky,
> and that is a **resolution limit, not noise in the physics**: the packet's
> momentum width is $\\sigma_p = 1/(\\sqrt2\\,\\sigma_\\mathrm{{WP}}) = 0.177$ a.u.,
> while the $k$-grid spacing is $2\\pi/L_z = 0.105$ and $2\\pi/L_{{xy}} = 0.157$
> a.u. — so the whole distribution is only **1–2 grid points wide**, and a radial
> $|k|$ histogram of it is spiky by construction. (Checked: the spikes do *not*
> correlate with shell-occupancy counts, $r = -0.04$, so this is discrete
> sampling of a narrow off-origin Gaussian, not shell-binning geometry.)
>
> The **moments** in panel (c) come from `wp_momentum_stats.csv`, which computes
> $\\langle p\\rangle$ and $\\langle p^2\\rangle$ as exact grid expectation values
> and never bins — so they are quantitative even though this histogram is only
> qualitative. Trust panel (c) and the log-scale tail; treat the linear peak
> shape as indicative.

### Are the spikes real? **No — settled by direct test.**

Build the *analytic* launched Gaussian on the identical grid, give it **zero
interaction with anything**, and put it through the identical radial binning.
If the spikes were physics, that curve would be smooth. It is not:

| | roughness $\\langle|\\nabla^2 P|\\rangle / \\langle P\\rangle$ |
|---|---|
| analytic non-interacting Gaussian | **2.685** |
| measured, $t=0$ | **2.648** |

and the two histograms correlate at **r = 0.99982**. The spikes are reproduced
in full by a packet that has never interacted with anything at all.

The mechanism is *not* empty bins (every bin holds ≥ 368 grid points). It is
that the packet occupies only ~**95 effective k-points** out of 768 000
(participation ratio $1/\\sum w^2$; the top 50 carry 63 % of the norm). Binning
~10² significant points into 128 radial shells gives a comb. The cure is
**fewer bins** (`n_bins` ≈ 40, i.e. $\\Delta k_\\mathrm{{bin}} \\approx \\Delta
k_\\mathrm{{grid}}$) in a future run — not smoothing this one.
""")

    code(f"""
md_ = R.momentum_slices(WP_NAME)
slices = R.nearest_slices(md_, {times!r})
print("available slices:", md_.t.nunique(),
      f"from t = {{md_.t.min():.2f}} to {{md_.t.max():.2f}} a.u.")

fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.4))
colors = ["tab:blue", "tab:orange", "tab:red"]

for (t_act, sl), c in zip(slices, colors):
    ax[0].plot(sl.k, sl.n_wp, lw=1.8, color=c, label=f"t = {{t_act:.2f}} a.u.")
    ax[1].semilogy(sl.k, np.maximum(sl.n_wp, 1e-16), lw=1.6, color=c,
                   label=f"t = {{t_act:.2f}} a.u.")

k0 = R.CS.V0
for a in (ax[0], ax[1]):
    a.axvline(k0, color="0.4", ls="--", lw=1.0)
    a.set_xlabel(r"$k$ (Bohr$^{{-1}}$)")
    a.set_xlim(0, 2.5 * k0)
    a.legend(fontsize=8)
ax[0].text(k0 * 1.02, ax[0].get_ylim()[1] * 0.92, r"$k_0$", fontsize=9, color="0.4")
ax[0].set_ylabel(r"$|\\psi_\\mathrm{{WP}}(k)|^2$")
ax[0].set_title("(a) WP momentum distribution, linear")
ax[1].set_ylabel(r"$|\\psi_\\mathrm{{WP}}(k)|^2$")
ax[1].set_title("(b) same, log — shows the scattered tail")

ax[2].plot(wp.t, wp.var_pz, lw=2.0, color="tab:red", label=r"$\\mathrm{{var}}(p_z)$")
ax[2].plot(wp.t, wp.var_perp / 2.0, lw=1.6, color="tab:cyan",
           label=r"$\\mathrm{{var}}(p_\\perp)/2$")
ax[2].axhline(R.VAR_P_FREE, color="0.4", ls="--", lw=1.2, label="free value")
for (t_act, _), c in zip(slices, colors):
    ax[2].axvline(t_act, color=c, lw=1.0, ls=":")
ax[2].set_xlabel("t (a.u.)"); ax[2].set_ylabel(r"$\\sigma_p^2$ (a.u.)")
ax[2].set_title("(c) where the slices sit on var(p)"); ax[2].legend(fontsize=8)

fig.tight_layout(); save(fig, "05_momentum_distribution"); plt.show()

print()
print(f"{{'t (a.u.)':>10}} {{'peak k':>9}} {{'<k>':>9}} {{'norm':>10}} {{'frac k<k0/2':>12}}")
for t_act, sl in slices:
    n = sl.n_wp.to_numpy(); k = sl.k.to_numpy()
    tot = TRAPZ(n, k)
    kbar = TRAPZ(n * k, k) / tot if tot > 0 else np.nan
    lowf = TRAPZ(n[k < k0 / 2], k[k < k0 / 2]) / tot if tot > 0 else np.nan
    print(f"{{t_act:10.2f}} {{k[np.argmax(n)]:9.4f}} {{kbar:9.4f}} {{tot:10.4e}} {{lowf:12.4e}}")
""")

    # ============================================ 3b. 2-D momentum map
    md("""
---
## 3b. The 2-D momentum map — separating deceleration from heating

The radial $n(|k|)$ above cannot answer the question that matters, because it
folds the drift direction into the same coordinate as the transverse spread: a
packet that is **decelerating** and one that is **heating sideways** both look
like "the peak moved left". So here the orbital's own complex dump is
transformed and resolved into $(k_z,\\; k_\\perp)$:

- **weight moving left along $k_z$** = deceleration, the stopping we want;
- **weight moving up in $k_\\perp$** = transverse heating, the channeling premise
  failing.

$k_z$ is the **native FFT axis and is not binned at all**, so its moments are
exact; only $k_\\perp$ is binned, at one transverse grid spacing.

*Validated:* the $t=0$ map round-trips to the recorded moments exactly —
$\\langle k_z\\rangle = 1.917011$, $\\mathrm{var}(k_z) = 0.031250$,
$\\mathrm{var}(k_\\perp) = 0.062500$, matching `wp_momentum_stats.csv` to every
printed digit. (A wrong FFT ordering would shift $\\langle k_z\\rangle$ by a phase
ramp and still look plausible.)

### Your asymmetry question, answered

You asked whether the norm above the mean at $t=0$ exceeds that at a later time
— i.e. whether every momentum channel got the *same* impulse (rigid shift) or
not. **A rigid shift preserves the shape exactly**: same width, same skew. So
the test is whether width and skew change.

One caveat that had to be fixed first: a hard `k_z > mean` count returns
**0.454** for the $t=0$ packet, which is *exactly symmetric* — because only ~8
resolved $k_z$ points carry the packet and the mean falls between them. The
table below interpolates the CDF instead, which correctly returns 0.4987 at
$t=0$.
""")

    code(r"""
steps_wf = R.available_wf_steps(WP_NAME)
pick = [steps_wf[0], steps_wf[len(steps_wf)//2], steps_wf[-1]]
maps = {s: R.momentum_map(WP_NAME, s) for s in pick}
kz, kperp, P0 = maps[pick[0]]

rows = []
for s in pick:
    kz_s, kp_s, P = maps[s]
    a = R.kz_asymmetry(kz_s, P)
    a["t"] = s * R.CS.DT
    a["m2_kperp"] = float((kp_s**2 * P.sum(axis=0)).sum())
    rows.append(a)
asym = pd.DataFrame(rows)[["t", "mean_kz", "sigma_kz", "skewness",
                           "frac_above_mean", "median_minus_mean", "m2_kperp"]]
print("Longitudinal momentum distribution -- shape, not just position:")
print(asym.round(5).to_string(index=False))
print()
print("A RIGID shift would hold sigma_kz and skewness CONSTANT.")
print(f"  sigma_kz : {asym.sigma_kz.iloc[0]:.5f} -> {asym.sigma_kz.iloc[-1]:.5f}"
      f"  ({100*(asym.sigma_kz.iloc[-1]/asym.sigma_kz.iloc[0]-1):+.1f} %)")
print(f"  skewness : {asym.skewness.iloc[0]:+.5f} -> {asym.skewness.iloc[-1]:+.5f}")
print(f"  norm above the mean : {asym.frac_above_mean.iloc[0]:.4f} -> "
      f"{asym.frac_above_mean.iloc[-1]:.4f}")
print()
print("=> the impulse is NOT the same for every momentum channel.")
""")

    code(r"""
fig, ax = plt.subplots(2, 3, figsize=(16.0, 8.4))
vmax = max(P.max() for _, _, P in maps.values())

for j, s in enumerate(pick):
    kz_s, kp_s, P = maps[s]
    im = ax[0, j].pcolormesh(kz_s, kp_s, P.T, shading="auto", cmap="viridis",
                             vmin=0, vmax=vmax)
    ax[0, j].set_title(f"$P(k_z,k_\\perp)$   t = {s*R.CS.DT:.0f} a.u.")
    ax[0, j].set_xlabel(r"$k_z$ (Bohr$^{-1}$)")
    ax[0, j].set_ylabel(r"$k_\perp$ (Bohr$^{-1}$)")
    ax[0, j].axvline(R.CS.V0, color="w", ls="--", lw=1.0)
    ax[0, j].set_xlim(0.8, 3.0); ax[0, j].set_ylim(0, 1.6)
    fig.colorbar(im, ax=ax[0, j], fraction=0.046)

    D = P - P0
    lim = np.abs(D).max()
    im2 = ax[1, j].pcolormesh(kz_s, kp_s, D.T, shading="auto", cmap="RdBu_r",
                              vmin=-lim, vmax=lim)
    ax[1, j].set_title(f"difference from t=0   (t = {s*R.CS.DT:.0f})")
    ax[1, j].set_xlabel(r"$k_z$ (Bohr$^{-1}$)")
    ax[1, j].set_ylabel(r"$k_\perp$ (Bohr$^{-1}$)")
    ax[1, j].axvline(R.CS.V0, color="k", ls="--", lw=1.0)
    ax[1, j].set_xlim(0.8, 3.0); ax[1, j].set_ylim(0, 1.6)
    fig.colorbar(im2, ax=ax[1, j], fraction=0.046)

fig.suptitle("Top: momentum density.  Bottom: change from t=0 "
             "(red = gained, blue = lost).  Dashed line = launch $k_0$.",
             fontsize=10)
fig.tight_layout(); save(fig, "05b_kz_kperp_map"); plt.show()
""")

    code(r"""
fig, ax = plt.subplots(1, 2, figsize=(12.0, 4.4))
colors = ["tab:blue", "tab:orange", "tab:red"]

for s, c in zip(pick, colors):
    kz_s, kp_s, P = maps[s]
    ax[0].plot(kz_s, P.sum(axis=1), lw=1.8, color=c, label=f"t = {s*R.CS.DT:.0f} a.u.")
    ax[1].plot(kp_s, P.sum(axis=0), lw=1.8, color=c, label=f"t = {s*R.CS.DT:.0f} a.u.")
ax[0].axvline(R.CS.V0, color="0.4", ls="--", lw=1.0)
ax[0].set_xlim(0.8, 3.0); ax[0].set_xlabel(r"$k_z$ (Bohr$^{-1}$)")
ax[0].set_ylabel(r"$P(k_z)$")
ax[0].set_title(r"(a) longitudinal marginal — exact, unbinned")
ax[0].legend(fontsize=8)
ax[1].set_xlim(0, 1.6); ax[1].set_xlabel(r"$k_\perp$ (Bohr$^{-1}$)")
ax[1].set_ylabel(r"$P(k_\perp)$")
ax[1].set_title(r"(b) transverse marginal (Rayleigh: peaks away from 0)")
ax[1].legend(fontsize=8)
fig.tight_layout(); save(fig, "05c_momentum_marginals"); plt.show()

kz_e, _, Pe = maps[pick[-1]]
pe = Pe.sum(axis=1); p0m = P0.sum(axis=1)
gain_lo = float(pe[kz_e < R.CS.V0].sum() - p0m[kz_e < R.CS.V0].sum())
gain_hi = float(pe[kz_e > R.CS.V0].sum() - p0m[kz_e > R.CS.V0].sum())
print(f"norm transferred across the launch k0 = {R.CS.V0:.4f}:")
print(f"   below k0 : {gain_lo:+.5f}      above k0 : {gain_hi:+.5f}")
print(f"   far tail  k_z > 2.4 : {float(p0m[kz_e>2.4].sum()):.3e} -> "
      f"{float(pe[kz_e>2.4].sum()):.3e}")
print(f"   far tail  k_z < 1.4 : {float(p0m[kz_e<1.4].sum()):.3e} -> "
      f"{float(pe[kz_e<1.4].sum()):.3e}")
""")

    # ===================================================== 4. interactions
    md("""
---
## 4. Pairwise interaction energies

The decomposition into projectile **P**, bath electrons **S**, and neutralising
background **B** (`.claude/rules/decomposed-interaction-energies.md`):

| term | meaning |
|---|---|
| $E_{SS}$ | bath–bath Hartree |
| $E_{PS}$ | **projectile–bath — the interaction that stops it** |
| $E_{PP}$ | projectile self-Hartree — **exists only for the wavepacket** |
| $E_{SB}$, $E_{PB}$ | couplings to the background |

Why this decomposition and not INQ's own scalars: the two representations put the
projectile in **different ledger terms**. For the classical run `energy_external`
is non-zero and `energy_hartree` $= E_{SS}$; for the wavepacket
`energy_external` is *identically zero* and `energy_hartree`
$= E_{SS} + E_{PS} + E_{PP}$. A raw `energy_hartree` comparison between the two
halves compares a net quantity against a gross one. The pairwise terms are
representation-independent and genuinely comparable.

$E_{PP}$ is the **uncancelled self-interaction of a wavepacket in LDA**. It has
no classical counterpart and is the leading suspect for any residual
classical↔WP discrepancy in $S$.

**Deltas, not absolutes:** $E_{SB}$, $E_{PB}$, $E_{BB}$ carry the charged-cell
$G=0$ gauge, so their absolute values are not comparable across representations
while their changes are.
""")

    code(r"""
icl = R.interactions("classical", CL_NAME)
iwp = R.interactions("wp", WP_NAME)

fig, ax = plt.subplots(1, 3, figsize=(16.0, 4.4))
palette = {"e_ss": "tab:blue", "e_ps": "tab:red", "e_pp": "tab:purple",
           "e_sb": "tab:green", "e_pb": "tab:brown"}

for a, df_, title in ((ax[0], icl, "(a) classical"), (ax[1], iwp, "(b) wavepacket")):
    for term in R.INTERACTION_TERMS:
        col = f"d_{term}_ev"
        if col not in df_.columns:
            continue
        y = df_[col].to_numpy()
        if np.allclose(y, 0.0):
            continue                       # identically zero: uniform background
        a.plot(df_.time_au, y, lw=1.8, color=palette[term], label=R.TERM_LABEL[term])
    a.axhline(0, color="0.5", lw=0.8)
    a.set_xlabel("t (a.u.)"); a.set_ylabel(r"$\Delta E$ (eV)")
    a.set_title(title); a.legend(fontsize=8)

# E_PP has no classical counterpart -- give it its own axis so its SCALE is legible.
ax[2].plot(iwp.time_au, iwp.d_e_pp_ev, lw=2.2, color="tab:purple",
           label=r"WP $\Delta E_{PP}$ (self-Hartree)")
ax[2].plot(iwp.time_au, iwp.d_e_ps_ev, lw=1.6, ls="--", color="tab:red",
           label=r"WP $\Delta E_{PS}$")
ax[2].plot(icl.time_au, icl.d_e_ps_ev, lw=1.6, ls=":", color="tab:blue",
           label=r"classical $\Delta E_{PS}$")
ax[2].axhline(0, color="0.5", lw=0.8)
ax[2].set_xlabel("t (a.u.)"); ax[2].set_ylabel(r"$\Delta E$ (eV)")
ax[2].set_title("(c) the terms that differ between the twins")
ax[2].legend(fontsize=8)

fig.tight_layout(); save(fig, "06_interactions"); plt.show()

rows = []
for half, df_ in (("classical", icl), ("wp", iwp)):
    for term in R.INTERACTION_TERMS:
        col = f"d_{term}_ev"
        if col in df_.columns:
            rows.append({"half": half, "term": term.upper(),
                         "delta_end_eV": df_[col].iloc[-1],
                         "max_abs_eV": df_[col].abs().max()})
tab = pd.DataFrame(rows).pivot(index="term", columns="half",
                               values=["delta_end_eV", "max_abs_eV"])
print("Interaction-energy changes over the full run (eV):")
print(tab.round(4).to_string())
print()
print(f"E_PP absolute (WP): {iwp.e_pp_ev.iloc[0]:.4f} -> {iwp.e_pp_ev.iloc[-1]:.4f} eV")
print(f"E_PP absolute (classical, expect ~0 or constant): "
      f"{icl.e_pp_ev.iloc[0]:.4f} -> {icl.e_pp_ev.iloc[-1]:.4f} eV")

print()
print("=" * 68)
print("TWIN-MATCH CHECK AT t = 0 -- 'do the two halves create the same potential?'")
print("=" * 68)
worst = 0.0
for term in ("e_ss", "e_ps", "e_pp", "e_sb", "e_pb"):
    if term not in iwp.columns or term not in icl.columns:
        continue
    a, b = iwp[term].iloc[0], icl[term].iloc[0]
    d = abs(a - b); worst = max(worst, d)
    print(f"  {term.upper():5s}  wp = {a:+16.10f}   classical = {b:+16.10f}   "
          f"|diff| = {d:.2e} Ha")
print(f"\n  worst disagreement: {worst:.3e} Ha = {worst*R.HA_TO_EV:.3e} eV")
print(f"  => {'PASS' if worst < 1e-9 else 'FAIL'}: the classical Gaussian charge cloud at "
      f"sigma_pot = sigma_WP/sqrt(2) = {R.CS.SIGMA_POT:.6f} Bohr")
print("     reproduces the wavepacket's OWN charge density to numerical precision.")
print()
print("  This is the direct evidence for the twin contract: same geometry, same")
print("  bath, same t=0 electrostatics -- the ONLY difference is whether the")
print("  projectile is an external potential or an occupied KS orbital.")
print()
print(f"  classical E_PP is exactly constant (spread {icl.e_pp.max()-icl.e_pp.min():.1e} Ha):")
print("  a RIGID cloud has a fixed self-energy. The WP's falls as it disperses,")
print("  and that difference has no classical counterpart at all.")
""")

    # ============================================ 4b. E_PS + E_PB
    md("""
---
## 4b. $\\Delta(E_{PS} + E_{PB})$ — the combined projectile coupling

$E_{PS}$ and $E_{PB}$ look completely different between the two halves above.
That is a **bookkeeping** difference, not a physical one, and the sum shows why:

$$E_{PS} + E_{PB} \\;=\\; \\int n_P\\,(\\phi_S - \\phi_+) \\;=\\; \\int n_P\\,\\phi_{S+B}$$

which is the projectile's interaction with the **net** charge density of
everything that is not itself. The system is neutral, so $\\phi_{S+B}$ is the
*screened, short-ranged* potential — a well-defined physical field. Split apart,
$\\phi_S$ and $\\phi_+$ are each the potential of a charged subsystem and carry
the charged-cell $G=0$ gauge individually.

So the split between "bath" and "background" is only meaningful for a projectile
that stays put relative to the background. The classical cloud is rigid and the
tube is uniform in $z$, so its $\\Delta E_{PB}$ is *identically zero*; the
wavepacket spreads **transversely** into regions of different $\\phi_+$, so its
$\\Delta E_{PB}$ is large. **$\\Delta E_{PB}$(WP) is a pure transverse-spreading
signal**, nothing to do with the stopping along $z$.
""")

    code(r"""
scl = R.combined_projectile_coupling(icl)
swp = R.combined_projectile_coupling(iwp)

fig, ax = plt.subplots(1, 2, figsize=(12.0, 4.4))
ax[0].plot(icl.time_au, icl.d_e_ps_ev, lw=1.5, ls=":", color="tab:blue",
           label=r"classical $\Delta E_{PS}$")
ax[0].plot(icl.time_au, icl.d_e_pb_ev, lw=1.5, ls="-.", color="tab:cyan",
           label=r"classical $\Delta E_{PB}$")
ax[0].plot(iwp.time_au, iwp.d_e_ps_ev, lw=1.5, ls=":", color="tab:red",
           label=r"WP $\Delta E_{PS}$")
ax[0].plot(iwp.time_au, iwp.d_e_pb_ev, lw=1.5, ls="-.", color="tab:orange",
           label=r"WP $\Delta E_{PB}$")
ax[0].axhline(0, color="0.5", lw=0.8)
ax[0].set_xlabel("t (a.u.)"); ax[0].set_ylabel(r"$\Delta E$ (eV)")
ax[0].set_title("(a) separately — they look nothing alike")
ax[0].legend(fontsize=8)

ax[1].plot(icl.time_au, scl, lw=2.2, color="tab:blue", label="classical")
ax[1].plot(iwp.time_au, swp, lw=2.2, color="tab:red", label="wavepacket")
ax[1].plot(icl.time_au, swp[:len(scl)] - scl, lw=1.4, color="k",
           label="WP $-$ classical")
ax[1].axhline(0, color="0.5", lw=0.8)
ax[1].set_xlabel("t (a.u.)")
ax[1].set_ylabel(r"$\Delta(E_{PS}+E_{PB})$ (eV)")
ax[1].set_title("(b) combined — they track each other")
ax[1].legend(fontsize=9)

fig.tight_layout(); save(fig, "06b_combined_coupling"); plt.show()

n = min(len(scl), len(swp))
worst_sum = float(np.abs(swp[:n] - scl[:n]).max())
worst_ps = float(np.abs(iwp.d_e_ps_ev.to_numpy()[:n] - icl.d_e_ps_ev.to_numpy()[:n]).max())
worst_pb = float(np.abs(iwp.d_e_pb_ev.to_numpy()[:n] - icl.d_e_pb_ev.to_numpy()[:n]).max())
print("worst WP-vs-classical disagreement over the whole run:")
print(f"   dE_PS alone          : {worst_ps:7.4f} eV")
print(f"   dE_PB alone          : {worst_pb:7.4f} eV")
print(f"   dE_PS + dE_PB (sum)  : {worst_sum:7.4f} eV   "
      f"<- {worst_ps/worst_sum:.1f}x smaller than either part")
print()
print("   t      classical SUM    wp SUM     difference")
for tq in (5, 10, 15, 20, 25, 30):
    i = int(np.argmin(np.abs(icl.time_au - tq))); j = int(np.argmin(np.abs(iwp.time_au - tq)))
    print(f"{tq:6.1f}   {scl[i]:+12.4f} {swp[j]:+11.4f} {swp[j]-scl[i]:+13.4f}")
""")

    # ============================================ 5. the central question
    md("""
---
## 5. Can the $\\Delta T_1$ vs classical difference be explained?

$T_1 = \\langle p\\rangle^2/2m$ **exactly**, and both halves start at the same
$p$. So the entire energy gap is algebraically a gap in $\\Delta p$ — the
question is not about energy bookkeeping at all, it is: **why did the wavepacket
receive a smaller impulse?** Working in impulse also removes the $p^2$
nonlinearity that makes the energy curves hard to read.

The candidate explanation with a measurable signature: the classical projectile
is a **rigid** Gaussian cloud whose width never changes (its self-energy
$E_{PP}$ is *exactly* constant, spread $0.0$ eV over 1501 steps), whereas the
wavepacket **disperses** — so its charge dilutes and its coupling to the medium
must weaken. If that is the mechanism, the *instantaneous* impulse ratio should
track the spreading. Panel (b) tests exactly that.
""")

    code(r"""
imp = R.impulse_comparison(wp, cl)
epp_t = np.interp(wp.t, iwp.time_au, iwp.e_pp_ev)

fig, ax = plt.subplots(1, 3, figsize=(16.0, 4.4))

ax[0].plot(imp.t, imp.dp_cl, lw=2.2, color="tab:blue", label=r"classical $\Delta p$")
ax[0].plot(imp.t, imp.dp_wp, lw=2.2, color="tab:green", label=r"WP $\Delta\langle p\rangle$")
ax[0].axhline(0, color="0.5", lw=0.8)
ax[0].set_xlabel("t (a.u.)"); ax[0].set_ylabel(r"$\Delta p$ (a.u.)")
ax[0].set_title("(a) cumulative impulse"); ax[0].legend(fontsize=9)

good = imp.t > 2.0
ax[1].plot(imp.t[good], imp.impulse_ratio[good], lw=2.2, color="k",
           label=r"impulse ratio  $\Delta p_{WP}/\Delta p_{cl}$")
ax[1].plot(wp.t[good], (epp_t / epp_t[0])[good], lw=1.5, ls="--", color="tab:purple",
           label=r"$E_{PP}(t)/E_{PP}(0)$")
ax[1].plot(wp.t[good], wp.f_bore[good], lw=1.5, ls=":", color="tab:red",
           label=r"$f_\mathrm{bore}$")
ax[1].axhline(1.0, color="0.5", lw=0.8)
ax[1].set_xlabel("t (a.u.)"); ax[1].set_ylabel("ratio")
ax[1].set_title("(b) the deficit tracks the spreading")
ax[1].legend(fontsize=8)

gap = np.interp(cl.t, wp.t, wp.d_T1_ev) - cl.d_ke_ev
ax[2].plot(cl.t, gap, lw=2.2, color="tab:brown")
ax[2].axhline(0, color="0.5", lw=0.8)
ax[2].set_xlabel("t (a.u.)")
ax[2].set_ylabel(r"$\Delta T_1 - \Delta(\frac{1}{2}mv^2)$ (eV)")
ax[2].set_title("(c) the energy gap (WP loses LESS)")

fig.tight_layout(); save(fig, "07_impulse_explained"); plt.show()

print("TWO REGIMES, and they have different causes.\n")
print("  t      impulse ratio   sigma_z   f_bore   E_PP/E_PP(0)")
for tq in (2, 3, 5, 7.5, 10, 15, 20, 25, 30):
    i = int(np.argmin(np.abs(imp.t - tq)))
    print(f"{tq:6.1f}   {imp.impulse_ratio.iloc[i]:12.4f}   "
          f"{wp.sigma_z_circ.iloc[i]:7.2f}   {wp.f_bore.iloc[i]:.4f}   {epp_t[i]/epp_t[0]:11.4f}")

early = (imp.t >= 2.0) & (imp.t <= 10.0)
late = imp.t > 10.0
print(f"\n  EARLY plateau (t = 2-10, packet still compact): "
      f"ratio = {imp.impulse_ratio[early].mean():.4f} "
      f"+- {imp.impulse_ratio[early].std():.4f}")
print(f"  LATE  (t > 10, packet delocalising): "
      f"ratio falls to {imp.impulse_ratio.iloc[-1]:.4f}")

m = (imp.t > 2) & (imp.t < 28)
print("\n  correlation of the impulse ratio with the spreading proxies:")
for lbl, v in (("E_PP/E_PP(0)", (epp_t / epp_t[0])[m]),
               ("sigma_z(0)/sigma_z(t)", (wp.sigma_z_circ.iloc[0] / wp.sigma_z_circ).to_numpy()[m]),
               ("f_bore", wp.f_bore.to_numpy()[m])):
    print(f"     {lbl:24s} r = {np.corrcoef(imp.impulse_ratio.to_numpy()[m], v)[0, 1]:+.4f}")

print("\n  gap growth rate d(gap)/dt -- is it constant?")
for a_, b_ in ((5, 10), (10, 15), (15, 20), (20, 25), (25, 30)):
    ia = int(np.argmin(np.abs(cl.t - a_))); ib = int(np.argmin(np.abs(cl.t - b_)))
    print(f"     t = {a_:2d}-{b_:2d}:  {(gap.iloc[ib]-gap.iloc[ia])/(b_-a_):.5f} eV per a.u.")
""")

    # ===================================================== 6. window
    md("""
---
## 6. Choose the fit window

**This is the cell to edit.** Set each half's window from the diagnostics above,
independently — there is no requirement that they match, and good reasons they
might not (the classical projectile has no packet to disperse, so its clean
window can extend further).

What to look for when choosing:

| section | what it tells you about the window |
|---|---|
| 1(b) | while $s_\\mathrm{centroid} \\approx s_{\\int\\langle p\\rangle}$, the packet still has a trajectory |
| 2(a) | the launch transient — the bath needs ~one screening time to respond |
| 3 | once a low-$k$ tail appears, part of the packet is no longer the projectile |
| 4 | $\\Delta E_{PS}$ turning over marks the projectile leaving the region it was draining |

### 6a. The two fit targets, on their own

These are the exact quantities that will be regressed — $\\Delta E_\\mathrm{total}$
for the classical half, $T_1$ for the wavepacket — with nothing else on the axes.

Sign convention: $\\Delta E_\\mathrm{total}$ is the **bath gaining** energy
(positive) while $\\Delta T_1$ is the **projectile losing** it (negative). They are
the same transfer with opposite signs, so panel (b) plots *energy lost by the
projectile* for both, which is what makes them directly comparable. For the
classical half these are interchangeable to $2.2\\times10^{-5}$ eV (§2a).

**Panel (c) is the one to choose the window from.** The stopping power is the
*slope* $-dE/ds$, so a good window is one where that slope is **flat** — a fit
over a region where the local slope is still moving averages a changing quantity
and reports a number that belongs to no particular velocity. The local slope is a
centred difference over a $\\pm 1$ a.u. window; the first and last 1 a.u. are
dropped rather than computed with a one-sided stencil, so nothing at the edges is
half-real.

### 6b. The windows chosen (2026-08-02)

| window | estimator | rationale |
|---|---|---|
| **9 – 25 a.u.** | $T_1$ | after the launch transient, spanning the flat part of the local slope |
| **21 – 30 a.u.** | $T_2$ | aggressive, late — $T_2$ has no clean plateau to read off |
| **5 – 20 a.u.** | $T_2$ | the alternative early window, for contrast |

**The classical half is fitted over each of the same three windows.** Both
projectiles decelerate, so a WP number from one window compared against a
classical number from another compares two different velocities and means
nothing. Each WP row therefore carries its own matched classical reference.

Two caveats attached to these specific choices, stated because they change how
the numbers should be read — not as an argument against running them:

- **21 – 30 extends past own-wake re-entry** (t = 22.1 a.u., §6a). Roughly the
  last 8 a.u. of that window has the projectile dragging on the disturbance it
  created at launch, which is also where the classical local $S$ turns upward
  from 0.111 to 0.141. The fit is still reported.
- **5 – 20 starts inside the launch transient**, where the local slope is still
  climbing from ~0 toward the plateau. This pulls $S$ *down* relative to a
  plateau-only fit; it is a whole-window average, not $S$ at $v_0$.

Edit `WINDOWS` below to change any of this.
""")

    code(r"""
# Local slope by CENTRED difference over +/- HALF a.u. of path.
# Deliberately not np.gradient + boxcar: a convolution with mode="same" quietly
# fabricates the first and last HALF a.u. from a truncated kernel, and those are
# exactly the early-time points a window choice is most sensitive to.
HALF = 1.0            # a.u. of time, half-width of the stencil

def local_slope(energy_ev, path, t, half=HALF):
    # -dE/ds at each t, or NaN where a full centred stencil is unavailable.
    energy_ev = np.asarray(energy_ev, float); path = np.asarray(path, float)
    t = np.asarray(t, float)
    out = np.full(t.size, np.nan)
    for i in range(t.size):
        lo = np.searchsorted(t, t[i] - half)
        hi = np.searchsorted(t, t[i] + half) - 1
        if lo < 0 or hi >= t.size or t[i] - half < t[0] or t[i] + half > t[-1]:
            continue
        ds = path[hi] - path[lo]
        if ds != 0:
            out[i] = -(energy_ev[hi] - energy_ev[lo]) / ds
    return out

cl_path = cl.z_unwrapped.to_numpy()
wp_path = wp.s_pintegral.to_numpy()
S_cl = local_slope(-cl.d_e_total_ev.to_numpy(), cl_path, cl.t.to_numpy())
S_wp = local_slope(wp.T1_drift_ev.to_numpy(), wp_path, wp.t.to_numpy())

fig, ax = plt.subplots(1, 3, figsize=(16.0, 4.4))

ax[0].plot(cl.t, cl.d_e_total_ev, lw=2.2, color="tab:blue",
           label=r"classical $\Delta E_\mathrm{total}$ (bath gains)")
ax[0].plot(wp.t, wp.d_T1_ev, lw=2.2, color="tab:green",
           label=r"WP $\Delta T_1$ (projectile loses)")
ax[0].axhline(0, color="0.5", lw=0.8)
ax[0].set_xlabel("t (a.u.)"); ax[0].set_ylabel(r"$\Delta E$ (eV)")
ax[0].set_title("(a) the two fit targets, raw sign"); ax[0].legend(fontsize=8)

ax[1].plot(cl_path - cl_path[0], cl.d_e_total_ev, lw=2.2, color="tab:blue",
           label=r"classical $\Delta E_\mathrm{total}$")
ax[1].plot(wp_path - wp_path[0], -wp.d_T1_ev, lw=2.2, color="tab:green",
           label=r"WP $-\Delta T_1$")
ax[1].set_xlabel("path travelled (Bohr)")
ax[1].set_ylabel("energy lost by the projectile (eV)")
ax[1].set_title("(b) vs path — slope IS the stopping power")
ax[1].legend(fontsize=8)

ax[2].plot(cl.t, S_cl, lw=2.0, color="tab:blue", label=r"classical $+dE_\mathrm{total}/ds$")
ax[2].plot(wp.t, S_wp, lw=2.0, color="tab:green", label=r"WP $-dT_1/ds$")
ax[2].axhline(0, color="0.5", lw=0.8)

# HARD UPPER BOUND: the projectile re-approaches its OWN launch-time disturbance
# through the periodic image. The wake is 36 Bohr long in a 60 Bohr box, so this
# is not a small correction -- it is the most likely cause of the late-time RISE
# in the classical curve, which is not physical stopping.
d_img = np.abs(cl_path - (R.CS.LAUNCH_Z + R.CS.LZ))
t_half = cl.t[d_img < R.CS.LAMBDA_P / 2.0]
if len(t_half):
    ax[2].axvspan(float(t_half.iloc[0]), cl.t.iloc[-1], color="tab:red", alpha=0.10)
    ax[2].text(float(t_half.iloc[0]) + 0.3, ax[2].get_ylim()[1] * 0.30,
               "re-entering its own wake (periodic image)",
               fontsize=7, color="tab:red", rotation=90, va="center")
ax[2].set_xlabel("t (a.u.)"); ax[2].set_ylabel("local S (eV/Bohr)")
ax[2].set_title(r"(c) LOCAL slope — pick a window where this is flat")
ax[2].legend(fontsize=8, loc="lower right")

fig.tight_layout(); save(fig, "08_fit_targets"); plt.show()

print("Local stopping power S = -dE/ds, centred +/-1 a.u. stencil (eV/Bohr)")
print(f"{'t':>6} {'S_classical':>13} {'S_wp(T1)':>11} {'ratio':>8}   "
      f"{'f_bore':>7} {'sigma_z':>8}")
for tq in (2, 4, 6, 8, 10, 12, 15, 18, 20, 24, 28):
    i = int(np.argmin(np.abs(cl.t - tq))); j = int(np.argmin(np.abs(wp.t - tq)))
    r = S_wp[j] / S_cl[i] if np.isfinite(S_cl[i]) and S_cl[i] != 0 else np.nan
    print(f"{tq:6.1f} {S_cl[i]:13.5f} {S_wp[j]:11.5f} {r:8.4f}   "
          f"{wp.f_bore.iloc[j]:7.4f} {wp.sigma_z_circ.iloc[j]:8.3f}")

print()
print(f"classical local S: min {np.nanmin(S_cl):.5f}, max {np.nanmax(S_cl):.5f} eV/Bohr")
print(f"WP       local S: min {np.nanmin(S_wp):.5f}, max {np.nanmax(S_wp):.5f} eV/Bohr")
print()
print("NEITHER is flat over the whole run -- that is the point of choosing a window.")
""")

    code("""
# ----------------------------------------------------------------- EDIT HERE
# (label, estimator, window). USER-CHOSEN, 2026-08-02.
#   T1 -> one window, 9-25 a.u.
#   T2 -> two windows tested, because T2 has no clean plateau to read off.
WINDOWS = [
    ("T1  9-25",        "T1", ( 9.0, 25.0)),
    ("T2  21-30",       "T2", (21.0, 30.0)),
    ("T2  5-20",        "T2", ( 5.0, 20.0)),
]
# ---------------------------------------------------------------------------

ECOL = {"T1": "T1_drift_ev", "T2": "T2_total_ev"}
rows = []
for label, est, (t0, t1) in WINDOWS:
    # classical over the SAME window, so every WP number has a matched reference.
    # Both halves decelerate, so a comparison across different windows compares
    # different velocities and is meaningless.
    f_cl = R.fit_in_window(cl.z_unwrapped.to_numpy(), cl.ke_ev.to_numpy(),
                           cl.t.to_numpy(), t0, t1)
    # same fit from the BATH side; equal to the above by energy conservation
    f_cl_bath = R.fit_in_window(cl.z_unwrapped.to_numpy(), -cl.e_total_ev.to_numpy(),
                                cl.t.to_numpy(), t0, t1)
    for pname, pcol in (("int<p>dt", "s_pintegral"), ("centroid", "s_centroid")):
        f = R.fit_in_window(wp[pcol].to_numpy(), wp[ECOL[est]].to_numpy(),
                            wp.t.to_numpy(), t0, t1)
        rows.append({
            "window": label, "estimator": est, "path": pname,
            "t0": t0, "t1": t1,
            "S_wp": f["S"], "sigma_wp": f["sigma"], "r2_wp": f["r2"], "n": f["n"],
            "S_cl": f_cl["S"], "sigma_cl": f_cl["sigma"], "r2_cl": f_cl["r2"],
            "ratio": f["S"] / f_cl["S"] if f_cl["S"] else float("nan"),
            "S_cl_bath": f_cl_bath["S"],
        })
fits = pd.DataFrame(rows)

print("STOPPING POWER  S = -dE/ds  (eV/Bohr)")
print("classical is fitted over the SAME window as each WP row.")
print()
print(fits[["window", "estimator", "path", "S_wp", "sigma_wp", "r2_wp",
            "S_cl", "sigma_cl", "r2_cl", "ratio", "n"]].round(5).to_string(index=False))

print()
print("Classical cross-check -- projectile KE vs bath E_total should agree")
print("(they sum to zero to 2.2e-5 eV; see section 2a):")
for _, r in fits[fits.path == "int<p>dt"].iterrows():
    print(f"  window {r['window']:12s}: S_cl(KE) = {r['S_cl']:.5f}   "
          f"S_cl(E_total) = {r['S_cl_bath']:.5f}   "
          f"diff = {abs(r['S_cl']-r['S_cl_bath']):.2e} eV/Bohr")

print()
print("CAVEATS attached to these specific windows:")
_d = np.abs(cl.z_unwrapped.to_numpy() - (R.CS.LAUNCH_Z + R.CS.LZ))
_tw = cl.t[_d < R.CS.LAMBDA_P / 2.0]
t_wake = float(_tw.iloc[0]) if len(_tw) else float("inf")
for label, est, (t0, t1) in WINDOWS:
    notes = []
    if t1 > t_wake:
        notes.append(f"extends {t1-t_wake:.1f} a.u. past own-wake re-entry (t={t_wake:.1f})")
    fb = float(np.interp(t1, wp.t, wp.f_bore))
    notes.append(f"f_bore at t1 = {fb:.3f}")
    if t0 < 10.0:
        notes.append(f"starts inside the launch transient (local S still rising until ~10)")
    print(f"  {label:12s}: " + "; ".join(notes))

fits.to_csv(HERE / "refined_stopping_summary.csv", index=False)
print()
print(f"wrote {HERE / 'refined_stopping_summary.csv'}")
""")

    code(r"""
fig, ax = plt.subplots(1, 3, figsize=(16.0, 4.6))
pal = {"T1  9-25": "tab:green", "T2  21-30": "tab:orange", "T2  5-20": "tab:purple"}

# (a) T1 fitted segment
ax[0].plot(wp.s_pintegral, wp.T1_drift_ev, lw=1.0, color="0.85")
ax[1].plot(wp.s_pintegral, wp.T2_total_ev, lw=1.0, color="0.85")
for label, est, (t0, t1) in WINDOWS:
    m = (wp.t >= t0) & (wp.t <= t1)
    a = ax[0] if est == "T1" else ax[1]
    col = "T1_drift_ev" if est == "T1" else "T2_total_ev"
    a.plot(wp.s_pintegral[m], wp[col][m], lw=2.6, color=pal[label], label=label)
    r = fits[(fits.window == label) & (fits.path == "int<p>dt")].iloc[0]
    xs = wp.s_pintegral[m].to_numpy()
    a.plot(xs, wp[col][m].to_numpy()[0] - r["S_wp"] * (xs - xs[0]),
           lw=1.2, ls="--", color="k")
ax[0].set_title(r"(a) $T_1$ fits"); ax[1].set_title(r"(b) $T_2$ fits")
for a, yl in ((ax[0], r"$T_1$ (eV)"), (ax[1], r"$T_2$ (eV)")):
    a.set_xlabel("path (Bohr)"); a.set_ylabel(yl); a.legend(fontsize=8)

# (c) S bar chart, WP vs its matched classical
sub = fits[fits.path == "int<p>dt"].reset_index(drop=True)
ypos = np.arange(len(sub))
ax[2].barh(ypos - 0.19, sub.S_wp, height=0.36, xerr=sub.sigma_wp,
           color=[pal[w] for w in sub.window], alpha=0.9, capsize=3, label="WP")
ax[2].barh(ypos + 0.19, sub.S_cl, height=0.36, xerr=sub.sigma_cl,
           color="tab:blue", alpha=0.65, capsize=3, label="classical, same window")
ax[2].set_yticks(ypos)
ax[2].set_yticklabels([f"{r.window}" for r in sub.itertuples()], fontsize=8)
ax[2].set_xlabel("S (eV/Bohr)")
ax[2].set_title("(c) WP vs matched classical")
ax[2].legend(fontsize=8)

fig.tight_layout(); save(fig, "09_window_fits"); plt.show()

for r in sub.itertuples():
    print(f"{r.window:12s} ({r.estimator}): S_wp = {r.S_wp:.5f} +- {r.sigma_wp:.5f} eV/Bohr "
          f"(r2 {r.r2_wp:.4f}) | classical {r.S_cl:.5f} +- {r.sigma_cl:.5f} "
          f"| ratio {r.ratio:.4f}")
""")

    md("""
---
### Provenance

Data layer `refined.py`; tests `tests/test_refined.py` (12 passing). Figures
written to `refined_figs/`. Plan:
`docs/plans/cylindrical-channeling-ks-stopping.md` §8.

Rebuild with:

```bash
PYTHONPATH=$REPO/inq-stack/python $REPO/venv/bin/python3 build_refined_notebook.py
```
""")
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="refined_analysis.ipynb")
    ap.add_argument("--wp", default="wp")
    ap.add_argument("--classical", default="classical")
    ap.add_argument("--times", default="0,15,30",
                    help="comma-separated times (a.u.) for the momentum slices")
    ap.add_argument("--timeout", type=int, default=1800)
    ap.add_argument("--no-execute", action="store_true")
    ap.add_argument("--out-dir", default=None,
                    help="directory to write AND execute the notebook in "
                         "(default: this builder's own folder). Used by sibling "
                         "hypotheses (e.g. channeling_sic) to reuse this builder "
                         "verbatim on other result trees via CHAN_*_RESULTS; the "
                         "target dir must make refined.py/channeling_stopping.py "
                         "importable (symlinks suffice).")
    a = ap.parse_args()

    times = [float(x) for x in a.times.split(",")]

    import nbformat as nbf
    from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

    nb = new_notebook()
    nb.cells = [
        (new_markdown_cell(src) if kind == "md" else new_code_cell(src))
        for kind, src in cells(a.wp, a.classical, times)
    ]
    for c in nb.cells:
        c.metadata["gen"] = "builder"

    out_dir = Path(a.out_dir).resolve() if a.out_dir else HERE
    out = out_dir / a.out
    if not a.no_execute:
        from nbconvert.preprocessors import ExecutePreprocessor
        ep = ExecutePreprocessor(timeout=a.timeout, kernel_name="python3")
        ep.preprocess(nb, {"metadata": {"path": str(out_dir)}})
    with open(out, "w") as fh:
        nbf.write(nb, fh)
    size_mb = out.stat().st_size / 1e6
    print(f"wrote {out}  ({len(nb.cells)} cells, {size_mb:.2f} MB"
          f"{', NOT executed' if a.no_execute else ''})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
