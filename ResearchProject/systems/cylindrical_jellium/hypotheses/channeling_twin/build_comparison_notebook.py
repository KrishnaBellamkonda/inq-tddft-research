#!/usr/bin/env python3
"""Build THE PHASE NOTEBOOK of the annular-tube channeling twin.

Plan: docs/plans/cylindrical-channeling-ks-stopping.md

This is the deliverable: one executed notebook that puts the classical and the
wavepacket halves side by side and lets a reader decide FROM THE PLOTS whether
the aim was met. The aim has three parts and each gets its own figure block:

  RESULT     does the KS-orbital stopping power land on the classical one?
             -> figure 3 (the S bar chart with uncertainties) and figure 4 (the
                fits themselves, with residuals — a good S from a bad fit is not
                a result)
  PREMISE    did the packet actually channel?
             -> figure 5 (f_bore / f_wall / <r_perp> against the bore radius)
  MECHANISM  did that freeze var(p), which is WHY the drift channel should equal
             the classical kinetic energy?
             -> figure 6 (var(p_z) and T1-T2 against their free-evolution values)

plus the kinematics both halves share (figure 2), the path-definition consistency
check (figure 7), the pairwise energy ledger including E_PP, the wavepacket
self-Hartree that has NO classical counterpart (figure 8), and the correctness
gates (figure 9). The density matrix GIF sits at the top as the visual intuition
(.claude/rules/notebook-density-gif.md).

The arithmetic is NOT done here — it is done in channeling_stopping.py, which has
its own known-case tests (tests/test_channeling_stopping.py). This file lays out
the narrative and the plots.

Usage:
    PYTHONPATH=<repo>/inq-stack/python <repo>/venv/bin/python3 build_comparison_notebook.py
        [--out channeling_twin_comparison.ipynb] [--wp NAME] [--classical NAME]
        [--timeout SECONDS]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]

sys.path.insert(0, str(HERE))
import channeling_stopping as CS  # noqa: E402  (import-time constants used below)


# ---------------------------------------------------------------------------
# Notebook cells. Code cells are written as source strings and EXECUTED, so
# every number and every figure in the shipped .ipynb is produced from the run
# data at build time — nothing is transcribed by hand.
# ---------------------------------------------------------------------------

def cells(wp_name: str, cl_name: str) -> list[tuple[str, str]]:
    """[(kind, source)] with kind in {'md', 'code'}."""
    out: list[tuple[str, str]] = []

    def md(s: str): out.append(("md", s.strip("\n")))
    def code(s: str): out.append(("code", s.strip("\n")))

    # ---------------------------------------------------------------- title
    md(f"""
# Channeling twin — does a Kohn–Sham wavepacket reproduce the classical stopping power?

**System.** A periodic annular jellium tube: positive background between
$R_\\mathrm{{in}} = {CS.R_IN:g}$ and $R_\\mathrm{{out}} = {CS.R_OUT:g}$ Bohr,
axis $\\parallel z$, hollow bore, $L_z = {CS.LZ:g}$ Bohr in a
${CS.LX:g}\\times{CS.LY:g}\\times{CS.LZ:g}$ Bohr fully periodic cell at
$\\mathrm{{d}}x = {CS.DX:g}$. $N = {CS.N_ELEC}$ bath electrons give
$n_0 = {CS.N0:.3e}\\,a_0^{{-3}}$, i.e. $r_s = {CS.RS:.3f}$,
$\\hbar\\omega_p = {CS.OMEGA_P*CS.HA_TO_EV:.2f}$ eV, $v_F = {CS.V_FERMI:.3f}$ a.u.

**The pair.** Two runs, identical in every physical parameter, differing ONLY in
how the projectile is represented:

| | classical | wavepacket |
|---|---|---|
| projectile | rigid Gaussian **charge**, $\\sigma_\\mathrm{{pot}} = \\sigma_\\mathrm{{WP}}/\\sqrt2 = {CS.SIGMA_POT:.4f}$ Bohr | occupied **KS orbital**, $\\sigma_\\mathrm{{WP}} = {CS.SIGMA_WP:g}$ Bohr |
| dynamics | velocity-Verlet Ehrenfest from its own Hellmann–Feynman force | unitary TDDFT propagation (no CAP) |
| enters the ledger as | external potential | one more electron |

Both are launched on-axis at $z_0 = {CS.LAUNCH_Z:g}$ with $v = k_0 = {CS.V0:.4f}$ a.u.
($E = 50$ eV, $v/v_F = {CS.V0/CS.V_FERMI:.2f}$) and propagated for
{CS.N_STEPS} steps of $\\mathrm{{d}}t = {CS.DT:g}$ = {CS.N_STEPS*CS.DT:g} a.u.

## The question, and what would answer it

In **bulk**, the wavepacket's kinetic energy was contaminated: the momentum-spread
term $\\mathrm{{var}}(p)$ grew through interaction with the bath while the drift
term $\\tfrac12\\langle p\\rangle^2$ stayed flat, so $-\\mathrm{{d}}T_1/\\mathrm{{d}}s$
was not a stopping power. $\\mathrm{{var}}(p)$ is **conserved under free evolution**,
so that growth was interaction and nothing else.

**Channeling** is the proposed fix: fly the packet down a vacuum bore so it couples
to the wall only through the smooth image force. If that works, three things must
be true together, and this notebook checks all three:

1. **Result** — $S_\\mathrm{{WP}}$ lands on $S_\\mathrm{{classical}}$.
2. **Premise** — the packet stayed in the bore ($f_\\mathrm{{bore}} \\approx 1$).
3. **Mechanism** — $\\mathrm{{var}}(p_z)$ stayed frozen.

Result *without* premise and mechanism is a coincidence, not a validation. That is
why the verdict at the bottom requires all three.

## Definitions ($S_{{ij}} = -\\,\\mathrm{{d}}T_i/\\mathrm{{d}}s_j$)

| symbol | meaning | source |
|---|---|---|
| $T_1 = \\langle p^2\\rangle/2m$ | full orbital kinetic energy | `wp_momentum_stats.csv:e_kin_ha` |
| $T_2 = \\langle p\\rangle^2/2m$ | **drift** kinetic energy — the classical analogue | `px_mean, py_mean, pz_mean` |
| $s_3$ | density centroid, **circular** estimator, unwrapped | `wp_real_space_stats.csv:z_mean_circ` |
| $s_4 = \\int\\langle p_z\\rangle\\,\\mathrm{{d}}t$ | path from the mean momentum | `pz_mean` |
| $f_\\mathrm{{bore}}$ | fraction of $|\\psi|^2$ inside $r_\\perp < R_\\mathrm{{in}}$ | `wp_radial_occupancy.csv` |
| $\\mathrm{{var}}(p_z)$ | momentum spread | `sigma_pz2` |

$T_1 - T_2 = \\mathrm{{var}}(p)/2m$ is the localisation energy and is **constant** at
$3/(4\\sigma_\\mathrm{{WP}}^2) = {CS.LOCALISATION_EV:.4f}$ eV under free evolution, so
its drift reads out the mechanism directly.

**$S_{{24}}$ is the headline.** It is built from $\\langle p_z\\rangle$ on both sides
(drift energy against drift path), so it is a stopping power whether or not
$\\mathrm{{var}}(p)$ is frozen. The other three are cross-checks; if they agree with
$S_{{24}}$ that is itself evidence the packet is behaving classically.

**No in-medium path correction is needed here, and that is a result rather than an
omission.** The slab study needed $s_5 = \\int f\\,v\\,\\mathrm{{d}}t$ because 25 of its
85 Bohr were vacuum. This tube is *uniform along z*: the medium fills every $z$ the
projectile visits, so the path IS the in-medium path.

## Sources

- run definitions: [`../../scripts/channeling_twin/classical/run.cpp`](../../scripts/channeling_twin/classical/run.cpp),
  [`../../scripts/channeling_twin/wp/run.cpp`](../../scripts/channeling_twin/wp/run.cpp),
  [`../../scripts/channeling_twin/gs/run.cpp`](../../scripts/channeling_twin/gs/run.cpp)
- locked geometry: [`../../shared/configs/channeling_tube_rs3.hpp`](../../shared/configs/channeling_tube_rs3.hpp)
- analysis engine: [`channeling_stopping.py`](channeling_stopping.py) (tested in [`tests/`](tests/))
- stopping definitions, imported unchanged: `systems/jellium/hypotheses/bulk_ks_stopping/ks_stopping.py`
- this builder: [`build_comparison_notebook.py`](build_comparison_notebook.py)
""")

    # ---------------------------------------------------------------- setup
    code(f"""
import sys, math, warnings
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt

HERE = Path.cwd()
sys.path.insert(0, str(HERE))
import channeling_stopping as CS

from inqview.visualisation import style
style.apply_theme()

FIGS = HERE / "comparison_figs"; FIGS.mkdir(exist_ok=True)
HA = CS.HA_TO_EV

wp = CS.load_wp({wp_name!r})
cl = CS.load_classical({cl_name!r})
cmp_ = CS.compare(wp, cl)
W = cmp_.window

print(f"WP        : {{wp.run_dir}}  ({{wp.steps_done}}/{{wp.steps_target}} steps, complete={{wp.complete}})")
print(f"classical : {{cl.run_dir}}  ({{cl.steps_done}}/{{cl.steps_target}} steps, complete={{cl.complete}})")
print(f"fit window: t = {{W[0]:.2f}} - {{W[1]:.2f}} a.u.  (derived from the MEASURED f_bore)")
if not (wp.complete and cl.complete):
    warnings.warn("At least one half is INCOMPLETE — every number below is partial.")

def save(fig, name):
    fig.savefig(FIGS / f"{{name}}.png", dpi=150, bbox_inches="tight")
    return fig
""")

    # ------------------------------------------------- density matrix GIF
    md("""
---
## 1. Visual intuition — the density matrix

The classical and wavepacket rows are on **one shared colour scale per column**, so
a visible difference between them is a real difference in the field and not a
rescaling. The third row is their difference: **that row is the quantum effect**,
directly. Each panel is LINEAR | LOG side by side, so the low-amplitude wake tail
is visible alongside the peak.

Columns: $n(x,z,t)$ · induced $\\Delta n = n(t)-n(0)$ · instantaneous $\\Delta n = n(t)-n(t-\\Delta t)$.

The dashed lines mark the periodic cell faces at $z=\\pm L_z/2$. The bore/wall
boundaries at $r_\\perp = R_\\mathrm{in}, R_\\mathrm{out}$ run *horizontally* across
these xz panels at $x = \\pm 10$ and $\\pm 14$ Bohr — the generic twin helper does not
draw them, so read them off the axis.
""")

    code("""
from inqview.visualisation import make_twin_density_matrix
from IPython.display import Image, display

gifs = make_twin_density_matrix(
    str(cl.run_dir), str(wp.run_dir), str(FIGS / "density_matrix"),
    dt=CS.DT, slab_face=CS.LZ / 2.0,
    total_subpath="raw/vti/density_total",
    frames_max=30, fps=8,
)
if not gifs:
    print("no density frames on one or both halves — re-run with CH_SAVE_EVERY > 0")
else:
    for row, col, path, title in gifs:
        print(title)
        display(Image(filename=path))     # base64-embedded: animates on reopen
""")

    # ----------------------------------------------------- 2. kinematics
    md("""
---
## 2. Kinematics — do the two projectiles even travel together?

Before comparing energy *slopes*, check the trajectories. Both projectiles are
light (mass $m_e$) and therefore **decelerate by design**
(`.claude/rules/light-projectile-stopping.md`); the question is whether they
decelerate the *same way*. The shaded band is the fit window.

If the two $z(t)$ curves separate visibly, the two halves are not sampling the same
part of the medium and any later agreement in $S$ would be luck.
""")

    code("""
fig, ax = plt.subplots(1, 3, figsize=(13.5, 3.8))

for a in ax:
    a.axvspan(W[0], W[1], color="0.85", zorder=0, label="_")

ax[0].plot(cl.t, cl.base.z, label="classical", lw=1.6)
ax[0].plot(wp.t, wp.base.s3, label=r"WP  $s_3$ (centroid)", lw=1.6)
ax[0].plot(wp.t, wp.base.s4, "--", label=r"WP  $s_4=\\int\\langle p_z\\rangle dt$", lw=1.2)
ax[0].set_xlabel("t (a.u.)"); ax[0].set_ylabel("z (Bohr, unwrapped)")
ax[0].set_title("trajectory"); ax[0].legend(frameon=False, fontsize=8)

ax[1].plot(cl.t, cl.base.vz, label="classical $v_z$", lw=1.6)
ax[1].plot(wp.t, wp.base.pz, label=r"WP $\\langle p_z\\rangle$", lw=1.6)
ax[1].axhline(CS.V0, color="0.5", ls=":", lw=1, label="$v_0$")
ax[1].axhline(0.85 * CS.V0, color="0.5", ls="--", lw=1, label=r"$0.85\\,v_0$")
ax[1].set_xlabel("t (a.u.)"); ax[1].set_ylabel("velocity (a.u.)")
ax[1].set_title("deceleration"); ax[1].legend(frameon=False, fontsize=8)

ax[2].plot(cl.t, cl.base.T * HA, label="classical $T$", lw=1.6)
ax[2].plot(wp.t, wp.base.T2 * HA, label="WP $T_2$ (drift)", lw=1.6)
ax[2].plot(wp.t, wp.base.T1 * HA, "--", label="WP $T_1$ (total)", lw=1.2)
ax[2].set_xlabel("t (a.u.)"); ax[2].set_ylabel("kinetic energy (eV)")
ax[2].set_title("energy loss"); ax[2].legend(frameon=False, fontsize=8)

fig.suptitle("Kinematics of the twin (shaded = fit window)", y=1.02)
save(fig, "02_kinematics"); plt.show()

print(f"classical: v_z/v_0 ends at {cl.v_fraction[-1]:.3f}; "
      f"WP: <p_z>/k_0 ends at {wp.base.pz[-1] / CS.V0:.3f}")
print(f"classical KE lost {(cl.base.T[0] - cl.base.T[-1]) * HA:.2f} eV; "
      f"WP drift KE lost {(wp.base.T2[0] - wp.base.T2[-1]) * HA:.2f} eV")
""")

    # ------------------------------------------------------ 3. the result
    md("""
---
## 3. THE RESULT — stopping powers side by side

Every wavepacket definition is fitted over the **same time window** as the
classical run, so the comparison is not confounded by the two projectiles being at
different velocities. Error bars are the OLS standard error and the
window-sensitivity systematic added in quadrature.

The classical bar `S_cl_same_window` is the one to compare against. The extra
`S_cl_initial_drag` bar is the convention
`.claude/rules/light-projectile-stopping.md` prescribes for a *standalone*
classical number ($v \\ge 0.85\\,v_0$); it is shown so the two conventions can be
seen not to disagree, not as a second reference.
""")

    code("""
tab = cmp_.table()
tab.to_csv(HERE / "stopping_summary.csv", index=False)

fig, ax = plt.subplots(figsize=(8.2, 4.2))
order = ["S_13", "S_14", "S_23", "S_24", "S_cl_same_window", "S_cl_initial_drag"]
t2 = tab.set_index("estimator").reindex(order).dropna(how="all")
colours = ["#7aa6c2" if h == "wp" else "#c2864a" for h in t2["half"]]
# The headline definition gets a distinct fill so it is not lost among the checks.
colours = [("#2b6b93" if i == CS.PRIMARY_ESTIMATOR else c)
           for i, c in zip(t2.index, colours)]
ax.bar(range(len(t2)), t2["S_ev_per_bohr"], yerr=t2["uncertainty"],
       capsize=4, color=colours, edgecolor="0.25", lw=0.6)
scl = float(t2.loc["S_cl_same_window", "S_ev_per_bohr"])
ax.axhline(scl, color="#c2864a", ls="--", lw=1.2,
           label=f"classical (same window) = {scl:.2f} eV/Bohr")
ax.axhspan(scl * 0.8, scl * 1.2, color="#c2864a", alpha=0.12,
           label=r"$\\pm$20 % agreement band")
ax.set_xticks(range(len(t2)))
ax.set_xticklabels(t2.index, rotation=20, ha="right")
ax.set_ylabel("S (eV / Bohr)")
ax.set_title(f"Stopping power, fitted over t = {W[0]:.1f}–{W[1]:.1f} a.u.")
ax.legend(frameon=False, fontsize=9)
save(fig, "03_stopping_bars"); plt.show()

print(tab.to_string(index=False, float_format=lambda x: f"{x:.4g}"))
""")

    # ------------------------------------------------------- 4. the fits
    md("""
---
## 4. Are those fits any good?

A slope is only a stopping power if the underlying $T(s)$ really is linear over the
window. Top row: the data and the fitted line. Bottom row: the residuals — a bowed
residual means $S$ is drifting with velocity and the single number is an average
over the window, not a value *at* $v_0$.
""")

    code("""
picks = [("classical", cmp_.cl_same_window), (CS.PRIMARY_ESTIMATOR, cmp_.wp_fits[CS.PRIMARY_ESTIMATOR]),
         ("S_13", cmp_.wp_fits["S_13"])]
fig, ax = plt.subplots(2, len(picks), figsize=(4.4 * len(picks), 6.2), sharex="col")
for j, (name, f) in enumerate(picks):
    ax[0, j].plot(f.s_fit, f.T_fit * HA, ".", ms=2.5, label="data")
    ax[0, j].plot(f.s_fit, f.T_model * HA, "-", lw=1.6,
                  label=f"fit: S = {f.S_ev_per_bohr:.3f} ± {f.uncertainty:.3f}")
    ax[0, j].set_title(f"{name}\\n$r^2$ = {f.r2:.5f}, n = {f.n_points}", fontsize=10)
    ax[0, j].set_ylabel("T (eV)" if j == 0 else "")
    ax[0, j].legend(frameon=False, fontsize=8)
    res = (f.T_fit - f.T_model) * HA
    ax[1, j].axhline(0, color="0.6", lw=0.8)
    ax[1, j].plot(f.s_fit, res, ".", ms=2.5)
    ax[1, j].set_xlabel("s (Bohr)")
    ax[1, j].set_ylabel("residual (eV)" if j == 0 else "")
    ax[1, j].set_title(f"max |residual| = {np.max(np.abs(res)):.4f} eV", fontsize=9)
fig.suptitle("The fits behind the bars", y=1.01)
save(fig, "04_fits"); plt.show()
""")

    # ---------------------------------------------------- 5. the premise
    md("""
---
## 5. THE PREMISE — did the packet actually channel?

$f_\\mathrm{bore}(t)$ is measured on the grid every step
(`inqkit::observables::radial_occupancy`), not modelled from a Gaussian ansatz —
which matters because the packet stops being Gaussian the moment it scatters.

The **fit window ends where $f_\\mathrm{bore}$ first drops below 0.95**. The
free-dispersion estimate ($2\\sigma_d = R_\\mathrm{in}$, dotted) is drawn for
comparison only: if the measured breach comes *earlier* than the formula, the
packet is being pushed into the wall rather than merely spreading into it, which
is itself a physical finding.

Right panel: $\\langle r_\\perp\\rangle(t) \\pm \\sigma_{r_\\perp}$ against the bore and
wall radii. A packet that stays well inside $R_\\mathrm{in}$ is channeling; one whose
band crosses it is not.
""")

    code("""
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.0))

ax[0].axvspan(W[0], W[1], color="0.85", zorder=0)
ax[0].plot(wp.t, wp.f_bore, lw=1.8, label=r"$f_\\mathrm{bore}$  ($r_\\perp<R_\\mathrm{in}$)")
ax[0].plot(wp.t, wp.f_wall, lw=1.4, label=r"$f_\\mathrm{wall}$  ($R_\\mathrm{in}\\leq r_\\perp<R_\\mathrm{out}$)")
ax[0].axhline(CS.F_BORE_MIN, color="0.4", ls="--", lw=1,
              label=f"channeling threshold {CS.F_BORE_MIN}")
ax[0].axvline(CS.T_2SIGMA_AT_WALL, color="0.4", ls=":", lw=1.2,
              label=rf"free dispersion $2\\sigma_d=R_{{in}}$ ({CS.T_2SIGMA_AT_WALL:.1f} a.u.)")
ax[0].set_xlabel("t (a.u.)"); ax[0].set_ylabel(r"fraction of $|\\psi|^2$")
ax[0].set_ylim(-0.02, 1.05)
ax[0].set_title("where the packet is, radially")
ax[0].legend(frameon=False, fontsize=8, loc="center left")

ax[1].axvspan(W[0], W[1], color="0.85", zorder=0)
ax[1].fill_between(wp.t, wp.r_mean - wp.sigma_r, wp.r_mean + wp.sigma_r,
                   alpha=0.25, label=r"$\\langle r_\\perp\\rangle\\pm\\sigma_{r_\\perp}$")
ax[1].plot(wp.t, wp.r_mean, lw=1.8, label=r"$\\langle r_\\perp\\rangle$")
ax[1].plot(wp.t, 2 * CS.sigma_d(wp.t), ":", lw=1.4,
           label=r"$2\\sigma_d(t)$ free dispersion")
ax[1].axhline(CS.R_IN, color="#c2864a", ls="--", lw=1.4, label=r"$R_\\mathrm{in}$ (bore wall)")
ax[1].axhline(CS.R_OUT, color="#c2864a", ls="-.", lw=1.2, label=r"$R_\\mathrm{out}$")
ax[1].set_xlabel("t (a.u.)"); ax[1].set_ylabel(r"$r_\\perp$ (Bohr)")
ax[1].set_title("radial extent vs the bore")
ax[1].legend(frameon=False, fontsize=8)

fig.suptitle("PREMISE: is the projectile channeling?", y=1.03)
save(fig, "05_channeling"); plt.show()

print(cmp_.channel.summary())
""")

    # -------------------------------------------------- 6. the mechanism
    md("""
---
## 6. THE MECHANISM — is $\\mathrm{var}(p_z)$ frozen?

This is the figure that separates *a validated method* from *a lucky number*.

$\\mathrm{var}(p)$ is conserved under free evolution. In the bulk study it grew by
+6.8 eV through interaction with the bath while the drift term stayed flat, which
is precisely why $-\\mathrm{d}T_1/\\mathrm{d}s$ was not a stopping power there. If
channeling works, the left panel should be flat at the free value
$1/(2\\sigma_\\mathrm{WP}^2)$ and the right panel flat at
$3/(4\\sigma_\\mathrm{WP}^2) = %.4f$ eV.

Any growth here is interaction that the bore was supposed to suppress; how much of
it survives sets how far $S_{13}/S_{14}$ can legitimately drift away from
$S_{23}/S_{24}$.
""" % CS.LOCALISATION_EV)

    code("""
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.0))

ax[0].axvspan(W[0], W[1], color="0.85", zorder=0)
ax[0].plot(wp.t, wp.var_pz, lw=1.8, label=r"$\\mathrm{var}(p_z)(t)$")
ax[0].axhline(CS.VAR_P_FREE, color="0.4", ls="--", lw=1.2,
              label=rf"free value $1/(2\\sigma^2)$ = {CS.VAR_P_FREE:.5f}")
ax[0].set_xlabel("t (a.u.)"); ax[0].set_ylabel(r"$\\mathrm{var}(p_z)$  (a.u.$^2$)")
ax[0].set_title("momentum spread"); ax[0].legend(frameon=False, fontsize=8)

ax[1].axvspan(W[0], W[1], color="0.85", zorder=0)
ax[1].plot(wp.t, wp.localisation_ev, lw=1.8, label=r"$T_1-T_2$")
ax[1].axhline(CS.LOCALISATION_EV, color="0.4", ls="--", lw=1.2,
              label=rf"free value $3/(4\\sigma^2)$ = {CS.LOCALISATION_EV:.4f} eV")
ax[1].set_xlabel("t (a.u.)"); ax[1].set_ylabel("$T_1 - T_2$  (eV)")
ax[1].set_title("localisation + scattering energy")
ax[1].legend(frameon=False, fontsize=8)

fig.suptitle("MECHANISM: does channeling freeze the momentum spread?", y=1.03)
save(fig, "06_var_p_freeze"); plt.show()

print(cmp_.freeze.summary())
""")

    # ----------------------------------------- 7. path/energy consistency
    md("""
---
## 7. Consistency of the two path and two energy definitions

$s_3$ (where the density centroid is) and $s_4$ ($\\int\\langle p_z\\rangle\\,
\\mathrm{d}t$) must agree for a packet moving coherently — that is the Ehrenfest
relation. Their difference is the cleanest single number for "is this still one
packet with a well-defined position?".

The centroid is the **circular** (Resta phase) estimator, which is not optional
here: the packet is launched 2 Bohr from the $-z$ face with a density std of
%.2f Bohr, so ~24 %% of it is on the far side of the periodic cell at $t=0$. The
naive $\\langle z\\rangle$ is shown for contrast — it slides smoothly to a wrong
answer rather than jumping, which is why it cannot be repaired after the fact.
""" % CS.SIGMA_POT)

    code("""
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.0))

ax[0].axvspan(W[0], W[1], color="0.85", zorder=0)
ax[0].plot(wp.t, wp.base.s3 - wp.base.s4, lw=1.6)
ax[0].axhline(0, color="0.6", lw=0.8)
ax[0].set_xlabel("t (a.u.)"); ax[0].set_ylabel("$s_3 - s_4$  (Bohr)")
ax[0].set_title(r"Ehrenfest residual (want $\\approx 0$)")

ax[1].plot(wp.t, wp.base.s3, lw=1.6, label="$s_3$ circular centroid (unwrapped)")
ax[1].plot(wp.t, wp.base.s3_naive, ":", lw=1.4, label=r"naive $\\langle z\\rangle$ (WRONG near a face)")
ax[1].plot(wp.t, wp.base.s4, "--", lw=1.2, label="$s_4$")
ax[1].set_xlabel("t (a.u.)"); ax[1].set_ylabel("z (Bohr)")
ax[1].set_title("why the circular estimator is mandatory")
ax[1].legend(frameon=False, fontsize=8)

save(fig, "07_path_consistency"); plt.show()

resid = np.max(np.abs(wp.base.s3 - wp.base.s4))
print(f"max |s3 - s4| = {resid:.4f} Bohr over the whole run "
      f"({100*resid/abs(wp.base.s4[-1]-wp.base.s4[0]):.2f} % of the path)")
""")

    # ---------------------------------------------- 8. pairwise ledger
    md("""
---
## 8. The pairwise energy ledger — where the two representations genuinely differ

INQ's own scalars **cannot** answer this question: the two projectiles sit in
different ledger terms (the classical one in `energy_external`, the wavepacket in
`energy_hartree`), so a raw comparison of those pits a net quantity against a gross
one. The pairwise decomposition into P (projectile), S (bath), B (background) is
representation-independent and *is* comparable
(`.claude/rules/decomposed-interaction-energies.md`).

Everything is plotted as a **change from $t=0$**. Absolute $E_{SB}$, $E_{PB}$,
$E_{BB}$ carry the charged-cell $G=0$ gauge and are not comparable across
representations; their changes are.

$E_{PP}$ — the projectile **self-Hartree** — is the one term with no classical
counterpart at all. It is the uncancelled self-interaction of a wavepacket in LDA,
and it is the leading suspect for any residual discrepancy in $S$. The
carried-over bulk reading (`docs/handovers/bulk-jellium-ks-stopping.md`) is that
this SIE is a property of the *orbital*, not of its environment — so channeling
should **not** remove it. Whether $E_{PP}$ stays flat while $T_1-T_2$ does or does
not is therefore a discriminating test, and it is on this figure.
""")

    code("""
ix_wp = CS.load_interactions(wp.run_dir)
ix_cl = CS.load_interactions(cl.run_dir)

fig, ax = plt.subplots(1, 3, figsize=(14.0, 4.0))
terms = [("d_e_ss_ev", r"$\\Delta E_{SS}$ bath-bath"),
         ("d_e_ps_ev", r"$\\Delta E_{PS}$ projectile-bath"),
         ("d_e_sb_ev", r"$\\Delta E_{SB}$ bath-background")]
for a, (col, ttl) in zip(ax, terms):
    a.axvspan(W[0], W[1], color="0.85", zorder=0)
    if col in ix_cl: a.plot(ix_cl["time_au"], ix_cl[col], lw=1.5, label="classical")
    if col in ix_wp: a.plot(ix_wp["time_au"], ix_wp[col], lw=1.5, label="wavepacket")
    a.axhline(0, color="0.6", lw=0.8)
    a.set_xlabel("t (a.u.)"); a.set_title(ttl, fontsize=10)
    a.legend(frameon=False, fontsize=8)
ax[0].set_ylabel("energy change from $t=0$ (eV)")
fig.suptitle("Pairwise Coulomb ledger, classical vs wavepacket", y=1.03)
save(fig, "08a_pairwise_ledger"); plt.show()

fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.0))
ax[0].axvspan(W[0], W[1], color="0.85", zorder=0)
ax[0].plot(ix_wp["time_au"], ix_wp["e_pp_ev"], lw=1.8, label="wavepacket $E_{PP}$")
if "e_pp_ev" in ix_cl:
    ax[0].plot(ix_cl["time_au"], ix_cl["e_pp_ev"], lw=1.4, ls="--",
               label="classical $E_{PP}$ (rigid Gaussian: constant by construction)")
ax[0].set_xlabel("t (a.u.)"); ax[0].set_ylabel("$E_{PP}$ (eV)")
ax[0].set_title("projectile self-Hartree — the WP-only term")
ax[0].legend(frameon=False, fontsize=8)

# Closure gate: the pairwise terms must sum back to INQ's own scalars.
obs_wp = CS.K._concat_segments(wp.run_dir / "raw" / "observables", "observables")
m = pd.merge(ix_wp[["step", "e_hartree_check", "e_external_check"]],
             obs_wp[["step", "energy_hartree", "energy_external"]], on="step")
ax[1].plot(m["step"], (m["e_hartree_check"] - m["energy_hartree"]) * HA,
           lw=1.4, label="$E_H$ closure residual")
ax[1].plot(m["step"], (m["e_external_check"] - m["energy_external"]) * HA,
           lw=1.4, label=r"$E_\\mathrm{ext}$ closure residual")
ax[1].set_xlabel("step"); ax[1].set_ylabel("residual (eV)")
ax[1].set_title("closure: pairwise terms vs INQ's own scalars")
ax[1].legend(frameon=False, fontsize=8)
save(fig, "08b_selfhartree_closure"); plt.show()

print(f"E_PP (WP): {ix_wp['e_pp_ev'].iloc[0]:.4f} -> {ix_wp['e_pp_ev'].iloc[-1]:.4f} eV "
      f"(drift {ix_wp['e_pp_ev'].iloc[-1]-ix_wp['e_pp_ev'].iloc[0]:+.4f} eV)")
print(f"max |E_H closure residual|   = {np.max(np.abs((m['e_hartree_check']-m['energy_hartree'])*HA)):.3e} eV")
print(f"max |E_ext closure residual| = {np.max(np.abs((m['e_external_check']-m['energy_external'])*HA)):.3e} eV")
""")

    # ------------------------------------------------ 9. correctness gates
    md("""
---
## 9. Correctness gates — is either run trustworthy at all?

Read this **before** believing anything above.

- **WP energy conservation.** There is no CAP, so the Hamiltonian is Hermitian and
  time-independent and $E_\\mathrm{total}$ must be constant. Any real drift means the
  propagation itself is untrustworthy and every $S$ above is meaningless. This is a
  far stronger gate than the norm monitoring a CAP'd run has to settle for.
- **Classical Ehrenfest conservation.** $E_\\mathrm{electronic} + \\mathrm{KE}_\\mathrm{proj}
  + U_\\mathrm{proj,bg}$ must be flat; a drift means the Hellmann–Feynman force and the
  perturbation potential disagree (typically a minimum-image mismatch).
- **Off-axis excursion and $F_x$.** The tube is axially symmetric on this grid, so
  the transverse force at $r_\\perp=0$ must vanish. The classical projectile is
  integrated in full 3-D and is free to leave the axis, so $x(t)\\approx 0$ is a
  *measured* statement about channeling stability rather than an imposed constraint.
- **WP norm.** Unitary propagation, no absorber $\\Rightarrow$ exactly conserved.
""")

    code("""
fig, ax = plt.subplots(1, 3, figsize=(14.0, 3.9))

if wp.e_total.size:
    ax[0].plot(wp.t[:wp.e_total.size], (wp.e_total - wp.e_total[0]) * HA, lw=1.5,
               label=f"WP $E_{{tot}}$  (drift {wp.energy_drift_ev:+.2e} eV)")
c = cl.conserved
if c.size:
    ax[0].plot(cl.t[:c.size], (c - c[0]) * HA, lw=1.5,
               label=f"classical conserved  (drift {cl.conserved_drift_ev:+.2e} eV)")
ax[0].axhline(0, color="0.6", lw=0.8)
ax[0].set_xlabel("t (a.u.)"); ax[0].set_ylabel("drift from $t=0$ (eV)")
ax[0].set_title("energy conservation"); ax[0].legend(frameon=False, fontsize=8)

ax[1].plot(cl.t, np.hypot(cl.x, cl.y), lw=1.5, label=r"classical $r_\\perp(t)$")
ax[1].set_xlabel("t (a.u.)"); ax[1].set_ylabel(r"$r_\\perp$ (Bohr)")
ax[1].set_title(f"off-axis excursion (max {cl.off_axis_max:.2e} Bohr)")
ax[1].legend(frameon=False, fontsize=8)
axb = ax[1].twinx()
axb.plot(cl.t, cl.force_x, lw=1.0, color="0.5", label="$F_x$")
axb.set_ylabel("$F_x$ (Ha/Bohr)", color="0.4")

ax[2].plot(wp.t, wp.base.norm, lw=1.5)
ax[2].set_xlabel("t (a.u.)"); ax[2].set_ylabel(r"$\\int|\\psi_\\mathrm{WP}|^2 dV$")
ax[2].set_title(f"WP norm (drift {wp.norm_drift:.2e})")

fig.suptitle("Correctness gates", y=1.03)
save(fig, "09_gates"); plt.show()

for label, value, want in [
    ("WP energy drift (eV)", wp.energy_drift_ev, "< 1e-3"),
    ("classical conserved drift (eV)", cl.conserved_drift_ev, "< 0.05"),
    ("classical max off-axis (Bohr)", cl.off_axis_max, "~ 0"),
    ("WP norm drift", wp.norm_drift, "< 1e-6"),
]:
    print(f"  {label:<34} {value:>12.4e}   (want {want})")
""")

    # ------------------------------------------------------- 10. verdict
    md("""
---
## 10. Verdict

The verdict below is **computed**, not asserted: `channeling_stopping.compare()`
requires all three of result, premise and mechanism, and says which one failed
when it does. Its branch logic is unit-tested in
[`tests/test_channeling_stopping.py`](tests/test_channeling_stopping.py) —
including the case that matters most scientifically, *clean channeling and frozen
$\\mathrm{var}(p)$ but the stopping powers still differ*, which is a real finding
and not a failure of the run.

Read it together with §9: a verdict from a run that failed its conservation gate
is not a verdict.
""")

    code("""
print(cmp_.verdict)
print()
print(f"aim_met = {cmp_.aim_met}")
print(f"wrote {HERE / 'stopping_summary.csv'} and {len(list(FIGS.glob('*.png')))} figures to {FIGS}")
""")

    md("""
### Scope boundary (stated so it is not over-read)

This validates the classical↔quantum bridge in the **fast, channeling** regime
($v/v_F \\approx 3$, projectile in a vacuum bore). It does **not** speak to the
near-peak physics at $v \\approx v_F$: a wavepacket slow enough to probe that
regime disperses before it can traverse, so the wavepacket method is intrinsically
a *high-velocity* probe. A single velocity point is also a single point — the aim
here is to show that the KS-orbital definition lands on the classical curve at
all, not to trace $S(v)$.

The Lindhard curves drawn in the per-run notebooks are **bulk** response functions
and the projectile here is in a vacuum bore, so they are an upper reference rather
than a prediction for this geometry.
""")

    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="channeling_twin_comparison.ipynb")
    ap.add_argument("--wp", default="wp")
    ap.add_argument("--classical", default="classical")
    ap.add_argument("--timeout", type=int, default=3600)
    ap.add_argument("--no-execute", action="store_true",
                    help="write the notebook without running it (for a quick layout check)")
    a = ap.parse_args()

    import nbformat as nbf
    from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

    nb = new_notebook()
    nb.cells = [
        (new_markdown_cell(src) if kind == "md" else new_code_cell(src))
        for kind, src in cells(a.wp, a.classical)
    ]
    for c in nb.cells:
        c.metadata["gen"] = "builder"

    out = HERE / a.out
    if not a.no_execute:
        from nbconvert.preprocessors import ExecutePreprocessor
        ep = ExecutePreprocessor(timeout=a.timeout, kernel_name="python3")
        ep.preprocess(nb, {"metadata": {"path": str(HERE)}})
    with open(out, "w") as fh:
        nbf.write(nb, fh)
    print(f"wrote {out}  ({len(nb.cells)} cells"
          f"{', NOT executed' if a.no_execute else ''})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
