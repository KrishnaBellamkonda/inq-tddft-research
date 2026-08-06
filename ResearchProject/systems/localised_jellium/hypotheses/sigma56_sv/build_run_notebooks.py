"""
Build + execute the notebooks for the sigma56_sv twin campaign.

    python build_run_notebooks.py                 # everything
    python build_run_notebooks.py 6.0:3.0         # one twin pair
    python build_run_notebooks.py synthesis       # campaign synthesis only

Produces, in this directory:
    run_<half>_s<sigma>_v<v>.ipynb    16 per-run notebooks (8 WP + 8 classical)
    twin_s<sigma>.ipynb               2 twin comparisons (classical vs WP)
    synthesis.ipynb                   the campaign: tables + both S figures

EVERY notebook opens with an animated density GIF (.claude/rules/notebook-density-
gif.md): the real-space n(x,z,t) matrix on the mid-y xz slice, kinds {density,
induced, instantaneous}, DISPLAYED inline via IPython.display.Image so the bytes
are embedded in the stored outputs and animate on reopen. A notebook that only
writes the GIF to disk and prints a path does not satisfy the rule.

VTIs are loaded in PHYSICAL order and are NEVER fftshifted
(.claude/rules/vti-coordinate-mapping.md); the battery functions handle this.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import nbformat as nbf
from nbclient import NotebookClient

HERE = Path(__file__).resolve().parent
SIGMAS = (5.0, 6.0)
VELOCITIES = (2.0, 2.5, 3.0, 3.5)
DT = 0.04
SLAB_FACE = 12.5
CAP_INNER = 40.0


def md(t: str) -> nbf.NotebookNode:
    return nbf.v4.new_markdown_cell(t)


def code(t: str) -> nbf.NotebookNode:
    return nbf.v4.new_code_cell(t)


PRELUDE = f"""\
import sys, warnings
from pathlib import Path
warnings.filterwarnings("ignore")
sys.path.insert(0, {str(HERE)!r})
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from IPython.display import Image, display
import s56_stopping as S
from inqview.visualisation import style
from inqview.visualisation.density_gifs import (make_density_gif_battery,
                                                make_twin_density_matrix)
style.apply_theme()
E_GS = S.e_gs_ha()
print(f"E_GS (L_z = 105, dx = 0.40) = {{E_GS:.9f}} Ha")
"""


# ---------------------------------------------------------------------------
# per-run notebook
# ---------------------------------------------------------------------------
def build_run(sigma: float, v: float, half: str) -> nbf.NotebookNode:
    name = f"{half} sigma_WP = {sigma:g}, v = {v}"
    is_wp = half == "wp"
    sd0 = sigma / 2**0.5

    # Built OUTSIDE the f-string: these contain apostrophes, and an f-string
    # expression may not contain a backslash before Python 3.12.
    if is_wp:
        projectile = (f"a Gaussian electron **wavepacket**, σ_WP = {sigma:g} Bohr "
                      f"(density std σ_WP/√2 = {sd0:.4f} Bohr at t = 0)")
        estimator = ("The wavepacket is part of the electronic ledger and the CAP "
                     "eventually removes it, so this is the medium's RETAINED "
                     "excitation. The norm correction E − T₁(1−norm) is applied: INQ "
                     "divides the orbital kinetic term by its CAP-decaying norm.")
    else:
        projectile = (f"a **classical** mass-1 electron as a moving **direct erf/r "
                      f"potential**, σ_pot = σ_WP/√2 = {sd0:.5f} Bohr — the width at "
                      f"which its charge cloud equals the wavepacket twin's t = 0 density")
        estimator = ("The projectile is an external potential and was never in the "
                     "electronic ledger, so no norm correction applies. With the CAP on, "
                     "this is the medium's retained excitation — the same estimator as "
                     "the wavepacket twin, which is why the absorber was added to this half.")

    nb = nbf.v4.new_notebook()
    nb.cells = [
        md(f"""\
# {name} — sigma56_sv

**Sweep:** `sigma56_sv` &nbsp;·&nbsp; **plan:** `docs/plans/sigma56-sv-twin.md`
&nbsp;·&nbsp; **half:** `{half}`

35 × 35 × **105** Bohr, periodicity(2), dx = 0.40, 25-Bohr jellium slab
(r_s = 4.183, N = 100), LDA/ALDA, dt = 0.04.
Launch **z = −27.5** (15 Bohr from the slab face, 12.5 Bohr from the CAP).
Absorbing bands η = −1 Ha, 12.5 Bohr per z face, inner edges at |z| = 40.

**Projectile:** {projectile}.

**The measurement:** S = [E_total(t_final) − E_GS] / 25 Bohr.
{estimator}
"""),
        code(PRELUDE),
        code(f"SIGMA, V, HALF = {sigma}, {v}, {half!r}\n"
             "RUN = S.run_dir(SIGMA, V, HALF)\n"
             "p = S.measure(SIGMA, V, HALF)\n"
             "print(RUN)\n"
             "for k, val in p.__dict__.items():\n"
             "    print(f'  {k:20s} {val}')\n"
             "assert p.complete, f'INCOMPLETE: {p.steps_done}/{p.steps_target} steps — "
             "do not quote S from this run; extend it with LJ_RESUME=1'"),

        md("""\
## 1. Density evolution (mandatory first read)

The real-space density on the mid-y **x–z** plane, animated. Three kinds:
`density` n(x,z,t), `induced` Δn = n(t) − n(0), `instantaneous` Δn = n(t) − n(t−Δt).
VTIs are in physical order — never fftshifted. Slab faces at z = ±12.5 are dashed;
the CAP inner edges sit at |z| = 40.

This is the most direct picture of what the projectile does: watch whether it
traverses, is reflected, disperses, or is captured."""),
        code(f"""\
gifs, vmax = make_density_gif_battery(
    str(RUN), str(Path.cwd() / f"run_{{HALF}}_s{sigma:g}_v{v}_figs"),
    run_label=HALF, dt={DT}, slab_face={SLAB_FACE}, cap_inner={CAP_INNER},
    frames_max=30, fps=10,
    run_title=f"{{'Wavepacket' if HALF=='wp' else 'Classical'}} sigma_WP={{SIGMA:g}}, v={{V}}")
print(f"{{len(gifs)}} GIFs, shared density vmax = {{vmax:.3e}}")
for cat, kind, path, cap in gifs:
    print(f"--- {{cat}} / {{kind}} — {{cap}}")
    display(Image(filename=path))"""),

        md("""\
## 2. The energy ledger and the deposit

`dE_raw` is INQ's `energy_total` relative to E_GS. `dE_corr` additionally removes
the norm artefact on the wavepacket half (identical to raw on the classical half).

**Read the plateau, not the endpoint alone.** S is the final value of `dE_corr`
divided by 25 Bohr; if the tail is still climbing, the run is too short and must be
extended (`LJ_RESUME=1` with a larger `LJ_N_STEPS`) rather than quoted."""),
        code("""\
tr = S.energy_trace(SIGMA, V, HALF)
fig, axes = plt.subplots(1, 2, figsize=(9, 3.2))
axes[0].plot(tr.t, tr.dE_raw, lw=1.0, label="raw")
axes[0].plot(tr.t, tr.dE_corr, lw=1.4, label="norm-corrected")
ti, to = S.transit_window(V)
for a in axes:
    a.axvspan(ti, to, color="0.85", zorder=0, label="in-slab transit")
    a.set_xlabel("t (a.u.)")
axes[0].set_ylabel(r"$E_{\\rm total}(t)-E_{\\rm GS}$ (eV)")
axes[0].legend(frameon=False, fontsize="small")
if HALF == "wp":
    axes[1].plot(tr.t, tr.norm, lw=1.2, color="#d62728")
    axes[1].set_ylabel("wavepacket norm")
else:
    axes[1].plot(tr.t, tr.dE_corr, lw=1.2)
    axes[1].set_ylabel(r"$\\Delta E$ (eV)")
plt.tight_layout(); plt.show()
print(f"S = {p.S_eV_per_Bohr:.4f} eV/Bohr   (raw {p.S_raw_eV_per_Bohr:.4f})")
print(f"E_absorbed = {p.E_absorbed_eV:.3f} eV over 25 Bohr")
print(f"plateau drift over the last 10% = {p.plateau_drift_eV:.4f} eV  "
      f"-> settled = {p.settled}")"""),

        md("""\
## 3. Pairwise interaction energies (P / S / B) and the closure gate

`interactions.csv` carries the electrostatic decomposition into projectile (P),
bath electrons (S) and neutralising background (B)
(`.claude/rules/decomposed-interaction-energies.md`). These terms are
representation-independent and are therefore the ONLY energies directly comparable
between the classical and wavepacket halves — INQ's own `energy_hartree` and
`energy_external` put the projectile in different places in the two cases.

E_PS is the term that stops the projectile. E_PP is the projectile self-Hartree —
a pure quantum residual with no classical counterpart."""),
        code("""\
ix = S._concat(S.run_dir(SIGMA, V, HALF) / "raw" / "observables", "interactions")
HA = 27.211386245988
fig, ax = plt.subplots(figsize=(6.2, 3.2))
for col, lab in [("e_ps", r"$E_{PS}$ (projectile–bath)"),
                 ("e_pp", r"$E_{PP}$ (projectile self-Hartree)"),
                 ("e_ss", r"$E_{SS}$ (bath–bath)")]:
    if col in ix:
        ax.plot(ix.time_au, (ix[col] - ix[col].iloc[0]) * HA, lw=1.2, label=lab)
ax.axvspan(ti, to, color="0.85", zorder=0)
ax.set_xlabel("t (a.u.)"); ax.set_ylabel(r"$\\Delta E$ from $t=0$ (eV)")
ax.legend(frameon=False, fontsize="small"); plt.tight_layout(); plt.show()

# CLOSURE GATE. The pairwise terms must sum back to INQ's own scalars.
# NOTE e_hartree_check / e_external_check are RECONSTRUCTED VALUES
# (= 0.5*int n_total*phi_total, ~233 Ha here), NOT residuals — printing them
# raw looks like a catastrophic failure when the ledger is in fact closing to
# ~1e-10 Ha. The residual is the DIFFERENCE against the observables scalar.
ob = S._concat(S.run_dir(SIGMA, V, HALF) / "raw" / "observables", "observables")
m = ix.merge(ob[["step", "energy_hartree", "energy_external"]], on="step")
if HALF == "wp":
    rH = (m.e_ss + m.e_ps + m.e_pp - m.energy_hartree).abs().max()
    rX = (m.e_sb + m.e_pb - m.energy_external).abs().max()
    print(f"closure  E_SS+E_PS+E_PP vs energy_hartree : max|resid| = {rH:.3e} Ha")
    print(f"closure  E_SB+E_PB      vs energy_external: max|resid| = {rX:.3e} Ha")
else:
    rH = (m.e_ss - m.energy_hartree).abs().max()
    rX = (m.e_sb + m.e_ps - m.energy_external).abs().max()
    print(f"closure  E_SS           vs energy_hartree : max|resid| = {rH:.3e} Ha")
    print(f"closure  E_SB+E_PS      vs energy_external: max|resid| = {rX:.3e} Ha")
assert rH < 1e-6 and rX < 1e-6, "pairwise ledger does NOT close against INQ's scalars"
# The bath-count column differs by half: the classical ledger writes norm_slab
# (the projectile is not in the density), the WP ledger writes norm_total and
# norm_wp (it is). Either way this is the CAP witness — loss is charge leaving.
_ncol = "norm_slab" if "norm_slab" in ix else "norm_total"
print(f"electron count [{_ncol}]: {ix[_ncol].iloc[0]:.4f} -> {ix[_ncol].iloc[-1]:.4f}"
      "   (the CAP witness: loss here is charge leaving the box)")
if "norm_wp" in ix:
    print(f"wavepacket norm:  {ix.norm_wp.iloc[0]:.4f} -> {ix.norm_wp.iloc[-1]:.4e}")"""),
    ]

    if is_wp:
        nb.cells.append(md("""\
## 4. Wavepacket momentum and width

T₁ = ⟨p²⟩/2m, T₂ = ⟨p⟩²/2m; their difference is the localisation energy
3/(4σ²) = 0.030 Ha (σ=5) / 0.021 Ha (σ=6) — two orders of magnitude below the
81.6 eV of the σ = 0.5 campaign, i.e. this packet is overwhelmingly drift energy.

The measured density width is compared against the free-dispersion law
σ_d(t) = √(σ²/2 + t²/2σ²). Departure from it inside the slab is interaction."""))
        nb.cells.append(code("""\
obs = S.run_dir(SIGMA, V, HALF) / "raw" / "observables"
mom, pos = S._concat(obs, "wp_momentum_stats"), S._concat(obs, "wp_real_space_stats")
fig, axes = plt.subplots(1, 3, figsize=(11, 3.0))
axes[0].plot(mom.time_au, mom.pz_mean, lw=1.2); axes[0].axhline(V, ls=":", color="0.5")
axes[0].set_ylabel(r"$\\langle p_z\\rangle$ (a.u.)")
T1 = mom.e_kin_ha * 27.211386245988
T2 = 0.5*(mom.px_mean**2 + mom.py_mean**2 + mom.pz_mean**2) * 27.211386245988
axes[1].plot(mom.time_au, T1, lw=1.2, label=r"$T_1=\\langle p^2\\rangle/2m$")
axes[1].plot(mom.time_au, T2, lw=1.2, label=r"$T_2=\\langle p\\rangle^2/2m$")
axes[1].set_ylabel("energy (eV)"); axes[1].legend(frameon=False, fontsize="small")
# inqkit writes VARIANCES (sigma_z2) plus a CIRCULAR std (sigma_z_circ) — there
# is no plain "z_std" column. Prefer the circular one: the naive second moment is
# discontinuous across the periodic z face once the packet reaches it.
if "sigma_z_circ" in pos:
    zwidth = pos.sigma_z_circ.to_numpy()
else:
    zwidth = np.sqrt(pos.sigma_z2.to_numpy())
axes[2].plot(pos.time_au, zwidth, lw=1.2, label="measured")
axes[2].plot(pos.time_au, S.sigma_d(pos.time_au.to_numpy(), SIGMA), ls="--",
             lw=1.0, label=r"free $\\sigma_d(t)$")
axes[2].set_ylabel(r"$\\sigma_d$ (Bohr)"); axes[2].legend(frameon=False, fontsize="small")
for a in axes: a.axvspan(ti, to, color="0.85", zorder=0); a.set_xlabel("t (a.u.)")
plt.tight_layout(); plt.show()
print(f"width at slab entry {p.sigma_d_entry:.3f} -> exit {p.sigma_d_exit:.3f} Bohr "
      f"(free-dispersion prediction; growth x{p.sigma_d_exit/p.sigma_d_entry:.2f})")
print(f"time-averaged equivalent label sigma_eq = {p.sigma_eq:.2f} Bohr (vs label {SIGMA:g})")"""))
    else:
        nb.cells.append(md("""\
## 4. Projectile trajectory

A mass-1 electron under free Ehrenfest **decelerates strongly** — that is physics,
not a fault (`.claude/rules/light-projectile-stopping.md`). Nothing gates on the
velocity drift. The direct erf/r representation means the projectile may leave the
box entirely; the in-box potential simply becomes the tail of its free-space
field, with no charge to clip and hence no exit transient."""))
        nb.cells.append(code("""\
pj = S._concat(S.run_dir(SIGMA, V, HALF) / "raw" / "observables", "projectile")
fig, axes = plt.subplots(1, 3, figsize=(11, 3.0))
axes[0].plot(pj.time_au, pj.proj_z, lw=1.2)
for f in (-12.5, 12.5): axes[0].axhline(f, ls="--", color="0.5", lw=0.8)
axes[0].axhline(40, ls=":", color="#d62728", lw=0.8)
axes[0].set_ylabel("projectile z (Bohr)")
axes[1].plot(pj.time_au, pj.proj_vz, lw=1.2); axes[1].axhline(V, ls=":", color="0.5")
axes[1].set_ylabel(r"$v_z$ (a.u.)")
axes[2].plot(pj.time_au, pj.force_z, lw=1.0)
axes[2].axhline(0, lw=0.6, color="0.5"); axes[2].set_ylabel(r"$F_z$ (a.u.)")
for a in axes: a.axvspan(ti, to, color="0.85", zorder=0); a.set_xlabel("t (a.u.)")
plt.tight_layout(); plt.show()
print(f"v: {pj.proj_vz.iloc[0]:.3f} -> {pj.proj_vz.iloc[-1]:.3f} a.u.  "
      f"(z: {pj.proj_z.iloc[0]:.1f} -> {pj.proj_z.iloc[-1]:.1f} Bohr)")"""))

    nb.cells.append(md("""\
## 5. Energy ledger — total, INQ components, pairwise decomposition

Three views of the same energy, coarse to fine.

**(a) Total.** ΔE_total(t) = E_total(t) − E_GS, the estimator's raw input. For the
WP the norm-corrected trace is also drawn: INQ divides the orbital kinetic term by
its CAP-decaying norm (`inq/src/hamiltonian/energy.hpp`), so the raw curve diverges
as norm → 0 and only the corrected one is meaningful.

**(b) INQ components** — kinetic / Hartree / xc / external, each as a change from
t = 0 so they share one scale.

**(c) Pairwise decomposition** into projectile (P), bath electrons (S) and
neutralising background (B), per `.claude/rules/decomposed-interaction-energies.md`:

| term | meaning |
|---|---|
| E_SS | bath–bath, ½∫n_S·φ_S |
| E_PP | projectile **self**-Hartree — the quantum residual; ≡ 0 for a classical point charge |
| E_PS | projectile–bath — **the interaction that stops it** |
| E_SB | bath–background |
| E_PB | projectile–background |

These are representation-INDEPENDENT; the INQ scalars are not. For the classical
half the projectile lives in `energy_external`; for the WP it is inside
`energy_hartree`. So a raw `energy_hartree` comparison across the twin compares a
net quantity with a gross one — only the pairwise terms are comparable."""))
    nb.cells.append(code("""\
obs = S.run_dir(SIGMA, V, HALF) / "raw" / "observables"
ob  = S._concat(obs, "observables")
ix  = S._concat(obs, "interactions")
tr  = S.energy_trace(SIGMA, V, HALF)
HA  = 27.211386

fig, ax = plt.subplots(1, 3, figsize=(13.5, 3.3))

ax[0].plot(tr.t, tr.dE_raw, lw=1.0, color="0.65", label="raw")
ax[0].plot(tr.t, tr.dE_corr, lw=1.5, color="#1f77b4",
           label="norm-corrected" if HALF == "wp" else "deposit")
ax[0].set_ylabel(r"$\\Delta E_{\\mathrm{total}}$ (eV)")
ax[0].set_title("(a) total energy"); ax[0].legend(fontsize="small", frameon=False)

for c, lab in [("energy_kinetic", "kinetic"), ("energy_hartree", "Hartree"),
               ("energy_xc", "xc"), ("energy_external", "external")]:
    if c in ob:
        ax[1].plot(ob.time_au, (ob[c] - ob[c].iloc[0]) * HA, lw=1.0, label=lab)
ax[1].set_ylabel(r"$\\Delta E$ from $t=0$ (eV)")
ax[1].set_title("(b) INQ components")
ax[1].legend(fontsize="small", frameon=False, ncol=2)

for c, lab in [("e_ss", r"$E_{SS}$"), ("e_pp", r"$E_{PP}$"), ("e_ps", r"$E_{PS}$"),
               ("e_sb", r"$E_{SB}$"), ("e_pb", r"$E_{PB}$")]:
    if c in ix:
        ax[2].plot(ix.time_au, (ix[c] - ix[c].iloc[0]) * HA, lw=1.0, label=lab)
ax[2].set_ylabel(r"$\\Delta E$ from $t=0$ (eV)")
ax[2].set_title("(c) pairwise P/S/B")
ax[2].legend(fontsize="small", frameon=False, ncol=2)

for a in ax:
    a.axvspan(ti, to, color="0.85", zorder=0)
    a.set_xlabel("t (a.u.)")
plt.tight_layout(); plt.show()"""))

    nb.cells.append(md("""\
## 6. The E_PS tail, and the closure gates

**Left — why the classical half needs a correction.** `S = [E_total(t_f) − E_GS]/L`
assumes the projectile–bath interaction has decayed to zero by t_final.

* **WP:** true. The CAP annihilates the packet (norm → ~1e−10), so E_PS(t_f) ≈ 1e−5 eV.
* **Classical:** never true. The projectile is a real moving charge that keeps
  going, and its monopole tail falls only as N_e/z — at t_f it sits at z ≈ 321 Bohr
  and still carries **≈ 8.5 eV**, which is 62–80 % of the raw classical "deposit".
  The dashed overlay is the bare monopole N_e/z; the two lie on top of each other.

`S_deposit_eV_per_Bohr` (printed in §1) removes it. Cross-checked against the
projectile's OWN kinetic-energy loss: the two agree to 2–5 %.

**Right — closure.** The pairwise terms must sum back to the INQ scalars, and this
is gated, not assumed:

* classical: `E_SS = E_hartree` and `E_SB + E_PS = E_external`
* WP: `E_SS + E_PS + E_PP = E_hartree` and `E_SB + E_PB = E_external`

Residuals are plotted on a log axis and asserted below 1e−6 Ha. A residual that
grows in time means the decomposition has drifted out of step with the propagator
and nothing else on this page can be trusted."""))
    nb.cells.append(code("""\
mg = ob[["step", "energy_hartree", "energy_external"]].merge(ix, on="step")
if HALF == "wp":
    rH = (mg.e_ss + mg.e_ps + mg.e_pp) - mg.energy_hartree
    rX = (mg.e_sb + mg.e_pb) - mg.energy_external
else:
    rH = mg.e_ss - mg.energy_hartree
    rX = (mg.e_sb + mg.e_ps) - mg.energy_external

fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.3))

ax[0].plot(ix.time_au, ix.e_ps * HA, lw=1.3, color="#d62728", label=r"$E_{PS}(t)$")
if HALF == "classical":
    zt = S.LAUNCH_Z + V * ix.time_au.to_numpy()
    ax[0].plot(ix.time_au, 100.0 / np.maximum(zt, 1e-9) * HA, ls="--", lw=1.0,
               color="0.35", label=r"monopole $N_e/z$")
ax[0].axhline(0, lw=0.6, color="0.5")
ax[0].axvspan(ti, to, color="0.85", zorder=0)
ax[0].set_xlabel("t (a.u.)"); ax[0].set_ylabel(r"$E_{PS}$ (eV)")
ax[0].set_title("(d) projectile-bath interaction")
ax[0].legend(fontsize="small", frameon=False)

ax[1].semilogy(mg.time_au, np.abs(rH) + 1e-20, lw=1.0, label="Hartree closure")
ax[1].semilogy(mg.time_au, np.abs(rX) + 1e-20, lw=1.0, label="external closure")
ax[1].axhline(1e-6, ls=":", color="#d62728", lw=0.9)
ax[1].set_xlabel("t (a.u.)"); ax[1].set_ylabel("|residual| (Ha)")
ax[1].set_title("(e) closure gates"); ax[1].legend(fontsize="small", frameon=False)
plt.tight_layout(); plt.show()

print(f"E_PS(t_final)       = {ix.e_ps.iloc[-1] * HA:10.4f} eV")
print(f"S raw               = {p.S_eV_per_Bohr:10.4f} eV/Bohr")
print(f"S with E_PS removed = {p.S_deposit_eV_per_Bohr:10.4f} eV/Bohr")
print(f"max |closure| H/X   = {np.abs(rH).max():.3e} / {np.abs(rX).max():.3e} Ha")
assert np.abs(rH).max() < 1e-6 and np.abs(rX).max() < 1e-6, \\
    "pairwise decomposition does not close against the INQ scalars"

if HALF == "classical":
    pj = S._concat(obs, "projectile")
    ke = pj.energy_proj_ke.to_numpy() * HA
    print(f"\\nindependent check - projectile KE loss over the whole track:")
    print(f"  S_KE = {(ke[0] - ke[-1]) / 25.0:.4f} eV/Bohr "
          f"(vs {p.S_deposit_eV_per_Bohr:.4f} from the field side)")"""))

    nb.cells.append(md(f"""\
## Takeaway

S(σ_WP = {sigma:g}, v = {v}) from the printed value above, by
S = [E_total(t_f) − E_GS − E_PS(t_f)]/25 Bohr — quote
`S_deposit_eV_per_Bohr`, **not** the raw `S_eV_per_Bohr` (§6). Its width-matched twin is the
`{'classical' if is_wp else 'wp'}` notebook at the same (σ, v) — identical in every
parameter except the projectile representation. The pair difference is the
quantum effect this campaign exists to measure."""))
    return nb


# ---------------------------------------------------------------------------
# twin notebook (one per sigma)
# ---------------------------------------------------------------------------
def build_twin(sigma: float) -> nbf.NotebookNode:
    nb = nbf.v4.new_notebook()
    nb.cells = [
        md(f"""\
# Twin comparison — σ_WP = {sigma:g} Bohr

Classical (direct erf/r potential, σ_pot = σ_WP/√2 = {sigma/2**0.5:.5f} Bohr) versus
the Gaussian electron wavepacket of the same σ_WP, at v = 2.0, 2.5, 3.0, 3.5.
Identical in every physical parameter except the projectile representation: same
ground state, box, dx, dt, step counts, CAP, launch z.

**Both halves carry the absorber**, so S = [E_total(t_f) − E_GS]/25 is the same
estimator on both — the change that makes this pair comparable at all. In every
previous campaign only the wavepacket half had a CAP, which is why its deposit
curve sat 3–5× below the classical one.
"""),
        code(PRELUDE),
        code(f"SIGMA = {sigma}\n"
             "rows = []\n"
             "for v in S.VELOCITIES:\n"
             "    for half in ('wp', 'classical'):\n"
             "        try:\n"
             "            p = S.measure(SIGMA, v, half)\n"
             "        except FileNotFoundError:\n"
             "            print(f'  MISSING {half} v={v}'); continue\n"
             "        if not p.complete:\n"
             "            print(f'  INCOMPLETE {p.run}: {p.steps_done}/{p.steps_target}')\n"
             "        rows.append(p.__dict__)\n"
             "d = pd.DataFrame(rows)\n"
             "d[d.complete][['half','v','S_eV_per_Bohr','E_absorbed_eV','norm_final',"
             "'settled','sigma_eq']].round(4)"),

        md("""\
## Density matrix: classical vs wavepacket vs their difference

Rows {classical, wavepacket, WP − classical} × columns {density, induced,
instantaneous}, on the mid-y x–z slice. The **difference row is the quantum
effect** made visible: everything the wavepacket does that a point-like classical
charge of the same width does not."""),
        code(f"""\
V_SHOW = 3.0    # representative; rerun with another velocity to compare
cl, wp = S.run_dir(SIGMA, V_SHOW, "classical"), S.run_dir(SIGMA, V_SHOW, "wp")
if cl.exists() and wp.exists():
    tiles = make_twin_density_matrix(
        str(cl), str(wp), str(Path.cwd() / f"twin_s{sigma:g}_figs"),
        dt={DT}, slab_face={SLAB_FACE}, cap_inner={CAP_INNER},
        frames_max=30, fps=10, total_subpath="raw/vti/density_total")
    for row, col, path, title in tiles:
        print(f"--- {{row}} / {{col}} — {{title}}")
        display(Image(filename=path))
else:
    print("twin pair not both present at v =", V_SHOW)"""),

        md("""\
## S(v): the pair

If the two curves coincide, the projectile's quantum nature no longer matters at
this width — which is the question the campaign was built to answer."""),
        code("""\
fig, ax = plt.subplots(figsize=(5.2, 3.4))
for half, ls, mfc in (("wp", "-", "#9467bd"), ("classical", "--", "none")):
    dd = d[(d.half == half) & d.complete].sort_values("v")
    if dd.empty: continue
    ax.plot(dd.v, dd.S_eV_per_Bohr, ls=ls, marker="o", color="#9467bd",
            mfc=mfc, ms=5, label=half)
ax.set_xlabel("$v$ (a.u.)")
ax.set_ylabel(style.axis_label("stopping_power", "$S$"))
ax.legend(frameon=False); plt.tight_layout(); plt.show()

w = d[(d.half=='wp') & d.complete].set_index('v').S_eV_per_Bohr
c = d[(d.half=='classical') & d.complete].set_index('v').S_eV_per_Bohr
common = w.index.intersection(c.index)
if len(common):
    rel = 100*(w[common]-c[common])/c[common]
    print("WP - classical, as % of classical:")
    print(rel.round(1).to_string())"""),
    ]
    return nb


# ---------------------------------------------------------------------------
# campaign synthesis
# ---------------------------------------------------------------------------
def build_synthesis() -> nbf.NotebookNode:
    nb = nbf.v4.new_notebook()
    nb.cells = [
        md("""\
# sigma56_sv — campaign synthesis

Every classical and wavepacket run of the σ_WP = 5 and 6 Bohr twin campaign, on
one S(v) axis, together with the existing σ_WP = 0.5 / 2 / 3 wavepacket traces.

**Plan:** `docs/plans/sigma56-sv-twin.md`.

### What is being tested

A wavepacket has no single width — σ_d(t) = √(σ²/2 + t²/2σ²) grows in flight — so
a classical projectile of fixed σ_pot is only a fair comparison when the packet's
label agrees with its time-average. At σ = 5/6 it does (σ_eq = 5.3–5.7 and
6.2–6.5 against labels of 5 and 6); at σ = 2 it does not (σ_eq spans 4.0–6.4
across the velocity grid). These two widths are therefore the first at which
"the classical and quantum results agree at width σ" is a well-posed statement.

### Read the caveats with the numbers

* σ = 5/6 ran at L_z = 105, launch z = −27.5; σ = 0.5/2/3 at L_z = 85, launch −24.
  Slab, r_s, dx, dt and CAP are identical.
* σ = 5/6 have CAP-on classical twins. The σ = 0.5 classical benchmark is
  **CAP-free** and therefore a different estimator (the medium's gain directly,
  not its retained excitation) — orientation only, drawn grey and dashed.
"""),
        code(PRELUDE),
        code("import build_sv_figure as B\n"
             "rc = B.main()\n"
             "print('figure builder exit code:', rc)"),
        md("## The figures"),
        code("for f in ('S_of_v_sigma56.png', 'S_of_sigma_eq.png'):\n"
             "    p = Path(f)\n"
             "    print('---', f)\n"
             "    display(Image(filename=str(p))) if p.exists() else print('  missing')"),
        md("""\
## The CAP cost on the classical half

The absorber was added to the classical runs so both halves measure the same
thing. This is what it cost: S with the CAP on minus S with η = 0, at v = 3.0,
from the same binary. A small difference means the symmetrisation was cheap; a
large one means emitted electrons carry a substantial share of the deposit and the
CAP-free and CAP-on estimators must not be mixed on one axis."""),
        code("cost = S.cap_cost()\n"
             "print(cost.round(4).to_string(index=False) if not cost.empty "
             "else 'no CAP-cost pairs yet')"),
        md("""\
## The collapse test

A σ = 6 packet at v = 2.0 has σ_eq = 6.45. The existing σ = 2, v = 2.0 run has
σ_eq = 6.35 — the same time-averaged width, reached completely differently (one
holds ~4.5 Bohr throughout, the other sweeps 2.5 → 6.6 Bohr). If their S agree,
time-averaged σ is a valid collapse variable and the pattern in the σ-sweep figure
is real. If they do not, the collapse is coincidental and each σ trace must be
read on its own."""),
        code("""\
new = S.table(); old = B.legacy_wp()
if not new.empty and not old.empty:
    a = new[(new.sigma_wp==6.0)&(new.v==2.0)&(new.half=='wp')&new.complete]
    b = old[(old.sigma_wp==2.0)&(old.v==2.0)]
    if len(a) and len(b):
        Sa, Sb = float(a.S_eV_per_Bohr.iloc[0]), float(b.S_eV_per_Bohr.iloc[0])
        print(f"sigma=6, v=2.0 : sigma_eq {float(a.sigma_eq.iloc[0]):.2f}  S = {Sa:.3f} eV/Bohr")
        print(f"sigma=2, v=2.0 : sigma_eq {float(b.sigma_eq.iloc[0]):.2f}  S = {Sb:.3f} eV/Bohr")
        print(f"difference     : {100*(Sa-Sb)/Sb:+.1f} % of the sigma=2 value")
    else:
        print("one of the two collapse-test points is missing or incomplete")"""),
    ]
    return nb


def execute(nb: nbf.NotebookNode, out: Path) -> int:
    """Execute in this directory (relative figure paths land beside the notebook)
    and count errors. allow_errors=True so one broken cell does not lose the rest;
    the error count is what the caller gates on."""
    try:
        NotebookClient(nb, timeout=3600, kernel_name="python3",
                       resources={"metadata": {"path": str(HERE)}},
                       allow_errors=True).execute()
    except Exception as exc:                                  # noqa: BLE001
        print(f"  EXECUTION FAILED {out.name}: {type(exc).__name__}: {exc}")
        nbf.write(nb, out)
        return 1
    n_err = sum(1 for c in nb.cells if c.get("cell_type") == "code"
                for o in c.get("outputs", []) if o.get("output_type") == "error")
    nbf.write(nb, out)
    print(f"  {out.name}: {n_err} error(s)")
    return n_err


def main() -> int:
    target = sys.argv[1] if len(sys.argv) > 1 else "all"
    errs = 0

    if target == "synthesis":
        return execute(build_synthesis(), HERE / "synthesis.ipynb")

    if ":" in target:
        s, v = target.split(":")
        sigma, vel = float(s), float(v)
        for half in ("wp", "classical"):
            errs += execute(build_run(sigma, vel, half),
                            HERE / f"run_{half}_s{sigma:g}_v{vel}.ipynb")
        return errs

    force = os.environ.get("S56_NB_FORCE", "0") == "1"
    for sigma in SIGMAS:
        built_any = False
        for v in VELOCITIES:
            for half in ("wp", "classical"):
                out = HERE / f"run_{half}_s{sigma:g}_v{v}.ipynb"
                why = skip_reason(sigma, v, half, out, force)
                if why:
                    print(f"  SKIP ({why}) s{sigma:g} v{v} {half}")
                    continue
                errs += execute(build_run(sigma, v, half), out)
                built_any = True
        # The twin notebook cross-plots both halves at every v, so it is only
        # meaningful once this sigma has complete runs to compare.
        if any(is_complete(sigma, v, h) for v in VELOCITIES
               for h in ("wp", "classical")):
            errs += execute(build_twin(sigma), HERE / f"twin_s{sigma:g}.ipynb")
        else:
            print(f"  SKIP (no complete runs) twin s{sigma:g}")
    errs += execute(build_synthesis(), HERE / "synthesis.ipynb")
    print(f"\ntotal errors: {errs}")
    return 0 if errs == 0 else 1


def S_exists(sigma: float, v: float, half: str) -> bool:
    import s56_stopping as _S
    return (_S.run_dir(sigma, v, half) / "raw" / "observables").exists()


def is_complete(sigma: float, v: float, half: str) -> bool:
    """True only if the run reached its full step target."""
    import s56_stopping as _S
    try:
        return bool(_S.measure(sigma, v, half, True).complete)
    except Exception:
        return False


def skip_reason(sigma: float, v: float, half: str, out: Path,
                force: bool) -> str | None:
    """Why this notebook should NOT be rebuilt, or None to build it.

    Two guards, both learned the hard way on 2026-08-03:

    * INCOMPLETE runs are skipped. A notebook off a half-finished run renders a
      density GIF that stops mid-flight and an S value from a trace that never
      plateaued -- it LOOKS finished, which is worse than absent.
    * ALREADY-FRESH notebooks are skipped. Each WP notebook is ~140 MB and takes
      ~13 min because it re-reads every VTI frame to embed the GIF. The build is
      routinely interrupted (session teardown, walltime), so it has to be
      resumable: rebuild only if the run's observables are NEWER than the
      notebook. S56_NB_FORCE=1 overrides.
    """
    if not S_exists(sigma, v, half):
        return "no run"
    if not is_complete(sigma, v, half):
        return "incomplete"
    if force or not out.exists():
        return None
    import s56_stopping as _S
    obs = _S.run_dir(sigma, v, half) / "raw" / "observables"
    newest = max((p.stat().st_mtime for p in obs.glob("*.csv")), default=0.0)
    # THIS builder counts as an input too. Adding a section (e.g. the energy
    # ledger, 2026-08-03) changes what the notebook should contain even though
    # the run data is untouched -- without this the guard would report every
    # notebook "up to date" and silently ship the old layout forever.
    newest = max(newest, Path(__file__).stat().st_mtime,
                 (HERE / "s56_stopping.py").stat().st_mtime)
    if out.stat().st_mtime < newest:
        return None
    return "up to date" if _is_executed(out) else None


def _is_executed(nb_path: Path) -> bool:
    """A notebook only counts as done if it PARSES and every code cell ran clean.

    mtime alone cannot tell a finished 136 MB notebook from one truncated by a
    walltime kill mid-write -- and a half-written notebook with a fresh mtime
    would be skipped forever. Parse failure, an unexecuted cell, or a stored
    error all mean rebuild.
    """
    try:
        nb = nbf.read(nb_path, as_version=4)
    except Exception:
        return False
    for c in nb.cells:
        if c.cell_type != "code":
            continue
        outs = c.get("outputs") or []
        if not outs or any(o.get("output_type") == "error" for o in outs):
            return False
    return True


if __name__ == "__main__":
    sys.path.insert(0, str(HERE))
    raise SystemExit(main())
