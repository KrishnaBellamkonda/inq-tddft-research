#!/usr/bin/env python3
"""Build the vacuum self-interaction notebook.

Plan: docs/plans/wp-self-interaction-correction.md

THE QUESTION THIS ANSWERS: how much extra spreading, and how much energy, does
the uncancelled self-Hartree term cost a sigma = 4 Bohr wavepacket?

The design is a difference measurement. One electron alone in vacuum has, exactly,
no self-interaction, so the exact answer is free-particle dispersion. Running the
SAME injected packet at three theory levels isolates the error without
implementing any correction:

    noninteracting  ->  the reference (and a check that the grid is converged)
    hartree         ->  + Hartree self-interaction
    lda             ->  + Hartree AND LDA xc self-interaction

Two further runs then apply the projected SIC kick and turn the measurement
into an intervention test (Tier V of the plan):

    sic_h           ->  LDA + Hartree-only correction (predicted to over-correct)
    sic_pzrun       ->  LDA + full run-consistent SIC (must match noninteracting)

Arithmetic lives in selfinteraction.py (tests: tests/test_selfinteraction.py).
This file lays out the narrative and the plots.

Usage:
    PYTHONPATH=<repo>/inq-stack/python <repo>/venv/bin/python3 \
        build_selfinteraction_notebook.py [--out selfinteraction.ipynb] [--suffix ""]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))


def cells(suffix: str) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []

    def md(s: str): out.append(("md", s.strip("\n")))
    def code(s: str): out.append(("code", s.strip("\n")))

    md(r"""
# How much does the self-Hartree term spread a wavepacket?

A single electron has, **exactly**, no self-interaction. So one electron alone in
a vacuum box must follow free-particle dispersion, which is known in closed form.
Running the *same* injected Gaussian at three theory levels therefore measures
the self-interaction **by difference**, with no self-interaction correction
implemented anywhere:

| run | contains | role |
|---|---|---|
| `noninteracting` | nothing | **the reference** — and the check that the grid is converged |
| `hartree` | Hartree self-interaction | |
| `lda` | Hartree **and** LDA xc self-interaction | |
| `sic_h` | LDA propagation **+ Hartree-only correction kick** | intervention — predicted to over-correct |
| `sic_pzrun` | LDA propagation **+ full run-consistent SIC kick** | intervention — must land back on the reference |

$$\text{lda} - \text{noninteracting} = \text{total self-interaction error}$$
$$\text{hartree} - \text{noninteracting} = \text{its Hartree part} \qquad
  \text{lda} - \text{hartree} = \text{its xc part}$$

The last split decides a design question in
`docs/plans/wp-self-interaction-correction.md`: whether removing the Hartree
self-term is sufficient, or whether the xc self-term has to go too.

The two `sic_*` runs close the loop: the correction
(`inqkit::SelfInteractionCorrection`, a projected phase kick applied to the
wavepacket orbital each step) is switched ON, so the notebook can verify —
rather than infer — what each subtracted term does. If the full correction is
exact, `sic_pzrun` must reproduce the closed-form free dispersion to the
accuracy of the grid itself.

**Why the initial state is identical across all three.** The ground state is
always computed non-interacting and is then *overwritten* by the injected packet
(`extra_states(0)` + `extra_electrons(1.0)` gives exactly one state). The runs
differ only in the theory used to **propagate**.

**Why the packet is stationary ($k_0 = 0$).** Spreading is frame-independent for
a free particle — the self-Hartree depends only on the packet's shape in its own
rest frame. A stationary packet needs no traversal length, so the box can be
small enough to afford a fine grid, and it cannot wrap in $z$.

### The closed-form reference (atomic units, $m_e = 1$)

For $\psi_0 \propto e^{-r^2/2\sigma^2}$ the **density** $|\psi|^2$ has per-axis
standard deviation

$$\sigma_\mathrm{dens}(t) = \sqrt{\tfrac{\sigma^2}{2} + \tfrac{t^2}{2\sigma^2}}$$

and the momentum distribution never changes at all:

$$\langle p_d\rangle = k_{0,d}, \qquad \mathrm{var}(p_d) = \frac{1}{2\sigma^2}
\quad \text{for all } t.$$

$\mathrm{var}(p)$ is the sharpest of the gates: it is *exactly* conserved under
free evolution, so a growing $\mathrm{var}(p)$ cannot be blamed on discretising a
spreading packet — nothing about it is supposed to change.
""")

    code(f"""
import sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt

HERE = Path.cwd()
sys.path.insert(0, str(HERE))
import selfinteraction as S

try:
    from inqview.visualisation import style
    style.apply_theme()
except Exception as exc:
    warnings.warn(f"inqview theme unavailable ({{exc}}); using matplotlib defaults")

FIGS = HERE / "si_figs"; FIGS.mkdir(exist_ok=True)
def save(fig, name):
    fig.savefig(FIGS / f"{{name}}.png", dpi=150, bbox_inches="tight")
    return fig

SUFFIX = {suffix!r}
runs = S.load_all(SUFFIX)
ref  = runs["noninteracting"]

for th, r in runs.items():
    print(f"{{th:16s}} {{len(r.t):5d}} steps, t = {{r.t[0]:.2f}} .. {{r.t[-1]:.2f}} a.u.  "
          f"complete={{r.complete}}  sigma_WP={{r.sigma_wp}}")
print()
print(f"analytic: sigma_dens(0)   = {{S.sigma_dens_free(0.0, ref.sigma_wp):.6f}} Bohr")
print(f"          sigma_dens(end) = {{S.sigma_dens_free(ref.t[-1], ref.sigma_wp):.6f}} Bohr")
print(f"          var(p_d)        = {{S.var_p_free(ref.sigma_wp):.6f}} (constant)")
""")

    # ---------------------------------------------------- 1. numerics gate
    md(r"""
---
## 1. Is the reference run actually free? (the gate that licenses everything else)

Before any self-interaction can be extracted by difference, the non-interacting
run has to reproduce the analytic solution. If it does not, the grid or the time
step is not converged, and the "self-interaction" measured below would be
contaminated by that instead.
""")

    code(r"""
g = S.numerics_gate(ref)
print("NUMERICS GATE on the non-interacting run")
print(f"  max |sigma_meas/sigma_exact - 1| = {g['max_rel_sigma_error']:.3e}   "
      f"({'PASS' if g['sigma_ok'] else 'FAIL'}, tol 5e-3)")
print(f"  max |var(p)/var_free - 1|        = {g['max_var_p_drift']:.3e}   "
      f"({'PASS' if g['var_p_ok'] else 'FAIL'}, tol 1e-3)")
print(f"  E_total drift                    = {g['e_total_drift_ev']:+.3e} eV")
print(f"  max wrap indicator               = {g['max_wrap_indicator']:.3e}  "
      f"(non-zero => density reached the box face)")
print(f"  => {'PASS' if g['passed'] else 'FAIL'}")
print()
print("CLOSURE: offline E_PP vs INQ's own energy_hartree")
print("(in vacuum the packet is the ONLY charge, so for an interacting theory")
print(" these must be the same number)")
for th, r in runs.items():
    c = S.closure(r)
    print(f"  {th:16s} max|residual| = {c['max_abs_residual_ha']:.3e} Ha   "
          f"gated={c['gated']}  {'PASS' if c['passed'] else 'FAIL'}")
""")

    # ---------------------------------------------------- 2. spreading
    md(r"""
---
## 2. The measured Gaussian width — the three theory levels

The direct answer to "how much spreading does the self-Hartree term cause": the
**measured** density width $\sigma_\mathrm{dens}(t)$ of the same injected packet
propagated at three levels of theory, against the closed-form free curve.

| curve | contains |
|---|---|
| non-interacting | nothing — must lie on the analytic curve |
| Hartree only | the bare self-repulsion, **uncancelled** |
| Hartree + LDA xc | the self-repulsion as an actual LDA calculation sees it |

Panel (a) is the width itself. Panel (b) is the same information as an
**absolute** excess $\sigma_\mathrm{meas}-\sigma_\mathrm{free}$ in Bohr, which is
the panel to read for the LDA curve: a ratio near 1.08 is easy to mistake for a
line thickness, whereas half a Bohr of extra width is a physical quantity you can
compare against the 4 Bohr bore radius.

The gap between the two interacting curves is the **LDA xc self-term cancelling
the Hartree self-term**. It is not a small correction — it is most of the effect.
""")

    code(r"""
THREE = ("noninteracting", "hartree", "lda")

fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))

ax[0].plot(ref.t, ref.sigma_free(), lw=3.2, color="0.78", zorder=0,
           label=r"analytic free $\sigma_\mathrm{dens}(t)$")
for th in THREE:
    r = runs[th]
    ax[0].plot(r.t, r.sigma_iso, lw=2.0, color=S.COLOR[th], label=S.LABEL[th])
    ax[0].annotate(f"{r.sigma_iso[-1]:.2f}",
                   xy=(r.t[-1], r.sigma_iso[-1]), xytext=(-2, 3),
                   textcoords="offset points", ha="right", fontsize=8,
                   color=S.COLOR[th])
ax[0].set_xlabel("t (a.u.)")
ax[0].set_ylabel(r"$\sigma_\mathrm{dens}$ (Bohr)")
ax[0].set_title("(a) measured Gaussian width")
ax[0].legend(fontsize=8, loc="upper left")

for th in THREE:
    r = runs[th]
    n = min(len(r.t), len(ref.t))
    ax[1].plot(r.t[:n], r.sigma_iso[:n] - r.sigma_free()[:n], lw=2.0,
               color=S.COLOR[th], label=S.LABEL[th])
ax[1].axhline(0.0, color="0.4", ls="--", lw=1.2)
ax[1].set_xlabel("t (a.u.)")
ax[1].set_ylabel(r"$\sigma_\mathrm{meas}-\sigma_\mathrm{free}$ (Bohr)")
ax[1].set_title("(b) the same, as an absolute excess")
ax[1].legend(fontsize=8, loc="upper left")

fig.tight_layout(); save(fig, "01a_measured_width"); plt.show()

print("MEASURED GAUSSIAN WIDTH sigma_dens (Bohr)")
print(f"{'t':>6} {'free':>9} {'non-int':>9} {'hartree':>9} {'lda':>9}   "
      f"{'H excess':>9} {'LDA excess':>11}")
for tq in (0, 5, 10, 15, 20, 25, 30):
    i = int(np.argmin(np.abs(ref.t - tq)))
    free = float(ref.sigma_free()[i])
    sn, sh, sl = (float(runs[t].sigma_iso[i]) for t in THREE)
    print(f"{ref.t[i]:6.1f} {free:9.4f} {sn:9.4f} {sh:9.4f} {sl:9.4f}   "
          f"{sh/sn:9.4f} {sl/sn:11.4f}")
print()
print(f"at t = {ref.t[-1]:.0f} a.u.:  Hartree alone widens the packet by "
      f"{100*(runs['hartree'].sigma_iso[-1]/runs['noninteracting'].sigma_iso[-1]-1):.1f} %,")
print(f"                 LDA (H + xc) by only "
      f"{100*(runs['lda'].sigma_iso[-1]/runs['noninteracting'].sigma_iso[-1]-1):.1f} % "
      f"-- the xc self-term cancels most of it.")
""")

    md(r"""
---
## 2b. All five levels, and the excess ratio

Panel (a) puts all five widths on the analytic curve. Panel (b) is the number
that matters: **measured width divided by the non-interacting run's width**, so
grid and propagator error (identical between them by construction) cancels.
The two `sic_*` curves read directly: `sic_pzrun` on the dashed line at 1.0 is
the correction working; `sic_h` **below** it is the over-correction — removing
only the self-repulsion leaves the attractive xc self-term, and the packet
self-binds.

The dashed line in (b) at 1.0 is "no self-interaction". The horizontal marker is
the excess measured in the *channeling* production run (**1.378**), on the SAME
3-D geometric-mean width definition used here — a definition mismatch is exactly
how this comparison went wrong once (see `selfinteraction.py`, the
`CHANNELING_*` block: the transverse $\langle r_\perp\rangle$ ratio for that run
is 1.467, and dividing a 3-D number into it understated the fraction).

It is shown for scale, and it is **not** a like-for-like situation — that packet
also met a bath and a tube wall which this vacuum run removes by construction.
What licenses transferring the vacuum number at all is the *separate* SIC-PZ
channeling run, which measures the same fraction directly; §4 compares the two.
""")

    code(r"""
fig, ax = plt.subplots(1, 3, figsize=(16.0, 4.4))

ax[0].plot(ref.t, ref.sigma_free(), lw=3.0, color="0.75",
           label=r"analytic $\sigma_\mathrm{dens}(t)$")
for th, r in runs.items():
    ax[0].plot(r.t, r.sigma_iso, lw=1.8, color=S.COLOR[th], label=S.LABEL[th])
ax[0].set_xlabel("t (a.u.)"); ax[0].set_ylabel(r"$\sigma_\mathrm{dens}$ (Bohr)")
ax[0].set_title("(a) packet width"); ax[0].legend(fontsize=8)

for th in S.INTERACTING:
    if th not in runs: continue
    e = S.effect(runs[th], ref)
    ax[1].plot(e.t, e.sigma_ratio, lw=2.2, color=S.COLOR[th], label=S.LABEL[th])
ax[1].axhline(1.0, color="0.4", ls="--", lw=1.2, label="no self-interaction")
ax[1].axhline(S.CHANNELING_EXCESS, color="k", ls=":", lw=1.2,
              label=f"channeling run ({S.CHANNELING_EXCESS:.3f}, NOT like-for-like)")
ax[1].set_xlabel("t (a.u.)")
ax[1].set_ylabel(r"$\sigma$ / $\sigma_\mathrm{non-interacting}$")
ax[1].set_title("(b) EXCESS spreading from self-interaction")
ax[1].legend(fontsize=8)

for th, r in runs.items():
    ax[2].plot(r.t, r.var_pz / S.var_p_free(r.sigma_wp), lw=2.0,
               color=S.COLOR[th], label=S.LABEL[th])
ax[2].axhline(1.0, color="0.4", ls="--", lw=1.2)
ax[2].set_xlabel("t (a.u.)")
ax[2].set_ylabel(r"$\mathrm{var}(p_z)$ / free value")
ax[2].set_title(r"(c) $\mathrm{var}(p)$ — exactly conserved when free")
ax[2].legend(fontsize=8)

fig.tight_layout(); save(fig, "01_spreading"); plt.show()

print("EXCESS WIDTH vs the non-interacting run:")
for tq in (5, 10, 15, 20, 25, 30):
    row = f"  t={tq:5.1f}: "
    for th in S.INTERACTING:
        if th not in runs: continue
        e = S.effect(runs[th], ref)
        i = int(np.argmin(np.abs(e.t - tq)))
        if e.t[i] <= e.t[-1]:
            row += f"{th} {e.sigma_ratio[i]:.4f}   "
    print(row)
""")

    # ---------------------------------------------------- 3. SIC verification
    md(r"""
---
## 3. The correction, verified against the closed form (Tier V)

`sic_pzrun` has to satisfy the **same numerics gate as the reference run** —
that is the Tier V acceptance criterion of
`docs/plans/wp-self-interaction-correction.md`: with the self-field removed the
packet is free again, so the analytic dispersion law must re-emerge to the
accuracy of the grid itself, not merely "get closer".

`sic_h` is the control that makes the xc question quantitative: removing only
the Hartree self-repulsion leaves the *attractive* xc self-term behind, and the
packet self-binds — visible in panel (a) as a distance from the closed form
*larger* than uncorrected LDA's. That asymmetry is the vacuum proof that the xc
self-term must be subtracted together with the Hartree term.

Under an active correction the KS total energy is **not** conserved, by design
(plan §0/D2 — the projected kick is the one-sided Lagrange form of TDSIC). The
conserved quantity is the *corrected* energy
$E_\mathrm{corr} = E_\mathrm{KS} - U_\mathrm{self} - E_{xc}[n_\mathrm{wp}]$,
written per step to `sic.csv` — panel (c).
""")

    code(r"""
print("TIER V GATE — the corrected run against the reference's own numerics gate:")
for th in ("noninteracting", "sic_pzrun"):
    if th not in runs: continue
    g = S.numerics_gate(runs[th])
    print(f"  {th:16s} max|sigma/analytic - 1| = {g['max_rel_sigma_error']:.3e}   "
          f"var(p) drift = {g['max_var_p_drift']:.3e}   "
          f"=> {'PASS' if g['passed'] else 'FAIL'}")

fig, ax = plt.subplots(1, 3, figsize=(16.0, 4.4))

for th in ("noninteracting", "lda", "sic_h", "sic_pzrun"):
    if th not in runs: continue
    r = runs[th]
    ax[0].semilogy(r.t, np.abs(r.sigma_iso / r.sigma_free() - 1.0) + 1e-18,
                   lw=1.8, color=S.COLOR[th], label=S.LABEL[th])
ax[0].set_xlabel("t (a.u.)")
ax[0].set_ylabel(r"$|\sigma / \sigma_\mathrm{analytic} - 1|$")
ax[0].set_title("(a) distance from the closed form")
ax[0].legend(fontsize=8)

for th in ("sic_h", "sic_pzrun"):
    if th not in runs or runs[th].u_self is None: continue
    r = runs[th]
    ax[1].plot(r.t, r.u_self * S.HA_TO_EV, lw=2.0, color=S.COLOR[th],
               label=f"{th}: $U_\\mathrm{{self}}$ (removed)")
    ax[1].plot(r.t, r.exc_self * S.HA_TO_EV, lw=2.0, ls="--", color=S.COLOR[th],
               label=(f"{th}: $E_{{xc}}[n_\\mathrm{{wp}}]$ "
                      + ("(measured, LEFT IN)" if th == "sic_h" else "(removed)")))
ax[1].axhline(0, color="0.5", lw=0.8)
ax[1].set_xlabel("t (a.u.)"); ax[1].set_ylabel("energy (eV)")
ax[1].set_title("(b) the self-terms the kick subtracts")
ax[1].legend(fontsize=8)

for th, r in runs.items():
    if r.e_corrected is not None:
        ax[2].plot(r.t, (r.e_corrected - r.e_corrected[0]) * S.HA_TO_EV,
                   lw=2.0, color=S.COLOR[th],
                   label=f"{th}: $E_\\mathrm{{corr}}$ drift")
    else:
        ax[2].plot(r.t, (r.e_total - r.e_total[0]) * S.HA_TO_EV,
                   lw=1.2, ls=":", color=S.COLOR[th],
                   label=f"{th}: $E_\\mathrm{{tot}}$ drift")
ax[2].set_xlabel("t (a.u.)"); ax[2].set_ylabel("drift (eV)")
ax[2].set_title("(c) the conserved quantity under SIC")
ax[2].legend(fontsize=7)

fig.tight_layout(); save(fig, "03_sic_verification"); plt.show()

for th in ("sic_h", "sic_pzrun"):
    if th not in runs or runs[th].cum_norm_removed is None: continue
    r = runs[th]
    print(f"{th:10s} cum. norm removed by Q-projection = {r.cum_norm_removed[-1]:.3e}"
          f"   (vacuum: no bath, so machine zero expected)")
""")

    # ---------------------------------------------------- 4. energy
    md(r"""
---
## 4. The energy cost

$E_{PP} = \tfrac12\int n_\mathrm{wp}\,\phi_\mathrm{wp}$ is the self-Hartree
energy — spurious in its entirety for a one-electron system. It **falls** as the
packet spreads (self-energy of a Gaussian $\propto 1/\sigma$), so it is a
*source*, not a store: the energy it releases goes into the packet's own internal
kinetic energy $\mathrm{var}(p)/2m$.

Note it is measured in the reference run too, where INQ reports
`energy_hartree = 0`. There it is a pure diagnostic of the packet's size — the
self-energy the packet *would* have had, which it never feels. The same holds
for `sic_pzrun`: its $E_{PP}$ curve must fall exactly along the reference's,
because the corrected packet has the reference's width at every time. Panel (c)
compares released against absorbed for the two *uncorrected* interacting runs
only — for a corrected run the released self-energy is taken out by the kick,
not stored in $\mathrm{var}(p)$.
""")

    code(r"""
fig, ax = plt.subplots(1, 3, figsize=(16.0, 4.4))

for th, r in runs.items():
    ax[0].plot(r.t, r.e_pp_ev, lw=2.0, color=S.COLOR[th], label=S.LABEL[th])
ax[0].set_xlabel("t (a.u.)"); ax[0].set_ylabel(r"$E_{PP}$ (eV)")
ax[0].set_title(r"(a) self-Hartree energy"); ax[0].legend(fontsize=8)

for th, r in runs.items():
    ax[1].plot(r.t, r.var_term_ev, lw=2.0, color=S.COLOR[th], label=S.LABEL[th])
ax[1].axhline(S.localisation_ev(ref.sigma_wp), color="0.4", ls="--", lw=1.2,
              label=r"free value $3/(4\sigma^2)$")
ax[1].set_xlabel("t (a.u.)"); ax[1].set_ylabel(r"$\mathrm{var}(p)/2m$ (eV)")
ax[1].set_title("(b) internal kinetic energy"); ax[1].legend(fontsize=8)

for th in ("hartree", "lda"):
    if th not in runs: continue
    e = S.effect(runs[th], ref)
    ax[2].plot(e.t, -e.d_e_pp_ev, lw=2.0, ls="--", color=S.COLOR[th],
               label=f"{th}: $E_{{PP}}$ released")
    ax[2].plot(e.t, e.d_var_term_ev, lw=2.0, color=S.COLOR[th],
               label=f"{th}: excess $\\mathrm{{var}}(p)/2m$")
ax[2].axhline(0, color="0.5", lw=0.8)
ax[2].set_xlabel("t (a.u.)"); ax[2].set_ylabel("energy (eV)")
ax[2].set_title("(c) released vs absorbed"); ax[2].legend(fontsize=7)

fig.tight_layout(); save(fig, "02_energy"); plt.show()

print("ENERGY LEDGER at the end of the run (eV):")
for th, r in runs.items():
    print(f"  {th:16s} E_PP {r.e_pp_ev[0]:8.4f} -> {r.e_pp_ev[-1]:8.4f}   "
          f"var(p)/2m {r.var_term_ev[0]:7.4f} -> {r.var_term_ev[-1]:7.4f}   "
          f"E_total drift {(r.e_total[-1]-r.e_total[0])*S.HA_TO_EV:+.3e}")
""")

    # ---------------------------------------------------- 5. the answer
    md(r"""
---
## 5. The answer

`excess_width_pct` is how much wider the packet is than it would be with no
self-interaction at all, at the end of the run. The `xc part` row is obtained as
`lda - hartree`, which is what decides whether a Hartree-only correction would
be sufficient. The two `sic_*` rows are the intervention outcomes: `sic_pzrun`
must sit at ~0 % excess, and `sic_h`'s **negative** excess is the over-correction
measured directly.
""")

    code(r"""
tab = S.summary_table(runs)
print("SELF-INTERACTION OF A sigma = %.1f BOHR WAVEPACKET, IN VACUUM" % ref.sigma_wp)
print(tab.round(5).to_string(index=False))
tab.to_csv(HERE / "selfinteraction_summary.csv", index=False)
print()
print(f"wrote {HERE / 'selfinteraction_summary.csv'}")

print()
print("SCALE AGAINST THE CHANNELING RUN (a ratio, NOT a decomposition --")
print("that packet also met a bath and a tube wall, which this run removes):")
cc = S.channeling_comparison(runs)
for k, v in cc.items():
    print(f"  {k:34s} {v:.4f}")
""")

    # ---------------------------------------------------- 6. density GIF
    md(r"""
---
## 6. The density, animated

`.claude/rules/notebook-density-gif.md`. In vacuum the wavepacket is the only
electron, so `density_total` **is** the wavepacket density — the five animations
below differ only in the theory used to propagate an identical initial packet
(and, for the `sic_*` runs, in the correction kick applied to it).
""")

    code(r"""
from IPython.display import Image, display
try:
    from inqview.visualisation import make_density_gif_battery
    made = False
    for th, r in runs.items():
        # vacuum box: no slab faces, no CAP -- empty cap_lines suppresses the
        # guide lines the battery would otherwise draw at +/-slab_face
        gifs, _vmax = make_density_gif_battery(
            str(r.run_dir), str(FIGS / f"density_{th}"),
            run_label=th, run_title=S.LABEL[th],
            dt=float(r.t[1] - r.t[0]) if len(r.t) > 1 else 0.02,
            slab_face=0.0, cap_inner=0.0, cap_lines=(),
            frames_max=12, fps=6)
        for _cat, _kind, path, _title in (gifs or []):
            print(f"{th}: {path}")
            display(Image(filename=path))
            made = True
    if not made:
        print("no density frames found -- re-run with WP_SAVE_EVERY > 0")
except Exception as exc:
    print(f"density GIF battery unavailable ({exc}); "
          f"VTI frames are on disk under each run's raw/vti/density_total")
""")

    md(r"""
---
### Provenance

Runs: `ResearchProject/systems/vacuum/scripts/wp_selfinteraction/results/`
(driver `shared/bin/run-wp-si.slurm`). Arithmetic: `selfinteraction.py`, tested
in `tests/test_selfinteraction.py`. Plan:
`docs/plans/wp-self-interaction-correction.md`.

The natural follow-up is the **$\sigma$ sweep**: $E_{PP}\propto 1/\sigma$, so
repeating this experiment across $\sigma_\mathrm{WP}$ maps how the error scales.
Same binary, one environment variable.
""")
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="selfinteraction.ipynb")
    ap.add_argument("--suffix", default="", help="run-name suffix, e.g. _smoke")
    ap.add_argument("--timeout", type=int, default=1800)
    ap.add_argument("--no-execute", action="store_true")
    a = ap.parse_args()

    import nbformat as nbf
    from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

    nb = new_notebook()
    nb.cells = [(new_markdown_cell(src) if k == "md" else new_code_cell(src))
                for k, src in cells(a.suffix)]
    for c in nb.cells:
        c.metadata["gen"] = "builder"

    out = HERE / a.out
    if not a.no_execute:
        from nbconvert.preprocessors import ExecutePreprocessor
        ExecutePreprocessor(timeout=a.timeout, kernel_name="python3").preprocess(
            nb, {"metadata": {"path": str(HERE)}})
    with open(out, "w") as fh:
        nbf.write(nb, fh)
    print(f"wrote {out}  ({len(nb.cells)} cells, "
          f"{out.stat().st_size/1e6:.2f} MB"
          f"{', NOT executed' if a.no_execute else ''})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
