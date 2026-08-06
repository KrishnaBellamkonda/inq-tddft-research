#!/usr/bin/env python3
"""Build the run notebooks and the phase notebooks for the NG mass ladder.

Plan: docs/plans/nazarov-gross-slab-mass-ladder.md  (step 15)

TWO KINDS OF NOTEBOOK, DELIBERATELY DIFFERENT
  RUN notebook   one per completed run. What THIS run did: setup, the density
                 GIF, its own deposit/width/energy curves, its own S. Answers
                 "is this run trustworthy and what did it measure?"
  PHASE notebook one per campaign phase. What the run SET means together. The
                 mass ladder's phase notebook is the one that adjudicates
                 Nazarov-Gross; the others gate it.

House narrative (.claude/skills/notebook-making): context -> formulas with every
term defined -> reconstructable setup -> linked source files -> results ->
takeaway. Density-matrix GIF near the TOP and DISPLAYED inline, not merely
written to disk (.claude/rules/notebook-density-gif.md).
"""
from __future__ import annotations

import sys
from pathlib import Path

import nbformat as nbf

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]
SCRIPTS = REPO / "ResearchProject/systems/localised_jellium/scripts/ng_mass_ladder"
FIGS = HERE / "figures"

LADDER = [("classical", "cl_inf"), ("classical", "cl_m1"),
          ("wp", "wp_m3"), ("wp", "wp_m1p2"), ("wp", "wp_m1"), ("wp", "wp_m0p5")]
SWEEP = [("wp", "wp_m1_s2p0"), ("wp", "wp_m1_s3p0"), ("wp", "wp_m1_s6p0")]
PILOT = [("classical", "pilot_cl_inf"), ("wp", "pilot_wp_m1"), ("wp", "pilot_wp_m0p5")]

PRELUDE = f'''import sys, warnings
from pathlib import Path
warnings.filterwarnings("ignore")
REPO = Path(r"{REPO}")
HYP  = REPO / "ResearchProject/systems/localised_jellium/hypotheses/ng_mass_ladder"
SCRIPTS = REPO / "ResearchProject/systems/localised_jellium/scripts/ng_mass_ladder"
sys.path.insert(0, str(HYP)); sys.path.insert(0, str(REPO / "inq-stack/python"))
import numpy as np, pandas as pd
import ng_analysis as NG
from IPython.display import Image, display, Markdown
pd.set_option("display.width", 160); pd.set_option("display.max_columns", 50)
print("v0 = %.6f a.u. = %.2f v_F   r_s = %.4f   E_F = %.2f eV" % (NG.V0, NG.V0/NG.KF, NG.RS, NG.EF_EV))'''

THEORY_MD = r"""## What is being tested

Nazarov & Gross (arXiv:2510.26222, 2025) treat the projectile fully quantum
mechanically via Exact Factorization, $\Psi(\mathbf R,\underline r,t)=
\chi(\mathbf R,t)\,\Phi_{\mathbf R}(\underline r,t)$, and obtain a friction
coefficient

$$Q=\int (\hat e\cdot\nabla)V_0^{(e)}(\mathbf r)\;
\partial_\omega\,\mathrm{Im}\,\chi_1^{(e)}(\mathbf r,\mathbf r',\omega)\big|_{\omega=0}\;
(\hat e\cdot\nabla')V_0^{(e)}(\mathbf r')\,\mathrm d\mathbf r\,\mathrm d\mathbf r' .$$

This is **identical to the classical result** except that the point-charge
Coulomb potential is replaced by $V_0^{(e)}$, the potential of the projectile's
own **ground-state density** $n_0^{(n)}$. The projectile's mass enters through
**one channel only** — the width of $\chi_0$, which solves
$-\tfrac{\hbar^2}{2M}\nabla^2\chi_0+V_0^{(n)}\chi_0=E_0\chi_0$ in the well its
screening cloud digs. Their conclusion, verbatim: the mass dependence is *"a
result of the differences in the wave packets' sizes of particles with different
masses."*

So there are exactly two testable statements:

| | statement | measured by |
|---|---|---|
| **the claim** | at fixed $Z$ and $v$, $S$ depends on $M$ | `does_S_depend_on_mass` |
| **the mechanism** | $M$ acts *only* through the width | `does_S_collapse_on_width` |

## The stopping power used here

$$S=\frac{\mathrm d E_{\mathrm{bath}}}{\mathrm d s}\quad\text{over } v\ge 0.85\,v_0,\ |z|\le L_z/2$$

- $E_{\mathrm{bath}}$ — energy of the electron liquid. For the classical half this
  is INQ's `energy_total` (the projectile is an external perturbation and owns no
  INQ energy). For the wavepacket half the packet *is* an occupied orbital, so
  $E_{\mathrm{bath}}=E_{\mathrm{total}}-(T_{\rm wp}+E_{PP}+E_{PS}+E_{PB})$ from the
  pairwise decomposition.
- $s$ — projectile path. Classical: the Verlet track. Wavepacket: $\langle z\rangle$.
- **Why the deposit and not $-\mathrm dT/\mathrm ds$:** the deposit is measured on
  the *medium*, so it means the same thing for both representations and for every
  mass. A wavepacket also absorbs energy into its own spreading,
  $\mathrm{var}(p)/2M$, which never reaches the bath — counting that as "stopping"
  would inflate $S$ for exactly the light rungs where the effect is being sought.
- **Why an early window:** the projectile decelerates, so a whole-run regression
  averages $S$ over every velocity from $v_0$ down, not $S$ *at* $v_0$.

## Regime, stated honestly

$v_0=1.0743$ a.u. $=1.40\,v_F$, which is $0.875$ of the Bragg peak — **below** it,
as required, but **not** the $v\to0$ friction limit in which $Q$ is defined. This
campaign therefore tests the *claim* and the *mechanism*, and does **not**
reproduce their $Q$.
"""


def md(t):
    return nbf.v4.new_markdown_cell(t)


def code(t):
    return nbf.v4.new_code_cell(t)


def _gif_cell(half: str, tag: str) -> list:
    return [
        md("## Density evolution (watch this first)\n\n"
           "The real-space density $n(\\mathbf r,t)$ in the propagation $x$–$z$ plane "
           "(mid-$y$ slice), as **total**, **induced** $\\Delta n=n(t)-n(0)$ and "
           "**instantaneous** $\\Delta n=n(t)-n(t-\\Delta t)$. A static carpet compresses "
           "time onto one axis and hides exactly the behaviour this campaign is about — "
           "how a light packet disperses as it crosses."),
        code(f'''from inqview.visualisation import make_density_gif_battery
res = SCRIPTS / "{half}" / "results" / "{tag}"
try:
    gifs, _ = make_density_gif_battery(str(res), str(HYP / "gifs" / "{tag}"),
                                       run_label="{tag}", dt=float(NG.parse_summary(res).get("dt", 0.02)),
                                       slab_face=NG.SLAB_HALF, cap_inner=NG.CAP_INNER, frames_max=24)
    for cat, kind, path, cap in gifs[:3]:
        display(Markdown(f"**{{cat}} / {{kind}}** — {{cap}}")); display(Image(filename=path))
except Exception as exc:
    print("no density frames for this run:", repr(exc))'''),
    ]


def build_run_notebook(half: str, tag: str) -> Path | None:
    res = SCRIPTS / half / "results" / tag
    if not (res / "run_summary.txt").exists():
        return None
    nb = nbf.v4.new_notebook()
    nb.cells = [
        md(f"# Run `{tag}` — Nazarov-Gross mass ladder\n\n"
           f"Half: **{half}**. Plan: `docs/plans/nazarov-gross-slab-mass-ladder.md`.\n\n"
           f"Source files: `{SCRIPTS/half}/run.cpp`, config "
           "`shared/configs/slab_n206_L30x30x120_rs2p5.hpp`, analysis "
           "`hypotheses/ng_mass_ladder/ng_analysis.py`."),
        code(PRELUDE),
        code(f'r = NG.load_run(SCRIPTS, "{half}", "{tag}")\n'
             'print("complete:", r.complete, " mass:", r.mass, " sigma_WP:", r.sigma_wp)'),
        *_gif_cell(half, tag),
        md("## Full configuration (verbatim `run_summary.txt`)\n\n"
           "Pasted unabridged so this run is reproducible from the notebook alone."),
        code('print((r.path / "run_summary.txt").read_text())'),
        md(THEORY_MD),
        md("## Kinematics — did the projectile actually cross at near-constant velocity?\n\n"
           "If it did not, `extract_S` has no valid window and every number below is void."),
        code('trk = NG.projectile_track(r)\n'
             'if not trk.empty:\n'
             '    print(f"z: {trk.z.iloc[0]:.2f} -> {trk.z.iloc[-1]:.2f} Bohr")\n'
             '    print(f"v: {trk.v.iloc[0]:.4f} -> {trk.v.iloc[-1]:.4f} a.u. '
             '(drop {100*(1-abs(trk.v.iloc[-1]/trk.v.iloc[0])):.1f}%)")\n'
             '    print(f"crossed slab: {trk.z.min() < -NG.SLAB_HALF and trk.z.max() > NG.SLAB_HALF}")\n'
             'else:\n    print("no track")'),
        md("## Stopping power (the primary number) and its conservation cross-check"),
        code('S = NG.extract_S(r)\n'
             'print(f"S = {S.S_ev_per_bohr:.5f} +- {S.stderr:.5f} eV/Bohr   r2={S.r2:.4f}  n={S.n_points}")\n'
             'print(f"window z = {S.window_z}, mean v = {S.mean_v:.4f}, v drop {100*S.v_drop_frac:.1f}%")\n'
             'print(f"note: {S.note or \'(none)\'}")\n'
             'print(f"[cross-check only, NEVER the headline] -dKE/ds = {NG.ke_cross_check(r):.5f} eV/Bohr")'),
        md("## Width — the mechanism variable\n\n"
           "$\\sigma_{\\rm iso}$ is the 3-D geometric mean of the per-axis density widths. "
           "Width definitions must never be mixed: a transverse-vs-3-D mismatch once moved "
           "a headline number in this project by 4 percentage points."),
        code('w = NG.wp_width(r)\n'
             'if not w.empty:\n'
             '    print(f"sigma_iso: {w.sigma_iso.iloc[0]:.3f} -> {w.sigma_iso.iloc[-1]:.3f} Bohr '
             '(x{w.sigma_iso.iloc[-1]/w.sigma_iso.iloc[0]:.2f})")\n'
             '    print(f"sigma_perp end {w.sigma_perp.iloc[-1]:.3f} Bohr; 4-sigma = '
             '{4*w.sigma_perp.iloc[-1]:.1f} vs 30 Bohr cell")\n'
             '    display(w.iloc[::max(1,len(w)//10)])\n'
             'else:\n    print("classical run — rigid cloud, width is frozen by construction")'),
        md("## Energy decomposition and the pairwise P/S/B interaction energies\n\n"
           "`observables.csv` is the gross INQ ledger; `interactions.csv` is the "
           "representation-independent decomposition. The two are **not** interchangeable: "
           "for a wavepacket `energy_external` is identically zero and `energy_hartree` "
           "silently contains $E_{SS}+E_{PS}+E_{PP}$, so a raw scalar comparison against "
           "the classical half compares a net quantity with a gross one."),
        code('cols = [c for c in r.obs.columns if c.startswith("energy")]\n'
             'display(r.obs[["time_au", *cols]].iloc[::max(1,len(r.obs)//8)])\n'
             'if not r.inter.empty:\n'
             '    display(r.inter.iloc[::max(1,len(r.inter)//8)])\n'
             '    if "e_hartree_check" in r.inter and "energy_hartree" in r.obs:\n'
             '        h = np.interp(r.inter.time_au, r.obs.time_au, r.obs.energy_hartree)\n'
             '        print(f"closure |E_SS+E_PS+E_PP - E_hartree| max = '
             '{np.abs(r.inter.e_hartree_check - h).max():.3e} Ha")'),
        md("## Kinetic channels — where the energy actually went"),
        code('k = NG.kinetic_channels(r)\n'
             'if not k.empty:\n'
             '    print(f"T1 drift loss      {k.T1_drift_ev.iloc[0]-k.T1_drift_ev.iloc[-1]:+.4f} eV")\n'
             '    print(f"var(p)/2M gain     {k.var_p_over_2m_ev.iloc[-1]-k.var_p_over_2m_ev.iloc[0]:+.4f} eV")\n'
             '    dep = NG.bath_energy_ev(r)\n'
             '    if dep is not None: print(f"bath deposit       {dep[-1]:+.4f} eV")\n'
             'else:\n    print("classical run — no orbital momentum moments")'),
        md("## Takeaway\n\n"
           "_Fill in after reading: did it cross cleanly, is the fit window honest, and "
           "how much of the drift loss reached the bath rather than the packet's own spreading?_"),
    ]
    out = HERE / f"run_{tag}.ipynb"
    nbf.write(nb, str(out))
    return out


def build_phase_notebook(name: str, title: str, specs: list, extra: list | None = None) -> Path:
    nb = nbf.v4.new_notebook()
    nb.cells = [
        md(f"# {title}\n\nPlan: `docs/plans/nazarov-gross-slab-mass-ladder.md`"),
        code(PRELUDE),
        md(THEORY_MD),
        md("## The run set"),
        code(f'SPECS = {specs!r}\n'
             'tbl = NG.ladder_table(SCRIPTS, SPECS)\n'
             'display(tbl)\n'
             'tbl.to_csv(HYP / "ng_ladder_table.csv", index=False)'),
        *(extra or []),
        md("## Figures"),
        code('import make_figures as MF\n'
             'made = MF.build_all(SCRIPTS, SPECS, HYP / "figures")\n'
             'for p in made: display(Markdown(f"**{p.name}**")); display(Image(filename=str(p)))'),
    ]
    out = HERE / f"{name}.ipynb"
    nbf.write(nb, str(out))
    return out


VERDICT_CELLS = [
    md("## Verdict 1 — THE CLAIM: does $S$ depend on mass?\n\n"
       "Quantum rungs only, at fixed charge, velocity and initial width. The classical "
       "run is the $M\\to\\infty$ reference and is excluded from the statistic, because "
       "a rigid cloud cannot depend on mass by construction."),
    code('import json\nv = NG.does_S_depend_on_mass(tbl)\nprint(json.dumps(v, indent=2, default=str))'),
    md("## Verdict 2 — THE MECHANISM: does $S$ collapse onto the width?\n\n"
       "If mass acts only through width, the mass ladder and the fixed-mass $\\sigma$ "
       "sweep must fall on ONE curve against the **measured** mid-transit width. A "
       "systematic offset between the two families is mass acting through some other "
       "channel — which would be a result against Nazarov-Gross's stated mechanism, "
       "not against their claim."),
    code('c = NG.does_S_collapse_on_width(tbl)\nprint(json.dumps(c, indent=2, default=str))'),
]


def main() -> int:
    HERE.mkdir(parents=True, exist_ok=True)
    made = []
    for half, tag in PILOT + LADDER + SWEEP:
        p = build_run_notebook(half, tag)
        if p:
            made.append(p)
            print(f"  run notebook: {p.name}")

    made.append(build_phase_notebook(
        "phase_P3_pilot", "Phase P3 — pilot: can the measurement work at all?", PILOT))
    made.append(build_phase_notebook(
        "phase_P5_mass_ladder", "Phase P5 — THE MASS LADDER (the Nazarov-Gross test)",
        LADDER, VERDICT_CELLS))
    made.append(build_phase_notebook(
        "phase_P6_sigma_sweep", "Phase P6 — sigma sweep: width decoupled from mass", SWEEP))
    made.append(build_phase_notebook(
        "phase_P9_synthesis", "Phase P9 — synthesis: the campaign's answer",
        LADDER + SWEEP, VERDICT_CELLS))
    print(f"{len(made)} notebooks written to {HERE}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
