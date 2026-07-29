#!/usr/bin/env python3
"""Build the MASTER campaign notebook: all phases of the localised-jellium
scattering study, systematically presented, with the headline figure per phase and
links to the per-phase study notebooks.

    PYTHONPATH=/local/data/public/skcb2/tddft/inq-stack/python \
    /local/data/public/skcb2/tddft/venv/bin/python3 build_campaign_report.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _nbreport import md, code, setup_cell, embed, set_outdir, build, SYS

HYP = os.path.join(SYS, "hypotheses")
OUT = os.path.join(HYP, "localised_jellium_campaign_study.ipynb")
set_outdir(os.path.dirname(OUT))

cells = []

# 1. Title + the question + phase map -----------------------------------------
cells += [md(
r"""# Localised jellium scattering — campaign master notebook

**All phases, end to end: a finite slab of electron gas, a projectile fired into it
from vacuum, and its stopping power.** This notebook is the systematic overview; each
phase has its own self-contained study notebook (linked below) with the full setup,
formulas, and per-run diagnostics.

## The overarching question

> Can we build a **localised** electron gas in INQ — a finite positive background
> confined to a slab, charge-neutral, injected as a *static perturbation with no
> engine edit* — fire a projectile into it from vacuum, absorb the projectile with a
> CAP, and read off a **stopping power** $S=\Delta E/x$?

## Phase map

| Phase | Goal | Result | Study notebook |
|---|---|---|---|
| **1** Implement | wrapper-only static background $v_\mathrm{bg}=-\mathrm{poisson}(n_+)$ | T0 ✓ (neutral, attractive) | [01](01_slab_validation/slab_validation_study.ipynb) |
| **2** Validate | GS interior = bulk HEG, clean surface, static well | T1 ✓ (2.0% flat), T3.4 ✓ (2e-8 Ha drift) | [01](01_slab_validation/slab_validation_study.ipynb) |
| **3** WP baseline | wave packet traverses, closed-box energy conserved | ✓ ($dE\approx-1.1$ mHa) | [02](02_projectile_slab/projectile_slab_study.ipynb) |
| **5** CAP stopping | two-sided sin² CAP; measure $S$ | **$S\approx0.71$ eV/Bohr** (PROVISIONAL) | [03](03_cap_stopping/cap_stopping_study.ipynb) |

**The full-suite re-run is complete** (`fullsuite_wp` p3wp/p5wp + `fullsuite_classical`
p5cl, all on inq-study), so the three-way **total/bath/WP** density decomposition,
momentum, GS-overlap excitations and per-component energetics are all available.

## Notebook map

| Notebook | What it holds |
|---|---|
| **this** campaign master | end-to-end overview, headline, links |
| [01 slab validation](01_slab_validation/slab_validation_study.ipynb) | GS = bulk HEG (T0/T1/T3.4) |
| [02 projectile (study)](02_projectile_slab/projectile_slab_study.ipynb) | Phase-3 WP baseline + 9-GIF decomposition |
| [03 CAP stopping (study)](03_cap_stopping/cap_stopping_study.ipynb) | Phase-5 WP+CAP & classical+CAP + decompositions + $S$ |
| [p3wp deep-dive](02_projectile_slab/p3wp_run_notebook.ipynb) | full per-run battery (WP, no CAP) |
| [p5wp deep-dive](03_cap_stopping/p5wp_run_notebook.ipynb) | full per-run battery (WP + CAP) |
| [p5cl deep-dive](03_cap_stopping/p5cl_run_notebook.ipynb) | full per-run battery (classical + CAP) |

Grounding: Lang–Kohn (jellium surfaces); GPAW charge-neutral background recipe.
Theory: `docs/notes/localised-jellium-theory.md`. Mechanism ADR:
`docs/adr/0008-localised-jellium-background-perturbation.md`. Campaign prompt:
`docs/campaigns/localised_jellium/localised_jellium_campaign.md`.
""")]

# 2. Conventions + headline parameters ----------------------------------------
cells += [md(
r"""## System at a glance

Atomic units (Hartree, Bohr); $1\,\mathrm{Ha}=27.2114\,\mathrm{eV}$.

| Quantity | Value |
|---|---|
| Cell | 50 Bohr cubic, periodic, $dx=0.5$ |
| Background | slab, 25 Bohr thick, full $x,y$, sharp edge |
| $N$ electrons | 234 → $n_0=3.744\times10^{-3}$, $r_s=3.995$ (Na metal) |
| XC / engine | LDA / stock `inq` (Ph 1–3), `inq-study` (Ph 5, complex CAP) |
| Projectile | $\sigma=0.5$ Bohr, $E=100$ eV, $k_0=2.711$ |
| CAP | two-sided sin², $\eta=-0.5$ Ha, 7.5 Bohr/side, in vacuum |
""")]

cells += [setup_cell()]
cells += [code(
"N=234; V=50.0*50.0*25.0; n0=N/V\n"
"rs=(3/(4*np.pi*n0))**(1/3); kF=(3*np.pi**2*n0)**(1/3)\n"
"k0=np.sqrt(2*100.0/27.21138625)\n"
"print(f'n0={n0:.4e} a0^-3   r_s={rs:.3f}   k_F={kF:.3f} a0^-1   projectile k0={k0:.3f} a0^-1')")]

# 3. Per-phase headline -------------------------------------------------------
cells += [md(
r"""## Phase 1+2 — the slab is a real electron gas

Planar-averaged GS density: flat at $n_0$ inside, peaks in the slab, ~0 in vacuum;
interior deviation 2.0%. Kinetic/electron within ~3% of the HEG value. The well is
static to machine precision. **Full detail:**
[`01_slab_validation/slab_validation_study.ipynb`](01_slab_validation/slab_validation_study.ipynb).""")]
cells += [embed(os.path.join(HYP, "01_slab_validation", "n_of_z.png"))]

cells += [md(
r"""## Phase 3 — wave-packet baseline (closed box)

A 100 eV WP traverses the slab; total energy conserved to ~1 mHa (no sink by
construction). This validates the projectile/observable machinery. **Full detail:**
[`02_projectile_slab/projectile_slab_study.ipynb`](02_projectile_slab/projectile_slab_study.ipynb).""")]
cells += [embed(os.path.join(HYP, "02_projectile_slab", "wp_response.png"))]

cells += [md(
r"""## Phase 5 — stopping power with the CAP

The classical Gaussian-electron ion decelerates crossing the slab; its kinetic-energy
loss across the faces gives the headline number. The CAP absorbs the projectile in
vacuum while the bath stays intact. **Full detail:**
[`03_cap_stopping/cap_stopping_study.ipynb`](03_cap_stopping/cap_stopping_study.ipynb)
and the [classical deep-dive](03_cap_stopping/p5cl_run_notebook.ipynb).""")]
cells += [embed(os.path.join(HYP, "03_cap_stopping", "classical_stopping.png"))]
cells += [code(
"# Headline S from the FULL-SUITE classical run (electron_track.csv, corrected KE).\n"
"c = np.genfromtxt(os.path.join(SYS,'scripts/fullsuite_classical/results/p5_classical/'\n"
"                  'raw/observables/electron_track.csv'), delimiter=',', names=True)\n"
"z, ke = c['z'], c['ke_ion_ha']; HALF=12.5\n"
"S = (np.interp(-HALF,z,ke) - np.interp(HALF,z,ke)) / (2*HALF)\n"
"print(f'HEADLINE  S = {S:.5f} Ha/Bohr = {S*27.21138625:.3f} eV/Bohr   (PROVISIONAL)')\n"
"print(f'work-energy check: -int(fz dz) over crossing matches dKE (see deep-dive)')")]

# 4. Bottom line --------------------------------------------------------------
cells += [md(
r"""## Bottom line

- **Localised jellium works in INQ** as a wrapper-only static perturbation — no
  `inq/` edit — and reproduces a bulk HEG interior with a clean surface (Phases 1–2).
- The **projectile + CAP machinery** is validated: closed-box energy conserved
  (Phase 3); the CAP absorbs the projectile in vacuum **without draining the bath**
  (Phase 5, Run A).
- **Headline:** $S \approx 0.71$ eV/Bohr for a 100 eV Gaussian electron through
  $r_s\approx4$ jellium (Phase 5, Run B) — **PROVISIONAL**.
- **Full-suite re-run delivered** the three-way **total/bath/WP** density
  decomposition (exact, via saved `density_wp`), momentum, GS-overlap excitations,
  and per-component energetics — see the three per-run deep-dives in the notebook
  map above. The classical stopping is corroborated by the work-energy theorem
  ($-\int f_z\,dz = \Delta\mathrm{KE}$).

### Open / next (not yet done)
- **T2 (Lang–Kohn)** absolute work-function & surface-energy benchmark +
  $dx{\times}\tfrac12$ grid convergence.
- **Convergence:** $dx\,0.5\to0.25$ and a multi-trajectory $S(v)$ sweep for a
  publishable stopping curve.
- **Loss function $L(q,\omega)$** — gated behind the Fourier-analysis methodology
  task (resolution caveat: $\Delta\omega\approx9$ eV at $\tau\approx18$ a.u.).
""")]

build(cells, OUT)
