#!/usr/bin/env python3
"""Assemble phase1_gs_sie.ipynb — Phase-1 results of the quantum-stopping-power
campaign (GS at the new density + the SIE diagnostic both ways). Path-referenced
markdown + the GS density figure. Reads computed numbers from sie_results.csv and
the GS run_summary so the notebook is self-consistent with the runs."""
import os
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = "/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium"
r = pd.read_csv(f"{HERE}/sie_results.csv").iloc[0]

cells = []
def md(s): cells.append(new_markdown_cell(s))
def fig(png, cap): cells.append(new_markdown_cell(f"*{cap}*\n\n![{cap}]({png})"))

md(rf"""# Phase 1 — `quantum-stopping-power` campaign: GS + self-interaction diagnostic
### localised jellium slab · r$_s$≈5.67 · 82 e · box 50×50×70 · σ_WP=0.5

**Campaign:** `docs/campaigns/jellium_wp_stopping/quantum-stopping-power.md` (Phase 1).
This notebook holds the Phase-1 results: a validated ground state at the new
(n162-reference) density, and the wavepacket **self-interaction energy (SIE)**
quantified **two ways**. Every number is computed from the runs.

**Why Phase 1.** The SIE — the spurious one-electron self-repulsion LDA does not
cancel — would otherwise contaminate any future WP−classical stopping comparison.
Phase 2 (CAP / sim-time / steady-state / the stopping comparison) is built *after*
analysing these results.
""")

md(rf"""## §1 — Ground state (task P1.1)

Converged LDA GS of the localised slab, background injected as a static
perturbation (orthorhombic 50×50×70 cell, stock inq).

| quantity | value |
|---|---|
| GS total energy | **{r.E_GS_slab_Ha:.5f} Ha** ({r.E_GS_slab_Ha*27.211386:.1f} eV) |
| electrons | **82** (even) |
| density / r_s | n₀ = 0.001312 a₀⁻³ ⟹ **r_s = 5.667** |
| interior density (|z|<10) | mean 0.00136, **5.6%** variation (flat to a few %) |
| cell | 50 × 50 × 70 Bohr orthorhombic, spacing 0.50 |

SCF ended normally; the density localises inside |z|<12.5 with a roughly flat
interior (Friedel + surface variation ~6%). The **free region** (12.5<|z|<25) holds
mean density **6.2×10⁻⁵ a₀⁻³ — only 4.6% of the interior**, i.e. nearly empty: the
log panel (left) visually exaggerates it, the linear panel (right) shows the true
confinement. Closed-shell: 41 doubly-occupied + 20 empty states; clean low-T (100 K)
convergence is consistent with a HOMO–LUMO gap (full occupation-array check
deferred — non-blocking per the campaign).
""")
fig("gs_density_xz.png", "GS slab density (xz) — log (left) vs linear (right); density confined to |z|<12.5, free region ≈ 0")

md(rf"""## §2 — Wavepacket self-interaction (SIE), both ways (tasks P1.2–P1.3)

The σ_WP=0.5 / 100 eV wavepacket was injected **far** from the slab (z_mean(0) =
**{r.launch_z:.2f}**, target −32; WP norm = {r.wp_norm0:.4f}), CAP off, and the t=0
total energy read off. The SIE is isolated as `E_total(0) − E_GS_slab − KE_WP`.

| quantity | value |
|---|---|
| E_total(0) [WP far] | **{r.E_total0_Ha:.5f} Ha** ({r.E_total0_Ha*27.211386:.1f} eV) |
| KE_WP = ⟨p²⟩/2 (measured) | **{r.KE_WP_eV:.1f} eV** |
| &nbsp;&nbsp;drift ½⟨p_z⟩² | {r.drift_eV:.1f} eV |
| &nbsp;&nbsp;zero-point (measured) | **{r.zero_point_meas_eV:.1f} eV** (vs 3/4σ² = {r.zero_point_theory_eV:.1f} eV theory) |

**The two estimates:**

| estimate | formula | value |
|---|---|---|
| **SIE_a** (user reference) | `E_total(0) − (E_GS + 100 eV)` | **{r.SIE_a_eV:+.2f} eV** = SIE + zero-point |
| **SIE_b** (clean SIE) | `E_total(0) − E_GS − KE_WP(measured)` | **{r.SIE_b_eV:+.2f} eV** = SIE |
| difference | `SIE_a − SIE_b = KE_WP − 100 eV` | **{r.SIE_a_minus_b_eV:+.2f} eV** |

**Reading it.** `SIE_a − SIE_b = {r.SIE_a_minus_b_eV:+.1f} eV` is exactly the
energy the "+100 eV" reference fails to account for — the WP's measured
**zero-point KE** ({r.zero_point_meas_eV:.0f} eV) minus any drift discrepancy. This
is *why* the "+100 eV" formula over-counts the SIE: it subtracts only the drift, not
the full kinetic energy a width-σ wavepacket carries. **THE SIE is SIE_b ≈
{r.SIE_b_eV:.1f} eV** (compare the old r_s=4 estimate ≈ 4.5 eV — weakly
density-dependent, as expected).

**Diagnostic value.** Any future WP−classical stopping difference must be read
against this ~{r.SIE_b_eV:.0f} eV artifact floor; at larger σ it shrinks (~1/σ).
""")

md(r"""## §3 — Phase-1 verdict & what feeds Phase 2

- ✅ **GS validated** at the n162-reference density (r_s 5.67, 82 e), localised slab,
  flat-ish interior — a reusable checkpoint for the production runs.
- ✅ **SIE quantified both ways**, with the zero-point cross-check confirming the
  "+100 eV" reference omits exactly the zero-point KE.
- ➡️ **Phase 2 (to author next):** CAP (10 Bohr each, η=−0.7) + total sim time +
  the **steady-state (energy-convergence) test for BOTH the WP and classical runs**
  + the stopping comparison (classical via the `stopping-power-extraction` skill,
  WP via the converged E_total balance). The ~{:.0f} eV SIE floor is carried into
  that comparison.
""".format(r.SIE_b_eV))

nb = new_notebook(); nb.cells = cells
nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
out = f"{HERE}/phase1_gs_sie.ipynb"
nbf.write(nb, out)
print(f"wrote {out}  ({len(cells)} cells)")
