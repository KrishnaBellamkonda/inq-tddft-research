#!/usr/bin/env python3
"""Assemble + execute the 4 phase study notebooks from the figures + summary CSV.

Run build_phase_notebooks.py first (produces figures/*.png). This assembles
nbformat notebooks (context -> runs -> figure -> table -> verdict) and executes
them to 0 errors, embedding outputs. Figures are path-referenced (KB-sized nb).
"""
from __future__ import annotations
from pathlib import Path
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbclient import NotebookClient

HERE = Path(__file__).resolve().parent
RESULTS = ("/local/data/public/skcb2/tddft/ResearchProject/systems/vacuum/"
           "scripts/wp_traversal_energy/results")

CONV = """## Conventions & symbols
| symbol | meaning |
|---|---|
| σ | wavepacket width σ_wp = 3 Bohr (charge std = σ/√2) |
| k0, E | launch momentum 5.421 au, E = ½k0² = 400 eV |
| E_reported | the kinetic energy INQ prints — the **per-particle mean** ⟨T⟩/⟨ψ\\|ψ⟩ (energy.hpp:50-55) |
| norm | physical orbital norm ⟨ψ\\|ψ⟩(t), from momentum-stats norm_check ratio |
| **E_ext** | **extensive** energy = E_reported · norm — the physically conserved quantity |
| η, W | CAP strength (Ha) and width (Bohr); one-sided band z∈[7.5, 22.5] |

Engine: **inq-study** (only fork that compiles a CAP; stock inq hard-errors on
`double += complex`). Vacuum, non-interacting, single WP electron — every non-kinetic
KS term is 0, so E_reported == energies.csv:kinetic == wp_momentum_stats:e_kin_ha.
"""

TABLE_CELL = f"""import pandas as pd
df = pd.read_csv("{RESULTS}/investigation_summary.csv")
def show(phase):
    d = df[df.phase.astype(str)==str(phase)][
        ["run","param","dE_rep_eV","dE_rep_pct","norm_T","E_ext_frac","n_steps"]]
    return d.round(3).reset_index(drop=True)
"""

PHASES = {
"phase0_baseline_and_proof": dict(
  title="Phase 0 + 6 — the phenomenon and the /norm proof",
  question=(
    "**Question.** Does the CAP energy really 'shoot up', and if so is it the "
    "reported (per-particle) kinetic or the extensive energy that rises? "
    "Phase 0 shows the phenomenon; Phase 6 proves INQ prints the norm-divided quantity."),
  runs="`nocap`, `cap` (1-sided η=−3.5), `nocap_long` — baselines.",
  figs=["figures/phase0_phenomenon.png", "figures/phase6_norm_proof.png"],
  phase_key="0",
  verdict=(
    "**Verdict.** The *reported* energy stays pinned/rises while the CAP drains the "
    "norm; the *extensive* E_ext = E_reported·norm decays smoothly to 0. Phase 6 proves "
    "the mechanism directly: `energies.csv:kinetic` lies exactly on "
    "`wp_momentum_stats:e_kin_ha` (the per-particle mean), NOT on `e_kin_ha·norm`. "
    "INQ reports ⟨T⟩/⟨ψ|ψ⟩ (energy.hpp:50-55). The no-CAP control shows ΔE≈0, norm≡1."),
),
"__phase0_worked_example__": dict(  # injected into phase0 only (see build())
  md_intro=(
    "## Worked example — why the *reported* energy rose (real `cap` data)\n\n"
    "The number INQ writes each step is `E_reported = ⟨ψ|T|ψ⟩ / ⟨ψ|ψ⟩` — the "
    "**per-particle mean** kinetic (energy.hpp `occ_sum` divides by the norm). "
    "The *extensive* (physical) kinetic is the numerator alone, "
    "`E_ext = ⟨ψ|T|ψ⟩ = E_reported · norm`, which I reconstruct in post-processing "
    "(`E_reported` from `energies.csv:kinetic`, `norm` from `wp_momentum_stats:norm_check` "
    "normalized to its t=0 value). Watch the two right-hand columns below:"),
  code=(
    "import pandas as pd, numpy as np\n"
    "HA = 27.211386\n"
    "run = '" + RESULTS + "/cap'\n"
    "en  = pd.read_csv(f'{run}/raw/observables/energies.csv', comment='#')\n"
    "mom = pd.read_csv(f'{run}/raw/observables/wp_momentum_stats.csv', comment='#')\n"
    "nc  = mom['norm_check'].to_numpy(float); norm = nc / nc[0]\n"
    "t   = en['time_au'].to_numpy(float)\n"
    "Ek  = en['kinetic'].to_numpy(float)          # E_reported = <T>/norm (from sim), Ha\n"
    "ni  = np.interp(t, mom['time_au'].to_numpy(float), norm)   # physical norm on E-grid\n"
    "ext = Ek * ni                                # extensive kinetic <T> = E_rep*norm\n"
    "rows = [0,300,400,500,600,799]\n"
    "tab = pd.DataFrame({\n"
    "    't_au':        t[rows].round(2),\n"
    "    'norm':        ni[rows].round(4),\n"
    "    'E_reported_eV': (Ek[rows]*HA).round(1),   # what INQ prints (per-particle mean)\n"
    "    'E_ext_eV':      (ext[rows]*HA).round(2),  # physical energy = E_rep*norm\n"
    "}).reset_index(drop=True)\n"
    "print('BEFORE: t=0   E_reported =', round(Ek[0]*HA,1), 'eV   norm =', round(ni[0],4),\n"
    "      '  E_ext =', round(ext[0]*HA,1), 'eV')\n"
    "print('AFTER : t=8   E_reported =', round(Ek[-1]*HA,1), 'eV   norm =', round(ni[-1],4),\n"
    "      '  E_ext =', round(ext[-1]*HA,2), 'eV')\n"
    "print()\n"
    "print('reported \"rose\"  402.1 ->', round(Ek[-1]*HA,1),\n"
    "      'eV  =  E_ext/norm =', round(ext[-1]*HA,3), '/', round(ni[-1],4),\n"
    "      '=', round(ext[-1]/ni[-1]*HA,1), 'eV')\n"
    "print('physical E_ext   402.1 ->', round(ext[-1]*HA,2), 'eV  (99.99% absorbed by the CAP)')\n"
    "tab"),
  md_toy=(
    "**The arithmetic in one line.** At the end the CAP has removed 99.99% of the "
    "wavepacket: `E_ext = 0.03 eV`, `norm = 1e-4`. INQ reports the ratio "
    "`0.03 / 1e-4 = 416 eV` — the mean kinetic of the last high-k sliver stuck at the "
    "CAP edge. It looks like the energy *rose* +14 eV, but nothing gained energy.\n\n"
    "**Minimal toy (two plane waves).** Let the packet be a slow mode (300 eV) and a "
    "fast mode (500 eV), each with weight 0.5 → mean = 400 eV, norm = 1. A slow mode "
    "dwells longer in the CAP, so it is absorbed first. Remove 98% of the slow weight "
    "and 60% of the fast: extensive `= 0.01·300 + 0.20·500 = 103 eV` (fell), "
    "norm `= 0.21` (fell more), reported `= 103/0.21 = 490 eV` (**rose**). Same "
    "mechanism: a ratio whose denominator collapses faster than its numerator.\n\n"
    "**Conclusion.** `E_reported` is a per-particle mean — a ratio — so it is only "
    "meaningful while `norm ≈ 1`. Under a CAP the physically correct energy is "
    "`E_ext = E_reported · norm`, which decays smoothly to zero as the packet is "
    "absorbed. The apparent 'total energy rise' is entirely the normalization artifact."),
  md_fix=(
    "## Post-processing fix — and when `× norm` is valid\n\n"
    "The ÷norm bias is analytic, so it is invertible **after the run, without touching "
    "the source**. Everything needed is already saved: `E_reported` (`energies.csv`), the "
    "WP orbital per-particle kinetic `T_WP` (`wp_momentum_stats.csv:e_kin_ha`), and the WP "
    "norm (`norm_check` ratio). But only the **kinetic** term is contaminated "
    "(energy.hpp:83); hartree / external / xc come from n(r) and are already extensive "
    "(self_consistency.hpp:191,196,245). So the correct general fix touches only the "
    "kinetic:\n\n"
    "> **General (correct):**  `E_corr = E_reported − T_WP · (1 − norm_WP)`\n\n"
    "This subtracts the inflated WP kinetic and restores the extensive one, leaving the "
    "density terms untouched.\n\n"
    "**Vacuum special case (this notebook).** Here the WP is the *only* electron and every "
    "non-kinetic term is 0, so `E_reported = T_WP`. Substituting collapses the fix to\n\n"
    "> **Vacuum:**  `E_corr = E_reported − E_reported·(1 − norm) = E_reported · norm`\n\n"
    "That is exactly the dashed `E_ext = E_reported · norm` curve plotted above — **exact**, "
    "not an approximation, because in vacuum total energy = kinetic.\n\n"
    "**Do NOT reuse `× norm` for the jellium WP run.** There the total is dominated by the "
    "extensive hartree/external/xc terms (which are already correct); a global `× norm` "
    "would wrongly shrink them too. Use the kinetic-only form "
    "`E_corr = E_reported − T_WP·(1 − norm_WP)` — this is what "
    "`wp_cap_energy_plateau/wp_kinetic_normalization_fix.py` applies (the 93.5 → 115.9 eV "
    "jellium plateau correction)."),
),
"phase1_geometry_independence": dict(
  title="Phase 1a — geometry independence (η-sweep): rise is absorption, not reflection",
  question=(
    "**Question.** Falsifier (b): if the residual scales with η it could be CAP "
    "reflection. Does the reported rise track η — or the absorbed fraction?"),
  runs="`exp1a_eta-{0.3, 0.7, 1.0, 2.0, 3.5}` — one-sided CAP, only η varies.",
  figs=["figures/phase1_eta_sweep.png"],
  phase_key="1",
  verdict=(
    "**Verdict.** ΔE_reported grows with |η| (1.3→14.0 eV) — but only because stronger "
    "CAPs absorb more (norm_T 0.44→0.00). Plotted against fraction-absorbed the points "
    "collapse to one curve, and E_ext/E0 == norm_T to 3 d.p. across the whole sweep. "
    "The rise is the normalization artifact of a shrinking norm, **not** reflection: "
    "falsifier (b) fails."),
),
"phase2_partial_absorption": dict(
  title="Phase 2 — partial-absorption ladder: E_ext/E0 is the identity line in norm",
  question=(
    "**Question.** The core prediction: E_ext = E_reported·norm should equal the "
    "surviving fraction. Tune weak CAPs so norm ends ~0.7/0.5/0.3 and test E_ext/E0 vs norm."),
  runs="`exp2_N0p{1,3,5}` (weak CAPs) + the Phase-1 η-sweep as extra points.",
  figs=["figures/phase2_identity_line.png"],
  phase_key="2",
  verdict=(
    "**Verdict.** Every run lands on y = x: E_ext/E0 = 0.708/0.549/0.317 at "
    "norm = 0.707/0.548/0.316. The extensive energy is exactly the absorbed norm — "
    "the cleanest confirmation of the hypothesis."),
),
"phase3_decisive_mask": dict(
  title="Phase 3 (DECISIVE) — mask ETRS vs CN, and the Crank-Nicolson confound",
  question=(
    "**Question.** Falsifier (c): a NORM-PRESERVING absorber that still shows the rise "
    "would kill the hypothesis. Same sin² spatial clip, toggle only the propagator: "
    "ETRS keeps the removal (norm-losing); Crank-Nicolson renormalizes each step "
    "(norm-preserving). Does the rise survive with norm held?"),
  runs="`exp3a_mask_etrs` (norm→0) vs `exp3b_mask_cn` (norm≡1).",
  figs=["figures/phase3_decisive.png"],
  phase_key="3",
  verdict=(
    "**Verdict — the surprise, resolved.** Both show ~+19 eV in E_reported, so on the "
    "*reported* energy alone falsifier (c) appears to fire. But E_ext separates them: "
    "ETRS → E_ext/E0 = 0.000 (pure normalization artifact, energy correctly absorbed); "
    "CN → E_ext/E0 = **1.048** — the extensive energy genuinely *rose* 4.8%. "
    "Crank-Nicolson does not merely preserve norm; it **renormalizes** the sin²-clipped "
    "orbital every step, and clip-then-renormalize injects real high-k content. So "
    "mask+CN is not a clean norm-preserving control — it changes the physics (energy "
    "pumping), a different mechanism from the CAP artifact. The real CAP runs use ETRS + "
    "`perturbations::absorbing` (norm-losing), where E_ext = E_reported·norm is the "
    "correct extensive energy. **Hypothesis stands; the red flag was a CN confound.**"),
),
}


def build(stem: str, spec: dict):
    nb = new_notebook()
    cells = [
        new_markdown_cell(f"# {spec['title']}\n\n{spec['question']}"),
        new_markdown_cell(CONV),
        new_markdown_cell(f"## What was run\n{spec['runs']}\n\n"
                          "See each run's `results/<run>/report/run_report.ipynb` for the "
                          "per-run deep dive (density GIFs, momentum, energetics)."),
    ]
    for f in spec["figs"]:
        cells.append(new_markdown_cell(f"![{Path(f).stem}]({f})"))
    if stem == "phase0_baseline_and_proof":
        wx = PHASES["__phase0_worked_example__"]
        cells.append(new_markdown_cell(wx["md_intro"]))
        cells.append(new_code_cell(wx["code"]))
        cells.append(new_markdown_cell(wx["md_toy"]))
        cells.append(new_markdown_cell(wx["md_fix"]))
    cells.append(new_markdown_cell("## Cross-run table"))
    cells.append(new_code_cell(TABLE_CELL + f"\nshow('{spec['phase_key']}')"))
    cells.append(new_markdown_cell(spec["verdict"]))
    # tag builder cells so reader annotations round-trip future rebuilds
    for c in cells:
        c.metadata["gen"] = "builder"
    nb.cells = cells
    out = HERE / f"{stem}.ipynb"
    NotebookClient(nb, timeout=300, kernel_name="python3").execute()
    nbf.write(nb, out)
    print(f"[nb] wrote {out.name}")


if __name__ == "__main__":
    for stem, spec in PHASES.items():
        if stem.startswith("__"):
            continue  # pseudo-entries (injected sub-sections), not standalone notebooks
        build(stem, spec)
    print("[nb] all phase notebooks built")
