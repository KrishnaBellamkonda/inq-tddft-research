#!/usr/bin/env python3
"""Assemble the investigation INDEX notebook — the one-glance tabulation.

Full cross-run table (investigation_summary.csv) + links to every per-run
run_report.ipynb and the 4 phase notebooks, grouped by plan phase. Executed to 0
errors so the table renders on open.
"""
from __future__ import annotations
from pathlib import Path
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbclient import NotebookClient

HERE = Path(__file__).resolve().parent
RES = ("/local/data/public/skcb2/tddft/ResearchProject/systems/vacuum/"
       "scripts/wp_traversal_energy/results")

INTRO = """# CAP energy-normalization investigation — index

**Hypothesis** (docs/plans/cap-energy-normalization-validation.md,
docs/notes/inq-energy-normalization-error.md): the CAP "energy shoots up" because
INQ reports the per-particle (norm-divided) kinetic energy (energy.hpp:50-55). The
extensive energy **E_ext = E_reported · norm** decays smoothly with the absorbed
norm. All runs: cheap vacuum WP sims on GPU 1 (inq-study), σ=3, E=400 eV, one-sided
CAP z∈[7.5, 22.5].

**Headline verdict:** confirmed for the CAP case (Phases 1, 2, 6 — E_ext/E0 == norm
to 3 d.p.). The mask+CN "surprise" (Phase 3) is a Crank-Nicolson renormalization
confound, not a refutation: E_ext rises there too (+4.8%), so CN pumps real energy
rather than cleanly preserving norm.
"""

PHASE_NB = """## Phase notebooks (study — run-SET verdicts)
| phase | notebook | question |
|---|---|---|
| 0 + 6 | [phase0_baseline_and_proof](phase0_baseline_and_proof.ipynb) | the phenomenon + INQ prints /norm |
| 1a | [phase1_geometry_independence](phase1_geometry_independence.ipynb) | rise is absorption, not reflection |
| 2 | [phase2_partial_absorption](phase2_partial_absorption.ipynb) | E_ext/E0 = norm identity line |
| 3 | [phase3_decisive_mask](phase3_decisive_mask.ipynb) | mask ETRS vs CN — the confound resolved |
"""

TABLE = f"""import pandas as pd
from pathlib import Path
RES = Path("{RES}")
df = pd.read_csv(RES/"investigation_summary.csv")
df["run_notebook"] = df.run.map(
    lambda r: f"results/{{r}}/report/run_report.ipynb"
    if (RES/r/"report"/"run_report.ipynb").exists() else "(none)")
cols = ["phase","run","param","dE_rep_eV","dE_rep_pct","norm_T","E_ext_frac","n_steps"]
df[cols].sort_values(["phase","run"]).round(3).reset_index(drop=True)
"""

LINKS = """## Per-run deep-dive notebooks
Each run's full battery (density GIFs, momentum, energetics) is at
`results/<run>/report/run_report.ipynb` under the dispatcher folder
`ResearchProject/systems/vacuum/scripts/wp_traversal_energy/`.
The `run_notebook` column above lists the relative path for each run.
"""


def main():
    nb = new_notebook()
    cells = [
        new_markdown_cell(INTRO),
        new_markdown_cell(PHASE_NB),
        new_markdown_cell("## Master cross-run table"),
        new_code_cell(TABLE),
        new_markdown_cell(LINKS),
    ]
    for c in cells:
        c.metadata["gen"] = "builder"
    nb.cells = cells
    NotebookClient(nb, timeout=180, kernel_name="python3").execute()
    out = HERE / "index.ipynb"
    nbf.write(nb, out)
    print(f"[index] wrote {out.name}")


if __name__ == "__main__":
    main()
