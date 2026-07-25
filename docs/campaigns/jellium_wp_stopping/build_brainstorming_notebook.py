#!/usr/bin/env python3
"""Assemble brainstorming-jellium-campaigns.ipynb — the record of the 2026-06-24
campaign brainstorming session (jellium WP-vs-classical stopping).

Path-referenced markdown + the one analysis figure produced by
make_sigma_lindhard_comparison.py. Output: brainstorming-jellium-campaigns.ipynb.
"""
import os
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell

HERE = os.path.dirname(os.path.abspath(__file__))
cells = []
def md(s): cells.append(new_markdown_cell(s))
def fig(png, cap): cells.append(new_markdown_cell(f"*{cap}*\n\n![{cap}]({png})"))

md(r"""# Brainstorming — jellium wavepacket-vs-classical stopping campaigns
### 2026-06-24 session · localised jellium slab (r$_s$≈4) · WP vs classical projectile

**Purpose.** Working surface for designing the next campaigns after the
`03_cap_stopping` baselines. Holds the data points/plots requested in the
brainstorm. Companion docs (same folder): `draft_campaigns.md` (rough draft of the
three campaigns), `notes_campaign1_sigma05_restrictions.md` (Campaign-1 known
limits). Literature anchor: `docs/sources/nazarov-gross-2025-quantum-projectile-stopping.md`
(Nazarov & Gross 2025 — quantum projectile stopping via Exact Factorization).

**Three campaigns (see draft for detail):**
1. **`quantum-stopping-power`** — small σ_WP=0.5 (point-like), quantum stopping vs matched classical vs Lindhard. *Formalise next via `/campaigns`.*
2. **large rigid σ** — σ_WP≈4 @ E≥300 eV, isolate quantum-vs-classical at fixed matched σ (no spreading).
3. **muon (future)** — needs an `inq-study` per-orbital-mass engine fork; classical-muon trivial.
""")

md(r"""## 1 — The σ decision for Campaign 1: is σ_WP=1.0 "good enough" vs point-Lindhard?

**Question (from the session).** σ_WP=0.5 was chosen because the classical projectile
converges to point-charge linear response (Lindhard). Would the *less strongly
spreading* σ_WP=1.0 still track Lindhard well enough to use instead, easing the
spreading problem?

**σ convention (the √2 trap — made explicit).** The localised campaign labels the
**wavepacket width** σ_WP; the classical Gaussian charge std is σ_pot = σ_WP/√2.
(The bulk jellium S(v) campaign instead labelled σ = σ_pot — so its "sigma0p5"
means σ_WP=0.707. Watch this when reusing bulk data.)

**Finding (validated, and it flips the expected approach):**
- ✅ **Point-charge Lindhard is trustworthy** — analytical `stopping_power_point`
  gives **0.717 eV/Bohr** at r_s=4, v=2.71, matching the **0.719** the baseline used.
- ✅ **σ_WP=0.5 classical ≈ point-Lindhard** — the localised baseline classical run
  gave **S = 0.706 eV/Bohr = 0.98× point**. Your premise is correct.
- ⚠️ **The analytical finite-σ Lindhard (`stopping_power_sigma`) is NOT reliable** —
  it predicts σ_WP=0.5 should be only **0.77× point**, contradicting the 0.98× run.
  It **over-suppresses**, so it **cannot** be used to judge σ_WP=1.0.
- ➡️ **Verdict: the σ=1 vs σ=0.5 choice is _undecidable from existing data_** — it
  needs one cheap σ_WP=1.0 classical S(v) run vs point-Lindhard.
- ✅ **Decision (2026-06-24): σ_WP = 0.5 LOCKED for the first test run.** σ=0.5 is
  the validated point-faithful choice; we accept its spreading limits (recorded in
  §2). We **revisit raising σ only if the first run's results motivate it** — the
  σ_WP=1.0 vs point-Lindhard check is **deferred** until then, not a prerequisite.
""")
fig("sigma_lindhard_comparison.png",
    "Classical stopping vs point-Lindhard — σ_WP=0.5 run sits on the POINT curve (0.98×); "
    "analytical finite-σ over-suppresses (dashed), so σ_WP=1.0 needs a dedicated run")

md(r"""## 2 — Campaign 1 known restrictions (recorded, NOT blockers)

Full detail in `notes_campaign1_sigma05_restrictions.md`. Summary: at σ_WP=0.5 the
E_total ledger is hard because (i) **72× free-Gaussian spreading**, (ii) the
**centroid stalls** (Δz for dE/dx undefined), (iii) **no-wrap vs full-absorption
are incompatible** in a 50-Bohr box, (iv) the ledger is contaminated by **zero-point
KE ≈ 82 eV** and **SIE ≈ 4.5 eV**. **Resolutions:** smoke-test a **force/work-integral
estimator** (local to the projectile, survives spreading) before committing; else
bigger box + longer run; always subtract zero-point and bound SIE (vacuum-WP control).
User decision: **proceed with Campaign 1**, carrying these as documented limits.
""")

md(r"""## 3 — Literature anchor (Nazarov & Gross 2025)

The premise is sound and current: a quantum projectile's **width/mass** changes the
stopping — a purely quantum effect; the classical (M→∞) limit recovers Lindhard.
They solve the WP-energy-partition problem with **Exact Factorization** (which INQ
does **not** implement — hence our reliance on the total-energy ledger / a force
estimator). Matching WP σ to the classical-potential σ is the right comparison.
See the source note for the equations and how it maps to our runs.
""")

nb = new_notebook()
nb.cells = cells
nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
out = os.path.join(HERE, "brainstorming-jellium-campaigns.ipynb")
nbf.write(nb, out)
print(f"wrote {out}  ({len(cells)} cells)")
