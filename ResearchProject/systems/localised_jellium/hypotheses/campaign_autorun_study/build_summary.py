#!/usr/bin/env python3
"""Build the campaign SUMMARY notebook: ideas, conclusions, problems, key plots,
and the expert-panel verdict. Re-runnable. (notebook-making house style.)

Run: PYTHONPATH=.../inq-stack/python <venv>/python3 build_summary.py
"""
from __future__ import annotations
import base64
from pathlib import Path
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell

LJ = "/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium"
CA = f"{LJ}/scripts/campaign_autorun"
OUT = Path(f"{LJ}/hypotheses/campaign_autorun_study"); OUT.mkdir(parents=True, exist_ok=True)
PNG = {"H0": "runs/h0/H0_base_difference.png", "H1": "runs/h1/H1_edge_model.png",
       "H2": "runs/h2/H2_gs_convergence.png", "H3": "runs/h3/H3_surface_energetics.png",
       "H4": "runs/h4/H4_wp_energetics.png", "H5": "runs/h5/H5_classical_subtraction.png"}

def md(s): return new_markdown_cell(s)
def img(key, cap):
    p = Path(CA) / PNG[key]
    if not p.exists(): return md(f"*(plot pending: {p.name})*")
    b = base64.b64encode(p.read_bytes()).decode()
    c = md(f"**{cap}**\n\n![{p.name}](attachment:{p.name})"); c.attachments = {p.name: {"image/png": b}}
    return c

VERDICT = r"""## Expert-panel verdict (provisional — user owns it)

*Convened via the `scientific-panel` skill: 4 independent opus/high experts
(TDDFT methodologist · condensed-matter generalist · jellium-surface specialist ·
run-data custodian) → openings → rebuttal → judge. They actively corrected each
other (E4's recompute killed E1's "divergent-drift" reading; E2's positivity bound
is independent of the slope argument).*

**Reduced-model ladder.** Isolated σ_WP=0.5 Gaussian electron (density std s=0.35):
Perdew–Zunger SIE = E_H[n]+E_xc[n]; self-Hartree 1/(2s√π)=22 eV, LDA-xc cancels ~80%
→ **net ≈ 4.4 eV, positive-definite, BC/box-independent.** Triangulates with
far-launch (4.55 eV) and PBC-at-r=40 (4.26 eV).

**Consensus.** **Periodicity-3 (PBC) is the physical branch; the open-z SIE must be
discarded for the charged WP.** A negative one-electron SIE is excluded by positivity
of E_H+E_xc. Measured excess = E_SIE (const) + E_cross(r) (real; BC-robust slope
≈ −0.046 eV/Bohr in both branches) + a BC-dependent charged-cell offset (~6 eV). The
offset is the net-charge × G=0 bookkeeping: PBC drops G=0 (jellium trick; residual
~0.8 eV in-plane Madelung); periodicity-2 sets G=0 = 0.5·rc² (rc=L_z, poisson.hpp:49)
— a monopole self-energy on the charged WP only, absent from the neutral GS. Absolute
GS energies (−161/−109/+60 Ha) are a **box/BC gauge, not physics** (interior n₀, N
correct).

**Best answer.** Use PBC; hand the stopping campaign **E_SIE ≈ 4.4 eV** (2 s.f.;
4.4 / 4.55 / 4.26 three-way), as a fixed-geometry difference, caveated that r→∞ may
drop ~0.3–0.8 eV. Discard open-z SIE for the charged WP. The open-z "more-negative
with r" trend is **not physics** — real E_cross(r) on a spurious offset. *Inference:*
the neutral GS under open-z (+60 Ha, L_z-flat) is itself a clean reference; open-z is
wrong only for the charged difference. **Branch verdict solid; number provisional.**

**Decisive next test (cheap, not yet run).** Isolated σ_WP=0.5 electron in an EMPTY
box, {periodicity 3,2} × {L_z=90,120}: PBC → ≈4.4 eV & L_z-independent (confirms SIE);
open-z → ∝ L_z² (proves the monopole). + in-plane face scan (L_xy=50/70/100) for the
Madelung tail.

**Open flags / unresolved.** (a) 3× monopole magnitude gap (envelope ~2 eV vs observed
~6 eV; unverified → INQ FFT normalisation). (b) Is 4.3 eV converged (r=40 still
drifting). (c) classical ghost omits ∫v_ghost·n₊ → WP-vs-classical comparison not yet
apples-to-apples.
"""

cells = [
    md("# Localised-jellium GS parameter study — CAMPAIGN SUMMARY\n\n"
       "*Ideas, conclusions, problems, key plots, and the expert-panel verdict for the "
       "autonomous ladder H0→H5. Built from the dense re-run data. Companion to the "
       "per-phase notebooks (`H0..H5_*.ipynb`) and `campaign_autorun_study.ipynb`.*"),
    md("## What this campaign did\n"
       "A localised jellium **slab** (50×50 face, 25 Bohr thick, N=82, r_s≈5.7, LDA) probed as an "
       "ordered ladder of falsifiable hypotheses, to pin the **energy reference + self-interaction "
       "error (E_SIE)** for a Gaussian projectile and to build analytical intuition. Each phase ran "
       "headless (Python orchestrator, idempotent), emailed a result, and has its own notebook.\n\n"
       "| Phase | Question | Key result |\n|---|---|---|\n"
       "| **H0** | Is the raw WP−classical t=0 gap the localisation energy? | **No** — WP excess ~stable ≈86 eV; classical ghost excess +188→+11 eV; raw gap −100→+75 eV, artifact-dominated |\n"
       "| **H1** | Does a finite edge width kill Gibbs while keeping Friedel? | edge ringing falls as w grows past the grid; interior n₀≈1.3e-3 ∀w |\n"
       "| **H2** | Is the interior box-independent; is open-z usable? | interior n₀ flat (L_z=50–150); open-z GS converges; **absolute E_GS is a box/BC gauge** |\n"
       "| **H3** | Liquid-drop σ_s + bulk limit? | E(N) **non-linear** (E_self-confounded) → σ_s RAW only; thin-slab n₀ rises (loses bulk) |\n"
       "| **H4** | E_SIE + PBC-vs-open-z verdict | **PBC E_SIE → 4.3 eV plateau (≈ known 4.5)**; open-z negative/drifting (artifact) |\n"
       "| **H5** | classical mirror + route-2 | classical excess strongly r-dependent; route-2 ghost-bg integral still needed |\n"),
    md("## Key plots"),
    img("H0", "H0 — base WP-vs-classical gap vs r (artifact-dominated, not localisation)"),
    img("H4", "H4 — WP E_SIE vs r: PBC plateau ≈4.3 eV (physical); open-z negative (charged-cell artifact)"),
    img("H2", "H2 — interior n₀ vs L_z (box-converged); open-z GS usable"),
    img("H3", "H3 — E(N) across thickness (non-linear: E_self confound — σ_s RAW)"),
    img("H1", "H1 — n(z) vs edge width (Gibbs vs Friedel)"),
    img("H5", "H5 — classical ghost excess vs r (PBC vs open-z)"),
    md(VERDICT),
    md("## Conclusions (campaign-level)\n"
       "1. **Energy reference for Campaign 1: E_SIE ≈ 4.4 eV (PBC), provisional.** Use fixed-geometry, "
       "fixed-BC differences only — absolute energies are a gauge.\n"
       "2. **Open-z is the right idea but needs a charged-cell correction.** It is clean for the neutral "
       "GS; the net-charged WP picks up the 0.5·rc² monopole term → discard open-z SIE until corrected.\n"
       "3. **The WP's base energy above GS is ≈ localisation (82 eV) + SIE (4.4 eV)**, distance-stable; "
       "the classical ghost's apparent gap is dominated by the omitted background term, not physics.\n"),
    md("## Problems / open follow-ups\n"
       "- **Open-z net-charge G=0 reference** (the 0.5·rc² monopole) — the headline methodology fix.\n"
       "- **3× monopole magnitude gap** (~2 eV envelope vs ~6 eV observed) — likely INQ FFT normalisation; unverified.\n"
       "- **H2 work function Φ** extractor; **H3 σ_s** E_self correction; **H5 ghost-background integral** — pre-gates still flagged.\n"
       "- **Convergence of the 4.3 eV plateau** — r=40 WP is only 7.5 Bohr from the box edge.\n"),
    md("## Decisive next step (panel-recommended, cheap)\n"
       "Run the **slab-free isolated σ_WP=0.5 Gaussian electron in an empty box**, "
       "{periodicity 3, 2} × {L_z=90, 120} (4 configs, seconds each) + an in-plane face scan. "
       "PBC→≈4.4 eV & L_z-flat confirms the SIE; open-z→∝L_z² proves the monopole and quantifies the "
       "correction. This locks the number before it is handed to the classical-projectile campaign.\n\n"
       "*Open questions for the user are listed in the panel verdict above (§7).*"),
]

n = new_notebook(); n.cells = cells
n.metadata.kernelspec = {"name": "python3", "display_name": "Python 3"}
p = OUT / "campaign_summary.ipynb"; nbf.write(n, str(p)); print("wrote", p)
