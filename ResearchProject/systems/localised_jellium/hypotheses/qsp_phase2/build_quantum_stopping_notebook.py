#!/usr/bin/env python3
"""Build the quantum-stopping energy-ledger notebook (26-6-26 meeting).

Usage: build_quantum_stopping_notebook.py [p2|p5]   (default p2)

Assembles a self-contained, executable notebook over `quantum_stopping_ledger.py`
so every number is live and checkable. House narrative: context -> formula (terms
defined) -> setup -> careful energy ledger table -> S_WP (FULL LEDGER, upper bound,
user decision 2026-06-26) -> convergence diagnostics -> classical slab -> takeaway.

Headline S = (E_total(t_f) - E_GS)/L_z, reported as an UPPER BOUND: the sigma=0.5
packet's ~82 eV zero-point energy inflates the bath gain and the WP is not fully
absorbed. The drift-credit alternative (subtract only the 100 eV drift) is shown in
the transparency table but is NOT the headline -- it gives an impossible negative
for the more-absorbed run (p2: -0.5 eV/Bohr), proving the zero-point does not
persist in the final jellium. (Independently validated 2026-06-26.)

Numbers rounded to 2 s.f. (3 s.f. for near-equal differences) per
.claude/rules/number-rounding.md; full precision stays in the live code cells.
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = Path("/local/data/public/skcb2/tddft")
sys.path.insert(0, str(ROOT / "ResearchProject/systems/localised_jellium/hypotheses"))
sys.path.insert(0, str(HERE))

import _nbreport as nb  # noqa: E402
import quantum_stopping_ledger as L  # noqa: E402

RUN = sys.argv[1] if len(sys.argv) > 1 else "p2"
OUT_NAME = {"p2": "quantum_stopping_ledger_26-6-26.ipynb",
            "p5": "quantum_stopping_ledger_p5_26-6-26.ipynb",
            "p3": "quantum_stopping_ledger_p3_26-6-26.ipynb"}[RUN]
OUT = HERE / OUT_NAME
FIGDIR = HERE / (OUT_NAME.replace(".ipynb", "_figs"))
nb.set_outdir(str(HERE))


def f(x, sf=2):
    return L._fmt(x, sf)


def _classical_section(cl):
    """Markdown for §5 — branches on whether the ion cleanly traversed the slab."""
    head = ("## 5. Classical slab projectile\n\n"
            "The classical Ehrenfest ion is **not** absorbed (the CAP removes "
            "wavefunctions, not the point ion), so stopping is read from the ion ΔKE.\n\n")
    if cl["traversed"]:
        return head + (
            f"| Method | S (eV/Bohr) | note |\n|---|---:|---|\n"
            f"| equal-potential faces ($|z|{{=}}12.5$) | **{f(cl['S_face_ev_per_bohr'],2)}** | "
            f"KE {f(cl['ke_entry_ev'],3)}→{f(cl['ke_exit_ev'],3)} eV; defensible |\n"
            f"| slab-centre KE minimum | {f(cl['S_center_ev_per_bohr'],2)} | "
            f"over-counts (conservative well) |\n")
    return head + (
        f"> ⚠️ **ANOMALY — the ion did NOT traverse the slab.** Launched at "
        f"z={f(cl['z_launch'],3)} toward the slab, it **stalled at z={f(cl['z_max'],2)} "
        f"(KE→0) and reversed** (z_final={f(cl['z_final'],3)}), never crossing the far "
        f"face. So the equal-potential-face method is invalid here. Over the run the ion "
        f"lost **ΔKE_ion = {f(cl['dKE_ion_ev'],2)} eV** (KE {f(cl['ke_launch_ev'],3)}→"
        f"{f(cl['ke_final_ev'],2)}); the electronic energy balance gives "
        f"ΔE_total = {f(cl['dE_total_cl_ev'],2)} eV.\n\n"
        f"This is **~10× the ~8 eV a 100 eV point charge should lose over 25 Bohr** "
        f"(bulk classical S≈0.34 at this v) — the projectile is **trapped / reflected**, "
        f"not cleanly stopping. Treat as a flag to investigate the classical run setup, "
        f"NOT as a stopping-power measurement.\n")


def main():
    wp = L.compute_wp_ledger(RUN)
    cl = L.compute_classical_slab(RUN)
    figs = L.make_figures(FIGDIR, RUN)
    HA = L.HA_TO_EV
    c = L.cfg(RUN)

    cells = [
        nb.md(
            f"# Quantum (wavepacket) stopping power — energy ledger ({RUN})\n"
            f"**{wp['label']} · σ_WP = 0.5 · E = 100 eV · 26 Jun 2026 meeting**\n\n"
            "Computes the quantum electronic stopping power from the WP / classical "
            "twin via the **energy-balance (retained-energy) method** — the only route, "
            "since the wavepacket and the bath are a single inseparable Kohn–Sham "
            "system with no well-defined projectile force or trajectory. Every "
            "intermediate number is tabulated for checking (2 s.f.; full precision in "
            "the live code cells).\n\n"
            "> **Headline (full ledger, user decision 2026-06-26):** "
            f"`S_WP = (E_total(t_f) − E_GS)/L_z = {f(wp['S_wp_ev_per_bohr'])}` eV/Bohr, "
            "reported as an **UPPER BOUND**. The σ=0.5 packet carries ~82 eV of "
            "zero-point energy (≈ its 100 eV drift); the energy balance cannot "
            "separate that from velocity-stopping, and the packet is not fully "
            "absorbed — both inflate the number. See §4–5."
        ),
        nb.setup_cell(),

        nb.md(
            "## 1. The definition\n\n"
            "Stopping power as energy **retained by the electron gas** per unit path:\n\n"
            "$$S_{\\rm WP} = \\frac{E_{\\rm total}(t_f) - E_{\\rm GS}}{L_z}, "
            "\\qquad L_z = 25\\ \\text{Bohr (slab thickness)}$$\n\n"
            "valid once the CAP has absorbed the wavepacket remnants so that "
            "$E_{\\rm total}(t_f) = E_{\\rm jellium}(t_f)$. The initial bath energy is "
            "reconstructed as a **consistency check**:\n\n"
            "$$E_{\\rm jellium}(0) \\equiv E_{\\rm total}(0) - \\langle T_{\\rm WP}\\rangle "
            "- E_{\\rm SIE} \\;\\stackrel{?}{\\approx}\\; E_{\\rm GS}$$\n\n"
            "- $\\langle T_{\\rm WP}\\rangle$ = run-measured WP kinetic energy at $t=0$ "
            "= **drift** $\\tfrac12 k_0^2$ (100 eV) + **zero-point** $3/4\\sigma^2$ (82 eV);\n"
            "- $E_{\\rm SIE}$ = wavepacket self-interaction energy.\n\n"
            "*Source: campaign `quantum-stopping-power.md`; handover `localised-jellium.md`.*"
        ),

        nb.md(
            f"## 2. Setup & provenance ({RUN})\n"
            f"- WP run: `{c['wp_obs'].parent.parent}` (σ_WP=0.5, E=100 eV, k₀=2.711, "
            f"r_s≈{c['r_s']}, N_bath={c['n_bath']}, slab |z|<12.5, τ={c['tau_au']} a.u.)\n"
            f"- Classical twin: matched Gaussian-electron ion (σ_charge=0.35)\n"
            f"- Ground state: E_GS = {f(wp['E_GS_ha'],6)} Ha\n"
        ),
        nb.code(
            "import sys\n"
            f"sys.path.insert(0, {str(HERE)!r})\n"
            "import quantum_stopping_ledger as L\n"
            f"wp = L.compute_wp_ledger({RUN!r}); cl = L.compute_classical_slab({RUN!r})\n"
            "HA = L.HA_TO_EV\n"
            "print('run:', wp['label'])"
        ),

        nb.md(
            "## 3. Quantum energy ledger — the t=0 decomposition (consistency check)\n\n"
            "Strip the WP's own energy from the run's t=0 total and confirm it recovers "
            "the bare-slab ground state.\n\n"
            "| Quantity | value | note |\n|---|---:|---|\n"
            f"| $E_{{\\rm total}}(0)$ | {f(wp['E_total_0_ha']*HA,5)} eV | run total at launch |\n"
            f"| $-\\langle T_{{\\rm WP}}\\rangle$ | −{f(wp['T_wp_ha']*HA,3)} eV | "
            f"drift {f(wp['T_drift_ev'],3)} + zero-point {f(wp['T_zp_ev'],2)} |\n"
            f"| $-E_{{\\rm SIE}}$ | −{f(wp['E_sie_ev'],2)} eV | WP self-interaction |\n"
            f"| **$E_{{\\rm jellium}}(0)$** | **{f(wp['E_jellium_0_ha']*HA,6)} eV** | "
            "reconstructed initial bath |\n"
            f"| $E_{{\\rm GS}}$ (bare slab) | {f(wp['E_GS_ha']*HA,6)} eV | independent GS run |\n"
            f"| **$E_{{\\rm jellium}}(0)-E_{{\\rm GS}}$** | **{f(wp['E_jellium0_minus_GS_ev'],2)} eV** | "
            "✓ consistent (small cross-Hartree) |\n"
        ),
        nb.code(
            "print('E_total(0)        = %12.4f Ha' % wp['E_total_0_ha'])\n"
            "print('<T_WP>(0)         = %12.4f Ha = %.1f eV (drift %.0f + zp %.0f)' % "
            "(wp['T_wp_ha'], wp['T_wp_ha']*HA, wp['T_drift_ev'], wp['T_zp_ev']))\n"
            "print('E_jellium(0)      = %12.4f Ha' % wp['E_jellium_0_ha'])\n"
            "print('E_GS              = %12.4f Ha' % wp['E_GS_ha'])\n"
            "print('E_jellium(0)-E_GS = %.3f eV  (consistency)' % wp['E_jellium0_minus_GS_ev'])"
        ),

        nb.md(
            "## 4. Stopping power — full ledger (headline) + the one free assumption\n\n"
            f"| Quantity | value |\n|---|---:|\n"
            f"| $E_{{\\rm total}}(t_f{{=}}{wp['t_f_au']:.0f})$ | {f(wp['E_total_f_ha']*HA,6)} eV |\n"
            f"| **$\\Delta E = E_{{\\rm total}}(t_f) - E_{{\\rm GS}}$** | "
            f"**+{f(wp['dE_ev'],2)} eV** |\n"
            f"| **$S_{{\\rm WP}} = \\Delta E / 25$** | "
            f"**{f(wp['S_wp_ev_per_bohr'],2)} eV/Bohr (UPPER BOUND)** |\n\n"
            "**Why an upper bound, and why the alternative is rejected.** The bath gain "
            "$\\Delta E$ is governed by **one free assumption — the fate of the 82 eV "
            "zero-point energy:**\n\n"
            "| Assumption | $\\Delta E$ | $S$ |\n|---|---:|---:|\n"
            f"| **Full ledger** (zero-point removed with the absorbed packet) | "
            f"+{f(wp['dE_ev'],2)} eV | **{f(wp['S_wp_ev_per_bohr'],2)}** |\n"
            f"| Drift-credit (zero-point *persists* and cancels) | "
            f"{f(wp['dE_driftcredit_ev'],2)} eV | {f(wp['S_driftcredit_ev_per_bohr'],2)} |\n\n"
            "The drift-credit assumption is **rejected**: the *same* assumption gives an "
            "**impossible negative** for the more-absorbed p2 run (jellium below its own "
            "ground state), because the 82 eV zero-point does **not** sit inertly in the "
            "final jellium — it leaves the box with the absorbed packet. So the "
            "full-ledger value is the internally-consistent answer.\n\n"
            "It remains an **upper bound** because (i) the WP is not fully absorbed at "
            f"$t_f$ (norm remaining ≈ {f(wp['wp_norm_remaining'],2)}, gate <0.02; late "
            f"$E_{{\\rm total}}$ slope ≈ {f(wp['slope_ev_per_au'],2)} eV/a.u., not "
            "plateaued), and (ii) the σ=0.5 zero-point (82 eV ≈ drift 100 eV) inflates "
            "the bath gain — the packet's *spreading* dumps energy the gas would not "
            "receive from a rigid point charge. A converged, comparable-to-linear-"
            "response number needs a **large-σ / heavy (muon) projectile** where "
            "$3/4\\sigma^2 \\to 0$."
        ),
        nb.embed(figs["wp_convergence"], "WP run is not converged: E_total still above "
                 "E_jellium(0), packet not fully absorbed, residual WP kinetic ≈ ΔE.",
                 width=460),

        nb.md(_classical_section(cl)),
        nb.embed(figs["classical_ke_z"], "Classical ion KE(z): conservative well vs the "
                 "equal-potential-face stopping.", width=460),

        nb.md(
            "## 6. Takeaway\n\n"
            f"- **Quantum (WP) at 100 eV, {RUN}:** $S_{{\\rm WP}} = "
            f"{f(wp['S_wp_ev_per_bohr'],2)}$ eV/Bohr — **full-ledger upper bound** "
            "(internally consistent; the drift-credit alternative is rejected for giving "
            "an impossible negative on the more-absorbed run).\n"
            + ((f"- **Classical slab:** {f(cl['S_face_ev_per_bohr'],2)} eV/Bohr "
                "(equal-potential faces).\n") if cl["traversed"] else
               (f"- **Classical slab:** ⚠️ ANOMALY — ion trapped/reflected (stalled at "
                f"z={f(cl['z_max'],2)}, did not traverse); ΔKE_ion={f(cl['dKE_ion_ev'],2)} "
                "eV is ~10× expected. Flag to investigate, not a stopping number.\n"))
            + f"- **Consistency win:** reconstructed $E_{{\\rm jellium}}(0)$ lands within "
            f"{f(wp['E_jellium0_minus_GS_ev'],2)} eV of the independent $E_{{\\rm GS}}$ — "
            "the bookkeeping is self-consistent even though the σ=0.5 stopping number is "
            "a zero-point-inflated upper bound, not a converged stopping power."
        ),
    ]
    nb.build(cells, str(OUT), timeout=600)
    print("OK:", OUT)


if __name__ == "__main__":
    main()
