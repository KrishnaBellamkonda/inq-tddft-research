#!/usr/bin/env python3
"""build_index.py — the guided-analysis INDEX notebook for the annular-tube S(v) sweep.

This is a *navigator*, not a re-analysis: it walks the reader through every analysis
notebook / figure set of the `cylindrical-jellium-projectile` campaign, one stop after
another, in a deliberate **trust-chain order** (validate the ground state → headline
synthesis → per-run evidence cleanest-first → quantum rung → open verdicts).

Every link is REAL: the builder resolves each target against this folder and only emits
a live markdown link when the file exists; a missing target is flagged `PENDING` in the
status cell instead of producing a dead link. Rebuild after adding/removing notebooks:

    /local/data/public/skcb2/tddft/venv/bin/python3 build_index.py

Writes ``annular_sv_index.ipynb`` beside this script (file-placement: run-tied analysis
lives in ``hypotheses/annular_sv/``).
"""
from __future__ import annotations
from pathlib import Path
import nbformat as nbf

HERE = Path(__file__).resolve().parent
OUT = HERE / "annular_sv_index.ipynb"


def exists(rel: str) -> bool:
    return (HERE / rel).exists()


def link(label: str, rel: str) -> str:
    """A live markdown link if the target resolves, else a struck-through PENDING tag."""
    return f"[{label}]({rel})" if exists(rel) else f"~~{label}~~ *(PENDING — not yet built)*"


# --------------------------------------------------------------- the reading order
# Each stop: (n, title, [ (label, relpath) ... ], method, what_to_look_for)
REPORT = "annular_sv_report.ipynb"
WPNB = "wp_rs6_v0p30_run_notebook.ipynb"


def per_run(label: str) -> str:
    return f"per_run_figs/{label}/"


STOPS = [
    (1, "Ground-state validation — trust the tube before the dynamics",
     [("radial n(r): r_s=6", "gs_validation/gs_radial_rs6.png"),
      ("radial n(r): r_s=4", "gs_validation/gs_radial_rs4.png"),
      ("radial n(r): r_s=2", "gs_validation/gs_radial_rs2.png"),
      ("xz/xy slices: r_s=6", "gs_validation/gs_slices_rs6.png"),
      ("xz/xy slices: r_s=4", "gs_validation/gs_slices_rs4.png"),
      ("xz/xy slices: r_s=2", "gs_validation/gs_slices_rs2.png")],
     "Converged LDA ground state of each empty tube: radial density profile n(r) and "
     "xz/xy density slices, per wall density r_s∈{6,4,2}.",
     "Interior density flat at n0; a clean annular wall between R_in=5 and R_out=13 Bohr; "
     "no spurious ripples or charge leaking into the bore; neutrality (N matches the "
     "design count 24/48/136). If the GS is wrong, every downstream number is suspect."),

    (2, "Synthesis §1 — stopping power S(v) and friction slope β(r_s)",
     [("report → §1 (S(v), β)", REPORT),
      ("β overview figure", "Sv_beta.png"),
      ("S(v) table (CSV)", "Sv_results.csv")],
     "S(v0) = initial drag = −d(KE_ion)/ds from a linear regression of ke_ion_ha vs path "
     "s over the early near-constant-velocity window vz≥0.85·v0 (light-projectile rule). "
     "β(r_s)=dS/dv is the low-velocity friction coefficient — the hypothesis quantity. "
     "NOTE: this synthesis panel uses the KINETIC (KE) channel.",
     "Is β(r_s) monotonic across r_s={6,4,2}? Are the three S(v) points per r_s resolved "
     "beyond their error bars? This is the campaign's headline claim — read it first, "
     "then let the per-run stops (6–8) tell you whether to trust each point."),

    (3, "Synthesis §2 — induced wall current (flow → current)",
     [("report → §2 (current_z)", REPORT)],
     "Raw axial current current_z(t) from each run's observables.csv, overlaid across "
     "(r_s, v). No transform — the time-domain hydrovoltaic signal.",
     "A coherent induced axial current building as the projectile glides through the "
     "bore; relative magnitude vs r_s and v. (The FFT of this signal lives per-run, "
     "stop 6–9, via the fourier-analysis 6-stage panel.)"),

    (4, "Synthesis §3 — wake structure",
     [("report → §3 (Δn wake)", REPORT)],
     "A mid-time induced-density frame Δn=n(t)−n(0), xz slice at y=0, shown linear and "
     "signed-log10 (canonical load_vti, physical order — no fftshift).",
     "A lagging wall-charge wake trailing the projectile (the retarding image field). "
     "Confirm it is an established interior wake, NOT a periodic-image / box-edge "
     "artifact — this bears on the L_z-adequacy verdict (stop 10)."),

    (5, "Synthesis §4 — quantum rung: wavepacket vs classical ghost",
     [("report → §4 (WP vs ghost)", REPORT)],
     "ΔE_system vs projectile path for the classical r_s=6, v=0.30 run (solid) and its "
     "matched WP run (dashed; centroid path from wp_real_space_stats). Slope = S.",
     "Does the quantum WP deposit energy at a different slope than the classical ghost? "
     "This is the leading-order quantum effect on S. Caveat (already in the report): the "
     "electron-as-cation proxy rests on charge-even S; Barkas is the charge-odd correction."),

    (6, "Per-run deep dives — r_s=6 (dilute wall, CLEAN channels)",
     [("report → run rs6_v0p15", REPORT), ("figs rs6_v0p15", per_run("rs6_v0p15")),
      ("report → run rs6_v0p30", REPORT), ("figs rs6_v0p30", per_run("rs6_v0p30")),
      ("report → run rs6_v0p45", REPORT), ("figs rs6_v0p45", per_run("rs6_v0p45"))],
     "Full per-run battery (read v=0.15→0.30→0.45): density-matrix GIFs with the moving "
     "projectile overlaid, z–t carpets, Method-A stopping (PRIMARY = ΔE_total slope; "
     "SANITY = −dKE_ion/dx), energy decomposition, current+dipole, and the 6-stage FFT.",
     "Start here to build intuition on the runs that AGREE: channel ratio (ΔE_total vs "
     "−dKE_ion) ≈ 1.01–1.05 with r²≥0.90. Note the projectile decelerating in v_z(t) and "
     "the early-window fit; confirm N(t)≈const (no absorption)."),

    (7, "Per-run deep dives — r_s=4 (mostly clean)",
     [("report → run rs4_v0p15", REPORT), ("figs rs4_v0p15", per_run("rs4_v0p15")),
      ("report → run rs4_v0p30", REPORT), ("figs rs4_v0p30", per_run("rs4_v0p30")),
      ("report → run rs4_v0p45", REPORT), ("figs rs4_v0p45", per_run("rs4_v0p45"))],
     "Same battery as stop 6, intermediate wall density.",
     "Mostly clean channel agreement; watch rs4_v0p15 — marginal (~12% divergence). Ask "
     "whether the S(v) point at v=0.15 is trustworthy or borderline."),

    (8, "Per-run deep dives — r_s=2 (dense wall, ⚠ ALL FLAGGED)",
     [("report → run rs2_v0p15", REPORT), ("figs rs2_v0p15", per_run("rs2_v0p15")),
      ("report → run rs2_v0p30", REPORT), ("figs rs2_v0p30", per_run("rs2_v0p30")),
      ("report → run rs2_v0p45", REPORT), ("figs rs2_v0p45", per_run("rs2_v0p45"))],
     "Same battery, densest wall — but the short L_z=10 tube gives a tiny traversal "
     "(rs2_v0p15 moves only ~2.5 Bohr).",
     "⚠ ALL THREE are FLAGGED: the two stopping channels DIVERGE — rs2_v0p15 ratio 2.04 "
     "r²=0.23; rs2_v0p30 ratio 1.18 r²=0.69; rs2_v0p45 ratio 1.11. The synthesis headline "
     "(stop 2) uses the KE channel (rs2_v0p15 S=0.0067) while Method-A PRIMARY gives 0.0136. "
     "This is the campaign's open question — carry your read into stop 10."),

    (9, "Quantum rung — WP run deep-dive (wp_rs6_v0p30)",
     [("standalone run-notebook", WPNB), ("tube-aware figs", per_run("wp_rs6_v0p30")),
      ("report → WP run section", REPORT)],
     "The full single-run WP battery: 9-GIF density matrix {n, Δn-vs-0, per-step Δn} × "
     "{total, wavepacket |ψ|², bath = total−WP}; 1D momentum n(k) incl/excl WP; KL "
     "divergence of the WP momentum vs launch; energy decomposition; 6-stage FFT.",
     "Momentum broadening / drift of the wavepacket (KL rising from 0); the bath (wall) "
     "anti-wake in the total−WP channel; how the quantum deposit compares to the classical "
     "ghost (stop 5). ⚠ CAVEAT: the standalone run-notebook is built by the shared, "
     "SLAB-oriented builder — its 'slab face ±12.5 Bohr' dotted lines and any ΔE/L_z "
     "slab-thickness stopping number are NOT physical for this periodic TUBE (the wall is "
     "radial at x=±5,±13; authoritative stopping is Method A in the report / per-run figs). "
     "The tube-aware figures are the per-run figs linked above."),

    (10, "Open verdicts — what to resolve (your call)",
     [],
     "The decisions this analysis must settle. Verdicts are yours (house rule): record "
     "them below; do not let a builder decide.",
     None),
]

VERDICTS = [
    "**r_s=2 stopping channel** — accept KE-channel S, accept Method-A ΔE_total S, or "
    "rerun r_s=2 with a longer L_z? (stops 2, 8)",
    "**β(r_s) robustness** — is β monotonic AND resolved once the r_s=2 point is decided? "
    "Is the r_s=6 small-gas (~24 e) finite-size/shell effect acceptable? (stops 2, 6)",
    "**Wake / L_z adequacy** — is the wake an established interior structure, or "
    "contaminated by periodic images / box edges, per run? (stops 4, 6–8)",
    "**Quantum effect on S** — does the WP deposit differ from the classical ghost beyond "
    "noise, and is the electron-as-cation (charge-even) proxy defensible here? (stops 5, 9)",
]


# ------------------------------------------------------------------------- assemble
cells: list = []

cells.append(nbf.v4.new_markdown_cell(
    "# Cylindrical (annular) jellium — GUIDED ANALYSIS INDEX\n\n"
    "*A navigator for the `cylindrical-jellium-projectile` campaign. It walks you "
    "through every analysis notebook and figure set **one stop after another** in a "
    "deliberate order, so you can analyse the results carefully without deciding the "
    "path yourself.*\n\n"
    "**How to use.** Go top to bottom. Each stop gives you: the **real links** to open, "
    "the **method** (how the number/figure was made), **what to look for**, and a blank "
    "**Your read:** line — write your observation there (it is your voice, not the "
    "builder's). The order is a **trust chain**: validate the ground state first, then "
    "the headline synthesis, then the per-run evidence *cleanest-first* (r_s=6 → 4 → 2), "
    "then the quantum rung, then the open verdicts.\n\n"
    "> Rebuild this index after adding/removing notebooks: "
    "`venv/bin/python3 build_index.py`. Links are validated at build time — a missing "
    "target shows as `PENDING`, never a dead link."))

# live manifest / link-integrity cell
cells.append(nbf.v4.new_code_cell(
    "from pathlib import Path\n"
    "HERE = Path.cwd()\n"
    "targets = {\n"
    "    'annular_sv_report.ipynb (run-SET report)': 'annular_sv_report.ipynb',\n"
    "    'wp_rs6_v0p30_run_notebook.ipynb (WP deep-dive)': 'wp_rs6_v0p30_run_notebook.ipynb',\n"
    "    'Sv_results.csv': 'Sv_results.csv',\n"
    "    'Sv_beta.png': 'Sv_beta.png',\n"
    "    'gs_validation/ (6 PNGs)': 'gs_validation',\n"
    "    'per_run_figs/ (10 run dirs)': 'per_run_figs',\n"
    "}\n"
    "print(f\"{'artefact':<48} {'status'}\")\n"
    "print('-'*64)\n"
    "for name, rel in targets.items():\n"
    "    print(f\"{name:<48} {'OK' if (HERE/rel).exists() else 'PENDING'}\")\n"))

cells.append(nbf.v4.new_markdown_cell(
    "## Stop 0 — north star & scope\n\n"
    "A charged projectile (classical electron: Gaussian UPF, mass mₑ, free Ehrenfest) "
    "glides on-axis down the hollow bore of a **periodic annular jellium tube** "
    "(R_in=5, R_out=13 Bohr wall). We measure its electronic stopping power S(v) and the "
    "low-velocity friction slope β(r_s)=dS/dv vs wall density r_s∈{6,4,2}. **North star:** "
    "nanotube hydrovoltaics / quantum friction — flow-induced electronic current, with the "
    "wall modelled as jellium of variable r_s. Sweep = 3 r_s × 3 v (9 classical) + 1 "
    "wavepacket run (`wp_rs6_v0p30`).\n\n"
    "*Provisional interpretation carried from the report/handover is tagged where it "
    "appears; the verdicts (stop 10) are left to you.*"))

for n, title, links, method, look in STOPS:
    parts = [f"## Stop {n} — {title}\n"]
    if links:
        parts.append("**Open:** " + "  ·  ".join(link(lab, rel) for lab, rel in links) + "\n")
    parts.append(f"**Method.** {method}\n")
    if look:
        parts.append(f"**What to look for.** {look}\n")
    if n == 10:
        parts.append("\n".join(f"- [ ] {v}" for v in VERDICTS) + "\n")
    parts.append("\n**Your read:** _______________________________________________")
    cells.append(nbf.v4.new_markdown_cell("\n".join(parts)))

nb = nbf.v4.new_notebook(cells=cells, metadata={
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"}})
nbf.write(nb, str(OUT))
print(f"wrote {OUT}  ({len(cells)} cells)")
n_ok = sum(exists(r) for _, _, ls, _, _ in STOPS for _, r in ls)
n_tot = sum(len(ls) for _, _, ls, _, _ in STOPS)
print(f"link targets resolved: {n_ok}/{n_tot}")
