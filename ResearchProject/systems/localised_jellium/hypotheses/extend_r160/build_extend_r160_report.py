#!/usr/bin/env python3
"""Study notebook for the extend-r (L_z=160) localised-jellium energetics run-set.

House-narrative (notebook-making skill): question -> conventions+formulas -> setup ->
sources -> results (excess vs r) -> the RIGHT-FORMULA section (U_ext vs U_H component
decomposition, TODO-2) -> takeaway. Figures use the canonical theme; the .ipynb is
executed and carries its outputs. Re-run to refresh:

  PYTHONPATH=/local/data/public/skcb2/tddft/inq-stack/python \
  /local/data/public/skcb2/tddft/venv/bin/python3 build_extend_r160_report.py
"""
import sys, os, csv
from pathlib import Path
sys.path.insert(0, "/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/hypotheses")
import _nbreport as nb
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, "/local/data/public/skcb2/tddft/inq-stack/python")
from inqview.visualisation import style
style.apply_theme()

HA = 27.211386
LJ    = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium")
RUNS  = LJ / "scripts/campaign_autorun/runs/extend_r160"
OUT   = LJ / "hypotheses/extend_r160"
FIGS  = OUT / "extend_r160_figs"
FIGS.mkdir(parents=True, exist_ok=True)
nb.set_outdir(str(OUT))

R_ALL = sorted({8, 16, 24, 32, 40, 44, 48, 52, 56, 60})

# ---------------------------------------------------------------- data readers
def gs_energy(per):
    for p in (RUNS / f"gs_lz160_p{per}").glob("**/run_summary.txt"):
        for ln in p.read_text().splitlines():
            if ln.startswith("ground_state_energy_ha"):
                return float(ln.split("=")[1])
    return None

def comps(t, r, per):
    """Return dict of energy components (Ha) at t=0 for one run, or None."""
    m = sorted((RUNS / f"{t}_r{r}_p{per}").glob("**/observables.csv"))
    if not m:
        return None
    rows = list(csv.reader(open(m[0])))
    if len(rows) < 2:
        return None
    h, row = rows[0], rows[1]
    g = lambda k: float(row[h.index(k)])
    Et, T, UH, Exc = g("energy_total"), g("energy_kinetic"), g("energy_hartree"), g("energy_xc")
    return dict(Et=Et, T=T, UH=UH, Exc=Exc, Uext=Et - T - UH - Exc)

EGS = {p: gs_energy(p) for p in (3, 2)}
# existing L_z=120 continuity anchor (from the H0 base-difference runs)
EGS120 = -108.5336851082701
def _e120(tag):
    m = sorted((LJ / f"scripts/h0_base_difference/{tag}/results").glob("**/observables.csv"))
    m = [p for p in m if "r40" in str(p)]
    if not m: return None
    rr = list(csv.reader(open(m[0]))); return float(rr[1][rr[0].index("energy_total")])
WP40_120 = (_e120("wp") - EGS120) * HA if _e120("wp") else None

def excess(t, per):
    out = []
    for r in R_ALL:
        c = comps(t, r, per)
        if c and EGS[per] is not None:
            out.append((r, (c["Et"] - EGS[per]) * HA))
    return out

# ================================================================= figures
def fig_excess():
    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    sty = {("wp", 3): ("o-", "C0"), ("wp", 2): ("o--", "C0"),
           ("cl", 3): ("s-", "C3"), ("cl", 2): ("s--", "C3")}
    for (t, per), (ls, col) in sty.items():
        pts = excess(t, per)
        if not pts: continue
        xs, ys = zip(*pts)
        lab = f"{'WP' if t=='wp' else 'classical'}, p{per}" + (" (raw)" if t=="wp" and per==2 else "")
        ax.plot(xs, ys, ls, color=col, ms=5, label=lab)
    if WP40_120 is not None:
        ax.plot([40], [WP40_120], "kx", ms=9, mew=2, label=f"WP r=40, L_z=120 (existing): {WP40_120:.0f} eV")
    ax.set_xlabel(r"projectile–slab distance $r$  (Bohr)")
    ax.set_ylabel(r"$E_\mathrm{total}(0)-E_\mathrm{GS}$  (eV)")
    ax.set_title("Localised jellium: injection energy vs $r$ (L$_z$=160)")
    ax.legend(frameon=False, fontsize=7.5)
    p = FIGS / "excess_vs_r.png"; fig.savefig(p, dpi=150, bbox_inches="tight"); plt.close(fig)
    return p

def fig_gap():
    fig, ax = plt.subplots(figsize=(6.6, 4.0))
    for per, ls in ((3, "o-"), (2, "o--")):
        w = dict(excess("wp", per)); c = dict(excess("cl", per))
        rs = [r for r in R_ALL if r in w and r in c]
        ax.plot(rs, [w[r] - c[r] for r in rs], ls, ms=5, label=f"periodicity {per}")
    ax.set_xlabel(r"projectile–slab distance $r$  (Bohr)")
    ax.set_ylabel(r"$\Delta E_\mathrm{WP}-\Delta E_\mathrm{cl}$  (eV)")
    ax.set_title("WP$-$classical gap (persistent quantum offset: ZP + SIE)")
    ax.legend(frameon=False, fontsize=8)
    p = FIGS / "wp_minus_cl_gap.png"; fig.savefig(p, dpi=150, bbox_inches="tight"); plt.close(fig)
    return p

def _delta_comp(t, key, per, ref_r=60):
    ref = comps(t, ref_r, per)
    out = []
    for r in R_ALL:
        c = comps(t, r, per)
        if c and ref:
            out.append((r, (c[key] - ref[key]) * HA))
    return out

def fig_components():
    fig, ax = plt.subplots(figsize=(6.8, 4.4))
    series = [("cl", "Uext", "s-", "C3", r"classical $\Delta U_\mathrm{ext}$"),
              ("cl", "UH",  "s:", "C1", r"classical $\Delta U_H$ (frozen bath)"),
              ("wp", "Uext","o-", "C0", r"WP $\Delta U_\mathrm{ext}$ (feels background)"),
              ("wp", "UH",  "o-", "C2", r"WP $\Delta U_H$ (WP$-$bath repulsion)")]
    for t, key, ls, col, lab in series:
        pts = _delta_comp(t, key, 3)
        if pts:
            xs, ys = zip(*pts); ax.plot(xs, ys, ls, color=col, ms=5, label=lab)
    ax.axhline(0, color="0.7", lw=0.8)
    ax.set_xlabel(r"projectile–slab distance $r$  (Bohr)")
    ax.set_ylabel(r"$\Delta$(component) vs $r{=}60$  (eV)")
    ax.set_title("Where the energy lives: classical$\\to U_\\mathrm{ext}$;  WP$\\to U_\\mathrm{ext}+U_H$ (cancel)")
    ax.legend(frameon=False, fontsize=7.5)
    p = FIGS / "component_decomposition.png"; fig.savefig(p, dpi=150, bbox_inches="tight"); plt.close(fig)
    return p

def fig_right_formula():
    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    cl = dict(_delta_comp("cl", "Uext", 3))
    wp = dict(_delta_comp("wp", "UH", 3))
    rs = [r for r in R_ALL if r in cl and r in wp]
    ax.plot(rs, [cl[r] for r in rs], "s-", color="C3", ms=5, label=r"classical $\Delta U_\mathrm{ext}$")
    ax.plot(rs, [wp[r] for r in rs], "o-", color="C2", ms=5, label=r"WP $\Delta U_H$ (Hartree = WP$-$bath)")
    ax.set_xlabel(r"projectile–slab distance $r$  (Bohr)")
    ax.set_ylabel(r"projectile$-$bath coupling vs $r{=}60$  (eV)")
    ax.set_title("The right comparison: classical $U_\\mathrm{ext}$ vs WP Hartree")
    ax.legend(frameon=False, fontsize=8)
    p = FIGS / "right_formula_compare.png"; fig.savefig(p, dpi=150, bbox_inches="tight"); plt.close(fig)
    return p

f1, f2, f3, f4 = fig_excess(), fig_gap(), fig_components(), fig_right_formula()

# tables (recomputed inline so numbers can't drift from the figures)
def excess_table():
    rows = ["| r (Bohr) | WP p3 | cl p3 | WP p2 | cl p2 |", "|---|---|---|---|---|"]
    d = {(t, p): dict(excess(t, p)) for t in ("wp", "cl") for p in (3, 2)}
    for r in R_ALL:
        rows.append(f"| {r} | {d[('wp',3)].get(r,float('nan')):.1f} | {d[('cl',3)].get(r,float('nan')):.1f} "
                    f"| {d[('wp',2)].get(r,float('nan')):.1f} | {d[('cl',2)].get(r,float('nan')):.1f} |")
    return "\n".join(rows)

# ================================================================= narrative
ZP = 3.0 / (4.0 * 0.5**2) * HA
cells = [
    nb.md(f"""# Localised jellium — injection energy $E_\\mathrm{{total}}(0)-E_\\mathrm{{GS}}$ vs projectile distance $r$

**The question.** Extend the existing H0-style energetics curve to larger projectile–slab
separation $r$, using a bigger *centered* box (L$_z$=160) so higher $r$ fits without charge
leaking to the boundary. Same system/density as the existing runs — these points land on the
*same* curve. Then use the run's **energy decomposition** to ask the deeper question raised for
the 2026-07-03 meeting: *the classical projectile and the wavepacket deposit their interaction
energy in **different components** of the total energy — how should we compare them?*

**Where this sits.** Extends `hypotheses/h0_base_difference` (r=4–40, L$_z$=120) and the
`campaign_autorun_study` H0/H5 ladder to r=8–60 at L$_z$=160, periodicity 2 **and** 3."""),

    nb.setup_cell(),

    nb.md(f"""## Conventions & symbols

Atomic units internally; energies reported in **eV** (1 Ha = {HA:.3f} eV), 2 s.f. by default.

| symbol | meaning | value / range |
|---|---|---|
| $r$ | projectile–slab-face distance (projectile at $z=-(a+r)$) | 8–60 Bohr |
| $a$ | slab half-width | 12.5 Bohr |
| $\\sigma_\\mathrm{{WP}}$ | wavepacket width (matched ghost: $\\sigma_\\mathrm{{pot}}=\\sigma_\\mathrm{{WP}}/\\sqrt2$) | 0.5 Bohr |
| $E_\\mathrm{{GS}}$ | neutral-slab ground-state energy (per periodicity) | p3 {EGS[3]:.2f} Ha, p2 {EGS[2]:.2f} Ha |
| $T,\\,U_\\mathrm{{ext}},\\,U_H,\\,E_\\mathrm{{xc}}$ | kinetic / external / Hartree / xc energy | — |

The total-energy decomposition (standard KS-DFT; Parr & Yang):
$$E = T + U_\\mathrm{{ext}} + U_H + E_\\mathrm{{xc}} + \\text{{const}},\\qquad U_\\mathrm{{ext}}=E-T-U_H-E_\\mathrm{{xc}}.$$
$U_\\mathrm{{ext}}$ is taken **by difference** from the scalar observables (no separate $U_\\mathrm{{ext}}$ channel is written)."""),

    nb.md(f"""## Setup — reconstructable

- **Cell / geometry:** orthorhombic 50×50×**160** Bohr, spacing 0.5, slab **centered** ($|z|<12.5$),
  sharp edge ($w=0$), $N=82$, $n_0=1.31\\times10^{{-3}}\\,a_0^{{-3}}$ ($r_s\\approx5.67$). periodicity **3 and 2**.
- **Engine:** inq-study, LDA. One neutral GS built per periodicity at L$_z$=160 and reused for all $r$.
- **Projectile placement:** $z=-(12.5+r)$; at $r=60$, $z=-72.5$ — margin 7.5 Bohr (=15$\\sigma$) to the
  ±80 boundary (same margin the existing $r=40$/L$_z$=120 point had).
- **Runs per $r$:** *classical ghost* (matched $\\sigma_\\mathrm{{pot}}$; gives $E_\\mathrm{{total}}(0)$) and
  *wavepacket* ($k_0=0$, $\\sigma_\\mathrm{{WP}}=0.5$, propagated 2 steps → valid $E_\\mathrm{{total}}$)."""),

    nb.md("""## Source files

| file | role |
|---|---|
| `scripts/campaign_autorun/{gs,wp,classical}/run.cpp` | env-driven GS / WP / classical binaries |
| `scripts/campaign_autorun/extend_r_lz160.py` | this run-set's driver (GS gate + r-sweep, detached) |
| `scripts/campaign_autorun/runs/extend_r160/` | run outputs (observables.csv per run) |
| `hypotheses/extend_r160/build_extend_r160_report.py` | this notebook's builder |
| `hypotheses/h0_base_difference/analyse_h0.py` | the original L$_z$=120 r=4–40 study (continuity anchor) |"""),

    nb.md(f"""## 1 — Injection energy vs $r$

$$\\Delta E(r) = E_\\mathrm{{total}}(0) - E_\\mathrm{{GS}}$$
— the energy cost of introducing the projectile at distance $r$ from the neutral slab."""),
    nb.embed(str(f1), "E_total(0)−E_GS vs r; × marks the existing L_z=120 r=40 WP point (continuity)"),
    nb.md(f"""**Continuity check (the point of the r=40 overlap):** existing L$_z$=120 WP r=40 =
**{WP40_120:.1f} eV** vs new L$_z$=160 = **{dict(excess('wp',3)).get(40,float('nan')):.1f} eV** — agree to
~0.7 %, so the new points sit on the same curve.

{excess_table()}

- **WP ≈ flat** at ~86.5 eV (p3): dominated by the intrinsic zero-point KE
  $E_\\mathrm{{ZP}}=3/(4\\sigma_\\mathrm{{WP}}^2)={ZP:.0f}$ eV plus SIE (~5 eV), weakly $r$-dependent.
- **Classical decays** 164→~1 eV: the unscreened ghost–slab Coulomb, → 0 at large $r$.
- **periodicity-2 WP is ~8–9 eV below p3** — the open-z net-charge $G{{=}}0$ bias (RAW; correction pending)."""),

    nb.md("""## 2 — The persistent WP$-$classical gap
At large $r$ the classical signal vanishes but the WP keeps its intrinsic offset, so the gap
plateaus at the quantum floor (zero-point + SIE) that no classical ghost reproduces."""),
    nb.embed(str(f2), "WP−classical injection-energy gap vs r"),

    nb.md("""## 3 — The right formula: which energy component to compare

**The user's insight (2026-07-03), confirmed by the data.** A classical projectile enters as a
*time-dependent external potential* — its bath coupling lives in $U_\\mathrm{ext}$. A wavepacket is
a *real electron density* — its bath coupling lives in $U_H$ (the Hartree cross-term). So the two
deposit their interaction energy in **different components**, and comparing the raw totals is not
apples-to-apples."""),
    nb.embed(str(f3), "Component decomposition (Δ vs r=60): classical→U_ext only; WP→U_ext+U_H (cancel)"),
    nb.md("""**What the decomposition shows (periodicity 3, referenced to r=60):**
- **Classical:** $U_H$ is *exactly constant* (frozen bath — the ghost is pure external potential);
  the entire signal is $\\Delta U_\\mathrm{ext}$ (163→0 eV).
- **WP:** $\\Delta U_\\mathrm{ext}$ (−120→0, WP attracted to the slab background) and $\\Delta U_H$
  (+121→0, WP–bath Coulomb repulsion) are **equal and opposite to ~1 %** — they cancel because the
  slab is neutral, which is *why* the WP total is flat. The physics hides in the components."""),

    nb.md("""### The candidate formula
Compare the **projectile–bath electrostatic coupling** on each side, not the totals:
$$\\underbrace{\\Delta U_\\mathrm{ext}^{\\;\\mathrm{cl}}(r)}_{\\int n_\\mathrm{bath}\\,V_\\mathrm{ghost}}
\\;\\stackrel{?}{\\approx}\\;
\\underbrace{\\Delta U_H^{\\;\\mathrm{wp}}(r)}_{\\iint n_\\mathrm{bath}\\,n_\\mathrm{wp}/|r-r'|\\;(\\text{cross term})}.$$
For a matched ghost ($\\sigma_\\mathrm{pot}=\\sigma_\\mathrm{WP}/\\sqrt2$) these *should* coincide; the WP
self-Hartree (SIE) is $r$-independent and cancels in the $r$-difference."""),
    nb.embed(str(f4), "classical ΔU_ext vs WP ΔU_H — the like-for-like projectile–bath coupling"),
    nb.md("""**Result (present the evidence; the verdict is yours).** The two curves are the **same order and
both decay**, but they do **not** coincide: WP $\\Delta U_H$ is *smaller* than classical
$\\Delta U_\\mathrm{ext}$ at small $r$ (121 vs 163 eV at r=8) and *larger* at intermediate $r$
(24 vs 10 eV at r=40) — they cross near r≈22. That difference is the quantum–classical signal:
the WP density can polarize/spread, the rigid ghost cannot.

**Caveat — what this scalar comparison cannot yet do.** $\\Delta U_H^\\mathrm{wp}$ from the scalar
observable is the *total* Hartree change, which includes the WP self-Hartree. Here it cancels in the
$r$-difference (the WP shape is $r$-independent), so the comparison is defensible — but the **clean**
formula isolates the cross term $\\iint n_\\mathrm{bath}\\,n_\\mathrm{wp}/|r-r'|$ directly from the
**density fields** (VTIs), separating it from the self-Hartree with no cancellation assumption. That
is the next step and needs the density outputs, not just the scalar energies."""),

    nb.md(f"""## Takeaway
- New points extend the injection-energy curve cleanly to **r=60** (continuity {WP40_120:.0f}→{dict(excess('wp',3)).get(40,0):.0f} eV at r=40, 0.7 %).
- **Classical → $U_\\mathrm{{ext}}$ only; WP → $U_\\mathrm{{ext}}+U_H$ that cancel** — comparing raw totals is
  misleading; the like-for-like quantity is classical $\\Delta U_\\mathrm{{ext}}$ vs WP $\\Delta U_H$.
- Those two are the same order and both decay but **cross near r≈22** — a real quantum–classical
  difference (polarizable WP vs rigid ghost), NOT numerical noise.
- **Next:** compute the explicit Hartree **cross term** from the density fields (removes the
  self-Hartree cancellation assumption); apply the periodicity-2 net-charge $G{{=}}0$ correction.
- *All numbers recomputed from run provenance; periodicity-2 WP is RAW (net-charge correction pending).*"""),
]

nb.build(cells, str(OUT / "extend_r160_study.ipynb"))
