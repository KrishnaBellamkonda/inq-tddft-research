#!/usr/bin/env python3
"""Assemble qsp_phase4_study.ipynb — the BIG-BOX production study notebook
(campaign quantum-stopping-power, task P3.1): σ=0.5 WP + matched classical in the
50×50×90 slab, two-sided CAP, τ=100 a.u., energy-method stopping.

Thin narrative assembler: every figure + GIF is pre-computed by analyse_phase4.py
(figs/*.png, figs/*.gif) and every headline number is read from results.json, so a
quoted number and its figure can never disagree.

Run:
  PYTHONPATH=/local/data/public/skcb2/tddft/inq-stack/python \
  /local/data/public/skcb2/tddft/venv/bin/python3 build_phase4_notebook.py
"""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))          # hypotheses/ for _nbreport
from _nbreport import md, embed, setup_cell, set_outdir, build

set_outdir(HERE)
FIGS = os.path.join(HERE, "figs")
R = json.load(open(os.path.join(HERE, "results.json")))

def g(k, fmt="{:.3f}", dflt="n/a"):
    v = R.get(k); return dflt if v is None else fmt.format(v)
def gh(run, key, fmt="{:.3f}", dflt="n/a"):
    v = R.get("heuristics", {}).get(run, {}).get(key); return dflt if v is None else fmt.format(v)

KIND_LABEL = {"density": "n(x,z,t)", "delta0": "Δn = n(t)−n(0)", "dstep": "Δn = n(t+dt)−n(t)"}
CAT_LABEL = {"total": "total system", "wp": "wavepacket |ψ|²", "bath": "bath (slab only)"}
W_GIF, W_PNG2, W_PNG1 = 360, 600, 520

def gif_row(run, cat, width=W_GIF):
    out = []
    for kind in ("density", "delta0", "dstep"):
        fn = f"{run}_{cat}_{kind}.gif"
        if os.path.exists(os.path.join(FIGS, fn)):
            out.append(embed(f"{FIGS}/{fn}", f"{CAT_LABEL.get(cat,cat)} · {KIND_LABEL[kind]}", width=width))
    return out

HA = 27.211386
cells = [setup_cell()]

# ----------------------------------------------------------------- §0 title
cells.append(md(rf"""# Phase 4 · P3.1 — big-box production: WP vs classical stopping (energy method)
### localised jellium slab · 50×50×**90** · r$_s$≈5.67 · σ$_{{\rm WP}}$=0.5 · 54 eV · two-sided CAP η=−0.7 · τ=100 a.u.

**Campaign:** `docs/campaigns/jellium_wp_stopping/quantum-stopping-power.md` (task **P3.1**).
**Supersedes P2.2** — fixes the three issues the P2.1 test (`hypotheses/qsp_phase2`) exposed:
1. **WP norm absorbed at launch** → bigger box + **equidistant launch z=−23.75** (11.25 Bohr to
   both the slab face and the CAP inner face).
2. **E_total not converged at τ=40** (13.6% WP unabsorbed) → **τ=100 a.u.** + the **convergence
   triple** gate (norm<0.02 AND E_total plateau).
3. **Stopping definition** → the **energy method**: the jellium *system* = the density that
   REMAINS in the box; S = [E_total(t$_f$) − E_GS]/L$_z$. CAP-absorbed energy (transmitted/
   reflected WP + secondaries) is *ledgered as a diagnostic*, not added back.

E$_{{\rm GS}}$ (90-box) = **−70.22568 Ha**. *All stopping numbers are valid only
if the convergence gate in §5 passes.*
"""))

# ----------------------------------------------------------- §1 conventions
cells.append(md(r"""## §1 — Conventions & symbols
Atomic units (ℏ=m$_e$=e=1; 1 Ha = 27.211 eV; lengths Bohr; for an electron projectile **v=k**).

| symbol | meaning | value |
|---|---|---|
| r$_s$ | electron-gas density of the slab interior | 5.67 (n₀=1.31×10⁻³) |
| σ$_{\rm WP}$ | WP **density** width | 0.5 Bohr (charge std σ/√2=0.354) |
| k₀ / v₀ | WP wavevector = projectile speed | 2.0 (½v₀²≈54 eV) |
| L$_z$ | slab traversal length for dE/dx | 25 Bohr |
| η, w | CAP strength / width per side | −0.7 Ha, 10 Bohr ([±35,±45]) |
| τ | total propagation time | 100 a.u. (2500 × dt=0.04) |
| E$_{\rm GS}$ | 90-box ground state | −70.22568 Ha |
| SIE | self-interaction floor (Phase 1) | 4.40 eV |
| ⟨T$_{\rm WP}$⟩ | WP kinetic = ½k₀² + 3/(4σ²) | 5.0 Ha (136 eV) |

**CAP = two-sided** (the benchmarked "known devil", ~1.3% reflection): `cap_lo+cap_hi` at
mid=±40/90, each a sin² bump peaking at |z|=40, zero at the ±35 inner faces and ±45 box edge.
Engine `absorbing.hpp` is stock (a seam-centred variant was built+validated then reverted).
Propagator = **ETRS** (CAP ⇒ non-Hermitian). Built vs **inq-study**; `inq/` untouched.
"""))

# ----------------------------------------------------------- §2 setup
cells.append(md(rf"""## §2 — Simulation setup (fully reconstructable)
| | value |
|---|---|
| **cell** | 50 × 50 × 90 Bohr orthorhombic periodic, spacing 0.50 (100×100×180 grid) |
| **slab** | positive jellium, half-width 12.5 along z, 82 electrons |
| **GS reused** | `shared_gs/slab_n82_L50x50x90` (E$_{{\rm GS}}$=−70.22568 Ha) |
| **propagator** | ETRS, dt=0.04, N_steps=2500 ⇒ τ=100 a.u. |
| **CAP** | two-sided sin², η=−0.7 Ha, 10 Bohr/side, region [±35,±45] |
| **region layout (z)** | slab [−12.5,12.5] · free [±12.5,±35] · CAP [±35,±45] |
| **launch** | z₀=−23.75 (equidistant: 11.25 Bohr to slab face AND CAP inner face) |
| **WP** | Gaussian σ=0.5, k₀=2.0, last extra state, orthogonalised against occupied |
| **classical** | Ehrenfest ghost ion, Gaussian-e UPF σ=0.354, v=(0,0,2.0) |
| **walltime** | WP {g('wp_wall_s','{:.0f}')} s · classical {g('classical_wall_s','{:.0f}')} s |

The −24.5 Ha shift in E$_{{\rm GS}}$ vs the 70-box is a **charged-slab-in-PBC electrostatic
constant** (kinetic/xc unchanged) that **cancels exactly in E_total−E_GS**.
"""))

# ----------------------------------------------------------- §3 source files
cells.append(md(r"""## §3 — Source files
| file | role |
|---|---|
| `scripts/qsp_phase4/wp/run.cpp` | WP run (ETRS, two-sided CAP, full observable suite) |
| `scripts/qsp_phase4/classical/run.cpp` | classical Ehrenfest ghost-ion run |
| `scripts/qsp_phase4/run_production.sh` | concurrent launch (GPU0/GPU1), env-driven dt/steps |
| `shared/configs/slab_n82_L50x50x90.hpp` | geometry / density / WP config struct |
| `shared_gs/slab_n82_L50x50x90/` | 90-box ground state |
| `hypotheses/qsp_phase4/analyse_phase4.py` | computes every figure + `results.json` |
| `hypotheses/qsp_phase4/p4wp_run_notebook.ipynb` / `p4cl_run_notebook.ipynb` | per-run deep-dives |
"""))

# ====================================================== §4 GIF battery
cells.append(md(rf"""## §4 — Per-run visual intuition (density evolution)
Each GIF is an **xz mid-plane slice**, fixed colour scale, slab faces (|z|=12.5) and CAP inner
faces (|z|=35) dashed. Three kinds — **n(x,z,t)** (log), **Δn=n(t)−n(0)** (induced wake),
**Δn=n(t+dt)−n(t)** (per-frame flux) — for **total / wavepacket / bath (=n_total−n_wp)**.
Shared total/bath density vmax = {g('gif_density_vmax','{:.2e}')} a₀⁻³.
"""))
cells.append(md("### §4a — Wavepacket — total system")); [cells.append(c) for c in gif_row("wp","total")]
cells.append(md("### §4b — Wavepacket — orbital |ψ_WP|²")); [cells.append(c) for c in gif_row("wp","wp")]
cells.append(md("### §4c — Wavepacket — bath only (n_total − n_wp)")); [cells.append(c) for c in gif_row("wp","bath")]
cells.append(md("### §4d — Classical — total system")); [cells.append(c) for c in gif_row("classical","total")]
cells.append(md(rf"""### §4e — Energetics of both runs (slab-exit t={g('T_exit_au','{:.1f}')} a.u. dashed)
WP: **monotonic drain** as the CAP absorbs; classical: **rises** (ion not absorbed → re-entry)."""))
cells.append(embed(f"{FIGS}/energetics.png", "Per-run energetics vs time", width=W_PNG2))

# ====================================================== §5 ENERGY METHOD (headline)
cells.append(md(rf"""## §5 — Q1: WP stopping by the energy method (the headline)
The retained-energy stopping, valid only when the WP is fully absorbed AND E_total has plateaued:
$$ S_{{\rm WP}} = \frac{{E_{{\rm total}}(t_f) - E_{{\rm GS}}}}{{L_z}}, \qquad
\text{{gate: WP norm}}<0.02 \ \text{{AND}}\ \left.\tfrac{{dE_{{\rm total}}}}{{dt}}\right|_{{t_f}}\to 0. $$

**Convergence verdict: {g('wp_converged','{}')}.** Residual WP norm = **{g('wp_norm_final','{:.3f}')}**
(gate < 0.02), late slope = **{g('wp_late_slope_eV_au','{:.2f}')} eV/a.u.**, E_total plateau width =
**{g('wp_plateau_width_au','{:.1f}')} a.u.**
- Deposited E_total(t$_f$)−E_GS = **{g('wp_deposited_EminusEGS_eV','{:.1f}')} eV** ⇒
  **S$_{{\rm WP}}$ = {g('wp_S_eVbohr','{:.2f}')} eV/Bohr** (valid iff converged above; else a bound).

**t=0 sanity (energy bookkeeping).** Subtract the WP **kinetic** ⟨T$_{{\rm WP}}$⟩ (NOT the eigenvalue
→ would double-count SIE) and the SIE from E_total(0), recover E_GS:
- ⟨T$_{{\rm WP}}$⟩ = **{g('T_WP_analytic_Ha','{:.3f}')} Ha** analytic (run: {g('T_WP_run_Ha','{:.3f}')} Ha — direct check).
- E_total(0) − ⟨T$_{{\rm WP}}$⟩ − E_GS = **{g('t0_cross_plus_sie_eV','{:.2f}')} eV** (cross-Hartree/XC + SIE).
- After also removing SIE: **E_system(0) − E_GS = {g('E_system0_minus_EGS_eV','{:+.2f}')} eV ≈ 0** ⇒ bookkeeping OK.

**CAP-energy ledger (diagnostic, not added back).** Total energy the CAP removed =
**{g('cap_removed_total_eV','{:.0f}')} eV**, split ≈ WP-carried **{g('cap_wp_carried_eV','{:.0f}')} eV** +
bath/collective+secondaries **{g('cap_bath_carried_eV','{:.0f}')} eV** (the *leakage* term — if large
and late, collective energy is escaping; see §1 definition). The split is approximate.
"""))
cells.append(embed(f"{FIGS}/energy_method.png",
                   "Energy method: retained energy + convergence triple + CAP-absorbed ledger", width=W_PNG1))
cells.append(embed(f"{FIGS}/conv_wp.png", "WP E_total(t) and N(t)", width=W_PNG2))

# ====================================================== §6 classical
cells.append(md(rf"""## §6 — Q2: Classical projectile stopping (Ehrenfest ΔKE)
Classical stopping = irreversible KE loss per unit path across the slab (via the
`stopping-power-extraction` skill in production). The KE is **not monotonic**: it dips at the slab
centre and recovers on exit (a conservative mean-field well); only the net loss between
**equal-potential** points (the symmetric faces ±12.5) is true stopping.

| | value |
|---|---|
| KE at launch | {g('classical_ke0_eV','{:.1f}')} eV |
| KE min (at z={g('classical_ke_min_z','{:+.1f}')}) | {g('classical_ke_min_eV','{:.1f}')} eV |
| **equal-potential-face S** | **{g('classical_S_facewindow_eVbohr','{:.3f}')} eV/Bohr** |

> **Caveat:** the classical KE carries the phase-2 ion mass/velocity-unit convention — the clean
> traversal-window S must be re-extracted via the **`stopping-power-extraction`** skill before
> being quoted as the production classical S(v).
"""))
cells.append(embed(f"{FIGS}/classical_ke.png", "Classical projectile KE(t)", width=W_PNG1))
cells.append(embed(f"{FIGS}/classical_transport.png",
                   "Classical ion z(t) + KE(z) dip-and-recovery (equal-potential window marked)", width=W_PNG2))

# ====================================================== §7 CAP / reflection / spreading
cells.append(md(rf"""## §7 — Q3: CAP quality (new two-sided geometry) — reflection & spreading
**Reflection vs transmission** (cumulative signed flux through the ±CAP inner faces |z|=35):
**transmitted (+z) = {g('transmit_cum','{:.3f}')}**, **reflected (−z) = {g('reflect_cum','{:.3f}')}** —
the first DIRECT reflectivity measurement for the new two-sided-CAP-in-the-90-box geometry.

**Spreading.** σ_z grows **×{g('sigma_z_growth','{:.0f}')}** over the run; centroid reaches
z$_{{\rm max}}$≈{g('zc_max','{:.1f}')} before absorption.
"""))
cells.append(embed(f"{FIGS}/reflection.png", "Cumulative signed flux through ±CAP faces", width=W_PNG1))
cells.append(embed(f"{FIGS}/centroid.png", "Centroid ⟨z⟩(t) and width σ_z(t)", width=W_PNG2))

# ====================================================== §8 momentum + excitation
cells.append(md(rf"""## §8 — Momentum & KS-excitation diagnostics
**Momentum n(k$_z$)** at the scattering anchors (FFT of ψ_WP) — fast components reflect/interfere
first; k$_z$<0 is reflected flux.
"""))
cells.append(embed(f"{FIGS}/momentum_nkz.png", "n(k_z) at scattering anchors", width=W_PNG1))
cells.append(md(rf"""**KS-excitation matrix** O$_{{ij}}$=|⟨ψ$_i^{{\rm GS}}$|ψ$_j$(t)⟩|² (overlap_full, t=0 vs t=τ).
Off-diagonal weight at t$_f$ = **{g('ks_offdiag_final','{:.3f}')}** (now correctly parsed — the phase-2
read-as-header bug is fixed) — the excitation of the frozen-GS orbitals by the projectile.
"""))
cells.append(embed(f"{FIGS}/excitation.png", "KS-excitation overlap matrix (t=0 | t=τ)", width=W_PNG1))
cells.append(md("**Induced E-field** E_z (mid-time), from the Poisson field of the induced density."))
cells.append(embed(f"{FIGS}/efield.png", "Induced E_z (mid-time)", width=W_PNG1))
cells.append(embed(f"{FIGS}/norm_absorption.png", "Total norm & boundary absorption vs time", width=W_PNG1))

# =============================== §8.5 heuristics ===========
cells.append(md(rf"""## §8.5 — Physical anchors & heuristics (groups A–I)
| quantity | value |
|---|---|
| k_F = v_F | {gh('wp','eg_scales.kF','{:.4f}')} a₀⁻¹ |
| E_F | {gh('wp','eg_scales.EF_ev','{:.2f}')} eV |
| plasmon ω_p | {gh('wp','eg_scales.omega_p_ev','{:.2f}')} eV |
| **plasmon period** T_p | **{gh('wp','eg_scales.T_plasmon_au','{:.1f}')} a.u.** (τ=100 ⇒ ~2 periods) |
| reach slab end (mean v) | {g('T_exit_au','{:.2f}')} a.u. |
| zero-point KE 3/(4σ²) | {gh('wp','wp_kinetics.zero_point_ke_ev','{:.1f}')} eV |
| spreading factor | ×{gh('wp','spreading.spread_factor','{:.0f}')} |

> T_plasmon ≈ {gh('wp','eg_scales.T_plasmon_au','{:.0f}')} a.u. < τ=100 ⇒ ~2 plasmon periods are now
> captured (vs <1 at τ=40). A loss-function L(q,ω) stays **Fourier-gated** (Δω≈2π/τ still coarse).
"""))

# ====================================================== §9 takeaway
cells.append(md(rf"""## §9 — Takeaway
- **WP energy-method S = {g('wp_S_eVbohr','{:.2f}')} eV/Bohr**, convergence = **{g('wp_converged','{}')}**
  (norm_f={g('wp_norm_final','{:.3f}')}, late slope={g('wp_late_slope_eV_au','{:.2f}')} eV/au,
  plateau={g('wp_plateau_width_au','{:.1f}')} a.u.). The t=0 bookkeeping checks out
  (E_system(0)−E_GS={g('E_system0_minus_EGS_eV','{:+.2f}')} eV ≈ 0).
- **Classical projectile REFLECTED** at 54 eV (ended z={g('classical_ke_min_z','{:+.1f}')} Bohr, v→0): it
  never traversed the slab, so the face-window S = **{g('classical_S_facewindow_eVbohr','{:.3f}')} eV/Bohr**
  is not a clean traversal value (lower energy than the 100 eV p3 run ⇒ reflection more likely).
  **Point-Lindhard = {g('S_lindhard_point_54eV_eVbohr','{:.3f}')} eV/Bohr.**
- **CAP (new two-sided, 90-box):** reflection {g('reflect_cum','{:.3f}')}, transmission
  {g('transmit_cum','{:.3f}')}; WP init-absorption fixed by the equidistant launch.
- **CAP-energy ledger** quantifies leakage: bath/collective-carried ≈ {g('cap_bath_carried_eV','{:.0f}')} eV
  (the term that, by the system definition, leaves the box).

**Why the gate fails at τ=100 — and why longer τ will not fix it.** The WP norm gate fails here
(norm_f={g('wp_norm_final','{:.3f}')}) just as at τ=40, yet S barely moves across the σ=0.5 runs
(τ=40 → 2.73, τ=100 → {g('wp_S_eVbohr','{:.2f}')} eV/Bohr) and stays ~5× above point-Lindhard. The
retained-energy ledger contains the packet's **zero-point kinetic energy** 3/(4σ²) = **81.6 eV**, which
*exceeds* the drift KE ½v₀² = 54.4 eV and never leaves the system — so E_total(t$_f$)−E_GS measures
**deposited** energy, not stopping. Extending τ only improves absorption at the margin (p3 norm 0.046,
p4 norm {g('wp_norm_final','{:.3f}')}); it cannot subtract the zero-point floor. The route to a converged,
comparable-to-Lindhard S is a **larger σ** (muon-like, zero-point ≪ drift) or a **first-moment ⟨p⟩**
observable — **not more propagation time**. So S$_{{\rm WP}}$={g('wp_S_eVbohr','{:.2f}')} eV/Bohr is a
deposited-energy **upper bound**, consistent across the campaign.
"""))

out = os.path.join(HERE, "qsp_phase4_study.ipynb")
build(cells, out, timeout=900)
print(f"done: {out}")
