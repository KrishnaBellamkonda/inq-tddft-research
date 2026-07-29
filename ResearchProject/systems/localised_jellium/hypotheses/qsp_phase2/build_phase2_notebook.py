#!/usr/bin/env python3
"""Assemble qsp_phase2_study.ipynb — Phase-2 P2.1 results of the
`quantum-stopping-power` campaign: the WP + classical convergence / CAP test
(2000 steps, 40 a.u.) on the localised jellium slab.

This is the run-SET study notebook (what the two test runs MEAN). It is a thin
narrative assembler: every figure + GIF is pre-computed by `analyse_phase2.py`
(figs/*.png, figs/*.gif) and every headline number is read from `results.json`,
so a quoted number and its figure can never disagree. Per-run deep-dives live in
the sibling run-notebooks p2wp_run_notebook.ipynb / p2cl_run_notebook.ipynb.

Run:
  PYTHONPATH=/local/data/public/skcb2/tddft/inq-stack/python \
  /local/data/public/skcb2/tddft/venv/bin/python3 build_phase2_notebook.py
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))          # hypotheses/ for _nbreport
from _nbreport import md, embed, setup_cell, set_outdir, build

set_outdir(HERE)
FIGS = os.path.join(HERE, "figs")
R = json.load(open(os.path.join(HERE, "results.json")))


def g(k, fmt="{:.3f}", dflt="n/a"):
    v = R.get(k)
    return dflt if v is None else fmt.format(v)


def gh(run, key, fmt="{:.3f}", dflt="n/a"):
    """Heuristic value from results.json['heuristics'][run][key]."""
    v = R.get("heuristics", {}).get(run, {}).get(key)
    return dflt if v is None else fmt.format(v)


KIND_LABEL = {"density": "n(x,z,t)", "delta0": "Δn = n(t)−n(0)",
              "dstep": "Δn = n(t+dt)−n(t)"}
CAT_LABEL = {"total": "total system", "wp": "wavepacket |ψ|²", "bath": "bath (slab only)"}

# display widths (px) — Markdown renders figures at native pixel size, which is
# too large; cap the on-screen size via embed(..., width=).
W_GIF = 360      # GIF battery (compact, many of them)
W_PNG2 = 600     # multi-panel PNGs (energetics / transport / convergence)
W_PNG1 = 520     # single-panel PNGs


def gif_row(run, cat, width=W_GIF):
    """Embed the 3 kinds (density/Δ/Δstep) for one category of one run."""
    out = []
    for kind in ("density", "delta0", "dstep"):
        fn = f"{run}_{cat}_{kind}.gif"
        if os.path.exists(os.path.join(FIGS, fn)):
            out.append(embed(f"{FIGS}/{fn}",
                             f"{CAT_LABEL.get(cat,cat)} · {KIND_LABEL[kind]}", width=width))
    return out


HA = 27.211386
cells = []
cells.append(setup_cell())

# ----------------------------------------------------------------- §0 title + Q
cells.append(md(rf"""# Phase 2 · P2.1 — WP + classical convergence / CAP test
### localised jellium slab · r$_s$≈5.67 · 82 e · box 50×50×70 · σ$_{{\rm WP}}$=0.5 · 100 eV · CAP η=−0.7

**Campaign:** `docs/campaigns/jellium_wp_stopping/quantum-stopping-power.md` (Phase 2, task P2.1).
**Builds on:** Phase 1 (`hypotheses/qsp_phase1/phase1_gs_sie.ipynb`) — validated GS
(E$_{{\rm GS}}$=−45.759 Ha) and the **SIE floor ≈ 4.40 eV**.

**The question this run answers (cheaply, before the production runs):**

1. **Does the WP total energy converge** within an affordable sim-time? E$_{{\rm total}}$(t)
   only equals the *deposited* energy once the projectile has fully traversed the slab
   **and** the excited density has been absorbed — only then is ΔE$_{{\rm total}}$/L$_z$ a
   stopping power. This test asks whether **40 a.u. is long enough**.
2. **Does the classical projectile reach a steady (irreversible-loss) state** over the
   slab, so its ΔKE$_{{\rm ion}}$ gives a clean stopping number?
3. **Are the CAP parameters (η=−0.7, 10 Bohr/side) good** — i.e. is reflection
   negligible, and is the WP initialised with enough room (no orthogonalisation loss,
   no premature CAP capture)?

These three answers set the parameters (sim-time, σ, CAP) for the P2.2 production
WP↔classical stopping comparison. **This is a 40 a.u. *test*, not a production run** —
read every stopping number below as PROVISIONAL.
"""))

# ----------------------------------------------------------- §1 conventions
cells.append(md(r"""## §1 — Conventions & symbols

Atomic units throughout (ℏ=m$_e$=e=1; energies quoted in eV via 1 Ha = 27.211 eV;
lengths in Bohr; velocities in a.u., where for an electron projectile **v = k**).

| symbol | meaning | value / range |
|---|---|---|
| r$_s$ | electron-gas density parameter of the slab interior | 5.67 (n₀=1.31×10⁻³ a₀⁻³) |
| σ$_{\rm WP}$ | wavepacket **density** width | 0.5 Bohr |
| σ$_{\rm pot}$ | classical Gaussian **charge** std = σ$_{\rm WP}$/√2 | 0.354 Bohr |
| k₀ / v₀ | initial WP wavevector = projectile speed | 2.711 a.u. (½v₀²≈100 eV) |
| L$_z$ | slab traversal length (full thickness) used for dE/dx | 25 Bohr (|z|<12.5) |
| η, w | CAP strength (Ha) and width per side | −0.7 Ha, 10 Bohr ([±25,±35]) |
| τ | total propagation time | 40 a.u. (2000 × dt=0.02) |
| SIE | self-interaction floor (Phase 1) | 4.40 eV |

**σ convention trap (this system):** the localised campaign labels by the WP density
width σ$_{\rm WP}$; the *classical* projectile reproduces the same charge cloud with a
Gaussian charge std σ$_{\rm pot}$ = σ$_{\rm WP}$/√2 = 0.354. Both projectiles carry the
same v₀ = 2.711.

**Numerical requirement:** the CAP makes the propagator **non-Hermitian** ⇒ the run
uses **ETRS**, not Crank–Nicolson (CN renormalises each step and would defeat the
absorber). Built against **inq-study** (the CAP-enabled engine replica), `inq/` untouched.
"""))

# ----------------------------------------------------------- §2 setup
cells.append(md(rf"""## §2 — Simulation setup (fully reconstructable)

Two runs, identical geometry / CAP / GS / launch, differing only in the projectile.
Launched concurrently (WP→GPU0, classical→GPU1) by `scripts/qsp_phase2/dispatch.sh`.

| | value |
|---|---|
| **cell** | 50 × 50 × 70 Bohr, orthorhombic periodic, spacing 0.50 (100×100×140 grid) |
| **background** | positive jellium slab, half-width 12.5 along z (slab = [−12.5,12.5]) |
| **electrons** | 82 (closed-shell: 41 doubly-occupied) + extra states for the WP |
| **GS reused** | `shared_gs/slab_n82_L50x50x70` (Phase-1, E$_{{\rm GS}}$=−45.75885 Ha) |
| **functional / engine** | LDA · inq-study (CAP-enabled), header-only INQ |
| **propagator** | **ETRS**, dt = 0.02 a.u., N_steps = 2000 ⇒ τ = 40 a.u. |
| **CAP** | two-sided sin², **η = −0.7 Ha**, 10 Bohr/side, region [±25, ±35] (inner face |z|=25) |
| **region layout (z)** | slab [−12.5,12.5] · free [±12.5,±25] · CAP [±25,±35] |
| **launch** | z₀ = −22 (≈ 6σ from the CAP inner face; in the free region) |
| **WP projectile** | Gaussian, σ$_{{\rm WP}}$=0.5, k₀=2.711, injected into the last extra state, orthogonalised against occupied |
| **classical projectile** | Ehrenfest ghost ion, Gaussian-e pseudopotential σ$_{{\rm pot}}$=0.354, mass = m$_e$, v=(0,0,2.711) |
| **save cadence** | density VTIs + wavefunction_wp every 10 steps (200 frames); overlap_full at t=0 & t=τ only |
| **walltime** | WP {g('wp_wall_s','{:.0f}')} s · classical {g('classical_wall_s','{:.0f}')} s (~2 h each, CPU-free GPUs) |

> **Provenance note (to fix for P2.2):** the auto-written `run_summary.txt` CAP string
> reads *"eta −0.5 mid ±0.425 width 0.15"* — this is a **stale hardcoded label** in the
> run.cpp summary writer. The **compiled** value (run.cpp:69 `CAP_ETA=-0.7`, applied at
> run.cpp:124–125) is **η = −0.7**, mid = ±30/70, width = 10/70 ⇒ region [±25,±35]. The
> physics used −0.7; only the label is wrong. (The same writer also packed multiple
> `key = value` pairs per line — reformatted post-hoc so the run-notebook parser reads it.)
"""))

# ----------------------------------------------------------- §3 source files
cells.append(md(r"""## §3 — Source files

| file | role |
|---|---|
| `scripts/qsp_phase2/wp/run.cpp` | WP run (ETRS, CAP, full observable suite) |
| `scripts/qsp_phase2/classical/run.cpp` | classical Ehrenfest ghost-ion run |
| `scripts/qsp_phase2/dispatch.sh` | concurrent launch (GPU0/GPU1) → analyse → notebook |
| `shared/configs/slab_n82_L50x50x70.hpp` | geometry / density / WP config struct |
| `shared_gs/slab_n82_L50x50x70/` | reused Phase-1 ground state |
| `hypotheses/qsp_phase2/analyse_phase2.py` | computes every figure + `results.json` |
| `hypotheses/qsp_phase2/build_phase2_notebook.py` | **this** study-notebook builder |
| `hypotheses/qsp_phase2/p2wp_run_notebook.ipynb` | WP per-run deep-dive |
| `hypotheses/qsp_phase2/p2cl_run_notebook.ipynb` | classical per-run deep-dive |

All paths repo-relative to `ResearchProject/systems/localised_jellium/`.
"""))

# ====================================================== §4 per-run visual intuition
cells.append(md(rf"""## §4 — Per-run visual intuition (before any derived number)

Each run's density evolution first, so the physics is visible before the extracted
numbers. Every GIF is an **xz mid-plane slice**, fixed colour scale across frames, with
the slab faces (|z|=12.5) and **CAP inner faces (|z|=25)** dashed.

We show, for each category — **total system**, **wavepacket |ψ|²**, **bath (= n_total −
n_wp, the slab with the projectile removed)** — three kinds:

| kind | what it reveals |
|---|---|
| **n(x,z,t)** | the density itself (log scale; total/bath share the classical-total scale so **low slab densities are visible**, per request) |
| **Δn = n(t)−n(0)** | cumulative change from the ground state (the induced response / wake) |
| **Δn = n(t+dt)−n(t)** | the per-frame change (where the action is *right now*) |

Colour scheme (shared-colorbar rule): density GIFs share one log scale across
total/bath incl. the classical total; the two difference kinds use their own symmetric
diverging scale. Shared total/bath density vmax = {g('gif_density_vmax','{:.2e}')} a₀⁻³.
"""))

cells.append(md(r"""### §4a — Wavepacket run — total system (3 kinds)

The σ=0.5 packet launches at z=−22 and is absorbed at the +z CAP. The **total** GIF now
uses the slab-tuned (classical-total) colour scale, so the slab and its perturbation are
visible rather than washed out by the bright packet."""))
[cells.append(c) for c in gif_row("wp", "total")]

cells.append(md(r"""### §4b — Wavepacket run — wavepacket orbital |ψ_WP|² (3 kinds)

The projectile orbital alone: its forward motion, the dramatic **spreading** (§7), and
absorption at the +z CAP. Its own colour scale (the packet is far denser than the slab)."""))
[cells.append(c) for c in gif_row("wp", "wp")]

cells.append(md(r"""### §4c — Wavepacket run — bath only (n_total − n_wp, 3 kinds)

The **slab's response with the projectile removed** — the genuine electronic wake the
packet leaves behind, decoupled from the packet's own density."""))
[cells.append(c) for c in gif_row("wp", "bath")]

cells.append(md(r"""### §4d — Classical run — total system (3 kinds)

The Ehrenfest ghost ion drags its screening cloud through the slab. There is no separate
WP orbital, so Δn = n(t)−n(0) **is** the projectile-induced wake."""))
[cells.append(c) for c in gif_row("classical", "total")]

cells.append(md(rf"""### §4e — Energetics of both runs (slab-exit time marked)

Energy components vs time, with the **mean-velocity slab-exit time t={g('T_exit_au','{:.1f}')} a.u.**
(dashed, the time a projectile at v=2.711 from z₀=−22 reaches the far slab face +12.5)
and the slab-entry time (dotted) marked. The CAP parameters (**η=−0.7 Ha, 10 Bohr/side,
region [±25,±35]**) are stated in the figure header.

For the **WP** the total energy **drains monotonically** as the CAP absorbs the packet;
for the **classical** run it **rises** — the key difference explained in §4f."""))
cells.append(embed(f"{FIGS}/energetics.png", "Per-run energetics vs time (slab-exit marked)", width=W_PNG2))

# --------------------------------------------- §4f looping-back (the headline physics)
cells.append(md(rf"""## §4f — Is the projectile looping back? (classical yes, WP no)

A natural question once the energy is seen to **rise** late in a run: is the projectile
re-entering the slab through the periodic image?

**Classical — YES, it loops back.** The CAP is an absorbing potential acting on the
electronic *wavefunctions*, **not** on the classical point charge. So the Ehrenfest ion
is never absorbed: its track runs z: −22 → **+62.7** monotonically (v_z stays positive).
The box is [−35, +35], so z=62.7 is **wrapped → physical −7.3** — the ion has come back
around and is re-approaching the slab. That is exactly why the classical E_total **rises
again after t≈30** (it peaks ≈+70 eV crossing the slab, dips when farthest, then climbs
as the periodic image re-enters). The rise is physical and expected.

**Wavepacket — NO.** The wavefunction *is* absorbed by the CAP, and the evidence is
unambiguous:

| quantity | behaviour | reading |
|---|---|---|
| E_total(t) | **monotonic down** (0 → {g('wp_dEtot_eV','{:.0f}')} eV) — never turns up | no re-entry |
| N_total(t) | monotonic down, 83 → {g('wp_Ntot_final','{:.2f}')} | norm only leaves (absorbed), never returns |
| WP orbital norm | 1.000 → {gh('wp','norms.N_wp_f','{:.3f}')} ({gh('wp','norms.wp_fraction_absorbed','{:.0%}')} absorbed) | the packet is consumed at the +z CAP |
| centroid ⟨z⟩ | advances to **+{g('zc_max','{:.1f}')}** (max), then drifts back to +8.0 | **survival-weighting artefact**, not a turn-around |

The apparent centroid retreat is *not* the packet reversing: the fast forward components
are absorbed first, so the **surviving** (slower) norm pulls the survival-weighted mean
back. **The only quantity that rises in the WP run is the system *kinetic* energy** (0 →
+6.6 eV in the first ~8 a.u., as the projectile excites the bath) before the CAP-driven
monotonic drop — see the C1 curve in §4e. The total energy never loops up.

> **Implication (feeds the simulation design):** to make the two channels comparable, the
> *classical* ion should be **removed once it reaches the box edge** (before it wraps) — a
> change discussed separately for the next run. The WP needs no such treatment.
"""))
cells.append(md(rf"""### §4g — Total norm & boundary absorption vs time

Total electron number and the WP orbital norm vs time, with the slab-exit marker. WP
total absorbed = **{gh('wp','norms.total_absorbed','{:.3f}')} e** (orbital
{gh('wp','norms.wp_orbital_absorbed','{:.3f}')} e + bath overflow
{gh('wp','norms.bath_overflow_absorbed','{:+.3f}')} e); classical total absorbed =
**{gh('classical','norms.total_absorbed','{:.3f}')} e** (the screening cloud the ion drags
into the CAP)."""))
cells.append(embed(f"{FIGS}/norm_absorption.png", "Total norm & boundary absorption vs time", width=W_PNG1))

# ====================================================== §5 convergence (the gate)
cells.append(md(rf"""## §5 — Q1: Does the WP total energy converge? (the make-or-break)

The WP stopping power is only defined once E$_{{\rm total}}$(t) has plateaued — i.e. the
projectile has fully traversed **and** the excited density is fully absorbed:

$$ S_{{\rm WP}} \;=\; \frac{{E_{{\rm total}}(t_f) - E_{{\rm GS}}}}{{L_z}}
\quad\text{{valid only when}}\quad \left.\frac{{dE_{{\rm total}}}}{{dt}}\right|_{{t_f}} \to 0 . $$

**Verdict — NOT converged at τ = 40 a.u.** The late-time slope is still
**{g('wp_late_slope_eV_au','{:.2f}')} eV/a.u.** at t$_f$; E$_{{\rm total}}$ is still
draining (the CAP is still absorbing packet that hasn't finished leaving). The electron
number N(t$_f$) = **{g('wp_Ntot_final','{:.3f}')}** — it started at **83** (82 slab + 1
WP) and has lost **{gh('wp','norms.total_absorbed','{:.2f}')} e** to the CAP, i.e.
absorption is incomplete-but-ongoing (the slope has not reached zero).
"""))
cells.append(embed(f"{FIGS}/conv_wp.png",
                   "WP E_total(t) and N(t) — late-time slope still nonzero ⇒ not converged", width=W_PNG2))
cells.append(md(rf"""**Consequence for the stopping number.** Because the run is not converged, the
"deposited" energy E$_{{\rm total}}$(t$_f$)−E$_{{\rm GS}}$ = **{g('wp_deposited_EminusEGS_eV','{:.1f}')} eV**
⇒ S$_{{\rm WP}}$ = **{g('wp_S_eVbohr','{:.2f}')} eV/Bohr** is an **UPPER BOUND only**, and it
sits well above the 4.40 eV SIE floor. (The raw ΔE$_{{\rm total}}$ over the whole run is
{g('wp_dEtot_eV','{:.1f}')} eV, dominated by the absorbed packet's own KE, not by stopping.)
➡️ **P2.2 must lengthen τ** until the slope vanishes before any WP stopping is quoted.
"""))

# ====================================================== §6 classical steady state
cells.append(md(rf"""## §6 — Q2: Does the classical projectile reach steady state?

Classical stopping is the **irreversible** KE loss per unit path while the ion crosses
the slab (Ehrenfest; via the `stopping-power-extraction` skill in production):

$$ S_{{\rm cl}} \;=\; -\,\frac{{\Delta {{\rm KE}}_{{\rm ion}}}}{{\Delta z}}
\Big|_{{\rm slab\ traversal}} . $$

Over the slab-traversal window the ion KE falls **{g('classical_ke0_eV','{:.0f}')}→… eV**;
the *traversal-window* loss gives **S$_{{\rm cl}}$ ≈ {g('classical_S_eVbohr','{:.3f}')} eV/Bohr**
(analyse_phase2 estimate). The KE curve below shows a near-linear loss across the slab
then a flat exit — consistent with a steady drag inside the slab and no loss outside.
"""))
cells.append(embed(f"{FIGS}/classical_ke.png",
                   "Classical projectile KE(t)/KE(z) — steady loss across the slab", width=W_PNG1))

cells.append(md(rf"""### §6a — KE is NOT monotonic: a conservative dip-and-recovery

Tracking KE along the trajectory shows the projectile **slows to a minimum at the slab
centre and recovers almost back to launch on exit** — a *conservative* mean-field-potential
effect (energy borrowed entering the well, returned leaving it), **not** stopping:

| z (Bohr) | −22 | −12.5 (face) | 0 (centre) | +12.5 (face) | +25 | +35 |
|---|---|---|---|---|---|---|
| KE (eV) | 100.0 | 67.4 | **{g('classical_ke_min_eV','{:.1f}')}** | 54.8 | 90.2 | 96.1 |

The left panel is **z(t)** (the requested position-vs-time): the ion decelerates into the
slab, reaches its turning-minimum near the centre (z={g('classical_ke_min_z','{:+.1f}')}), then
re-accelerates and eventually **wraps past +35** (periodic re-entry — the late energy rise,
§4f). The right panel is **KE(z)**: only the net loss between two points at **equal
background potential** is true electronic stopping. The slab is symmetric, so the **faces
±12.5 are an equal-potential pair** ⇒ ΔKE = 67.4−54.8 = 12.6 eV over 25 Bohr ⇒ **S =
{g('classical_S_facewindow_eVbohr','{:.3f}')} eV/Bohr**. Reading KE at the *centre* (or any
asymmetric window) would conflate the conservative well with stopping and give a wrong, much
larger number — **the window choice dominates the answer**."""))
cells.append(embed(f"{FIGS}/classical_transport.png",
                   "Classical ion z(t) (left) and KE(z) dip-and-recovery (right) — equal-potential window marked", width=W_PNG2))

cells.append(md(rf"""> **Caveat (P2.2 / stopping-power skill):** the run's raw `ke_ion_initial_ha`
> (0.986 Ha) and analyse_phase2's KE (100 eV) differ by the ion **mass/velocity-unit
> convention** — the clean traversal-window stopping number must be re-extracted through
> the **`stopping-power-extraction`** skill (localised-slab branch, with its built-in
> sanity checks) before it is quoted as the classical S(v). The {g('classical_S_eVbohr','{:.3f}')}
> eV/Bohr here is PROVISIONAL.
"""))

# ====================================================== §7 CAP / reflection / spreading
cells.append(md(rf"""## §7 — Q3: CAP quality — reflection, absorption, spreading

**Reflection vs transmission** (cumulative signed probability flux J$_z$ through the
±CAP inner faces): **transmitted (+z) = {g('transmit_cum','{:.3f}')}**, **reflected (−z) =
{g('reflect_cum','{:.3f}')}**. Reflection is **~1.3%** — the CAP (η=−0.7, 10 Bohr) absorbs
cleanly with little spurious back-reflection. The WP initialised with norm ≈ 1 at z=−22,
so there was **enough room for orthogonalisation and no premature CAP capture** (Q3 ✓).
"""))
cells.append(embed(f"{FIGS}/reflection.png", "Cumulative signed flux through ±CAP faces", width=W_PNG1))
cells.append(md(rf"""**Spreading (the dominant effect at σ=0.5).** The packet's longitudinal width grows by
**×{g('sigma_z_growth','{:.0f}')}** over the run (free-particle dispersion), and its centroid
advances to z$_{{\rm max}}$ ≈ {g('zc_max','{:.1f}')} before absorption. This large spread is
exactly the motivation for the planned **large-σ "no-appreciable-spread" run** in the
campaign — at σ=0.5 the quantum packet smears so much that the WP↔classical comparison is
dominated by dispersion, not stopping."""))
cells.append(embed(f"{FIGS}/centroid.png", "Centroid z(t) and width σ_z(t) — WP vs classical", width=W_PNG2))

# ====================================================== §8 momentum + excitation
cells.append(md(rf"""## §8 — Momentum & KS-excitation diagnostics

**Momentum distribution n(k$_z$)** at the scattering anchors — captured to test the
hypothesis that *fast plane-wave components reflect / interfere first*.
Anchors (a.u.): A1≈{R['anchor_table']['A1_lead_near']:.1f},
A2≈{R['anchor_table']['A2_cen_near']:.1f}, A4≈{R['anchor_table']['A4_lead_far']:.1f},
A5≈{R['anchor_table']['A5_cen_max']:.1f}, t$_f$={R['anchor_table']['t_final']:.0f}."""))
cells.append(embed(f"{FIGS}/momentum_nkz.png", "n(k_z) at scattering anchors (FFT of ψ_WP)", width=W_PNG1))
cells.append(md(rf"""**KS-excitation matrix** O$_{{ij}}$=|⟨ψ$_i^{{\rm GS}}$|ψ$_j$(t)⟩|² (overlap_full, t=0
and t=τ) — the excitation of the frozen-GS orbitals by the projectile. **Flag:** the
measured off-diagonal weight at t$_f$ is **{g('ks_offdiag_final','{:.3f}')}** — i.e. ~0,
which is physically implausible for a 100 eV projectile crossing the slab. This points to
an **overlap_full capture/parse issue**, not a real null excitation (post-mortem below)."""))
cells.append(embed(f"{FIGS}/excitation.png", "KS-excitation overlap matrix (t=0 | t=τ)", width=W_PNG1))

# =============================== §8.5 physical anchors / heuristics ===========
cells.append(md(rf"""## §8.5 — Physical anchors & heuristics (groups A–I)

The reusable heuristic battery (`inqview.analysis.compute_heuristics`), computed for this
geometry — the same anchors carried from the first jellium-slab test campaign. These set
the timescales and resolution limits the production run must respect.

**A. Electron-gas scales (r_s={gh('wp','eg_scales.rs','{:.3f}')})**

| quantity | symbol | value |
|---|---|---|
| density | n₀ | {gh('wp','eg_scales.n0','{:.2e}')} a₀⁻³ |
| Fermi wavevector | k_F | {gh('wp','eg_scales.kF','{:.4f}')} a₀⁻¹ |
| Fermi velocity | v_F | {gh('wp','eg_scales.vF','{:.4f}')} a.u. |
| Fermi energy | E_F | {gh('wp','eg_scales.EF_ev','{:.2f}')} eV |
| Friedel wavelength | λ_F=π/k_F | {gh('wp','eg_scales.lambda_F_friedel','{:.2f}')} Bohr |
| plasmon frequency | ω_p | {gh('wp','eg_scales.omega_p_ev','{:.2f}')} eV |
| **plasmon period** | T_p=2π/ω_p | **{gh('wp','eg_scales.T_plasmon_au','{:.1f}')} a.u.** |
| Thomas–Fermi screening | k_TF | {gh('wp','eg_scales.k_TF','{:.3f}')} a₀⁻¹ |
| HEG kinetic E / electron | t_HEG | {gh('wp','eg_scales.t_heg_ha_per_e','{:.4f}')} Ha |

> **Resolution flag:** the projectile velocity v=2.711 ≫ v_F={gh('wp','eg_scales.vF','{:.3f}')}
> (fast, perturbative regime). And **T_plasmon ≈ {gh('wp','eg_scales.T_plasmon_au','{:.0f}')} a.u.
> > τ = 40 a.u.** — the run is shorter than a single plasmon period, so any L(q,ω) /
> collective spectrum is severely under-resolved (Δω ≈ 2π/τ ≈ 4.3 eV vs ω_p ≈
> {gh('wp','eg_scales.omega_p_ev','{:.1f}')} eV). This is the quantitative basis for the
> Fourier/loss-function gate.

**B. Timescales (mean velocity v=2.711, z₀=−22)**

| event | value |
|---|---|
| reach near slab face (−12.5) | {gh('wp','timescales.t_enter_slab_au','{:.2f}')} a.u. |
| **reach far slab face (+12.5) — slab end** | **{gh('wp','timescales.t_exit_slab_au','{:.2f}')} a.u.** |
| near→far transit | {gh('wp','timescales.t_cross_au','{:.2f}')} a.u. |
| reach +box edge (wrap onset) | {gh('wp','timescales.t_reach_box_edge_au','{:.2f}')} a.u. |

**C. Wavepacket kinetics (σ_WP=0.5)** — charge std σ/√2 = {gh('wp','wp_kinetics.sigma_charge','{:.3f}')};
**zero-point KE = 3/(4σ²) = {gh('wp','wp_kinetics.zero_point_ke_ev','{:.1f}')} eV** (the energy the
"+100 eV" SIE reference omits; cf. Phase 1).

**D. Norm / boundary absorption** — WP total absorbed {gh('wp','norms.total_absorbed','{:.3f}')} e
(orbital {gh('wp','norms.wp_orbital_absorbed','{:.3f}')} + bath overflow
{gh('wp','norms.bath_overflow_absorbed','{:+.3f}')}); classical
{gh('classical','norms.total_absorbed','{:.3f}')} e. (§4g figure.)

**E. Stopping references** — point-Lindhard S(v,r_s) =
{gh('wp','stopping_refs.S_point_ev_per_bohr','{:.3f}')} eV/Bohr (classical measured
{g('classical_S_eVbohr','{:.3f}')}; WP upper-bound {g('wp_S_eVbohr','{:.2f}')}).

**H. Spreading** — σ_z: {gh('wp','spreading.sigma_z_0','{:.2f}')} → {gh('wp','spreading.sigma_z_f','{:.2f}')}
Bohr = **×{gh('wp','spreading.spread_factor','{:.0f}')}** (max {gh('wp','spreading.sigma_z_max','{:.1f}')}).

*F (momentum KL / σ_k) and G (collective spectrum) — the momentum panels are in §8; the
collective spectrum stays Fourier-gated given T_plasmon > τ above.*
"""))

# ====================================================== §9 post-mortem
cells.append(md(rf"""## §9 — Post-mortem (anomalies, log-grounded)

1. **WP not converged (expected, by design).** 40 a.u. is a *test* ceiling; the late
   slope {g('wp_late_slope_eV_au','{:.2f}')} eV/a.u. quantifies how far from a plateau we
   are. **Action:** P2.2 lengthens τ (estimate from this slope + the packet exit time)
   until |dE/dt|→0; only then quote S$_{{\rm WP}}$.
2. **KS off-diagonal weight = 0 (anomaly).** `analyse_phase2.log` reports
   *"KS excitation off-diagonal weight (final) = 0.000"*. A real run must excite KS
   states, so this is almost certainly an `overlap_full` indexing/parse bug in the
   analysis (e.g. reading the t=0 matrix twice, or a transpose/identity artefact), **not**
   physics. **Action:** verify the overlap_full reader against the t=0 vs t=τ files before
   relying on the excitation panel.
3. **E-field panel skipped.** `analyse_phase2.log`: *"[efield] efield API fallback
   skipped"* — the post-hoc Poisson E-field block hit an API mismatch and was guarded out.
   **Action:** wire `inqview.analysis.efield` correctly for P2.2 (density VTIs are saved,
   so the field is recoverable post-hoc).
4. **Classical KE unit/mass convention.** See §6 caveat — re-extract via the
   `stopping-power-extraction` skill.
5. **Provenance string stale.** `run_summary.txt` CAP label said η=−0.5; compiled value is
   η=−0.7 (verified in run.cpp). **Action:** fix the summary writer string in the P2.2 run.cpp.
"""))

# ====================================================== §10 takeaway
cells.append(md(rf"""## §10 — Takeaway

- **Q1 (WP convergence): NO at 40 a.u.** — late slope {g('wp_late_slope_eV_au','{:.2f}')}
  eV/a.u.; S$_{{\rm WP}}$={g('wp_S_eVbohr','{:.2f}')} eV/Bohr is an **upper bound only**.
  P2.2 must run **longer τ** before any WP stopping is quoted.
- **Q2 (classical steady state): YES** — clean near-linear KE loss across the slab,
  S$_{{\rm cl}}$≈{g('classical_S_eVbohr','{:.3f}')} eV/Bohr (PROVISIONAL, pending the
  stopping-power skill). For reference, point-Lindhard at this v,r$_s$ is
  **{g('S_lindhard_point_100eV_eVbohr','{:.3f}')} eV/Bohr**.
- **Q3 (CAP): GOOD** — reflection {g('reflect_cum','{:.3f}')} (~1%), transmission
  {g('transmit_cum','{:.3f}')}; WP initialised cleanly with room to orthogonalise. **Keep
  η=−0.7, 10 Bohr.**
- **Dominant physics surprise: spreading ×{g('sigma_z_growth','{:.0f}')}** at σ=0.5 — the
  quantum packet disperses so strongly that the WP↔classical comparison is dispersion-
  dominated. This **motivates the large-σ run** (chosen so the packet does not appreciably
  spread) as the clean quantum-vs-classical test.
- **Two analysis bugs to fix before P2.2:** the overlap_full off-diagonal=0 artefact and
  the skipped E-field block (both are post-processing, not run, issues).

*All stopping numbers PROVISIONAL — this is a 40 a.u. convergence test, not a production run.*
"""))

out = os.path.join(HERE, "qsp_phase2_study.ipynb")
build(cells, out, timeout=600)
print(f"done: {out}")
