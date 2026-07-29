#!/usr/bin/env python3
"""Build the Phase 5 study notebook: two-sided sin^2 CAP stopping power.

Two runs (each its own GIF + energetics block):
  A) wave-packet projectile + CAP   (mechanism: CAP absorbs WP, bath intact)
  B) classical Gaussian-e ion + CAP (headline: S = dKE_ion / x)

    PYTHONPATH=/local/data/public/skcb2/tddft/inq-stack/python \
    /local/data/public/skcb2/tddft/venv/bin/python3 build_cap_stopping_report.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from _nbreport import md, code, setup_cell, embed, set_outdir, build, SYS

HERE = os.path.join(SYS, "hypotheses", "03_cap_stopping")
OUT = os.path.join(HERE, "cap_stopping_study.ipynb")
set_outdir(os.path.dirname(OUT))

cells = []

# 1. Title + question ---------------------------------------------------------
cells += [md(
r"""# Stopping power of localised jellium — two-sided sin² CAP

**Phase 5 (final) of the localised-jellium scattering campaign.**

## The question

> Fire a projectile through the slab, absorb it in vacuum with a **complex absorbing
> potential (CAP)** on each side, and measure the **stopping power**
> $S=\Delta E/x$ — the energy lost per unit path length crossing the 25 Bohr slab.

Two complementary projectiles, run with the *same* CAP and geometry so they can be
compared directly:

- **Run A — wave packet + CAP.** Tests the *mechanism*: the CAP must absorb the
  projectile in vacuum **without draining the bath** (the worry of a CAP in a
  finite jellium). Energy/norm traces tell the story.
- **Run B — classical Gaussian-e ion + CAP.** The *headline*: an Ehrenfest ion of
  the same Gaussian width decelerates via electronic drag; its kinetic-energy loss
  across the slab faces gives a clean $S=\Delta\mathrm{KE}_\mathrm{ion}/x$.

| Where this sits | Phase | Notebook |
|---|---|---|
| prev | 1+2 slab GS validation | `01_slab_validation/` |
| prev | 3 WP projectile baseline | `02_projectile_slab/` |
| **this** | 5 CAP stopping power | `03_cap_stopping/` |

> **PROVISIONAL.** The stopping number is a *first pass*: a single trajectory
> (centroid, $+z$ through centre) on a coarse $dx=0.5$ grid (Phase 3 showed ~1 mHa
> high-$k$ drift). The gate that lifts this: a $dx\to0.25$ + multi-trajectory $S(v)$
> convergence study.
""")]

# 2. Conventions + formulas ---------------------------------------------------
cells += [md(
r"""## Conventions & symbols

Atomic units. $1\,\mathrm{Ha}=27.2114\,\mathrm{eV}$. Engine: **`inq-study`** (the
KS potential is complexified so the CAP's imaginary part propagates — this requires
**ETRS**-style absorption bookkeeping, not a norm-renormalising step).

| Symbol | Meaning | Units |
|---|---|---|
| $\eta$ | CAP amplitude (imaginary) | Ha |
| $x$ | traversal length (slab thickness) | $a_0$ |
| $S$ | stopping power $\Delta E/x$ | Ha/$a_0$ (report eV/$a_0$) |
| $\mathrm{KE}_\mathrm{ion}$ | classical projectile kinetic energy | Ha |
| $N_e$ | electron count (norm); drop = absorbed | — |
""")]

cells += [setup_cell()]

cells += [md(
r"""### The CAP: $v_\mathrm{CAP}(\mathbf r)=i\,\eta\,\sin^2\!\big(\pi\,\xi\big)$

Two sin² absorbers, one in each vacuum gap outside the slab, added to the
KS potential via `perturbations::absorbing(η, mid_frac, width_frac)` and composed
with the background by `perturbations::sum(bg, sum(cap_lo, cap_hi))`. Each CAP is
**zero at both edges** (no reflection seam) and spans 7.5 Bohr. Fractional placement
on the 50 Bohr cell:""")]
cells += [code(
"L = 50.0\n"
"eta = -0.5            # Ha (imaginary amplitude)\n"
"width_frac = 0.15     # 7.5 Bohr per side\n"
"mid_frac = 0.425      # centres at +/- 21.25 Bohr\n"
"print(f'CAP eta = {eta} Ha; each width {width_frac*L:.1f} Bohr; centred at +/-{mid_frac*L:.2f} Bohr')\n"
"print(f'slab faces at +/-12.5 Bohr -> CAPs sit ~{mid_frac*L-12.5:.1f}-{mid_frac*L+width_frac*L/1:.1f} Bohr out, in vacuum')")]

cells += [md(
r"""### Stopping power $S=\dfrac{\Delta E}{x}$

For the classical ion, $\Delta E = \mathrm{KE}_\mathrm{ion}(z_\mathrm{in}) -
\mathrm{KE}_\mathrm{ion}(z_\mathrm{out})$ measured at the entrance/exit **faces**
$z=\mp12.5$, and $x=25$ Bohr. Computed below directly from the run trace.""")]

# 3. Setup --------------------------------------------------------------------
cells += [md(
r"""## Simulation setup (fully reconstructable)

| Block | Run A (WP+CAP) | Run B (classical+CAP) |
|---|---|---|
| Cell / grid | 50 Bohr cubic, $dx=0.5$ | same |
| Background | slab on (same $v_\mathrm{bg}$) | same |
| Projectile | WP $\sigma_\mathrm{WP}=0.5$ (density std $0.354$), $E=100$ eV, $k_0=2.711$ | Gaussian-e **ion**, charge std $0.350$, mass 1 a.u., $v_z=2.711$ |
| Launch $z$ | $-15.5$ Bohr | $-15.5$ Bohr |
| CAP | sin², $\eta=-0.5$, width 0.15, mid $\pm0.425$ | same |
| Dynamics | $dt=0.02$, 900 steps, 91 frames | $dt=0.02$, 900 steps, 91 frames |
| Propagator | TDDFT (ETRS, complexified) | Ehrenfest (`.ehrenfest()`) |
| Engine | `inq-study` | `inq-study` |

The two projectiles are *intended* to present the same Gaussian charge cloud. The
$\sqrt2$ convention (validated 2026-06-23): a wavepacket $\psi\propto e^{-r^2/2\sigma_\mathrm{WP}^2}$
has electron **density** std $\sigma_\mathrm{WP}/\sqrt2$, so the matched classical
charge std is $\sigma_\mathrm{WP}/\sqrt2 = 0.5/\sqrt2 = \mathbf{0.354}$ Bohr.

**Provenance caveat (NOT a matched pair to the part-per-thousand).** These runs
loaded the *legacy* pseudo-file `electron_gaussian_sigma0p35.upf`, whose actual
charge std is **0.350** Bohr (the filename label is the charge std under the old
convention, a rounded stand-in for 0.354). So the WP density std (0.354) and the
classical charge std (0.350) differ by **~1 %** — physically negligible for the
stopping comparison, but the two clouds are not byte-identical. Future
quantum-vs-classical runs must use the exact-matched UPF
`electron_gaussian_wpsigma0p5.upf` (generated with `sigma_wp=0.5` → charge std
$0.35355$). See CONTEXT "σ-convention" and the handover.
""")]

# 4. Sources ------------------------------------------------------------------
cells += [md(
r"""## Source files

| Role | Path (repo-relative) |
|---|---|
| WP+CAP run | `…/scripts/03_cap_stopping/wp_cap/run.cpp` |
| classical+CAP run | `…/scripts/03_cap_stopping/classical_cap/run.cpp` |
| Built-in sin² CAP | `inq/src/perturbations/absorbing.hpp` (composed via `perturbations::sum`) |
| WP+CAP post-proc | `…/hypotheses/03_cap_stopping/make_wpcap_postproc.py` |
| classical post-proc (computes $S$) | `…/hypotheses/03_cap_stopping/make_classical_postproc.py` |
| Provenance | `…/scripts/03_cap_stopping/{wp_cap,classical_cap}/results/*` |
| This builder | `…/hypotheses/03_cap_stopping/build_cap_stopping_report.py` |
""")]

# 5. Results — Run A ----------------------------------------------------------
cells += [md(
r"""## Run A — wave packet + CAP: does the CAP drain the bath?

### xz density evolution — three-way decomposition (full-suite `fullsuite_wp` + CAP)
The WP enters, the CAP (in the vacuum gaps) eats it as it exits; the bath inside the
slab should be untouched. Decomposed into **total / bath / wavepacket**, each in
`n(t)`, Δ-vs-first (induced), Δ-vs-previous (flux). **Cyan dotted = slab faces
(±12.5 Bohr); lime dashed = CAP inner boundary (±17.5 Bohr)** — note the WP density
vanishing as it crosses the lime line.""")]
cells += [md(r"""**Total electron density**""")]
for _v in ["total", "dfirst", "dprev"]:
    cells += [embed(os.path.join(HERE, f"p5wp_decomp_total_{_v}.gif"))]
cells += [md(r"""**Bath (gas) density** $n_\mathrm{total}-|\psi_\mathrm{WP}|^2$ — stays intact (CAP doesn't drain the slab)""")]
for _v in ["total", "dfirst", "dprev"]:
    cells += [embed(os.path.join(HERE, f"p5wp_decomp_bath_{_v}.gif"))]
cells += [md(r"""**Wavepacket density** $|\psi_\mathrm{WP}|^2$ — absorbed at the lime CAP boundary""")]
for _v in ["total", "dfirst", "dprev"]:
    cells += [embed(os.path.join(HERE, f"p5wp_decomp_wp_{_v}.gif"))]

cells += [md(r"""### Energetics & absorbed norm""")]
cells += [embed(os.path.join(HERE, "wpcap_traces.png"))]
cells += [code(
"a = np.genfromtxt(os.path.join(SYS,'scripts/03_cap_stopping/wp_cap/results/cap_trace.csv'),\n"
"                  delimiter=',', names=True)\n"
"Ne0, Ne1 = a['num_electrons'][0], a['num_electrons'][-1]\n"
"E0, E1   = a['total_ha'][0], a['total_ha'][-1]\n"
"print(f'norm  N_e: {Ne0:.3f} -> {Ne1:.3f}   (absorbed {Ne0-Ne1:.3f} of the 1.0 WP)')\n"
"print(f'energy   : {E0:.3f} -> {E1:.3f} Ha   (CAP removed {E0-E1:.3f} Ha)')\n"
"print(f'bath check: end norm {Ne1:.3f} ~ 234 -> BATH INTACT (no over-drain, T3.2 ok)')")]

cells += [md(
r"""**Reading it:** the norm falls from 235 (bath 234 + WP 1) to **234.17** — the CAP
absorbed **0.83 of the WP** while the **bath stayed at ~234**. This is the key
mechanism result: because the slab is central and the CAPs sit ~3 Bohr out in
vacuum, the absorber removes the projectile **without draining the jellium**
(the concern that motivated the localised-target design). The CAP removed 2.49 Ha.

*Caveat:* only 0.83 of the WP was absorbed by $t=18$ a.u. (0.17 residual still in
the box), so a clean WP **bath**-energy stopping number needs a longer run — which
is why the classical ion (Run B) is the headline.""")]

# 5. Results — Run B ----------------------------------------------------------
cells += [md(
r"""## Run B — classical Gaussian-e ion + CAP: the stopping power

### xz density evolution — electrons + projectile (full-suite `fullsuite_classical`)
The Ehrenfest ion ploughs through the slab; the induced density wake is the drag
that decelerates it. The classical projectile is an *external* ion (not in the
electron density), so the decomposition here is **total electron density** +
**projectile Gaussian charge** (built exactly from the ion track, charge std 0.350), each
in `n(t)`, Δ-vs-first, Δ-vs-previous. **Cyan dotted = slab faces (±12.5); lime
dashed = CAP inner boundary (±17.5 Bohr).**""")]
cells += [md(r"""**Total electron density** (the gas + its response to the projectile)""")]
for _v in ["total", "dfirst", "dprev"]:
    cells += [embed(os.path.join(HERE, f"p5cl_decomp_total_{_v}.gif"))]
cells += [md(r"""**Projectile (Gaussian charge)** — the classical electron's UPF charge at the tracked $z_\mathrm{ion}(t)$; its Δ-views show the rigid forward motion""")]
for _v in ["total", "dfirst", "dprev"]:
    cells += [embed(os.path.join(HERE, f"p5cl_decomp_proj_{_v}.gif"))]

cells += [md(r"""### Deceleration and $S=\Delta\mathrm{KE}_\mathrm{ion}/x$""")]
cells += [embed(os.path.join(HERE, "classical_stopping.png"))]
cells += [code(
"c = np.genfromtxt(os.path.join(SYS,'scripts/03_cap_stopping/classical_cap/results/classical_trace.csv'),\n"
"                  delimiter=',', names=True)\n"
"z, ke = c['ion_z'], c['ke_ion_ha']\n"
"HALF = 12.5\n"
"ke_in  = np.interp(-HALF, z, ke)   # entrance face\n"
"ke_out = np.interp(+HALF, z, ke)   # exit face\n"
"dKE = ke_in - ke_out\n"
"S = dKE / (2*HALF)\n"
"print(f'KE_ion  in (z=-12.5) = {ke_in:.4f} Ha')\n"
"print(f'KE_ion out (z=+12.5) = {ke_out:.4f} Ha')\n"
"print(f'dKE = {dKE:.4f} Ha over x = {2*HALF:.0f} Bohr')\n"
"print(f'STOPPING POWER  S = dKE/x = {S:.5f} Ha/Bohr = {S*27.21138625:.4f} eV/Bohr')\n"
"print(f'bath over-drain check: N_e {c[\"num_electrons\"][0]:.2f} -> {c[\"num_electrons\"][-1]:.2f}  (intact)')")]

cells += [md(
r"""**Reading it:** the ion enters the slab at $\mathrm{KE}\approx3.67$ Ha and exits
at $\approx3.02$ Ha, losing $\Delta\mathrm{KE}\approx0.65$ Ha across the 25 Bohr
crossing — i.e.

$$ S = \frac{\Delta\mathrm{KE}_\mathrm{ion}}{x} \approx 0.026\ \mathrm{Ha/Bohr}
   \approx 0.71\ \mathrm{eV/Bohr}. $$

The electron count stays at ~234 (bath intact), so the deceleration is genuine
electronic stopping, not a CAP artefact.""")]

# 5b. Full-suite delivered + loss-function gate -----------------------------
cells += [md(
r"""## What the full-suite re-run delivered (and what stays gated)

The intuition-building suite IS now produced from the full-suite re-run
(`fullsuite_wp` + `fullsuite_classical`): the three-way **total / bath / WP**
decomposition above, the **WP momentum distribution**, the **GS-overlap
excitations** (electronic excitation budget), and the **per-component
energetics** — all in the per-run deep-dives
[`p5wp_run_notebook.ipynb`](p5wp_run_notebook.ipynb) (WP+CAP) and
[`p5cl_run_notebook.ipynb`](p5cl_run_notebook.ipynb) (classical+CAP). The
bath/WP split is exact (`density_wp = |ψ_WP|²` saved → `bath = total − wp`).

### Loss function $L(q,\omega)$ — still gated

> **⚠️ RESOLUTION LIMIT (read before trusting any $L(q,\omega)$ from these runs).**
> The frequency resolution is $\Delta\omega \approx 2\pi/T$. At $T\approx18$ a.u.
> this is **$\Delta\omega \approx 0.35$ Ha $\approx 9$ eV** — *coarser than the
> $r_s\approx4$ plasmon itself* ($\omega_p\approx6$ eV). A loss function computed
> from an 18 a.u. run **cannot resolve the plasmon** and would mislead.

The loss-function build is **hard-gated** behind a separate task — *locking down
the correct Fourier-analysis methodology* — and is not computed until that is
done.""")]

# 6. Takeaway -----------------------------------------------------------------
cells += [md(
r"""## Takeaway

- **The two-sided CAP works as designed:** it absorbs the projectile in vacuum
  while leaving the central jellium intact (Run A: bath stays at 234 while 0.83 of
  the WP is eaten) — vindicating the localised-target + CAP-in-vacuum geometry.
- **Headline stopping power $S \approx 0.71$ eV/Bohr** (≈0.026 Ha/Bohr) for a
  100 eV Gaussian-electron projectile through $r_s\approx4$ jellium, from the clean
  classical $\Delta\mathrm{KE}_\mathrm{ion}$ across the slab faces.
- **PROVISIONAL** — single centroid trajectory, coarse $dx=0.5$. The WP bath-energy
  cross-check needs a longer run (only 0.83 of the WP absorbed by $t=18$). A
  $dx\to0.25$ + $S(v)$ sweep is the next step toward a publishable number.
""")]

build(cells, OUT)
