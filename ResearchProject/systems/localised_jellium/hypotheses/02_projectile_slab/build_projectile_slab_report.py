#!/usr/bin/env python3
"""Build the Phase 3 study notebook: wave-packet projectile baseline (no CAP).

    PYTHONPATH=/local/data/public/skcb2/tddft/inq-stack/python \
    /local/data/public/skcb2/tddft/venv/bin/python3 build_projectile_slab_report.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from _nbreport import md, code, setup_cell, embed, set_outdir, build, SYS

HERE = os.path.join(SYS, "hypotheses", "02_projectile_slab")
OUT = os.path.join(HERE, "projectile_slab_study.ipynb")
set_outdir(os.path.dirname(OUT))

cells = []

# 1. Title + question ---------------------------------------------------------
cells += [md(
r"""# Wave-packet projectile through the localised jellium slab — baseline (no CAP)

**Phase 3 of the localised-jellium scattering campaign.**

## The question

With the slab background validated (Phase 1+2), inject an electronic **wave packet**
from vacuum and let it traverse the slab in a **closed box (no absorber)**:

> Does the WP cross the slab and excite the electron gas while **total energy is
> conserved**? A closed, Hermitian system *must* conserve energy — this run is the
> baseline that proves the projectile machinery (wave-packet injection, screens,
> observables) is sound before we add the absorbing CAP in Phase 5.

Energy *loss* cannot be read off here (it is a closed system: the energy the bath
gains, the WP loses). The point is the **conservation check** and the qualitative
wake, which set the reference for the CAP run.

| Where this sits | Phase | Notebook |
|---|---|---|
| prev | 1+2 slab GS validation | `01_slab_validation/` |
| **this** | 3 WP projectile baseline | `02_projectile_slab/` |
| next | 5 CAP stopping power | `03_cap_stopping/` |
""")]

# 2. Conventions + formulas ---------------------------------------------------
cells += [md(
r"""## Conventions & symbols

Atomic units. $1\,\mathrm{Ha}=27.2114\,\mathrm{eV}$.

| Symbol | Meaning | Units |
|---|---|---|
| $E$ | projectile kinetic energy | eV |
| $k_0$ | central wavevector | $a_0^{-1}$ |
| $\sigma$ | WP spatial width (Gaussian) | $a_0$ |
| $d_z$ | dipole moment along $z$ | $a_0$ (a.u.) |
| $\tau$ | total propagated time | a.u. |
""")]

cells += [setup_cell()]

cells += [md(r"""### Projectile wavevector $k_0=\sqrt{2E}$ (a.u.)

A free electron of kinetic energy $E$ has $k_0=\sqrt{2E}$; its group velocity is
$v=k_0$. Convert $E=100$ eV to Hartree first.""")]
cells += [code(
"E_eV = 100.0\n"
"E_ha = E_eV / 27.21138625\n"
"k0 = np.sqrt(2.0 * E_ha)\n"
"print(f'E = {E_eV} eV = {E_ha:.4f} Ha  ->  k0 = sqrt(2E) = {k0:.4f} a0^-1  (v = {k0:.4f} a.u.)')")]

cells += [md(r"""### Time to cross the slab $\Delta t = L_\mathrm{slab}/v$

At $v=k_0$ the WP centroid crosses the 25 Bohr slab in $\approx L_\mathrm{slab}/v$;
the launch point $z=-23$ adds ~10.5 Bohr of run-up to the $z=-12.5$ face.""")]
cells += [code(
"L_slab = 25.0\n"
"t_cross = L_slab / k0\n"
"runup = 23.0 - 12.5\n"
"print(f'slab crossing time ~ {t_cross:.2f} a.u.;  run-up {runup:.1f} Bohr -> face at t~{runup/k0:.2f} a.u.')")]

# 3. Setup --------------------------------------------------------------------
cells += [md(
r"""## Simulation setup (fully reconstructable)

| Block | Value |
|---|---|
| Cell / grid | 50 Bohr cubic, $dx=0.5$, periodic |
| Background | slab on (half-width 12.5, axis $z$), same static $v_\mathrm{bg}$ as Phase 1+2 |
| Projectile | electronic WP, $\sigma=0.5$ Bohr, $E=100$ eV, $k_0=2.711$, launched $+z$ from $z=-23$ (4$\sigma$ off the $-z$ wall) |
| Injection | orthogonalised against occupied states, placed in an extra (empty) state, norm 1 |
| Dynamics | $dt=0.02$, $N_\mathrm{steps}=880$, $\tau=17.6$ a.u. |
| Diagnostics | 20 feature-aligned plane screens ($z=-24+k\cdot48/19$), 89 density frames, dipole, current, energy |
| Engine | stock `inq` (closed, Hermitian — no CAP) |

The WP is injected with `inqkit::WavePacket{}.center(0,0,-23).sigma(0.5).k0(0,0,k0)`,
orthogonalised against the occupied manifold and dropped into the last extra state
so it carries unit norm without disturbing the bath occupation.
""")]

# 4. Sources ------------------------------------------------------------------
cells += [md(
r"""## Source files

| Role | Path (repo-relative) |
|---|---|
| WP run | `…/scripts/02_projectile_slab/wp_slab/run.cpp` |
| Wave-packet injector | `inq-stack/include/inqkit/wavepacket/wavepacket.hpp` |
| Background perturbation | `inq-stack/include/inqkit/jellium/background_perturbation.hpp` |
| Post-processing | `…/hypotheses/02_projectile_slab/make_wp_postproc.py` |
| Provenance | `…/scripts/02_projectile_slab/wp_slab/results/{run_summary.txt,energy_dipole_vs_time.csv}` |
| This builder | `…/hypotheses/02_projectile_slab/build_projectile_slab_report.py` |
""")]

# 5. Results ------------------------------------------------------------------
cells += [md(
r"""## Results

### xz density evolution — the WP crosses the slab

The visual sanity check: a compact packet enters from below ($z=-23$), traverses the
slab (dotted faces at $z=\pm12.5$), and leaves the top. A negative WP drives a
**depletion anti-wake** in its trail (charge-conjugate of the textbook positive-ion
wake).

All views below are from the **full-suite re-run** (`fullsuite_wp`, no CAP),
loaded through `inqview.load_vti` (physical order, **no fftshift**). The density is
decomposed **three ways** — *total* (bath + WP), *bath* (the gas alone,
$n_\mathrm{bath}=n_\mathrm{total}-|\psi_\mathrm{WP}|^2$), and *wavepacket*
($|\psi_\mathrm{WP}|^2$) — each in three views: absolute `n(t)`, induced
`Δn=n(t)−n(0)`, and flux `Δn=n(t)−n(t−Δt)` ($\propto-\nabla\cdot j$). Cyan dotted =
slab faces (±12.5 Bohr).""")]
cells += [md(r"""**Total electron density** (bath + projectile)""")]
for _v in ["total", "dfirst", "dprev"]:
    cells += [embed(os.path.join(HERE, f"p3wp_decomp_total_{_v}.gif"))]
cells += [md(r"""**Bath (gas) density** $n_\mathrm{bath}=n_\mathrm{total}-|\psi_\mathrm{WP}|^2$ — the gas response in isolation""")]
for _v in ["total", "dfirst", "dprev"]:
    cells += [embed(os.path.join(HERE, f"p3wp_decomp_bath_{_v}.gif"))]
cells += [md(r"""**Wavepacket density** $|\psi_\mathrm{WP}|^2$ — the projectile alone (a clean compact packet, the depletion it leaves shows up in the bath row above)""")]
for _v in ["total", "dfirst", "dprev"]:
    cells += [embed(os.path.join(HERE, f"p3wp_decomp_wp_{_v}.gif"))]

cells += [md(
r"""> **Full per-run deep-dive** (momentum distribution, GS-overlap excitations,
> per-component energetics, carpets) is in
> [`p3wp_run_notebook.ipynb`](p3wp_run_notebook.ipynb). The bath/WP separation
> above is exact: `density_wp = |ψ_WP|²` is now saved, so
> `bath = total − wp` is a subtraction, not a proxy.""")]

cells += [md(r"""### Provenance""")]
cells += [code(
"print(open(os.path.join(SYS,'scripts/02_projectile_slab/wp_slab/results/run_summary.txt')).read())")]

cells += [md(
r"""### Response: dipole drift and energy conservation

`wp_response.png` shows the $z$-dipole tracking the WP centroid across the cell and
the total energy versus time.""")]
cells += [embed(os.path.join(HERE, "wp_response.png"))]
cells += [code(
"d = np.genfromtxt(os.path.join(SYS,'scripts/02_projectile_slab/wp_slab/results/energy_dipole_vs_time.csv'),\n"
"                  delimiter=',', names=True)\n"
"dE = d['total_ha'][-1] - d['total_ha'][0]\n"
"print(f'E(0)    = {d[\"total_ha\"][0]:.6f} Ha')\n"
"print(f'E(tau)  = {d[\"total_ha\"][-1]:.6f} Ha')\n"
"print(f'dE      = {dE*1e3:.2f} mHa over tau=17.6 a.u.  (closed system: expected ~0)')\n"
"print(f'dipole_z: {d[\"dipole_z\"][0]:.2f} -> {d[\"dipole_z\"][-1]:.2f} a.u.  (centroid swept -23 -> +)')" )]

cells += [md(
r"""**Reading it:** the total energy drifts by only **$-1.1$ mHa** over the whole
traversal — small and consistent with a closed Hermitian system (the residual is
grid/time-step discretisation at high $k$, flagged as the convergence caveat). The
dipole sweeps monotonically as the WP crosses, confirming clean transit. There is
**no energy sink** here by construction; extracting a stopping number requires the
absorbing CAP of Phase 5.""")]

# 6. Takeaway -----------------------------------------------------------------
cells += [md(
r"""## Takeaway

- The WP **injection + propagation + screen/observable machinery works**: a
  $\sigma=0.5$, 100 eV packet cleanly traverses the slab and drives a visible
  depletion anti-wake.
- **Total energy is conserved to $\sim1$ mHa** over 17.6 a.u. — the closed-system
  baseline. The small residual is the high-$k$ grid/time-step signal that motivates
  the $dx{\to}0.25$ convergence caveat carried into Phase 5.
- This run **cannot** give a stopping power (closed box, no sink). It is the
  reference against which the Phase-5 CAP run is read.
""")]

build(cells, OUT)
