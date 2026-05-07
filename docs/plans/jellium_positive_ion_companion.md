# Plan: positive-ion (proton) companion run for the WP-jellium scattering project

**Status:** scaffolded; awaiting user OK to launch the GS + propagation.
**Linked entries:**
- `docs/journals/researchproject/plasmons-and-stopping-power.md` §6
  (charge-conjugate prediction table)
- `docs/journals/researchproject/2026-05-05_run_base_n162_L50_E1p5.md`
  §4 (anti-wake interpretation)
- `docs/sources/correa-2018-electronic-stopping-power.md` Eq. (10)
  (stopping power from energy slope)

## Goal

Run a single H⁺ (proton) projectile through the same L=50 cubic
jellium bath at the matched velocity $v = v_F \approx 0.337$ a.u.
($\hbar k_F \approx 0.337\,\hbar/\mathrm{Bohr}$, $E_\text{kin}^\text{proton} =
\tfrac{1}{2} m_p v^2 \approx 1037$ eV — note $m_p = 1836\, m_e$, so the
"proton at $v_F$" carries far more kinetic energy than the WP at the
same velocity, but the *electronic stopping power* depends on $v$ not
on $E_\text{kin}$).

The companion serves three purposes:

1. **Confirm the charge-conjugate symmetry of the wake.** Predicted:
   density *accumulation* behind the proton (textbook electron wake),
   versus the *depletion* (anti-wake) we see behind the negative WP.
2. **Read out the stopping power $S(v)$ directly** via Eq. (10) of
   Correa 2018: $S = \langle dE/dt\rangle / v_\text{proj}$. Compare to
   Eq. (3) Lindhard prediction for a homogeneous electron gas at our
   $r_s$.
3. **Tie our jellium WP project to the textbook electronic-stopping
   literature.** This is the bridge from "novel WP-jellium scattering"
   to "standard ion-jellium electronic stopping" — the same physics
   transparently expressed.

## Method

### Cell + bath

- L = 50 Bohr cubic, periodic boundaries (same as `run_base_n162_L50_E1p5`).
- N = 162 electrons (closed shell). With the proton added the system
  is effectively `N_e = 163` once the proton's `Z = 1` is counted, but
  here we keep the **bath** at N=162 and treat the proton as an
  external moving Coulomb attractor (its 1 valence electron is
  *not* added to the electron pool — the proton enters as a bare H⁺
  pseudopotential nucleus). This is the standard
  electronic-stopping setup (Correa 2018 §5: "isolating the electronic
  stopping from a simulation with moving host atoms is possible but
  it is not discussed here. The atoms of the host are held fixed
  during the simulation").

  Implementation tweak: INQ wants a non-zero electron count to match
  the pseudopotential's expected valence, so we may need to *replace*
  the 1 most-loosely-bound bath electron with the proton's electron.
  This is a one-line accounting change; document it explicitly.

### GS preparation

Two viable initial conditions, with a clean physical interpretation
each:

**Option G1 — Reuse pure-jellium GS, inject moving proton at $t=0$.**
- Use the existing checkpoint `save_gs/gs_L50_cubic_N162_dx1p0/`.
- At $t=0^+$, add a proton at $\mathbf r_0 = (0, 0, -L/4) = (0, 0, -12.5)$
  Bohr (so it has room to propagate forward to $+L/4$ before periodic
  reentry), with velocity $\mathbf v_0 = (0, 0, +v_F)$.
- The wavefunctions are not consistent with the new external potential
  (proton's Coulomb attractor wasn't there at GS time), so the first
  few steps are a violent **injection transient** — analogous to the
  WP injection transient.
- Same transient-exclusion preprocessing rule (`t_start_au` cutoff)
  applies to extracted spectra and stopping-power slopes.

**Option G2 — Fresh GS with proton at rest, then accelerate at $t=0$.**
- Run a new `save_gs_L50_proton_static/` calculation: jellium + proton
  at rest at the launch position. The GS will have a localised
  electron cloud screening the proton.
- At $t=0^+$, suddenly impose $\mathbf v_0 = (0, 0, +v_F)$ on the
  proton. The localised cloud now lags as the proton moves out — the
  textbook *adiabatic-cloud-spreading-into-wake* picture.
- This is the standard QBall protocol (Correa 2018 §5, Listing 1:
  `set_atom_v +H 2.0 0.0 0.0` after a GS converged with the proton at
  rest).

**Recommendation: Option G2** for the canonical comparison — it
matches Correa 2018 / QBall convention and produces a cleaner
steady-state plateau in $\Delta E(t)$. Option G1 is faster (no fresh
GS needed) but has a noisier transient region.

### Propagation parameters (matched to `run_base_n162_L50_E1p5`)

| Parameter | Value | Note |
|---|---|---|
| dt_au | 0.005 | Smaller than the WP run's 0.020; needed because the proton's nuclear motion couples to the GHz-scale electronic dynamics and ETRS stability requires a finer step. *Verify this is the right value with a 100-step convergence test before the full run.* |
| N_STEPS | 6000 | Total time 30 a.u. ≈ 0.726 fs (matches WP run). |
| WRITE_EVERY | 8 | Same total-density GIF cadence (~750 frames). |
| ion_dynamics | impulsive | Fixes velocity at $v_F$; matches QBall pattern. |
| Propagator | ETRS | Default; matches WP run. |

### Observables

Identical writer set to `run_base_n162_L50_E1p5`, plus:

- **Proton position** in `observables.csv` extra columns
  `proton_x, proton_y, proton_z` and velocity `proton_vx, …` (from
  `ions.positions()[0]` and `ions.velocities()[0]` at every step).
- **$\langle F_z^\text{ion}\rangle(t)$** — the z-force on the proton,
  extracted from `ham.forces()` (or `observables::forces_stress`) at
  every step. The cumulative integral $-\int_0^t F_z\,dt = \Delta P_e$
  is the electronic momentum imparted by the proton (direct stopping
  diagnostic, à la QBall `analyse.py`).
- **`density_rt_total/`** VTI series — for the density-wake snapshot
  comparison in §3 below.
- **Standard energy components** — exactly the same four columns as
  the WP run, so the bookkeeping is directly comparable.

## Expected results (predictions before launch — falsifiable)

| Observable | Negative WP (this project, **measured**) | Positive proton (this run, **predicted**) |
|---|---|---|
| Sign of trailing-density Δn | **negative** (depletion / anti-wake) | **positive** (accumulation / wake) |
| Sign of Δ`energy_hartree` | **−0.598 eV** (down) | **negative** (down) — same sign |
| Sign of Δ`energy_kinetic` (system) | **+0.502 eV** (up) | **positive** (up) — same sign |
| Sign of Δ`energy_xc` | **+0.0955 eV** (up) | **positive** (up) — same sign |
| Cumulative $-\int F_z\,dt$ | (cod_z slope analogue: WP slows from 0.373 → 0.103 a.u.) | **negative** (force on proton points $-\hat z$, slowing it down) |
| Effective stopping power $S = \langle dE/dt\rangle/v$ | (computed via WP KE drop / Δt) | direct from the energy-slope diagnostic |

The Lindhard reference value at our $r_s = 5.72$ Bohr (N=162),
$v = v_F \approx 0.337$ a.u., from Eq. (3) of Correa 2018 (or its
graphical form in Fig. 4): $S^\text{Lindhard}(v_F) \sim 0.2 \cdot E_h /
a_0 \approx$ a few eV/Bohr. We will compare numerically once the
slope is extracted; if the proton run reproduces a value in this
ballpark, the WP run's similar slope (after appropriate sign-flip and
charge normalisation) will validate the WP-as-charge-conjugate-projectile
interpretation.

## Verdict criteria

A **report** under `docs/reports/positive-ion-vs-wp-verdict.md` will
be written using the report-writing skill, with the following
sections:

1. **Methods** (matched-velocity comparison, transient excluded).
2. **Wake-sign verification** — side-by-side `system_yz.gif` /
   `delta_z_profile.gif` snapshots for the two runs, at the matched
   propagation distance.
3. **Energy-component table** — 4 components × {WP, proton} =
   8 numbers; expected: same signs across the row.
4. **Stopping-power numbers** — $S$ extracted from each run, compared
   to the Lindhard prediction at $r_s = 5.72$.
5. **Verdict on the charge-conjugate hypothesis**:
   - **YES (charge-conjugate confirmed)** if Δn signs flip between
     the two while all energy-component signs match. Implies the WP
     slowdown is the textbook electronic-stopping wake mechanism with
     a sign-flipped density signature, and the "hole-as-attractor"
     interpretation is the right physics.
   - **NO (different mechanism)** if the energy-component signs
     differ, or if the proton run does not show a clean wake. In
     that case the WP physics is more exotic and we need a separate
     interpretation.

## Implementation scaffold

Created in this plan run (not launched):

- `ResearchProject/systems/jellium/shared/configs/positive_ion_L50_v0p33.hpp`
  — config struct `Positive_Ion_L50_v0p33` with the parameters above.
- `ResearchProject/systems/jellium/run_positive_ion_L50_v0p33/run.cpp`
  — draft run.cpp following the QBall-INQ Li pattern
  (`QuantumKickExtension/inq-codebase/Li/run_propagate_v0p0123_extensive/run.cpp`)
  adapted to a single H ion + L=50 jellium cell. **Build pending.**

## Next actions (when launched)

1. Build with `inq-run` (auto-detects). Expect ~1 minute compile.
2. **Smoke test**: 50 steps at dt=0.005, check `observables.csv`
   exists, energy drift < 1 mHa, proton position drifts forward by
   $\sim 0.084$ Bohr per step (= dt × v_F, sanity).
3. **GS**: if Option G2 chosen, run a fresh GS with proton at rest
   (`save_gs/gs_L50_proton_static_N162/`); ~5 min wall.
4. **Full propagation**: ~25 min wall (6000 steps at dt=0.005, A30
   GPU).
5. **Postprocess**: standard pipeline + new positive-ion-specific
   diagnostics (proton position, $-\int F_z$).
6. **Write the verdict report** as scoped in §"Verdict criteria".

## Risks / open questions

- **N_e accounting**: does INQ accept a non-self-consistent
  electron-count when an H pseudopotential is added to a pure-jellium
  cell? Needs a quick test with `electrons.states().num_electrons()`
  before the full run.
- **dt convergence at proton velocity**: the proton at $v_F = 0.337$
  a.u. moves $\Delta x = 0.337 \cdot 0.020 = 0.0067$ Bohr per step at
  the WP run's dt; it moves $\Delta x = 0.337 \cdot 0.005 = 0.0017$
  Bohr per step at the proposed dt = 0.005 — well under the grid
  spacing dx = 1.0 Bohr. dt = 0.020 might in fact be fine; verify
  with a smoke convergence test.
- **Proton-electron Coulomb singularity at the grid origin**: the
  pseudopotential should regularise this on a few-Bohr scale, but
  with grid spacing dx = 1.0 Bohr we're sampling the pseudopotential
  coarsely. May need to drop to dx = 0.5 Bohr near the proton — but
  that requires re-doing the GS at a finer grid. *Probable issue*.
