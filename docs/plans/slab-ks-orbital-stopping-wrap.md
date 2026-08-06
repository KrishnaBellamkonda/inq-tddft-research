# Plan: KS-orbital stopping power on the jellium SLAB, CAP-free with wrap-around

**Created:** 2026-07-31
**Branch:** `quantum-stopping-power`
**Machine:** CSD3, `ampere` partition (A100), account `mphil-nikiforakis-skcb2-sl2-gpu`
**Status:** design locked with user (2026-07-31, 4 decisions); implementation in progress

Parent work:
- `docs/plans/bulk-jellium-ks-stopping.md` §4 — the four KS-orbital definitions this
  study transplants onto a slab.
- `docs/handovers/bulk-jellium-ks-stopping.md` — the bulk results at r_s = 5.70 and
  3.99 that these slab runs are compared against.
- `docs/campaigns/localised_jellium/classical-highdensity-sv-benchmark.md` — the
  classical S(v) curve whose geometry is reproduced exactly here.
- `docs/handovers/wavepacket-highdensity-sv-twin.md` — the CAP'd σ = 0.5/2/3 WP twin
  this supersedes for the stopping measurement.

---

## 1. Objective and the problem being solved

The KS-orbital-dependent stopping definitions (T = ⟨p²⟩/2m or ⟨p⟩²/2m, s = circular
centroid or ∫⟨p_z⟩dt) work on **bulk** jellium. On the **slab** they have so far only
been extractable over a very short time window — the CAP'd WP twin could fit only
~4 a.u. at σ_WP = 0.5 and ~16 a.u. at σ_WP = 2, because the packet either reached the
absorbing bands or its transverse periodic images overlapped. A stopping power fitted
over such a window is not scientifically defensible.

**This study removes the CAP and lets the wavepacket wrap.** The packet crosses the
slab, exits the +z face, re-enters at −z, and crosses again — ~14 slab crossings over
362 Bohr of path. The fit window becomes the whole run.

**Hypothesis (H):** the bulk KS-orbital definitions transfer to a slab geometry once
the measurement is made over a multi-pass, CAP-free, wrap-around trajectory; the
resulting S(v) is comparable with the classical S(v) curve for the same slab.

### Why no boundary-condition change is needed — the key engine fact

`periodicity(2)` is consulted **only** by the Poisson solver, ionic replicas and the
kick gauge (`inq/src/solvers/poisson.hpp:189,206`). The wavefunction basis and the
kinetic operator are a plain 3-D FFT, periodic in **all three** directions
(`inq/src/basis/fourier_space.hpp:60-151`,
`inq/src/hamiltonian/ks_hamiltonian.hpp:200-204`). A KS orbital travelling +z
therefore **already wraps** and re-enters at −z — which is precisely why the CAP had
to be added to the σ = 0.5/2/3 campaign in the first place
(`scripts/wp_highdensity_sv/wp/run.cpp:24-35`).

Turning the CAP off does not *introduce* the wrap; it *restores* it. So the cell,
the boundary condition, the ground state and every bath parameter stay byte-identical
to the classical benchmark. `periodicity(2)` is retained: the electrostatics stay
z-open, so the slab has **no** spurious periodic images along z.

### Two bonuses that fall out of the long runs

1. **The plasma period finally fits inside the run.** The bulk study's stated caveat
   was that 2π/ω_p = 49 a.u. exceeded its 26 a.u. run, so no wake could form. Here at
   r_s = 4.18 (T_plasmon = 31.0 a.u.) a 181 a.u. run spans **5.8 plasma periods** —
   at or above the 5-period wake criterion for the first time in this line of work.
   At r_s = 5.67 (T_plasmon = 49.0 a.u.) it spans 3.7.
2. **The deceleration sweeps a velocity range.** Each run is its own S(v) scan over
   v₀ → v_final, per `.claude/rules/light-projectile-stopping.md`.

---

## 2. Locked configuration (user decisions, 2026-07-31)

| Quantity | Value | Source |
|---|---|---|
| Cell | 35 × 35 × 85 Bohr, orthorhombic, **periodicity(2)** | classical benchmark, unchanged |
| Slab | half-width 12.5 (25 Bohr thick), centred z = 0, erfc edge 1.0 Bohr | unchanged |
| Grid spacing | **dx = 0.40 Bohr** | user |
| Time step | dt = 0.04 a.u. | unchanged |
| Theory | LDA ground state, ALDA real time; ETRS propagator | unchanged |
| Engine | `inq-study` | unchanged |
| **CAP** | **NONE** (`LJ_CAP_ETA=0`) | user |
| σ_WP | **2.0 Bohr** both halves (classical σ_pot = 2/√2 = 1.41421) | user; matches the bulk KS study exactly |
| Launch z₀ | −24.0 Bohr | classical benchmark, unchanged |
| Projectile mass / charge | 1 / −1, free Ehrenfest (classical half) | unchanged |
| Velocities | v₀ ∈ {2.0, 2.5, 3.0, 3.5} | user (4.0/4.5 dropped) |
| Run length | 1.5 × the classical CAP-free step count | user |

### Densities — the only variable between the two systems

Both use the SAME box, slab thickness, grid, dt and σ. Only N changes.

| | N | n₀ (a.u.⁻³) | r_s | k_F = v_F | E_F | ħω_p | T_plasmon | states |
|---|---|---|---|---|---|---|---|---|
| **D1** | 100 | 3.2653e-3 | **4.183** | 0.4590 | 2.87 eV | 5.51 eV | 31.0 a.u. | 74 |
| **D2** | 40 | 1.3061e-3 | **5.674** | 0.3383 | 1.56 eV | 3.49 eV | 49.0 a.u. | 44 |

Density ratio 2.50×. D1 is the exact system of the ongoing classical S(v) curve
(GS already computed). D2 lands on the project's long-standing jellium reference
density and on the bulk KS study's low-density point (r_s = 5.702), so both ends are
anchored to existing results.

### Velocity grid and run lengths

`N_wp = 1.5 × N_classical(no-CAP)`. Every velocity covers the same 362.3 Bohr of
path = 4.26 box lengths ≈ 14.5 slab crossings.

| idx | v₀ | KE₀ | N_cl | **N_steps** | t (a.u.) | plasma periods (D1 / D2) | save/ | wf/ | ckpt/ |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 2.0 | 54 eV | 3019 | **4529** | 181.2 | 5.8 / 3.7 | 15 | 45 | 906 |
| 1 | 2.5 | 85 eV | 2415 | **3623** | 144.9 | 4.7 / 3.0 | 12 | 36 | 725 |
| 2 | 3.0 | 122 eV | 2013 | **3020** | 120.8 | 3.9 / 2.5 | 10 | 30 | 604 |
| 3 | 3.5 | 167 eV | 1725 | **2588** | 103.5 | 3.3 / 2.1 | 9 | 27 | 518 |

v = 4.0 and 4.5 are dropped by user instruction. NOTE for the record: their original
justification (z-momentum aliasing) **no longer applies at σ_WP = 2** — σ_p =
1/(√2σ) = 0.354 against k_Nyq = π/0.4 = 7.85, and the σ = 2 smoke gate measured
momentum weight past Nyquist at 7e-60 %. They remain runnable as extra array indices
if the user later wants the curve extended.

---

## 3. The dispersion problem, stated honestly, and how S is extracted

At σ_WP = 2 the packet spreads as σ_d(t) = √(σ²/2 + t²/(2σ²)) = √(2 + t²/8):

| t (a.u.) | 0 | 16 | 34.6 | 60 | 120 | 181 |
|---|---|---|---|---|---|---|
| σ_d (Bohr) | 1.41 | 5.83 | 12.5 | 21.2 | 42.5 | 64.0 |

- transverse periodic images overlap (6σ_d = L_xy = 35) at **t ≈ 16 a.u.**
- packet is wider than the slab half-thickness at **t ≈ 34.6 a.u.**
- packet is wider than the box at **t ≈ 120 a.u.**

This is not fixable by choosing a better σ: the minimum width any Gaussian can have at
time T is √T, so nothing stays under ~13 Bohr for 180 a.u. The packet is
projectile-like for roughly the first pass and delocalised thereafter, **independent
of velocity** (dispersion is set by σ, not v).

That is acceptable because ⟨p⟩, ⟨p²⟩ and ∫⟨p_z⟩dt remain exact expectation values of a
unitarily-evolved orbital no matter how spread it is, and because a spread orbital
samples slab and vacuum in a *measurable* proportion. Two windows are therefore
reported side by side (user decision):

### Window A — first-pass drag (comparable to the classical single-pass number)

Fit over t ∈ [t_transient, min(t_wrap, 34.6 a.u.)], where t_wrap = (42.5 − z₀)/v is the
first +z face crossing (33.3 / 26.6 / 22.2 / 19.0 a.u. for the four velocities) and
t_transient = 20 % of the window per the `stopping-power-extraction` skill. This is
the localised-projectile regime and is directly comparable with the classical
benchmark's single-pass S.

### Window B — whole run, slab-overlap-weighted path

Define the **in-slab occupancy**

    f(t) = ∫_{|z| ≤ 12.5} |ψ_wp(r,t)|² dV        (WP norm inside the slab)

and the **in-slab path**

    s_slab(t) = ∫₀ᵗ f(t') · ⟨p_z⟩(t')/m dt'

Then, since only the in-slab fraction of the orbital feels drag,

    dT/dt = −F · v · f      and      ds_slab/dt = f · v      ⇒     −dT/ds_slab = F = S_slab

so S_slab = −dT_i/ds_slab is the stopping power **per Bohr of path travelled inside
the slab**, valid in both the localised and the delocalised regime, and reducing to
the ordinary −dT/ds when the packet is entirely inside the slab (f = 1). When the
packet is uniformly spread, f → 25/85 = 0.294 and the estimator automatically applies
the geometric factor that a naive whole-run fit would miss.

`f(t)` is measured **exactly, per step**, by a new observable (§5.1) — not modelled
from a Gaussian ansatz. A wrapped-Gaussian model built from the already-recorded
`z_mean_circ` and `sigma_z_circ` is computed in post-processing as an independent
cross-check, and the two are reported together.

### The four definitions (unchanged from the bulk study)

| | T | source |
|---|---|---|
| T₁ | ⟨p²⟩/2m | `wp_momentum_stats.csv`, `e_kin_ha` |
| T₂ | ⟨p⟩²/2m | ½(px_mean² + py_mean² + pz_mean²) |

| | s | source |
|---|---|---|
| s₃ | circular density centroid, unwrapped to a monotone path | `wp_real_space_stats.csv`, `z_mean_circ` + `np.unwrap` |
| s₄ | ∫⟨p_z⟩dt | cumulative trapezoid of `pz_mean` |
| s₅ | **in-slab path** (new, window B) | ∫ f·⟨p_z⟩/m dt |

T₁ − T₂ = 3/(4σ²) = 0.1875 Ha = **5.10 eV** at t = 0. The bulk study established that
its *growth* is a fixed ≈0.043 eV/Bohr offset independent of bath density — most
likely self-interaction error, since the WP is an occupied KS orbital whose own charge
enters the Hartree potential. **S₂ (drift) is the defensible stopping power; S₁ must
be reported as "S₂ minus a spreading term", never as an independent measurement.**
That conclusion is carried over here and re-tested, not re-derived.

s₃ vs s₄ is a **validation check, not independent physics**: with no ions and no CAP
the KS Hamiltonian is local and Ehrenfest gives d⟨z⟩/dt = ⟨p_z⟩/m exactly, so the two
must agree to numerical precision. Any deviation localises to wrap handling or norm
leakage.

---

## 4. Classical twin (user decision: wrapped, σ_WP = 2)

The published classical curve is σ_WP = 0.5 and single-pass — a 4× width mismatch and
a different number of slab crossings, so it cannot serve as the reference here. New
classical twins are run at σ_WP = 2 with the projectile **wrapped in z**, same step
counts, same densities.

Two code changes are required for the twin to be honest:

1. **Projectile position wrap.** `inqkit::dynamics::Projectile` has no position
   setter; add one (`set_R`) and wrap R_z into [−L_z/2, +L_z/2) each step behind an
   opt-in `LJ_WRAP_Z` flag, so existing runs are unaffected.
2. **Minimum-image Gaussian charge.** `inqkit::jellium::gaussian_density`
   (`projectile_background_energy.hpp:56-70`) uses a plain Cartesian distance, so a
   Gaussian at the box face is **clipped**, whereas the wavepacket wraps smoothly on
   the FFT grid. Without a minimum-image variant the two twins would differ exactly at
   the boundary we are deliberately introducing. Added as a NEW function, leaving the
   existing one untouched so the published campaign binaries stay reproducible.

**Expected behaviour, stated in advance so it is not mistaken for a bug:** a classical
mass-1 electron at v = 2.0 has KE = 54 eV and the benchmark measured 27 eV deposited
per slab crossing, so it will **stop after ~2 crossings** and then sit. At v = 3.5
(KE = 167 eV, 12.7 eV per crossing) it survives the whole run. This is real physics
under free Ehrenfest (`.claude/rules/light-projectile-stopping.md`) and is why S is
extracted as the initial drag, not a full-run regression.

---

## 5. Implementation

### 5.1 Library (inqkit) — two additive changes plus one new observable

| File | Change | Test |
|---|---|---|
| `inqkit/dynamics/projectile.hpp` | `set_R(Vec3)` | extend `inq-stack/tests/include/inqkit/dynamics/test_projectile.cpp` — set/get round-trip and a wrap sequence |
| `inqkit/jellium/projectile_background_energy.hpp` | `gaussian_density_minimum_image(basis, center, sigma, wrap_axes)` | engine test: identical to `gaussian_density` for a centred blob; ∫ = 1 for a blob straddling the face (the existing one loses norm) |
| `inqkit/observables/slab_occupancy.hpp` (new) | `slab_occupancy(field, axis, center, half_width)` → ∫ over the slab band | engine test: uniform field → 2h/L; narrow Gaussian at centre → 1 |

### 5.2 Run machinery — new sweep folder (ADR 0007 layout)

`ResearchProject/systems/localised_jellium/scripts/slab_ks_wrap/`

- `wp/run.cpp` — fork of `scripts/wp_highdensity_sv/wp/run.cpp` with: CAP default 0,
  σ default 2.0, N/n₀/GS taken from env so ONE binary serves both densities, and the
  per-step in-slab occupancy column appended to `wp_real_space_stats.csv`.
- `classical/run.cpp` — fork of `scripts/classical_highdensity_sv/dyn/run.cpp` with
  `LJ_WRAP_Z` and the minimum-image Gaussian; emits `proj_z` (wrapped) plus
  `proj_z_unwrapped` (cumulative path).
- `gs/run.cpp` — ground state for D2 (N = 40), reusing the existing D1 GS.

### 5.3 Config header

`shared/configs/slab_n40_L35x35x85.hpp` — identical to `slab_n100_L35x35x85.hpp`
except `N_ELECTRONS = 40` and `SPACING_BOHR = 0.40`.

### 5.4 Ground states

| | path | status |
|---|---|---|
| D1 | `shared_gs/slab_n100_L35x35x85_dx0p4_per2` | **exists** |
| D2 | `shared_gs/slab_n40_L35x35x85_dx0p4_per2` | to compute (~10 min) |

### 5.5 SLURM (parallel, user instruction)

- `shared/bin/run-slab-ks-gs.slurm` — D2 ground state.
- `shared/bin/run-slab-ks-wp.slurm` — `smoke` stage (builds + t=0 gates) then
  `--array=0-7` (index = density × 4 + velocity), all 8 concurrent.
- `shared/bin/run-slab-ks-classical.slurm` — same shape, 8 concurrent.
- `shared/bin/submit-slab-ks-wrap.sh` — the chain: D2 GS → wp smoke → wp array →
  classical smoke → classical array → notebooks.

Every run keeps the interior + final checkpointing already in the forked binaries
(`.claude/rules/final-timestep-checkpoint.md`, `checkpoint-dont-block.md`): 5 retained
numbered checkpoints plus the rolling one, `LJ_RESUME=1` extends.

### 5.6 Cost

At the measured 2.75 s/step (σ = 2, dx = 0.4, 74 states): D1 WP ≈ 3.5 / 2.8 / 2.3 /
2.0 h; D2 (44 states) ≈ 0.6×; classical similar. **≈ 29 GPU-hours over 16 runs**,
~3.5 h wall per wave if 8 run concurrently.

---

## 6. Validation

### Correctness gates (block the sweep)

- **WP, t = 0** (already in the forked binary): T₁ − T₂ = 3/(4σ²) = 5.10 eV; density
  std = σ/√2 = 1.414; ⟨p_z⟩ = k₀; max overlap with the occupied manifold < 1e-3.
- **WP, whole run:** with no CAP the Hamiltonian is time-independent ⇒ `energy_total`
  conserved (target: drift < 1e-3 eV) and WP norm = 1 ± 1e-6. This is a *stronger*
  gate than the CAP'd campaign could use, and it is the main reason the CAP-free
  design is scientifically cleaner.
- **s₃ ≡ s₄** to numerical precision (the Ehrenfest identity above).
- **Classical:** E_electronic + KE_proj + U_proj_bg flat; wrap events visible as
  single L_z jumps in `proj_z` and absent from `proj_z_unwrapped`.
- **D2 ground state:** ∫n dV = 40; SCF converged; E_GS recorded.

### Sanity checks (report, never gate)

- S(v) monotone decreasing over v ∈ [2, 3.5] (Bethe tail, above the Lindhard peak).
- S_slab(D1)/S_slab(D2) against the 2.50× density ratio and against the bulk study's
  measured scaling (S₂ scaled 2.73× for a 2.92× density ratio).
- Window A vs Window B agreement in the overlap region — if they disagree badly, the
  overlap-weighting is doing something wrong and must be diagnosed before the numbers
  are used.
- WP S₂ vs classical S at the same σ and velocity — the headline comparison.

---

## 7. Deliverables

- `hypotheses/slab_ks_wrap/` — `slab_ks_stopping.py` (loader + the five definitions +
  both windows), per-run notebooks with the mandatory density-matrix GIF
  (`.claude/rules/notebook-density-gif.md`), a synthesis notebook carrying S(v) for
  WP and classical at both densities, and `S_summary.csv`.
- `docs/handovers/slab-ks-orbital-stopping-wrap.md` — rolling handover.
- Gmail per phase (`email-notifications` skill), four-part structure, plot attached.

---

## 8. Open / accepted limitations

1. The packet is projectile-like only for the first pass; windows A and B split the
   measurement accordingly. This is inherent to a Gaussian over 180 a.u., not a
   choice that can be tuned away.
2. Transverse periodic images overlap from t ≈ 16 a.u. The transverse box is 35 Bohr
   because the classical benchmark's density demanded it; widening it changes r_s.
   Consequence: from t ≈ 16 a.u. the object being dragged is effectively a periodic
   array of packets, not one packet. Stated in the notebooks; not correctable within
   a like-for-like comparison.
3. The T₁ − T₂ growth (SIE, per the bulk study) is present here too and is why S₁ is
   never quoted alone.
4. Classical projectiles at v = 2.0 stop mid-run (§4); their S comes from the initial
   drag only.
