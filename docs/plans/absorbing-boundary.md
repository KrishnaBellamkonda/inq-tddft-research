# Plan — Absorbing-boundary benchmarking & parameter search (overnight, unattended)

**Source prompt:** `docs/prompts/absorbing_boundary/benchmarks_and_parameter_search.md`
**Paper:** De Giovannini, Larsen & Rubio, *Modeling electron dynamics coupled to
continuum states in finite volumes*, arXiv:1409.1689 (2014).
PDF: `ResearchProject/literature/tddft-quantum-projectile/resources/modeling-electron-dynamics-coupled to-continuum-states-in-finite-volumes.pdf`
**Glossary:** see the "Absorbing boundaries" section of `/local/data/public/skcb2/tddft/CONTEXT.md`.
**Crystallised by grill session 2026-06-13.** All decisions below are LOCKED.

---

## 0. Goal & two tasks

1. **Task 1 — CAP feasibility (analysis only).** Evaluate whether a complex
   absorbing potential (sin² CAP) can be implemented in **inq-study** (never
   `inq/`). Deliver a source-grounded go/no-go + a NumPy 1D toy demo. No
   inq-study modification, no INQ runs.
2. **Task 2 — Mask function absorber (implement + benchmark + parameter study).**
   Implement the sin² **mask absorber** (paper Eq. 13) entirely in the
   **inq-stack** wrapper, validate it, then sweep ε(E,L) and produce a notebook.

Deliverables: two executed notebooks in `docs/reports/absorbing-boundary/`.

---

## 1. LOCKED decisions (grill outcomes)

| # | Decision | Rationale (verified) |
|---|---|---|
| D1 | Gated pipeline: Task 1 (analysis) autonomous; Task 2 sweep **blocked** unless gate-1 + gate-2 pass | A masking bug would waste 76 runs |
| D2 | **Mask mechanism = in-callback mutation** (no restart, no engine edit) | Verified GPU run: fidelity PASS (`|ΔN|=0`, `|Δz|=0`), feedthrough PASS (norm 0.9999→0.0190). See §3.1 |
| D3 | ε = **inner-region surviving norm** `ε = ∫_{z<z_abs0}|ψ_WP(τ)|²` (paper Eq. 7) | At τ the ideal free packet has exited the box → inner region holds only the reflected wave → no per-point reference run |
| D4 | No-absorber **anchor runs** supply the ε≈1 asymptote | Avoids the expensive E→0 corner |
| D5 | **Quasi-1D in 3D**, minimal transverse box | Free-particle factorization: ψ=φ_x·φ_y·φ_z, mask ⟂-independent, ⟂ norm=1 → 3D inner integral = 1D ε exactly |
| D6 | Grid: **L∈{5,10,20,30,40,50}** × **12 k₀ (E≈0.5–490 eV)** + ~4 anchors ≈ **76 runs** | k₀≤6 keeps dx=0.1 below Nyquist (dx_max=0.38) — single resolution, no dual-grid |
| D7 | Runs in **`ResearchProject/systems/vacuum/`**; run-type = `free_wp` | per user |
| D8 | Observable set = lean+momentum (see §4); reduced cadence; density+wavefunction only on ~6 showcase runs | overnight-cheap, still analysable |
| D9 | Task 1 depth = **Level 2** (analysis + NumPy 1D toy CAP demo) | per user |
| D10 | dt=0.01, dx=0.1 (match paper); fully **periodic** cell | absorber M=0 at right cell edge ⇒ no wraparound; stop at τ before left wall matters |

---

## 2. Physics & coordinate mapping (bake into config generator)

Per energy point, with `k₀` (Bohr⁻¹):
- `σ = 4√2/k₀`  (inqkit `WavePacket.sigma` argument = this σ; density std = σ/√2)
- `E = k₀²/2` Ha `= 13.6 k₀²` eV  (paper's `5k₀²/4` prefactor is an x-axis relabel only)
- cell length on propagation axis (z): `Lcell = 6σ + L`
- absorber region: `z ∈ [z_abs0, z_abs0+L]`, `z_abs0 = (6σ − L)/2`
- WP launch: `z0 = −L/2` (= `z_abs0 − 3σ`); left wall at `z0 − 3σ = −Lcell/2`
- mask: `M(z)=1` for `z ≤ z_abs0`; `1 − sin²((z−z_abs0)π/(2L))` for `z_abs0<z<z_abs0+L`; `0` beyond
- propagation time: `τ = 2(3σ + L)/k₀`; `N_STEPS = round(τ/dt)`, `dt=0.01`
- **ε measured at the final step** = `∫_{z<z_abs0} |ψ_WP|² dV`
- transverse box `Lperp` minimal (target `N_perp ≈ 8–16`); pilot verifies ε is
  insensitive to `Lperp` (empirical check of D5)
- cutoff: `ec = ½(π/dx)²` Ha (as in `test_free_wp_engine.cpp`)
- Nyquist (PASS for all points): `dx_max = π/(k₀+3σ_k)`, `σ_k=1/(σ√2)`; worst k₀=6 → 0.38 ≫ 0.1

Free-WP recipe (from `inq-stack/tests/include/inqkit/wavepacket/test_free_wp_engine.cpp`):
empty `ions`, `options::theory{}.non_interacting()`, ghost occupied via
`extra_electrons(2.0)` (kept; non-interacting so it never touches the WP), WP
injected into the single extra state with occupation 1.0.

---

## 3. Task 2 — implementation

### 3.1 Mask mechanism (LOCKED, already verified)
The per-step `real_time::propagate` callback captures the **outer non-const
`electrons`** (the same object propagate holds by reference); the callback fires
**after** each ETRS step. Multiplying `electrons.kpin()` by `M(z)` in the callback
is exactly Eq. 12 (`ψ(t+dt)=M·U·ψ(t)`) and feeds into the next step.
INQ's `viewables` observer is const and is NOT used for mutation.
`inq/` and `inq-study` stay byte-identical for Task 2.

Verification artefact (relocate to `systems/vacuum/tests/mask_mechanism_check/`):
`ResearchProject/systems/absorbing-boundary/mask_mechanism_check/run.cpp` —
exit 0, FIDELITY PASS, FEEDTHROUGH PASS (norm 0.999999→0.018999, drop 0.981).

### 3.2 New inqkit code
- `inq-stack/include/inqkit/absorbers/mask_absorber.hpp`
  - builder: axis (z), `z_abs0`, `L`, functional form (sin²); `.apply(electrons&)`
    multiplies the WP orbital in place on GPU (mirror the injection loop in
    `wavepacket.hpp:246–257`, using `basis.point_op().rvector_cartesian`).
- ε reducer: `inner_region_norm(electrons, axis, z_abs0)` → GPU reduce of
  `|ψ_WP|²` over `z<z_abs0` (mirror the norm reduce `wavepacket.hpp:213–221`).
  **Formula-bearing** → run the **formula-validation agent** on `ε=∫|ψ|²·𝟙[z<z_abs0]`
  before locking.

### 3.3 Tests (gate-1 — must pass before any sweep)
1. **mechanism** (already passing): fidelity (M≡1 ⇒ bit-identical) + feedthrough.
2. **mask shape** (pure unit test): `M(z)` returns 1 / `1−sin²` / 0 in the three
   regions; `M(z_abs0)=1`, `M(z_abs0+L)=0`, `M` monotone decreasing — vs analytic
   values, tolerance 1e-12.
3. **ε reducer** (engine unit test): on a hand-placed orbital fully left of
   `z_abs0`, ε=‖ψ‖²; fully right, ε=0; half/half, ε=½ — known-case, pre-accepted.
   Audited by the **test-validation agent**.
Record rows in `docs/validation/test-catalogue.md`. Catalogue gate = `FIDELITY PASS && FEEDTHROUGH PASS && shape PASS && eps PASS`.

### 3.4 Config & run generation (`systems/vacuum/`)
- `systems/vacuum/shared/configs/mfa_mask_sin2.hpp` — base config (geometry
  formulas in §2, the mask, the observable selection).
- generator script emits one `run_mfa_E{e}_L{l}/` per (E,L) with `run.cpp`
  (config include + RUN_NAME) + the manifest (`RunType::free_wp`).
- anchors: `run_anchor_E{e}/` = same but mask disabled (M≡1) ⇒ expect ε≈1.

### 3.5 Gate-2 — pilot
One mid-point (e.g. k₀≈1.5, L=20): check WP norm∈[0.95,1.05], |ΔE|<1 mHa,
ε∈(0,1) and sane, **ε insensitive to `Lperp`** (run twice at 2 transverse sizes,
ε agree to ~1%). Pass ⇒ dispatch batch.

---

## 4. Observable set (`free_wp`, reduced cadence)

**All 76 runs:**
- `epsilon.txt` — ε at final step (the result)
- `inner_norm_vs_time.csv` — inner-region norm (~N/200 samples)
- `observables.csv` — E_total, E_kinetic, current, dipole (~N/100)
- `wp_real_space_stats.csv` — N, centroid_z, σ(t) (~N/100); validate vs analytic free law
- `wp_momentum_stats.csv` — ⟨p_z⟩, σ_p, E_kin (~N/100, plus exact t=0 & t=τ)
- `momentum_distribution` at **t=0 and t=τ** — `|ψ̃_WP(k)|²` (the vacuum excitation metric)
- final WP occupation / whole-box norm
- `observables_manifest.json`

**~6 showcase (E,L) runs add:** density VTI (total/system/wp; note total≈ghost+wp,
system≈flat ghost, wp = physics) ~60 frames; final WP **wavefunction** (complex
field, `ComplexField3DWriter`) for plane-wave/GS projection in post.

Showcase picks (provisional): one clearly-reflecting low-E (k₀≈0.3) and one
clearly-absorbed high-E (k₀≈4) at L∈{10,50}, plus 2 mid-points.

---

## 5. Dispatch (2 GPUs, NVML broken)

`nvidia-smi`/NVML is broken on this box (driver mismatch) — **do not poll it**.
- Before launch: `cudaMemGetInfo` probe per device (0,1). If a device has little
  free memory (occupied by another user — lm2153/fb638 seen at grill time), drop
  to single-GPU and **WARN** (email + handover).
- 2-worker queue: worker i → `CUDA_VISIBLE_DEVICES=i`, pulls next run; subprocess
  env per skill §6b (PATH/INQ_SHARE_PATH/PSEUDOPOD_SHARE_PATH).
- After each run: `analyse.py` (venv python only) → update run catalogue
  (`tddft-run-catalogue` `scan_runs.py --run`).
- Cost: ~1–2 h total; low-k₀ corner (k₀=0.2,L=50) ~few min, dominates.

---

## 6. Post-processing & notebooks (`docs/reports/absorbing-boundary/`)

- `mfa_reflectivity_study.ipynb` (Task 2): the ε(E,L) curves (one line per L,
  ε vs E log-x — paper Fig. 3 analogue), anchor ε≈1 markers, momentum
  before/after for showcase runs, density frames, method/validation narrative,
  canonical theme (`inqview.visualisation.style`). Figures `.png`.
- `feasibility_cap.ipynb` (Task 1, Level 2): §7 below.

## 7. Task 1 — CAP feasibility (analysis + NumPy toy)

Walk each step of the prompt's CAP workflow; verdict + challenge + mitigation,
citing inq-study source. Key grounded findings (verified at grill):
- `ks_hamiltonian` is **templated on `PotentialType`**; `scalar_potential_` is a
  **real-space** `field_set` (`inq/src/hamiltonian/ks_hamiltonian.hpp:47`) → a
  **spatially-restricted CAP is trivial** (real-space field nonzero on `[0,L]`);
  kinetic via FFT laplacian is untouched.
- Complex V ⇒ **non-Hermitian H** ⇒ breaks **ETRS** (assumes unitary). Use
  **Crank–Nicolson** (already in INQ; a linear solve tolerant of non-Hermitian H).
  Paper's `U_CAP=exp(−i(H₀+V_CAP)Δt)` is non-unitary.
- Main work = expose `PotentialType=complex` path + CN propagator in inq-study.
- ε via same inner-region survival logic.
Plus a **pure-NumPy 1D split-operator** demo: sin² CAP absorbing a Gaussian,
showing the expected ε(E,L) trend (no INQ). Ends with go/no-go + impl sketch.

---

## 8. Notifications (chiddukanna@gmail.com)
- email on **gate-1/gate-2 failure** (sweep halted) + diagnostic.
- **summary email** at completion: ε(E,L) PNG, run count, wall-time, failures.
- no per-run emails.

---

## 9. Risks / watch-items
- ε formula must clear the formula-validation agent before the sweep (§3.2).
- transverse-insensitivity (D5) is empirical — gate-2 confirms it; if it FAILS,
  stop and reconsider (would mean the factorization assumption is violated).
- another user occupying a GPU → single-GPU fallback (longer wall-time, still
  overnight-OK).
- low-k₀ runs are the cost driver (~1/k₀³); k₀_min=0.2 is the floor — do not lower.
- inq-study and inq/ MUST remain byte-identical for Task 2 (mechanism needs no edit).

## 10. Execution order (tonight)
1. relocate mechanism-check → `systems/vacuum/tests/`.
2. write `mask_absorber.hpp` + ε reducer + unit tests; formula/test agents; **gate-1**.
3. config + generator + pilot; **gate-2** (norm/drift/ε/transverse).
4. gates pass ⇒ dispatch 76 runs on 2 GPUs; catalogue each.
5. `mfa_reflectivity_study.ipynb` (Task 2).
6. `feasibility_cap.ipynb` (Task 1, Level 2).
7. summary email; update handover.
