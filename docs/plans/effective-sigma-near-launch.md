# Plan — Effective-σ hypothesis: near-launch σ=0.5 wavepacket sweep (localised jellium)

Branch `quantum-stopping-power`. Machine **CSD3**, `ampere`, account
`mphil-nikiforakis-skcb2-sl2-gpu`. Started 2026-08-01.

Parent campaign: `docs/campaigns/localised_jellium/classical-highdensity-sv-benchmark.md`
(id `classical-highdensity-sv`).
Parent handover: `docs/handovers/wavepacket-highdensity-sv-twin.md`.

---

## 1. Hypothesis

**It is not the *launch* σ_WP that sets the wavepacket–slab interaction — it is the
packet's width when it *arrives* at the slab.**

The existing σ campaigns cannot separate the two, because launch width and arrival
width are locked together by the flight path. At σ_WP = 0.5 launched from
z = −24 the packet disperses at 1/(√2σ) = 1.414 Bohr/a.u. and is **4.7–8.1 Bohr
wide on arrival** (parent handover, 2026-07-31) — it is simply not a σ = 0.5
projectile by the time it does the physics.

**The test.** Hold σ_WP = 0.5 and move the launch point from z = −24 to
**z = −14** (1.5 Bohr outside the slab face at −12.5). The packet then arrives
essentially undispersed, at density std σ/√2 = 0.354 Bohr.

- If S(v) is unchanged → the launch σ label is what matters; dispersion en route
  is irrelevant.
- If S(v) changes → the launch σ label is **not** the controlling parameter, and
  every S(σ_WP) curve in the parent campaign is really an S(arrival width) curve
  mislabelled.

Either outcome is informative; the second is the hypothesis.

---

## 2. Geometry and the locked decisions (user, 2026-08-01)

Cfg `localised_jellium::config::SlabN100_L35x35x85`
(`ResearchProject/systems/localised_jellium/shared/configs/slab_n100_L35x35x85.hpp`):
slab faces at z = ±12.5, `EDGE_WIDTH_BOHR = 1.0` (erfc-softened), Lz = 85.

| Decision | Value | Rationale |
|---|---|---|
| Launch z | **−14.0** (1.5 Bohr beyond the face) | Just outside the 1 Bohr softening |
| Scan floor | **−14.0** — retreat only | 1.5 Bohr is the floor, never go closer |
| Scan step | **0.5 Bohr** outward (−14.0, −14.5, −15.0, …) | vs 1.23 Bohr spill-out decay length |
| Accept criterion | **orthogonalisation-removed weight < 3 %** | User-set |
| Scan velocity | **v = 2.0 only**, result applied to all four | Worst-case overlap; one common geometry |
| Gaussianity | **reported, never vetoes** | 3 % alone decides |
| σ | **0.5 only**, 4 velocities (2.0, 2.5, 3.0, 3.5) | Decisive single-variable test |
| N_STEPS | **identical** to far-launch (3623/2898/2415/2070) | 10 Bohr shorter path ⇒ more post-exit plateau; comparable time budget |
| Vacuum CAP controls | 4 new, at the same launch z | Baseline is only subtractable if launch matches |
| CAPs | ON, η = −1 Ha, 12.5 Bohr/face | Already the campaign configuration |

Everything else (GS, dx = 0.40, dt = 0.04, cadences, LDA, CAP geometry) is held
at the parent campaign's values so the only changed variable is launch z.

---

## 3. Why a scan is needed at all (measured + inferred)

Bath spill-out, planar-averaged from the converged GS
(`.../wp/results/v2p0/raw/vti/density_gs_system/density_gs_system.vti`) — MEASURED:

| z | dist. from face | n(z)/n₀ |
|---|---|---|
| −12.5 | 0 | 0.39 |
| −13.5 | 1.0 | 0.18 |
| **−14.0** | **1.5** | **0.135** |
| −15.2 | 2.7 | 0.052 |
| −16.0 | 3.5 | 0.026 |
| −18.0 | 5.5 | 0.004 |
| −24.0 | 11.5 | 1e-5 |

Decay length ≈ **1.23 Bohr** ⇒ removed weight should fall ~2.25× per Bohr of
retreat, so the scan converges in a few trials.

**Inference (order-of-magnitude, NOT a measurement).** Scaling the one recorded
datum — `max_overlap = 3.7e-4` at z = −24 (`v2p0/raw/observables/wp_config.txt`)
— by the 1e4 density ratio suggests removed weight ~1 % at z = −14, i.e. the
scan likely accepts −14.0 on the first trial. Ignores the k-space matching
factor. Used only to bracket the scan, never as a result.

k_F = (3π²n₀)^⅓ = 0.457 Bohr⁻¹ at n₀ = 3.2653e-3, versus k₀ = 2.0–3.5 and
σ_p = 1/(√2σ) = 1.414 — so the packet sits well outside the occupied Fermi
sphere but has real low-k weight. v = 2.0 is the worst case: (2.0−0.457)/1.414 =
1.09 σ_p, against 2.15 σ_p at v = 3.5.

---

## 4. Engine facts established (verified this session, file:line)

- **Orthonormalisation is a t = 0 injection-only operation.** Two-pass modified
  Gram–Schmidt against all occupied states, then renormalise —
  `inq-stack/include/inqkit/wavepacket/wavepacket.hpp:299-392`. These runs use
  the **ETRS** propagator (INQ default; `wp/run.cpp:492` overrides nothing), and
  `inq-study/src/real_time/etrs.hpp` contains no `orthogonalize` call. Only
  `inq-study/src/real_time/crank_nicolson.hpp:139,162` orthogonalises per step.
  ⇒ the deformation risk is entirely at injection and needs no propagation to
  measure.
- **`InjectionReport` cannot currently express the deformation.**
  `norm_after` is measured AFTER renormalisation (`wavepacket.hpp:394-405`) so it
  is ≈1 by construction; `max_overlap` is only the single largest ⟨ψᵢ|ψ_WP⟩.
  The quantity wanted — Σᵢ|⟨ψᵢ|ψ_WP⟩|² — is computed and discarded. **This is the
  one library change this plan requires.**
- **`LJ_LAUNCH_Z` is already an env knob** (`wp/run.cpp:177`, default −24.0), so
  moving the launch point needs **no C++ change to the production run**.
- Complex WP orbitals are written as VTI with separate `wavefunction_real` /
  `wavefunction_imag` DataArrays; `inqview.load_vti` reads `GetArray(0)` only
  (`visualisation/field_io.py:86`) so it needs an array-selection parameter.

---

## 5. Work items

### 5.1 inqkit — expose the orthogonalisation loss (library change)

`inq-stack/include/inqkit/wavepacket/injection_report.hpp`:
add `norm_pre_ortho`, `norm_pre_renorm`, `removed_weight`, `sum_overlap_sq`.

`inq-stack/include/inqkit/wavepacket/wavepacket.hpp`: measure ‖ψ‖² immediately
after Gaussian construction and again after GS but *before* renormalisation;
accumulate Σᵢ|⟨ψᵢ|ψ_WP⟩|² on pass 0.

    removed_weight = 1 − (norm_pre_renorm / norm_pre_ortho)²

Purely additive — existing fields and behaviour unchanged, so the parent
campaign's binaries stay valid.

**Internal consistency check (free, exact in theory):** in exact arithmetic
`sum_overlap_sq == norm_pre_ortho² − norm_pre_renorm²`, because the KS states are
mutually orthonormal. Both routes are computed and compared — a built-in known-case
test on every real injection.

### 5.2 Test (validation-gates rule)

`inq-stack/tests/include/inqkit/wavepacket/test_wp_ortho_loss_engine.cpp`,
registered in `inq-stack/tests/include/engine/CMakeLists.txt`. Known cases:

1. WP far from any occupied state ⇒ `removed_weight` ≈ 0, `norm_pre_renorm` ≈ 1.
2. WP deliberately overlapped with an occupied manifold ⇒ `removed_weight` > 0
   and the two independent routes of §5.1 agree to ~1e-10.
3. `norm_after` ≈ 1 regardless (renormalisation still applied) — guards the
   back-compat contract.

### 5.3 Injection scan program

`ResearchProject/systems/localised_jellium/scripts/wp_highdensity_sv/inject_scan/run.cpp`
— loads the production GS, injects at `LJ_LAUNCH_Z`, **propagates nothing**, and
writes per-trial: removed weight, Σ|overlap|², per-occupied-state overlap
spectrum, real-space moments (centroid, density std) and momentum moments,
plus the t = 0 WP orbital VTI for the k_z profile. Minutes per trial.

### 5.4 Gaussianity post-processing (reported, non-vetoing)

`inqview.load_vti(..., array=...)` extension + a scan analysis script producing,
per trial: |ψ̃(k_z)|² marginal vs the analytic Gaussian N(k₀, σ_p²), fit R²,
residual, skew/kurtosis, and the real-space |ψ(z)| profile against σ/√2.

### 5.5 Dispatch

- `shared/bin/run-wp-hd-scan.slurm` — the 0.5 Bohr outward scan at v = 2.0,
  stopping at the first trial below 3 %.
- `run-wp-hd-wp.slurm` / `run-wp-hd-vac.slurm` — accept `LJ_LAUNCH_Z` /
  `WPC_LAUNCH_Z` and a `nl_` (near-launch) run-name prefix, leaving every
  existing name resolving unchanged.
- `shared/bin/submit-wp-hd-nearlaunch.sh` — scan → sweep(0–3) → vac → notebooks.

### 5.6 Analysis

Extend `hypotheses/wp_highdensity_sv/wp_hd_stopping.py` for the near-launch
campaign and add the decisive figure: **S(v) for far-launch vs near-launch at
identical σ_WP = 0.5**, with the σ = 2 / 3 far-launch traces for context.
Per-run notebooks carry the mandatory density-matrix GIF
(`.claude/rules/notebook-density-gif.md`).

---

## 6. Validation status

| Item | Status |
|---|---|
| Slab geometry / softening = 1.0 Bohr | VERIFIED (`slab_n100_L35x35x85.hpp:47,49`) |
| ETRS does not re-orthogonalise | VERIFIED (grep over `inq-study/src/real_time/`) |
| `LJ_LAUNCH_Z` already a knob | VERIFIED (`wp/run.cpp:177`) |
| Bath spill-out profile | MEASURED (GS VTI, table §3) |
| Removed weight ~1 % at z = −14 | **INFERENCE ONLY** — the scan measures it |
| `removed_weight` two-route agreement | NOT YET — §5.2 |
| Near-launch physics | NOT YET |

---

## 7. Open / deliberately excluded

- No classical twin at the near-launch point. The classical benchmark projectile
  does not disperse, so its σ_pot = 0.354 is already its arrival width — the
  existing classical curve **is** the correct narrow-arrival reference. Stated,
  not re-run.
- σ = 2 and 3 near-launch twins not run (user: σ = 0.5 only). Would upgrade the
  falsification test to a full 2-D (σ × launch distance) collapse test.
