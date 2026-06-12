# Plan: Overnight Gaussian classical projectile in jellium — S(v) mapping

Status: LOCKED via grill-with-docs session 2026-06-12. Ready for unsupervised
execution. Supersedes the pasted external prompt wherever they conflict (all
conflicts resolved below). One rolling handover will track execution at
`docs/handovers/overnight-gaussian-classical-jellium.md`.

## Goal

Implement a controlled erf/σ-smoothed **classical** electron projectile in the
existing rt-TDDFT jellium box, validate it, measure the q-resolved loss function
via a finite-q kick, and map electronic stopping power S(v) with Ehrenfest
dynamics. Secondary deliverable: **Subtask-4 (finetune) dogfooding data** for the
ecosystem rejuvenation — log where the codified skills/rules/hooks help or
collide with an external plan.

## Locked decisions (grilling outcomes)

1. **Density: r_s = 5.74 is authoritative** (Section 1 of the prompt). L = 50
   Bohr, N = 162, dx = 0.40, ~100 K smearing. **Reuse the existing GS** at
   `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/checkpoints/gs_L50_cubic_N162_dx0p40`
   — NO rebuild. Every r_s=2 / L=17.575 number in the prompt's Stages 1.2/4/5
   and the "≈35 Bohr" two-traversal figure is a **stale leftover** and is
   rewritten to r_s=5.74. Authoritative derived constants (log in report header,
   use in ALL checks):
   - k_F = 0.334, v_F = 0.334 a.u., E_F = 0.0558 Ha, ω_p = 0.1262 Ha
   - Fermi λ = 18.8 Bohr, Friedel λ = π/k_F = 9.4 Bohr, plasmon period = 49.8 a.u.
   - Landau onset q_c ≈ ω_p/v_F ≈ 0.38 a.u.
   - q_cut = π/0.40 = 7.85 a.u.; σ·q_cut = 3.9 (σ=0.5, resolvable)
   - Plasmon peak prediction at q₁: ω_p(1+3q²v_F²/(10ω_p²)) ≈ **0.13 Ha** (NOT 0.67)
   - Caveat (state plainly): at r_s=5.74 ALDA is less reliable and the box fits
     only a few screening lengths — trust trends/slopes/internal consistency,
     not absolute RPA agreement.
   - Path length = **ONE traversal = 50 Bohr** at L=50 (user decision; NOT two).

2. **Branch:** merge `rejuvenation/claude-ecosystem` → `main` FIRST (it is
   COMPLETE per its plan; user-gated, now authorised), then create
   `overnight-gaussian-classical` off `main`. All overnight research lives there.
   Handle stray uncommitted changes first: `inqview.egg-info` diffs are
   regenerated build metadata (discard/ignore); inspect and commit
   `.claude/settings.json` properly.

3. **Commit messages:** use the `commit-messages` rule format
   `action(scope): description` (NOT the prompt's `stage-N:`); scope = `jellium`
   or `inqkit`/`inqview`; no forbidden words; no Co-Authored-By trailer; split
   research vs internal-config commits.

4. **Stage-4 finite-q kick:** implement in **inqkit C++** as a new header
   `inq-stack/include/inqkit/perturbations/cosine_kick.hpp`, modeled on
   `inq/src/perturbations/kick.hpp`. INQ's kick is uniform (q=0) only, but
   orbitals are mutable in real space via `electrons.kpin()[ik].hypercubic()
   [ix][iy][iz][ist]`. Multiply every occupied orbital by exp(i·η·cos(q·r)) once
   at t=0, box-commensurate q = 2πn/L along z. NOT a Python `kick_q.py`.
   - **INQ duplicate is pre-authorised** as a genuine last resort: only if the
     inqkit-only path truly fails do we make a renamed physical INQ duplicate in
     `tddft/`, edit only the duplicate (original byte-for-byte untouched),
     repoint build config, and document every change in a report section "INQ
     duplicate modifications". Attempt inqkit-only first.

5. **Provenance/reporting → full ecosystem mapping** (NOT the prompt's bespoke
   artifacts):
   - `runs.json` → **tddft-run-catalogue** skill (CSV) — upsert each run.
   - per-run `analyse.py` + `REPORT.md` (full inqview pipeline) per
     `feedback_per_run_analyse_py`.
   - report → `docs/reports/overnight-gaussian-classical-jellium/` markdown
     (report-writing) with PNG figures via the canonical theme
     `inqview.visualisation.style` (ADR 0004, report-figures).
   - **PLUS an executed top-to-bottom `.ipynb`** in that same reports dir (user
     wants to move run-reports toward notebooks; PNGs still go through the
     theme). NOT at repo root.
   - rolling handover (handover-update) instead of root `MORNING_REPORT.md`;
     the one-screen morning summary becomes the handover's top section.
   - journal entries (journal-writing) with verbatim `run_summary.txt`.
   - `docs/sources/` notes (literature-review) for Lindhard 1954, Lindhard &
     Winther 1964, Echenique et al. 1981/1986, Correa 2018.

6. **Launch geometry:** velocity along **+z, on-axis** (x=y=0), matching the
   existing classical runs — reuses `analyse.py` / center-of-density / density
   slicing wholesale. NO transverse offset (symmetric-path risk accepted; watch
   via the piecewise-slope self-wake check). The prompt's x-axis +
   (0.13,0.27,0.41)·L is dropped. The `boundary_rule` (4σ/1σ) does NOT apply to
   a σ=0.5 classical projectile (it is a WP-envelope rule); use a fixed launch
   offset near the −z face like the existing E100 classical run.

7. **GPU:** default GPU; pair Stage-6 runs two-at-a-time via `build-run`
   Multi-GPU section (`CUDA_VISIBLE_DEVICES=0/1`). Pick lowest-util AND
   lowest-committed-mem GPUs; never steal a parked job. If `nvidia-smi` throws
   the NVML driver mismatch, that is monitoring-only — verify with the cuda
   probe and run on GPU anyway; never fall back to CPU
   (`reference_gpu_driver_mismatch`, `feedback_gpu_default_expectation`).

8. **Projectile mass & dynamics (CORRECTED after user grilling — the plan hid a
   mass change; this overrides the prompt):**
   - **Mass = m_e = 1.0 a.u.** (the existing classical-run value). The prompt's
     M = 1836 is **struck** — it was a fictitious mass to fake constant velocity.
   - **Dynamics = free `ehrenfest`** (`options::real_time{}.ehrenfest()`), the
     existing classical convention. The −1 charge **genuinely decelerates**
     ("electronic forces decelerate ion"). NOT constant-velocity. (INQ's
     `ionic::propagator::impulsive` would give true constant-v mass-free, but
     the user chose free dynamics — the physically literal real-electron case.)
   - **S(v) extraction (replaces the prompt's linear-fit method):** each run
     sweeps a velocity RANGE as it decelerates, so S(v) = **local**
     −d(KE_proj)/dx evaluated pointwise and **binned by instantaneous v(t)**,
     cross-checked against +d(E_electrons)/dx; total-energy drift = integrator
     health. NO single linear slope, NO fixed "2-traversal" assumption.
   - **Path/N_STEPS per run:** cap at **ONE traversal = 50 Bohr** (user decision)
     OR until the projectile stops, whichever first. N_STEPS computed per run
     from v₀ and the 50-Bohr cap, not fixed. (Removes any self-wake re-entry on
     the second pass entirely.)
   - **M = 18360 duplicate DROPPED** — a mass-independence check only meaningful
     in the constant-v limit; meaningless under free m_e dynamics.

8b. **Charge sign — VERIFIED electron, not proton.** Existing psp `PP_LOCAL` is
   sign-flipped and **positive** (+9.06 Ha core) ⇒ repulsive to bath electrons ⇒
   −1 charge. New erf psp must be **positive** `+erf(r/(σ√2))/r`, V(0)=+1.60 Ha.

8c. **Projectile spatial model — σ = 0.5 Bohr Gaussian charge (Option B).**
   Replaces the existing near-point ONCV psp (V(0)≈9 Ha) with a controlled
   softened charge (V(0)=1.60 Ha). σ=0.4 is the sensitivity duplicate. ALL
   references carry the form factor exp(−q²σ²). Existing near-point psp retired
   from production.

8d. **Inherited modeling facts (not changed tonight, recorded for clarity):**
   projectile = INQ species **"H"** (label only) + custom UPF + m_e override;
   **z_valence = 0** ⇒ pure external potential, contributes NO charge to
   SCF/Hartree/Poisson and no electron to the count (N stays 162); cell stays
   neutral (162 e⁻ + jellium background); the −1 charge is felt by the bath only
   through V_loc. This matches all existing classical runs.

9. **Run-menu approved as scoped** for unsupervised launch (validation gate
   satisfied): validation runs (static σ=0.5 impurity GS, dt check) → Stage-4
   kicks (2 linearity + 4 production q₁–q₄) → Stage-6 ladder paired on 2 GPUs →
   Stage-7 k-points (non-blocking). Drop-order if behind: Stage 7 → highest-v
   runs → σ=0.4 duplicate → kicks 4→2. Stage 8 report = final 45 min.

   **Stage-6 ladder (decel-aware, free Ehrenfest, m_e):** initial velocities
   v₀ (a.u.) = **{3.0, 2.0, 1.3, 0.8, 0.6, 0.4}** (decelerating sweeps tile the
   v axis; high v₀ first-traversal-only by finite size; v₀≤0.4 for the friction
   onset) **+ σ=0.4 duplicate at v₀=1.0** = **7 production runs**. v₀=2.5, 3.0
   capped/flagged for wake≫L; v₀=0.2 dropped (stops instantly, sub-v_F). All
   read instantaneous v(t); S binned by v.

## File placement (resolved against local rules)

| Artifact | Location |
|---|---|
| Generated psp σ=0.5, σ=0.4 | `ResearchProject/systems/jellium/shared/pseudopotentials/electron_gaussian_sigma{σ}.upf` |
| psp generator | `inq-stack/python/inqview/io/` (numpy-only util) |
| Lindhard refs, S_LR(v;σ) | `inq-stack/python/inqview/analysis/` (numpy/scipy) |
| Stopping extraction | `inq-stack/python/inqview/analysis/` |
| Finite-q cosine kick | `inq-stack/include/inqkit/perturbations/cosine_kick.hpp` |
| Run configs (.cpp + Cfg) | `ResearchProject/systems/jellium/run_<...>/` + `shared/configs/*.hpp` |
| Report (md + ipynb + PNG) | `docs/reports/overnight-gaussian-classical-jellium/` |
| Handover | `docs/handovers/overnight-gaussian-classical-jellium.md` |
| Source notes | `docs/sources/` |
| Run catalogue | via tddft-run-catalogue skill (CSV) |

## Existing assets to reuse (do not reinvent)

- GS checkpoint: `.../jellium/checkpoints/gs_L50_cubic_N162_dx0p40`
- Existing sign-flipped bare-local Coulomb electron psp template:
  `.../jellium/shared/pseudopotentials/electron-ONCV-1.2.upf` (z_valence=0,
  number_of_proj=0, PP_LOCAL sign-flipped + extended to r=50; ONCV-smoothed —
  Stage 2 replaces its PP_LOCAL with a *controlled* erf(r/(σ√2))/r form).
- Twin Cfg pattern: `shared/configs/electron_proj_E100_L50_cubic.hpp` etc.
- Existing classical run template: `run_classical_n162_L50_E100/run.cpp`.

## Stage map (r_s=5.74 throughout)

S1 setup/timing → S2 erf psp + static impurity validation → S3 dt check →
S4 cosine kick + ELF + f-sum (∫ω·Im[−1/ε]dω = (π/2)ω_p²) → S5 Lindhard refs
(with exp(−q²σ²) Gaussian form factor — Option B) → S6 S(v) ladder + extraction
→ S7 k-points (non-blocking) → S8 report (reserved 45 min).

Every new module ships a known-case test BEFORE GPU use (`code-test`). No
correctness claim without a completed validation. Failed runs: capture log, mark
FAILED in catalogue, never delete outputs, continue to the next independent
stage.

## Subtask-4 (finetune) observations to log

Seed list (conflicts the grilling already surfaced — each is a data point that
the ecosystem correctly overrode the external plan):
1. Python modules mislabeled "inqkit/*.py" (inqkit is C++; Python is inqview).
2. `stage-N:` commit format vs `action(scope):` rule.
3. Root `MORNING_REPORT.md` vs `docs/reports/`/handover.
4. `runs.json` vs tddft-run-catalogue.
5. `.ipynb` report vs project .md/LaTeX norm (user chose to pioneer ipynb).
6. x-axis launch vs z-axis pipeline; `boundary_rule` misapplied to a classical
   projectile.
7. Pervasive r_s=2 contamination vs the FIXED r_s=5.74 box.
8. **(Most significant)** The plan hid a **fictitious projectile mass M=1836**
   that silently redefined the S(v) measurement (faked constant velocity). User
   caught it; reverted to m_e + free Ehrenfest. Lesson for Subtask-4: external
   plans can smuggle physics changes inside an innocuous-looking parameter — a
   full parameter ledger (real-value vs plan-value vs meaning) before launch is
   the mitigation, and should arguably be a standing step in `tddft-simulations`.
Append any new friction/smoothness observations during execution.
