# Handover: Jellium folder reorganisation

## Current status

**Pending compile gate (Step 8)** — code generation phase complete.
- ✅ Step 1: plan written (`docs/plans/jellium_reorg.md`)
- ✅ Step 2: this handover
- ✅ Step 3: legacy material moved to `Tutorial/jellium-legacy/`; `jellium-wp-rt/hypotheses/` promoted to `jellium/hypotheses/`
- ✅ Step 4: skeleton dirs created (`shared/{configs,cpp}`, `scripts`, `save_gs`, `checkpoints`, `configurations/jellium_wp_rt_base`)
- ✅ Step 5: 4 shared C++ headers ported (`results_paths.hpp`, `eigenvalues_writer.hpp`, `leed_screen_layout.hpp`, `run_template.hpp`)
- ✅ Step 6: `shared/configs/base.hpp` + `save_gs/gs_L40_cubic_N38/run.cpp`
- ✅ Step 7: 6 variant headers + 7 run.cpp wrappers + `save_gs/gs_L40_cubic_N40/run.cpp`
- 🟡 Step 8: compile gate — awaiting user decision on validation menu

Auto-mode is active; the assistant is executing Steps 2–7 without
intermediate confirmation, then stopping at Step 8 (compile gate) per
`.claude/rules/testing.md` ("Suggest test options to the user before
running expensive simulations").

## What changed

**Documents created:**
- `docs/plans/jellium_reorg.md` — full redesign plan.
- `docs/handovers/jellium_reorg.md` — this handover.

**Legacy moved out of `ResearchProject/systems/jellium/`:**
- `01_ground_state/` → `Tutorial/jellium-legacy/01_ground_state/`
- `02_ground_state_convergence/` → `Tutorial/jellium-legacy/02_ground_state_convergence/`
- `03_free_gaussian_wp_propagation/` → `Tutorial/jellium-legacy/03_free_gaussian_wp_propagation/`
- `jellium-analytical/` → `Tutorial/jellium-legacy/jellium-analytical/`
- `jellium-wp-rt/run_0{1..7}_*/` → `Tutorial/jellium-legacy/jellium-wp-rt/`
- `jellium-wp-rt/{compare_observables.py, jellium_hypotheses.py, jellium_spectra.py, jellium-wp-rt.log, run_all_wp_rt.sh, results}` → `Tutorial/jellium-legacy/jellium-wp-rt/`
- `jellium-wp-rt/hypotheses/` → **promoted** to `ResearchProject/systems/jellium/hypotheses/` (legacy figures still in place; will be regenerated from new runs).

**New canonical tree at `ResearchProject/systems/jellium/`:**
- 4 shared C++ headers under `shared/cpp/`
- 7 variant config headers under `shared/configs/` (base + 6 variants)
- 7 production run dirs (`run_base`, `run_E50_s0p53`, `run_E400_s0p53`, `run_E200_s0p53_tilt45`, `run_E200_s2p0`, `run_E200_s0p265`, `run_E200_s0p53_N40`)
- 2 GS save dirs (`save_gs/gs_L40_cubic_N38`, `save_gs/gs_L40_cubic_N40`)
- Empty skeleton: `checkpoints/`, `configurations/jellium_wp_rt_base/`, `scripts/`

**Outstanding (post-compile-gate):**
- Step 9: port `coronene_postprocess.py` → `scripts/jellium_postprocess.py`, including the new overlap-matrix visualisations.
- Steps 10–15 per `docs/plans/jellium_reorg.md` §11.

## Files touched

**Plans / handovers:**
- `/local/data/public/skcb2/tddft/docs/plans/jellium_reorg.md`
- `/local/data/public/skcb2/tddft/docs/handovers/jellium_reorg.md`

**Shared C++ library (new):**
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/cpp/results_paths.hpp`
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/cpp/eigenvalues_writer.hpp`
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/cpp/leed_screen_layout.hpp`
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/cpp/run_template.hpp`

**Shared variant configs (new):**
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/configs/base.hpp`
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/configs/E50_s0p53.hpp`
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/configs/E400_s0p53.hpp`
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/configs/E200_s0p53_tilt45.hpp`
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/configs/E200_s2p0.hpp`
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/configs/E200_s0p265.hpp`
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/configs/E200_s0p53_N40.hpp`

**save_gs (new):**
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/save_gs/gs_L40_cubic_N38/run.cpp`
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/save_gs/gs_L40_cubic_N40/run.cpp`

**run_*/run.cpp (new):**
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/run_base/run.cpp`
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/run_E50_s0p53/run.cpp`
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/run_E400_s0p53/run.cpp`
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/run_E200_s0p53_tilt45/run.cpp`
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/run_E200_s2p0/run.cpp`
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/run_E200_s0p265/run.cpp`
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/run_E200_s0p53_N40/run.cpp`

## Commands run

- `mv` for the legacy moves (Step 3) — plain `mv`, not `git mv`, because the
  legacy on-disk files were untracked. The git index already shows old
  paths as `D` (deleted) from a prior reorganisation. The full set of
  staging changes (untracked-add for the new tree + delete for stale tracked
  paths) is deferred to a post-compile-gate commit.
- `mkdir -p` for skeleton dirs (Step 4).

## Tests and validation

No tests run yet. Validation menu is in `docs/plans/jellium_reorg.md` §10.
The hard gate is **Step 8**: compile-only (`inq-run --build-only` or
equivalent) for `save_gs/gs_L40_cubic_N38/run.cpp` and `run_base/run.cpp`,
plus header unit checks for `jellium::results::*` paths and
`compute_screen_window`. The assistant will pause here and ask the user
which expensive validations (#4 GS sanity, #6 short prop, #7 full base, #8
full suite) to authorise.

## Trusted sources used

- Coronene canonical layout: `ResearchProject/systems/coronene/` (in-repo).
- Coronene `run_template.hpp`, `results_paths.hpp`, `leed_screen_layout.hpp`,
  `eigenvalues_writer.hpp`, `tsubonoya_2014_base.hpp` — used as port targets.
- Existing jellium numerics frozen as the new Base from
  `ResearchProject/systems/jellium/jellium-wp-rt/run_01_base/run.cpp`
  (parameters: L=40 bohr, N=38, σ=0.53 Å, E=200 eV, +z, dt=0.020 a.u.,
  N_steps=417). **Internal-reference, not literature** (no published paper
  on file pinning this exact configuration).
- `inqkit::observables::OrbitalOverlapMatrix::snapshot()` — implementation in
  `inq-stack/include/inqkit/observables/orbital_overlap.hpp`. Computes the
  full n_ref × (n_ref+1) overlap matrix per step, header-documented.

## Attribution notes

- New shared library is a port of coronene's. Each ported header has a top
  comment crediting the coronene source headers it was adapted from.
- Tsubonoya, Hu & Watanabe, *Phys. Rev. B* **90**, 035416 (2014) is **not**
  used by jellium runs. The paper-window mechanism is preserved only for
  template compatibility; jellium's `T1_AU = 0`, `T2_AU = DT_AU * N_STEPS`.
- Per repo memory rule, no Claude attribution in commit messages.

## Known issues / blockers

- None at this point.
- Pre-existing `Tutorial/jellium-wp-rt/` and `Tutorial/angelo-jellium/` are
  separate user-created folders and will be left untouched. The new legacy
  destination is `Tutorial/jellium-legacy/` (distinct path).
- After Step 3, `Tutorial/jellium-legacy/jellium-wp-rt/run.cpp` files
  diverge from `Tutorial/jellium-wp-rt/` (already confirmed via `diff -rq`).

## Assumptions still in play

(All from `docs/plans/jellium_reorg.md` §12; all defaults are taken since
the user authorised "start the execution".)

- A. Open-shell variant uses N=40 (matches legacy `run_07_open_shell`).
- B. Closed-shell variants 02–06 share `checkpoints/gs_L40_cubic_N38`;
  variant 07 uses `checkpoints/gs_L40_cubic_N40`.
- C. Tilted variant keeps L=40 cubic, only `WP_KX/WP_KZ` change.
- D. Cell stays `.periodic()` Γ-only (physics — homogeneous gas).
- E. `hypotheses/` regenerated from new runs; legacy figures not carried over.
- N_STEPS values are inherited from the legacy variants without re-derivation
  (loop-back margin already hand-tuned for cubic-periodic). Future work:
  derive at compile time via `compute_n_steps()` adapted to periodic.
- VTI metadata for `[0, L]` cubic-periodic: assumed correct in
  `inqkit::io::RealField3DWriter`. To be verified during Step 9
  (post-processor visualisation review).
- `T1_AU/T2_AU` for jellium set to full-propagation; the per-screen physics
  window mechanism in `compute_screen_window` is otherwise preserved as in
  coronene. Re-derivation for periodic boundaries is a follow-up.

## Exact next steps

1. Mark Step 2 complete; start Step 3 (legacy moves).
2. `git mv` legacy folders into `Tutorial/jellium-legacy/` per
   `docs/plans/jellium_reorg.md` §3 mapping table. Specifically:
   ```
   ResearchProject/systems/jellium/01_ground_state               → Tutorial/jellium-legacy/01_ground_state
   ResearchProject/systems/jellium/02_ground_state_convergence   → Tutorial/jellium-legacy/02_ground_state_convergence
   ResearchProject/systems/jellium/03_free_gaussian_wp_propagation → Tutorial/jellium-legacy/03_free_gaussian_wp_propagation
   ResearchProject/systems/jellium/jellium-analytical            → Tutorial/jellium-legacy/jellium-analytical
   ResearchProject/systems/jellium/jellium-wp-rt/run_0{1..7}_*   → Tutorial/jellium-legacy/jellium-wp-rt/run_0{1..7}_*
   ResearchProject/systems/jellium/jellium-wp-rt/{*.py,*.sh,*.log,results} → Tutorial/jellium-legacy/jellium-wp-rt/
   ```
3. Promote `jellium-wp-rt/hypotheses/` up one level to
   `ResearchProject/systems/jellium/hypotheses/`. Per Assumption E, delete
   stale figures (PNGs) and keep only README.md + metadata.csv if any.
4. Create skeleton dirs (Step 4).
5. Steps 5–7 then proceed in order; pause at Step 8 for user input.

## Recovery if interrupted

To resume: re-read this handover and `docs/plans/jellium_reorg.md`,
re-run `git status` to see which moves landed, then continue from the
highest unchecked item in §"Current status".
