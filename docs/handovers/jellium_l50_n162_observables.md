# Handover: jellium WP-jellium scattering — L=50, N=162 closed-shell base run

Updated 2026-05-05 ~17:34 BST. Branch
`features/jellium-ks-energy-observables` (with the merged
`features/occupations-vs-time`).

## Current status

**Phase 2 (project-base run at the *true* closed shell N=162) is
complete in numerics and partly in postprocessing.** A journal entry
for the N=162 run has not yet been written. The pre-N=162 L=50 N=138
partial-shell run is fully postprocessed and journaled (now flagged as
"partial-shell, do not use as production base"). The earlier L=30 run
journal contains a corrective addendum noting that the "−9.69 eV ΔE_WP
unaccounted" framing was a band-structure-vs-total bookkeeping artefact,
not missing energy.

What is done:
- **OccupationsWriter** (`inq-stack/include/inqkit/observables/occupations_writer.hpp`)
  added on a feature branch, smoke-compiled, smoke-tested, merged into
  `features/jellium-ks-energy-observables`. Per-step dump of f_i to
  `results/raw/observables/occupations_vs_time.csv` at 5×WRITE_EVERY
  cadence; verified to capture the static (frozen-by-construction) f_i
  with std-over-time at machine epsilon.
- **`occupations` postprocess phase** (`inq-stack/python/inqview/postprocess/occupations.py`)
  reads the CSV and emits two animated bar GIFs
  (`occupations_absolute.gif`, `occupations_delta.gif`) with HOMO
  dashed line. Registered in `pipeline.py`.
- **`gs` phase** extended in `ground_state.py` to emit
  `analysis/ground_state/gs_occupations.png` (static GS occupations
  bar chart with HOMO line).
- **`state_energies` phase** extended to draw HOMO dashed lines on the
  four bar GIFs.
- **`fourier.py`** got `zero_pad` (default 4) and `smooth_sigma_bins`
  options for QBall-style smoothness; **`plots.py::plot_spectrum`**
  now uses `ScalarFormatter(useOffset=False)` per the no-offset
  styling rule.
- **Closed-shell magic-N audit:** previous claim "N=138 = closed
  shell at |G|²≤6" was wrong — the true closed shell at |G|²≤6 is
  **N=162** (24 spatial states from |G|²=6 × 2 spin = 48, on top of
  cumulative 114 at |G|²≤5). N=138 is a *partial-shell* fill.
  `docs/sources/free-electron-gas-magic-numbers.md` already had the
  correct table; the bug was our run-targeting picking 138.
- **Two L=50 runs completed:**
  1. `run_base_n138_L50_E1p5/` — partial-shell, fully postprocessed.
     Cleanly forward-monotonic cod_z trajectory (no PBC revival),
     6.3% velocity retardation, +2.0 eV bath kinetic gain matching
     1.64 eV WP kinetic loss. Journal entry written but flagged as
     non-canonical (partial shell).
  2. `run_base_n162_L50_E1p5/` — true closed shell. GS converged to
     E=−11.906 Ha with **81 states at f=2.000000 exactly + 20 states
     at f≈4×10⁻¹¹** (no smearing — the closed-shell signature).
     Propagation done (715 s wall, drift 1.4×10⁻⁷ eV). ΔE_WP
     ⟨H⟩ = −2.56 eV; bath kinetic gain ≈ +2.25 eV from the explicit
     Δenergy_kinetic component. WP velocity 0.373 → 0.103 Bohr/a.u.
     ⇒ KE drop 1.75 eV.
- **Postprocess for the N=162 run is in flight** — currently in the
  density phase (64/70 GIFs at 17:34). Expect screens phase next
  (~20 min), then overlap (~1 min). Total ETA ~50-60 min.

What is partially done:
- **Journal entry for the N=162 run** is queued (task 49). Will use
  the kinetic+Hartree+xc decomposition rather than the misleading
  band-structure framing of the L=30 entry.

What is not done:
- `energy_balance` postprocess phase (task 41). The numbers are easy
  to extract by hand (see Section "Tests and validation" below); the
  phase formalises this into a single auto-generated PNG.
- Read about WP revival dynamics (task 42).
- Confirm the partial-fill interpretation of the L=30 hole-behind-WP
  observation against the closed-shell N=162 result.

## What changed (since the last handover)

- Renamed and corrected the closed-shell convention:
  N=138 → N=162. The Cfg `Base_N162_L50_E1p5` derives from
  `Base_N138_L50_E1p5` and only overrides `N_ELECTRONS = 162`.
- New observable: per-step occupations dump.
- Postprocess: HOMO line on KS-energy GIFs; static GS occupations
  chart; new occupations animations; FFT smoothing; no-offset axis.
- Two corrective notes appended to the L=30 journal entry (Section 3
  energy bookkeeping clarification).
- Two L=50 runs (N=138 and N=162) executed end-to-end.

## Files touched

C++:
- `/local/data/public/skcb2/tddft/inq-stack/include/inqkit/observables/occupations_writer.hpp` (new)
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/cpp/run_template.hpp` (wire in occupations_wr)
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/cpp/results_paths.hpp` (add `occupations_csv_path()`)
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/configs/base_n138_L50_E1p5.hpp` (new)
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/configs/base_n162_L50_E1p5.hpp` (new)
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/save_gs/gs_L50_cubic_N138_dx1p0/run.cpp` (new)
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/save_gs/gs_L50_cubic_N162_dx1p0/run.cpp` (new)
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/run_base_n138_L50_E1p5/run.cpp` (new)
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/run_base_n162_L50_E1p5/run.cpp` (new)

Python (inqview):
- `/local/data/public/skcb2/tddft/inq-stack/python/inqview/postprocess/occupations.py` (new)
- `/local/data/public/skcb2/tddft/inq-stack/python/inqview/postprocess/state_energies.py` (HOMO line + `_read_homo_index`)
- `/local/data/public/skcb2/tddft/inq-stack/python/inqview/postprocess/ground_state.py` (gs_occupations chart)
- `/local/data/public/skcb2/tddft/inq-stack/python/inqview/postprocess/pipeline.py` (register `occupations` phase)
- `/local/data/public/skcb2/tddft/inq-stack/python/inqview/postprocess/eigenvalues_gs.py` (schema-tolerance adapter — earlier turn)
- `/local/data/public/skcb2/tddft/inq-stack/python/inqview/postprocess/gamma_transitions.py` (same adapter)
- `/local/data/public/skcb2/tddft/inq-stack/python/inqview/postprocess/layout.py` (parse `<L>^3` cell)
- `/local/data/public/skcb2/tddft/inq-stack/python/inqview/fourier.py` (zero_pad, smooth_sigma_bins)
- `/local/data/public/skcb2/tddft/inq-stack/python/inqview/plots.py` (ScalarFormatter)

Docs:
- `/local/data/public/skcb2/tddft/docs/observables_reference.md` (§13 jellium addendum incl. occupations dump rules)
- `/local/data/public/skcb2/tddft/docs/journals/researchproject/2026-05-05_run_base_n138_L30_E5.md` (corrective addendum)
- `/local/data/public/skcb2/tddft/docs/journals/researchproject/index.md`
- `/local/data/public/skcb2/tddft/.claude/rules/jellium-base-run-spec.md` (project-base rule, earlier turn)

Run output trees (preserved):
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/run_base_n138_L50_E1p5/results/`
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/run_base_n162_L50_E1p5/results/` (postprocess in flight)

Lost (do not exist anymore):
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/run_base_n138_L30_E5/results/`
  was destroyed in a `mv` race during the smoke-test of the new
  observable. The journal attachments are preserved at
  `docs/journals/researchproject/attachments/2026-05-05_run_base_n138_L30_E5/`.
  Re-running takes ~3 min propagation + ~70 min postprocess.

## Commands run

GS + propagation (chained):
```
cd /local/data/public/skcb2/tddft/ResearchProject/systems/jellium/save_gs/gs_L50_cubic_N162_dx1p0
nohup bash -c '
  CUDA_VISIBLE_DEVICES=0 timeout 1800 inq-run > gs_run.stdout 2>&1
  if grep -q "Done." gs_run.stdout; then
    cd /local/data/public/skcb2/tddft/ResearchProject/systems/jellium/run_base_n162_L50_E1p5
    CUDA_VISIBLE_DEVICES=0 timeout 14400 inq-run > full_run.stdout 2>&1
  fi
' &
```

Postprocess:
```
source /local/data/public/skcb2/tddft/venv/bin/activate
python3 -c "
from inqview.postprocess import pipeline
from pathlib import Path
res = pipeline.run(
  Path('ResearchProject/systems/jellium/run_base_n162_L50_E1p5/results'),
  run_name='run_base_n162_L50_E1p5',
  skip_paraview=True)
"
```

Verification (one-off check used to derive the energy-balance numbers
quoted below):
```
python3 ResearchProject/systems/jellium/scripts/verify_smoke_outputs.py \
  ResearchProject/systems/jellium/run_base_n162_L50_E1p5/results
```

Branch admin:
```
git checkout -b features/occupations-vs-time
# implement OccupationsWriter, smoke compile, smoke test
git add ...
git commit -m "add occupations-vs-time observable; new L=50 E=1.5 base config"
git checkout features/jellium-ks-energy-observables
git merge --no-ff features/occupations-vs-time -m "merge: ..."
```

## Tests and validation

**Run-base verification on N=162 (L=50)** — 10/10 standard checks pass:
- `run_completed = true`, final t=30 a.u.
- cod_z monotonically forward (no PBC revival).
- density_l2 starts at 0, grows.
- 89 unique state energies (now 101 at N=162).
- WP momentum peak at t=0 close to k₀.
- All raw VTI series populated.

**Energy decomposition (N=162 L=50, occ-weighted, all in eV):**
- Δenergy_total = 1.4 × 10⁻⁷ (conservation ✓)
- Δenergy_kinetic = +0.50
- Δenergy_hartree = −0.60
- Δenergy_xc = +0.10
- ΔE_WP ⟨H⟩ = −2.56 (state_energies, single state f_WP=1)
- Σf_iΔε_i bath = +1.99
- WP velocity drop: 0.373 → 0.103 Bohr/a.u. ⇒ KE drop 1.75 eV
- Bath kinetic gain (system Δkinetic minus WP kinetic drop) = +2.25 eV

**OccupationsWriter audit:** max std over time = 9 × 10⁻¹⁶ (machine
epsilon); 65/101 states bit-identical across all timesteps. INQ TDDFT
propagates with frozen f_i by construction; the dump captures this
correctly.

**GS occupations at N=162:** **81 states at f = 2.000000 exactly,
20 at f ≈ 4 × 10⁻¹¹** — the closed-shell signature. (Compare with
N=138: states 60-80 had occupations 0.95-1.0 from band-edge smearing.)

## Trusted sources used

- Ashcroft & Mermin, *Solid State Physics*, Ch. 2 (free-electron gas
  shell structure).
- Robinett, "Quantum wave packet revivals," *Phys. Rep.* 392, 1
  (2004) — to be read; required for interpreting the L=30 cod_z
  reversal.
- INQ source: `inq/src/observables/`, `inq/src/operations/`, the
  `viewables.hpp` extension we ported (`state_energy_writer.hpp`).
- `docs/sources/free-electron-gas-magic-numbers.md` (already correct).

## Attribution notes

- `state_energy_writer.hpp` ports the algorithm from the professor's
  `viewables.hpp` (`projected_occupation_array` + state-energy methods).
  Header comment cites the original.
- The `occupations_writer.hpp` is project-original (no upstream
  reference — INQ stores f_i but doesn't expose a per-step writer).

## Known issues / blockers

- **N=162 postprocess is still running** (in density phase, 64/70
  GIFs at 17:34). Will probably complete by ~18:30 BST.
- **WP injection `max_overlap` for the L=50 runs is 0.084** (i.e. 8 %
  WP-bath overlap before orthogonalisation). The wide WP (σ=5)
  has substantial overlap with bath plane-wave states, and the
  one-shot Gram-Schmidt absorbs it. Not a bug — but should be noted
  in the journal entry as a caveat for the energy balance.
- **`gamma_transitions` postprocess output is informative but
  occupation thresholds for "occupied / unoccupied" need tuning at
  the N=138 partial-shell case.** At N=162 the thresholds work fine
  out of the box (f=2 vs f=0).
- **Lost run data:** `run_base_n138_L30_E5/results/` was destroyed.
  Journal + attachments preserved.

## Assumptions still in play

- INQ's TDDFT keeps occupations f_i FROZEN. Confirmed by audit. If
  this changes (e.g. a future propagation method that updates f_i),
  the OccupationsWriter would automatically capture the time
  evolution.
- The "WP slowdown" in the L=50 runs is genuinely inelastic, NOT a
  PBC revival. The cleanly monotonic cod_z trajectory (vs the
  reversal at L=30) supports this — but a formal kinematic
  comparison vs the Robinett-revival timescale is not yet done.
- The N=162 closed-shell GS is uniform on the integrated grid —
  uniformity verified empirically in the GS density slice for N=138
  L=50 (was already very flat); should be even cleaner at N=162.
- The energy-balance bookkeeping holds for the WP slot the same way
  as for bath orbitals — i.e. the WP's f_WP=1.0 is treated identically
  to a half-filled occupied state in `Σ f_i Δε_i`. We have not
  formally verified this against INQ's own internal accounting.

## Exact next steps

(Updated 2026-05-05 evening — journal entries and topical analysis
landed.)

### Done in this session

- ✓ Journal entry written for `run_base_n138_L50_E1p5` →
  `docs/journals/researchproject/2026-05-05_run_base_n138_L50_E1p5.md`.
- ✓ Journal entry written for `run_base_n162_L50_E1p5` →
  `docs/journals/researchproject/2026-05-05_run_base_n162_L50_E1p5.md`.
- ✓ Topical journal entry on plasmons / e-h regime classification →
  `docs/journals/researchproject/plasmons-and-stopping-power.md`. Substitutes
  L=30 / L=50/N=138 / L=50/N=162 parameters into Correa 2018 Eq. (3)-style
  thresholds; concludes WP velocity at v_F is sub-threshold for plasmon
  excitation and ~14× above the e-h threshold; suggests bath
  polarisation-cloud mechanism, not "hole-as-attractor", as the slowdown
  driver.
- ✓ `docs/journals/researchproject/index.md` updated.
- ✓ Source note `docs/sources/correa-2018-electronic-stopping-power.md`
  written, cross-referenced from the topical entry.
- ✓ Plan `docs/plans/jellium_orthonormalisation_rerun.md` for the
  professor-PDF orthonormalisation rerun (Options A / B / C, three
  verdict criteria, a companion per-orbital cod tracking phase).
- ✓ Report `docs/reports/qball-spectra-comparison.md` documenting the
  diff between QBall `td_kicks` and `inqview/fourier.py`. Headline:
  inqview defaults to linear-detrend whereas QBall uses plateau-mean;
  the latter is what Correa 2018 §6 means by "exclude the transient".
  Recommends `t_start_au` parameter + `detrend_strategy="plateau_mean"`
  default for jellium runs.
- ✓ New rule `.claude/rules/journal-entries.md`: every run-based
  journal entry must paste `run_summary.txt` as a markdown table.
- ✓ `docs/observables_reference.md` extended with §13.5 (per-component
  energy graphs, mandatory) and §13.6 (transient exclusion, mandatory).

### Outstanding (next session)

1. **`energy_balance` postprocess phase** (task 41). Numerical
   evaluation has now been done by hand; the standalone phase remains
   to be packaged. Numbers found:
   - **N=162 L=50:** Δ(Σf_iε_i) = −0.5687 eV, ΔE_WP = −2.5553 eV,
     Σf_iΔε_i bath = +1.9866 eV, **Δ(∫v_xc·n) = +0.1243 eV** (closure),
     verifying the "+0.12 eV prediction" recorded in
     `plasmons-and-stopping-power.md §7` to within 10⁻⁴ eV.
   - **N=138 L=50:** Δ(Σf_iε_i) = −0.4289 eV, ΔE_WP = −2.6167 eV,
     Σf_iΔε_i bath = +2.1878 eV, Δ(∫v_xc·n) = +0.1278 eV.
   - **Phase work remaining:** wire the closure-equation into a
     postprocess script that auto-emits the table; OR (preferred) add
     a one-shot ∫v_xc·n GPU reduction to the run template so the
     measurement is direct, not inferred.
2. **Implement transient-exclusion FFT** per
   `docs/reports/qball-spectra-comparison.md §4.1`–§4.2. After landing,
   regenerate the spectra in all three open journal entries and update
   their §"Spectra — caveat" sections.
3. **Run the orthonormalisation rerun** per
   `docs/plans/jellium_orthonormalisation_rerun.md`. Start with
   Option A (CGS2 / twice-is-enough). Verdict report:
   `docs/reports/orthonormalisation-rerun-verdict.md`.
4. **Per-orbital cod-vs-time postprocess** (the jellium analogue of
   Correa Eq. (10)) per the orthonormalisation plan §4. Emits
   `analysis/observables/orbital_cod_z_vs_time.png`. Resolves
   "is cod_z slope = WP velocity?".
5. Read the WP-revival papers (Robinett 2004; Aronstein–Stroud 1997)
   and decide whether the L=30 reversal needs a separate analysis.
6. Decide whether to re-run the L=30 E=5 propagation to recover the
   destroyed `results/` tree.
7. Plan the next physics step (vary WP energy / σ / box size? move
   to a multi-WP run? do a v > v_th^plasmon ≈ 0.94 a.u. run that can
   actually excite the bulk plasmon kinematically?).
8. **Positive-ion classical comparison run.** Scaffolded:
   - Plan: `docs/plans/jellium_positive_ion_companion.md` (Methods,
     Option G1 vs G2 GS-prep, predicted observables table, verdict
     criteria, build/launch checklist).
   - Cfg: `ResearchProject/systems/jellium/shared/configs/positive_ion_L50_v0p33.hpp`
     (struct `Positive_Ion_L50_v0p33` deriving from `Base_N162_L50_E1p5`,
     no WP, single H ion at (0,0,−L/4) at v_z = v_F = 0.332 a.u.).
   - Run.cpp draft: `ResearchProject/systems/jellium/run_positive_ion_L50_v0p33/run.cpp`
     adapted from QBall-INQ Li `run_propagate_v0p0123_extensive/run.cpp`.
   - **Build pending.** Open API questions documented in the plan
     (N_e accounting, pseudopotential grid, dt convergence). User OK
     needed before launching the GS + propagation (~30 min wall total
     for Option G2; ~25 min for Option G1).
