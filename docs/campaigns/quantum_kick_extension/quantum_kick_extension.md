---
id: quantum-kick-extension
area: quantum_kick_extension
title: "Li multi-k kick: omega_peak(v) drift verification"
status: ready
hypothesis: "The non-monotonic drift of Li's peak response frequency omega_peak(v) with impulsive kick velocity, seen in the single-k Gamma-only runs, is a genuine physical effect reproduced by the multi-k near-0K (2x2x2 shifted, 400 K) configuration across the full velocity sweep -- not an artefact of Gamma-only sampling or 1000 K smearing."
handover: docs/handovers/quantum-kick-extension.md
tasks:
  - { name: "GS load check: reuse checkpoint li_54_2x2x2_T200_xyz; verify n_electrons=162", done: true }
  - { name: "Add 5 per-velocity Cfg structs (v0p0375/100/175/250/375) to base_li_54.hpp + build-once binary", done: true }
  - { name: "Pilot smoke test @ v=0.175 (~100 steps): build OK, GS loads, energy real, no NaN", done: true }
  - { name: "5 production runs (0.0375, 0.100, 0.175, 0.250, 0.375), 15500 steps each", done: false }
  - { name: "[UNGATED 2026-06-25] per-run FFT energy+dipole_x via fourier-analysis skill (mean->Hann->4x->coherent-gain, peak-in-band) -> omega_peak; gamma_transitions + state_energy_spectra attribution", done: false }
  - { name: "[UNGATED 2026-06-25] combined omega_peak(v) figure (9 multi-k) + single-k overlay + study notebook + handover + frontmatter flips", done: false }
blocked_reason: "UNBLOCKED 2026-06-25: check-stopping-power is done — the fourier-analysis skill now exists (.claude/skills/fourier-analysis/, locked FFT pipeline + peak-attribution rule + loss-locator caveat). The FFT/peak-extraction method is trusted, so the gated tasks 5-6 (omega_peak FFT + combined figure) are runnable. Tasks 4-6 still pending; ~150 GPU-h of production runs not yet launched."
---

# Li multi-k kick: omega_peak(v) drift verification

<gate>
**BLOCKED — do not launch any runs.** This campaign cannot start until the
**`fourier-analysis` skill** exists (task 5 of campaign `check-stopping-power`,
`docs/campaigns/check_logic/check_stopping_power_calculation.md`). The whole
campaign is gated, not just the analysis: the user's decision (2026-06-22) is that
the ~150 GPU-h of production runs must not be spent before the omega_peak(v) FFT /
peak-attribution method is locked and trusted, because the loss-function /
spectral methodology was gotten wrong once before. When `check-stopping-power`
flips to `status: done`, a session may set this campaign to `ready`/`running` and
proceed from task 0.
</gate>

<identity>
You are a scientific computing researcher working on first-principles
TDDFT simulations. You understand the solid-state impulsive-kick paradigm
(Santervas-Arranz, Stengel, Artacho, Phys. Rev. Research 7, 033292 (2025)), write
scientific-standard code, and adhere to this repository's rules and workflows.
Note: the QuantumKickExtension tree is its OWN git repository
(`/local/data/public/skcb2/tddft/QuantumKickExtension/`); run git from inside it.
</identity>

<description>
**The observation.** In the 17 single-k (Gamma-only, 1000 K) Li kick runs, the
peak of the electronic response (FFT of the per-unit-cell excess energy ΔE/N_uc,
and of dipole_x) was seen to **drift non-monotonically** as the impulsive kick
velocity v increases: low-v peaks sit near the plasmon (~6.5 eV) and high-v peaks
fall to ~2.6 eV (e-h / crossover regime), but the path between them is not a clean
monotone. It was not understood whether this drift is physics or an artefact of
the crude single-k sampling + hot 1000 K smearing.

**The plan.** A multi-k configuration much closer to 0 K (2x2x2 shifted MP grid,
400 K Fermi smearing) was started; 4 velocities are already done
(v = 0.0123, 0.0626, 0.300, 0.450). This campaign **replicates that multi-k config
at 5 more single-k velocities** so that omega_peak(v) can be recorded across all
three regimes and the single-k drift independently verified.

**The decision it informs.** Whether the non-monotonic omega_peak(v) is a real,
config-robust feature of Li's kick response (worth interpreting physically and
writing up) or a sampling artefact (to be discarded).

**Success / failure (falsifiable).**
- *Confirmed* if the multi-k omega_peak(v) reproduces the single-k non-monotonic
  shape (rise in low-v, turnover through medium-v, ~2.6 eV plateau in high-v)
  within FFT resolution.
- *Falsified* if multi-k omega_peak(v) is flat/monotonic, or differs
  systematically from single-k beyond resolution — i.e. the drift was an artefact.
</description>

<observables_set>
**Reuse the existing locked observable bundle in
`QuantumKickExtension/inq-codebase/Li/shared/cpp/run_template.hpp` verbatim. NO
new observables, NO new kernels → no code-test pre-gate required.**
- `observables.csv` every step: total/component energies, current, dipole,
  density_l2, centre_of_density (`.observables_current().observables_dipole()`).
- StateEnergyWriter every 10 steps (per-(k,state) ε(t); ~1.6 eV Nyquist).
- OccupationsWriter every 10 steps.
- density VTI every 100 steps (155 frames over the run).
Primary signals for the hypothesis: ΔE/N_uc(t) (from energy columns) and
dipole_x(t) → their FFT peak = omega_peak.
</observables_set>

<resolved_decisions>
All grounded in `inq-codebase/Li/shared/configs/base_li_54.hpp` +
`run_propagate_v0p0626_xyz/results/run_summary.txt` (2026-06-22 read):

- **Geometry** = 54-atom (3x3x3) BCC Li supercell, A_super = 10.53 Å,
  `../shared/li_54_3x3x3.xyz` (.xyz centred convention), 162 electrons.
- **Electronic structure** = ONCV-PBE, ecut 74 Ry, 20 extra states, **2x2x2
  shifted MP k-grid**, **Fermi smearing 400 K**.
  ⚠️ *Caveat:* the checkpoint is named `...T200...` but smearing was bumped from
  the 200 K target to **400 K** after a Broyden divergence (see base_li_54.hpp
  header). Use 400 K; the name is historical.
- **GS source** = REUSE checkpoint
  `inq-codebase/Li/checkpoints/li_54_2x2x2_T200_xyz` (E_GS ≈ −389.314 Ha,
  73 SCF iters). No new GS run. Task 0 = load-check (n_electrons → 162).
- **Propagator** = default INQ propagator (ETRS) via `real_time::propagate`,
  `options::theory{}.pbe()` (adiabatic-PBE in propagation), `.impulsive()` ion
  kick. Reuse `run_template.hpp` unchanged.
- **Duration + I/O** = dt = 0.04 a.u., N_steps = 15500 → T = 14.997 fs (~620 a.u.);
  FFT resolution Δω = 2π/T ≈ 0.28 eV (8x zero-pad interpolates; adequate for a
  6.5→2.6 eV drift). VTI/100, state-energy & occupations/10, observables.csv/step.
- **Kick** = impulsive, velocity gauge, +x direction (`KICK_DIRECTION_X = 1.0`),
  per-velocity override of `KICK_VELOCITY_AU` only.
- **Velocity matrix** (locked 2026-06-22) — 5 NEW runs at EXACT single-k values,
  3 points per regime once combined with the 4 existing multi-k runs:
  - Low: 0.0123✓, **0.0375**, 0.0626✓
  - Medium (crossover, previously unsampled in multi-k): **0.100**, **0.175**, **0.250**
  - High: 0.300✓, **0.375**, 0.450✓
  New Cfg structs to add: `Li_54_v0p0375`, `Li_54_v0p100`, `Li_54_v0p175`,
  `Li_54_v0p250`, `Li_54_v0p375` (subclass `Base_Li_54`, override
  `KICK_VELOCITY_AU` + `RUN_NAME = run_propagate_v<vel>` only).
- **Analysis scope** = omega_peak(v) curve + standard diagnostics: per-run
  plateau-detrend FFT of ΔE/N_uc and dipole_x → omega_peak; reuse the existing
  `gamma_transitions` and `state_energy_spectra` per-run attribution; combine all
  9 multi-k into THE omega_peak(v) figure with the single-k sweep overlaid. The
  Galilean/lab-frame density-spectra deep dive
  (`docs/plans/li_v0p0626_plasmon_density_analysis.md`) is OUT of scope here.
- **File placement** = runs in
  `QuantumKickExtension/inq-codebase/Li/run_propagate_v<vel>/` (existing pattern);
  combined study notebook + omega_peak(v) figure in that tree (e.g.
  `inq-codebase/Li/omega_peak_sweep/`). QuantumKickExtension is a separate git
  repo — commit there, not in main.
</resolved_decisions>

<guard_rails>
- **Hard gate (see `<gate>`):** launch nothing until the `fourier-analysis` skill
  exists. If the dependency is not done, STOP and surface it.
- **Pilot-first:** before the 5 production runs, a ~100-step smoke test at
  **v = 0.175** must pass: binary builds, GS checkpoint loads, n_electrons → 162,
  energy stays real, no NaN/complex. Only then launch the 5.
- **Abort conditions:** NaN or complex total energy; GS load failure; GPU occupied
  by another user (warn) — use the `cudaMemGetInfo` probe (NVML is broken; GPU is
  the default; do not CPU-fall-back on an nvidia-smi error).
- **Boundary/cadence note:** the 4σ/1σ launch-stop and jellium-WP VTI rules do NOT
  apply — this is a periodic solid Li crystal under an impulsive ion kick, not a
  jellium wavepacket. VTI/100 (155 frames) is the established cadence here.
- **Resolution caveat (state in the notebook):** Δω ≈ 0.28 eV; omega_peak
  differences smaller than this are not resolved — temper any "non-monotonic"
  claim against the bin size.
- **Reuse-or-surface:** the single-k overlay reuses the existing 17 single-k
  tddft.dat outputs; if a single-k peak cannot be re-extracted consistently,
  surface it rather than re-running single-k.
</guard_rails>

<tasks>
The executing agent flips the matching frontmatter `done` flag and updates the
handover as each completes. Tasks 0–3 are runnable once unblocked; 4–5 are
analysis, additionally requiring the `fourier-analysis` skill.

0. **GS load check.** Confirm `checkpoints/li_54_2x2x2_T200_xyz` loads and
   integrates to 162 electrons. *Done = load succeeds, n_electrons → 162.*
   Skill: `tddft-simulations`.
1. **Add 5 Cfg structs + build-once.** Append the 5 per-velocity structs to
   `base_li_54.hpp`; create the 5 `run.cpp` wrappers (mirroring
   `run_propagate_v0p0626_xyz/run.cpp`); build the binary once via `inq-run`.
   *Done = all 5 compile against the shared template.* Skill: `build-run`.
2. **Pilot smoke test @ v=0.175.** ~100 steps; verify the `<guard_rails>` pilot
   criteria. *Done = pilot passes all numeric checks.* Skill:
   `simulation-validation`.
3. **5 production runs.** Dispatch 0.0375, 0.100, 0.175, 0.250, 0.375 (15500 steps
   each) on GPU with per-run Gmail notification; ~30 h wall each (~150 GPU-h,
   ~3 days on 2 GPUs). *Done = 5 run_summary.txt show status complete,
   n_electrons → 162, energy real throughout.* Skill: `tddft-simulations`.
4. **[GATED] Per-run omega_peak + attribution.** Using the locked `fourier-analysis`
   skill: plateau-detrend FFT of ΔE/N_uc and dipole_x → omega_peak for each of the
   5; reuse `gamma_transitions` + `state_energy_spectra` for attribution.
   *Done = each run has a recorded omega_peak (energy + dipole_x) with attribution.*
   Skills: `fourier-analysis` (gating dep), `notebook-making`.
5. **[GATED] Combined omega_peak(v) + verdict.** Build THE omega_peak(v) figure
   over all 9 multi-k velocities with the single-k sweep overlaid; write the study
   notebook (per-run rows → combined curve → confirmed/falsified verdict against
   the hypothesis, tempered by the Δω ≈ 0.28 eV resolution). Update handover; flip
   frontmatter `done`/`status`. *Done = notebook executed, figure produced,
   verdict stated.* Skills: `notebook-making`, `fourier-analysis`.
</tasks>

<rules>
- **NEVER** launch a run while `status: blocked` (the Fourier gate).
- **ALWAYS** use the EXACT single-k velocity values (1:1 overlay depends on it).
- **NEVER** add or alter observables/config beyond the 5 velocity overrides — the
  whole point is config parity with the 4 existing multi-k runs.
- **ALWAYS** run git inside `QuantumKickExtension/` (separate repo).
- **NEVER** fftshift a VTI loaded for any density plot (physical-order rule).
- State the Δω ≈ 0.28 eV resolution limit in bold in the notebook before any
  non-monotonicity claim.
</rules>

<preflight>
Re-verify from this prompt alone BEFORE burning GPU (after the gate lifts):
- [ ] Gate lifted: `check-stopping-power` is `done` and the `fourier-analysis`
      skill exists. (If not — STOP.)
- [ ] Intent self-contained: hypothesis + confirmed/falsified criteria; every task
      has a done-criterion.
- [ ] Setup reproducible: 54-atom BCC Li / 162 e- / A=10.53 Å; GS = reused
      checkpoint li_54_2x2x2_T200_xyz; ETRS + adiabatic-PBE + dt 0.04 a.u. x 15500
      = 14.997 fs; 2x2x2 shifted, 400 K; observables bundle reused; file placement
      in QuantumKickExtension/inq-codebase/Li/.
- [ ] No new code: 5 Cfg velocity overrides only; no new observable/kernel.
- [ ] Validation & guard rails: pilot-first @ v=0.175; abort on NaN/complex/GS-fail;
      resolution caveat noted.
- [ ] Autonomous mechanics: GPU via cudaMemGetInfo probe (warn if occupied);
      per-run Gmail; study notebook auto-built; handover present; agent flips
      done/status.
- [ ] Grounding: paper-cited (Santervas-Arranz PRR 2025); config values carry
      file refs (base_li_54.hpp / run_summary.txt).
</preflight>
