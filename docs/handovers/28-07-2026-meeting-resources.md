# Handover — 28-07-2026 meeting resources (formulae, plots, tables)

Task: produce meeting-slide resources (theory write-ups, plots, workflow
diagrams, tables, run confirmations) for the 28-07-2026 meeting on isolating
the wavepacket stopping power. The user prepares the slides themselves; this
task only produces the resources. All resources live under
`/local/data/public/skcb2/tddft/docs/reports/28-07-2026-meeting/`.

## 2026-07-28 — theory document (DONE)

- Created `/local/data/public/skcb2/tddft/docs/reports/28-07-2026-meeting/formulae/wp_kinetic_energy_decomposition.md`:
  - Split A (momentum space): E_kin = ⟨p⟩²/2m + Var(p)/2m = T_drift + T_loc;
    Gaussian numbers T_loc(0) = 3/(8σ_r²) → 82 eV (σ_WP=0.5), 20 eV (σ_WP=1).
  - Split B (real space, Madelung): ψ = √n e^{iS} → T_W (von Weizsäcker,
    density-only) + T_v (flow, needs j); drift extraction T_drift = |∫j|²/(2mN).
  - Mapping between the splits (coincide at t=0; dispersion converts T_W→T_v),
    S_drift definition, Scheme 1/Scheme 2 correspondence, slide-ready shortlist.
- Content grounded in `/local/data/public/skcb2/tddft/docs/notes/refined-stopping-two-schemes.md`
  (2026-07-28) — all numbers taken from there, none invented. Plain-text
  equations by design; user will LaTeX them for slides.
- Verified: nothing to run; document is a write-up of already-validated theory.

## 2026-07-28 — equation PNGs (DONE)

- Generator: `/local/data/public/skcb2/tddft/docs/reports/28-07-2026-meeting/formulae/gen_equation_pngs.py`
  (matplotlib mathtext, Computer Modern, 300 dpi, white bg, grey one-line
  caption per equation; re-run with the venv python to regenerate).
- 12 PNGs in the same folder, verified written (nonzero size; not previewed,
  per the no-image-preview convention):
  - Split A: `eq_A1_ekin_momentum_split`, `eq_A2_variance_def`,
    `eq_A3_gaussian_variance`, `eq_A4_tloc`, `eq_A5_tloc_values`.
  - Split B: `eq_B1_polar_form` (caption flags polar/Madelung form of the WP
    orbital, per-orbital validity), `eq_B2_tw`, `eq_B3_tv` (captions note the
    polar-representation origin), `eq_B4_current`, `eq_B5_ekin_realspace_sum`.
  - Drift/stopping: `eq_C1_drift_extraction`, `eq_C2_sdrift`.

## 2026-07-28 — corridor-definition setup plots (DONE)

- Generator: `/local/data/public/skcb2/tddft/docs/reports/28-07-2026-meeting/formulae/gen_setup_plots.py`
  (venv python; inqview.load_vti physical-order loader, canonical theme,
  inferno). Outputs in the same folder:
  - `setup_t0_total_density.png` — qsp_phase5/p5_wp_v1p3 total density
    n(x,z,t=0), mid-y slice, linear+log panels; dashed slab faces ±12.5,
    dashed CAP inner faces ±35 + hatched CAP regions ±35..±45; vacuum
    corridors labelled with arrows; WP annotated at z=−23.75 (σ_WP=0.5,
    k0=1.3, 23 eV).
  - `setup_tlate_total_density.png` — same slice/scales at t*=28.8 a.u.
    (step 720, 0.70 fs): reflected lobe 0.11 e at z≈−27, transmitted lobe
    0.46 e at z≈+19, slab excess ≈0; caption defines projectile ≡ n−n_gs in
    the corridors → T_v, T_drift; coverage N(t*)=0.61 e (0.39 e CAP-absorbed).
- Shared colour scales across the pair (directly compared): linear vmax
  1.6e-3 (slab peak; WP core saturates, stated on-figure), fixed log
  [1e-7, 1.2] (per-frame scaling would hide dispersion/absorption).
- Verified from data (probe scripts, /tmp): density_total INCLUDES the WP in
  this run (N_total(0)=83.0, N_wp=1.0) — the campaign_autorun WP-exclusion
  trap does NOT apply here; frame choice from GS-subtracted corridor
  occupancy scan (step 720 optimal: both lobes in corridors, slab empty).
- PNGs not visually previewed (user previews, per convention).

## 2026-07-28 — quantum kick: omega_peak(v) figure + peaks table (DONE)

- New subfolder `/local/data/public/skcb2/tddft/docs/reports/28-07-2026-meeting/quantum_kick/`
  (user-requested "quantum kick" showcase from the quantum-kick-extension
  campaign, `docs/campaigns/quantum_kick_extension/quantum_kick_extension.md`).
- **Campaign data status** (checked on disk): 4 multi-k (2x2x2 shifted, 400 K)
  runs COMPLETE — v = 0.0123, 0.0626, 0.300, 0.450 (15500 steps, T = 620 a.u.,
  `QuantumKickExtension/inq-codebase/Li/run_propagate_v*/`); the 5 new runs
  (0.0375, 0.100, 0.175, 0.250, 0.375) were NEVER LAUNCHED (empty dirs — the
  "incomplete" part). Previous results = 17 single-k (Gamma-only, 1000 K) runs,
  T = 160 a.u., `QuantumKickExtension/systems/li/td_kick_v*/results/tddft.dat`.
- Generator: `.../quantum_kick/extract_omega_peaks.py` (venv python; locked
  fourier-analysis pipeline: mean baseline → Hann → 4× zero-pad →
  coherent-gain rfft, angular eV axis; peak searched in the physical band
  [1.8, 9.0] eV, never global argmax).
- Outputs (all in the subfolder):
  - `omega_peak_v.png` — headline presentation figure (title on canvas,
    `style.save_presentation`, 600 dpi): energy-channel ħω_peak vs kick
    velocity; 17 single-k points (grey, line) + 4 multi-k points (red
    squares); error bars = ±Δω/2 bin half-width (0.28 eV multi-k / 1.07 eV
    single-k).
  - `omega_peaks.csv` — 25 rows: dataset, v, channel, ω_peak (mean baseline),
    detrend comparison, baseline shift, T, Δω, provenance path.
  - `diagnostics/` — mandatory 3×2 FFT pipeline panel per reported spectrum
    (25 panels, 200 dpi).
- **Verified:** known-answer anchor reproduced exactly — multi-k v=0.0626
  dipole_x in [5.5, 8] eV → 6.480 eV (journal 2026-05-06 locked value);
  energy-channel peaks baseline-invariant (mean vs detrend) at all 21 points.
- **Headline numbers (2 s.f.):** single-k drift 5.6 → 6.7 eV (max near
  v ≈ 0.04) → falls to ~2.1 eV (v ≈ 0.3) → rises again to 3.2 eV (v = 0.45) —
  the non-monotonic shape. Multi-k: 5.7 (0.0123), 3.4 (0.0626), 2.6 (0.30),
  2.6 (0.45) eV. Multi-k confirms the ~2.6 eV high-v plateau; at v = 0.0626
  the multi-k energy peak (3.4 eV) sits well BELOW single-k (~6.1 eV) — the
  3.4 eV feature is cross-channel real (dipole_x dominant peak 3.24 eV, plasmon
  6.48 eV secondary), a genuine multi-k vs single-k difference to discuss.
- **Data-quality flags:** dipole_x EXCLUDED from the headline figure — at
  v = 0.0123/0.3/0.45 its band peak is a drift skirt (monotone low-ω comb;
  baseline-dependent at high v). Kept in CSV + diagnostics only. Single-k
  Δω ≈ 1.07 eV is coarse — sub-eV peak differences unresolved (stated as
  error bars, informational not a gate).
- Caption (slide spec, user's wording to refine): "Peak frequency of the
  electronic energy response of BCC Li after an impulsive ion kick, vs kick
  velocity. Grey: previous single-k (Gamma-only, 1000 K) sweep, 17 runs.
  Red: multi-k (2x2x2 shifted, 400 K) replication, 4 completed runs. Error
  bars are FFT bin half-widths (record lengths 160 vs 620 a.u.)."
- PNGs not previewed (user previews, per convention).

## 2026-07-28 — quantum kick: per-point hw_peak walkthrough notebooks (DONE)

- User request: show, in Jupyter notebooks, how the ħω_peak of each plotted
  multi-k point was calculated, processing step by step.
- Builder: `/local/data/public/skcb2/tddft/docs/reports/28-07-2026-meeting/quantum_kick/build_peak_walkthrough.py`
  (nbformat + nbclient; venv python with PYTHONPATH=inq-stack/python; rerun to
  rebuild; harvest-before-rebuild keeps reader-added markdown cells).
- Four EXECUTED notebooks in the same folder (0 errors, outputs embedded):
  `peak_walkthrough_v0p0123.ipynb`, `peak_walkthrough_v0p0626.ipynb`,
  `peak_walkthrough_v0p3.ipynb`, `peak_walkthrough_v0p45.ipynb`.
- Structure per notebook (house narrative): title/question → conventions +
  symbol table → setup printed VERBATIM from the run's `run_summary.txt` →
  source-file table → 10 steps, each formula immediately above its code cell:
  (1) raw E_total(t), (2) ΔE/N_uc response (N_uc=27), (3) mean baseline,
  (4) Hann window, (5) ×4 zero-pad, (6) coherent-gain rfft (/Σw, interior ×2),
  (7) angular axis ħω=2πf·27.211 eV + Δω=2π/T, (8) peak in band [1.8, 9] eV
  with linear+log spectrum, (9) robustness (detrend + Hamming/Blackman all
  equal reference), (10) production 3×2 fft_pipeline_panel + assert equality
  with `omega_peaks.csv` — all four print MATCH.
- **Verified in outputs:** 5.722 / 3.378 / 2.620 / 2.620 eV for v = 0.0123 /
  0.0626 / 0.30 / 0.45; every peak invariant across baselines and windows.
- Deliberate deviation, surfaced to user: NO density-evolution GIF in these
  notebooks — they are method walkthroughs of the peak extraction for the
  meeting, not per-run study notebooks; the runs' 155 VTI frames stay in
  their run dirs (user may request GIF-bearing run notebooks separately).

## 2026-07-28 — quantum kick: case-study 4-panel + spectra 2×2 + S(v) (DONE)

- Builder: `/local/data/public/skcb2/tddft/docs/reports/28-07-2026-meeting/quantum_kick/build_case_study_figures.py`
  (venv python, PYTHONPATH=inq-stack/python; presentation mode,
  `style.save_presentation`). Outputs in the quantum_kick folder:
  - `case_study_v0p45_fourier_steps.png` — 2×2 four-panel Fourier walkthrough
    of multi-k v=0.45: (1) ΔE/N_uc(t), (2) mean-removed + Hann (envelope
    dashed), (3) zero-pad ×4, (4) spectrum with band [1.8, 9] eV shaded and
    ħω_peak = 2.62 eV marked.
  - `spectra_2x2.png` — energy-response spectra of all 4 multi-k runs,
    normalised A/A_peak (shared axes 0–12 eV), band + peak marked per panel
    with the peak value in the panel title.
  - `S_v.png` + `stopping_sv.csv` + `stopping_fits/` (21 per-run fit
    diagnostics, S annotated on each) — per-ion electronic stopping power
    S(v), **Method A of the stopping-power-extraction skill** (Correa 2018
    Eq. 10): fixed 20%-of-TIME transient cut, free-intercept slope fit of
    ΔE_total vs common ion displacement x = v·t, normalised per ion
    (multi-k ÷54; single-k ÷2 — run.cpp confirms 2-atom BCC cell, 3×3×3
    Γ-centred ≡ 54-atom Γ-only).
- **Headline S numbers (multi-k, per ion):** 0.00018 (v=0.0123, ~zero within
  noise), 0.0023 (0.0626), 0.0033 (0.30), 0.0069 (0.45) eV/Bohr —
  monotonically increasing, small regression errors.
- **Caveats (surfaced to user):**
  - Single-k S(v) points are NOISE-DOMINATED (r² ≈ 0.01, e.g. v=0.30 gives
    −0.0071 ± 0.0021 eV/Bohr, unphysical sign): the 160 a.u. records are too
    short for a slope against the oscillation; error bars shown are
    regression stderr and understate the true (autocorrelated) uncertainty.
    Multi-k points are the trustworthy S(v) datapoints.
  - KE sanity channel (−dKE_ion/dx) UNAVAILABLE — neither output records ion
    kinetic energy (ions move at constant v; x_atom0 linear to machine
    precision). Stated per stopping skill §5. No absorber → N-drainage guard
    not applicable (unitary propagation).
- ΔE_total is the electronic total energy (excludes ionic KE) — verified to
  rise continuously (multi-k v=0.45: +31 Ha over the run), so the Method A
  slope channel is valid.

## 2026-07-28 — quantum kick: ALL-runs spectral map (DONE)

- User asked why not all multi-k simulations were used. Re-verified on disk
  (both tree mirrors, all QKE git branches, qball-codebase): exactly 4
  completed multi-k INQ runs exist; the 5 campaign runs are empty dirs (never
  launched); QBall td_kicks holds input templates only, no outputs. Nothing
  completed was omitted. HOWEVER the 17 "single-k" runs are literally
  multi-k too — `systems/li/td_kick_v*/run.cpp` uses a 2-atom BCC cell with
  a 3×3×3 Gamma-centred k-grid (the "single-k" campaign label means
  54-atom-supercell-EQUIVALENT Gamma sampling). Counting those, "everything"
  = 21 runs.
- New figure: `/local/data/public/skcb2/tddft/docs/reports/28-07-2026-meeting/quantum_kick/spectra_all_runs.png`
  (builder `build_all_runs_spectrum_map.py`, venv + PYTHONPATH): spectral map
  of ALL 21 completed runs — one row per run (energy-response spectrum,
  locked pipeline, normalised to its band peak), two stacked panels
  (17-row previous sweep / 4-row new config) with SHARED colour scale
  (inferno, 0–1), one colorbar outside at panel height, white open circles
  marking each row's band peak, x = ħω 0–12 eV, y ticks = kick velocities.
- Verified: all 21 per-row peaks reproduce omega_peaks.csv values exactly.

## 2026-07-28 — setup-plot label edits + strict-vacuum S recompute (DONE)

- `gen_setup_plots.py` edited per user: removed the "vacuum corridor
  (entrance)/(exit)" arrow labels and the "WP core saturates linear scale"
  note; WP-characteristics annotation kept. Both setup PNGs re-rendered.
- **Strict-vacuum criterion measured** (from `density_gs_system.vti`,
  planar-mean == planar-max verdict): n_gs < 0.1% of bulk (1.5e-6 vs bulk
  1.5e-3 e/Bohr³) for **|z| >= 21.0 Bohr** → strict vacuum = 21 <= |z| <= 35.
  Original exit plane (+22) already compliant; original entrance plane
  (−15.5) sits at ~5–7% of bulk.
- Recompute: `/local/data/public/skcb2/tddft/docs/reports/28-07-2026-meeting/formulae/recompute_strict_vacuum_S.py`
  — pipeline functions replicated verbatim from
  `hypotheses/qsp_phase5/build_momentum_stopping_report.py`, parametrised
  detector planes, null test re-run per plane set. Outputs beside it:
  `strict_vacuum_S_comparison.png`, `strict_vacuum_stopping.json`.
- **Baseline reproduction check PASSED** (planes −15.5/+22): S_drift matches
  published momentum_stopping_summary.json to 2 d.p. for all five runs.
- **Strict planes (±21) results:** null systematic grows ±0.41 → ±1.43
  eV/Bohr (entrance plane only 2.75 Bohr from launch z=−23.75 → almost no
  TOF flight distance → velocity resolution collapses; the null test
  captures it). v1p3: S = −0.19 ± 1.50 (was 0.30 ± 0.42); v3p0: 0.55 ± 1.88
  (was 1.34 ± 0.56). Strict and baseline AGREE within combined errors →
  baseline S is robust to the vacuum definition; the strict variant is a
  consistency check, not a better number (resolution-limited by geometry).

## 2026-07-28 — energy_oscillation subfolder: equations, repro figures, norm-corrected stopping (DONE)

All in `/local/data/public/skcb2/tddft/docs/reports/28-07-2026-meeting/energy_oscillation/`
(user-designated folder for this session's resources; underscore name).

- **Equation PNGs** (`build_equations.py`, mathtext, 600 dpi, transparent):
  `eq1_line83_kinetic_accumulation.png` (kinetic_ += occ_sum, defs line),
  `eq2_line55_occ_sum.png` (occ_sum = Σ occ·Re(arr)/Re(nor)),
  `eq3_inq_kinetic_energy.png` (E_kin = Σ occ⟨ψ|T|ψ⟩/⟨ψ|ψ⟩). Verified against
  `inq/src/hamiltonian/energy.hpp:55,83` before rendering.
- **Vacuum minimal-repro figures** (`build_vacuum_repro.py`): `setup_vacuum_cap.png`
  (t=0 total density of vacuum wp_traversal_energy/results/cap, dashed CAP band
  7.5–22.5, "CAP" only annotation), `vacuum_energy_vs_time.png` (CAP vs no-CAP,
  402.1 eV flat vs 360–420 swing), `vacuum_energy_channels.png` (all channels;
  only kinetic nonzero — non-interacting), `jellium_pairwise_channels.png`
  (ΔE_PP/PS/PB/SS/SB from wp_cap_energy_plateau cap interactions.csv — vacuum has
  no pairwise channels; interactions.csv confirmed already in eV),
  `vacuum_norm_proof.png` (reported == ⟨T⟩/⟨ψ|ψ⟩, ≠ corrected ⟨T⟩·norm; user
  renamed labels: "INQ reported" / "⟨T⟩ (corrected)").
- **Norm-corrected stopping** (`build_norm_corrected_stopping.py` +
  `norm_corrected_stopping.md`): E_corr = E_rep − T_pp·(1−N);
  E_dep = (E_corr − T_pp·N) − baseline. Outputs `replica_corrected_energy.png`,
  `plateau_corrected_energy.png`, `deposited_energy.png`.
  - In-slab traversal slope agrees across both runs: **S ≈ 0.24 eV/Bohr**
    (replica 0.244, plateau 0.234; norm ≈ 1.0000 in window → correction-free).
  - **CONVENTION FIX (later same day):** cap_wp_kinetic.csv:T_wp_ha is
    EXTENSIVE (correction −T·(1/N−1), verified vs its correction_eV column);
    wp_momentum_stats:e_kin_ha is PER-PARTICLE. First pass under-corrected the
    plateau absorption phase by ~22 eV. Corrected: plateau E_dep 4.4 eV →
    **S = 0.18 eV/Bohr** (NOT 1.1 as first computed); corrected plateau tail
    −115.8 eV matches the independent 115.9 eV corrected-gap result.
    All three deposition estimates now agree ~0.2 eV/Bohr (replica partial
    0.17 at t=57, N=0.085) — factor ~3 below σ_WP=1 linear response (0.57–0.70),
    ~6 below point-charge LR (1.2), r_s=3.32 (interior r_s verified 3.31 from
    GS z-profile). Density-scaled classical benchmark (r_s 4.18→3.32, LR ratio
    1.7–1.9): 1.5–1.7 eV/Bohr. σ-arg convention of stopping_power_sigma is
    ambiguous across call sites — pin down before quoting the finite-σ LR.
  - Rerun build_norm_corrected_stopping.py when the replica finishes (reads
    live CSVs). replica CAP run in flight (~step 2860/8000); no-CAP not started.
  - Also produced: `localised_highdens_energy_vs_time.png` (raw reported
    ΔE_total, no-CAP flat 0.00 eV vs CAP plateau −86 eV ± 5, 100 a.u.) and
    `localised_highdens_naive_correction.png` (naive −T0·(1−N) from streamed
    norm vs full correction: naive tail −206 eV OVER-corrects by ~90 eV vs full
    −116 eV — absorbed weight carried ~22 eV, not its initial 120.4 eV;
    builder `build_localised_highdens_energy.py`). High-density BULK jellium
    E(t) also built (`build_highdens_energy.py`, n162_L30 runs: conserved to
    0.44/0.13 meV at dt 0.02/0.01).

## 2026-07-28 — per-definition method stories (equations + plots) (DONE)

- User confirmed scope: ALL THREE v1p3 method notebooks (Methods #2 snapshot
  kinematics, #3 Ehrenfest drag, #4 momentum-space), each as a self-contained
  story subfolder under
  `/local/data/public/skcb2/tddft/docs/reports/28-07-2026-meeting/formulae/`:
  `snapshot_kinematics/`, `ehrenfest_drag/`, `momentum_space/`.
- Each subfolder: `story.md` (setup → equations → quick validation → energy/
  result plots, with headline numbers), setup-plot copies, method equation
  PNGs (13 total, mathtext; formulae transcribed from the notebooks' own
  equations cells; headline numbers from the summary JSONs), the notebooks'
  executed figures extracted as PNGs (8 total, ~1000–1600 px as executed),
  and `methods_comparison_v1p3.png`.
- Generators (all in `formulae/`, venv python): `gen_method_eqs.py`
  (equations → subfolders), `extract_method_figs.py` (notebook figures →
  subfolders), `gen_methods_comparison.py` (synthesis figure from the three
  summary JSONs; copies to root + all subfolders),
  `gen_snapshot_story_figs.py` (NEW figures for #2, whose notebook had only
  tables: `snapshot_val_tw_dispersion_gate.png` — lobe T_W vs free-dispersion
  law + gate numbers; `snapshot_result_S.png` — T_drift/N and S_snap vs t*
  with coverage, TOF band, deposit line).
- Headline numbers (2 s.f., from summary JSONs): #2 S_snap 0.16–0.53 (t*
  scan, brackets TOF 0.30 ± 0.42); #3 S(v0) = 1.2 at v̄ = 1.4, interior
  3.0 ± 0.55 at v̄ = 0.47; #4 S_ms = 0.38 ± 0.27 at u ≈ 1.7; all vs
  S_deposit = 2.4.
- mathtext gotchas fixed in generators: no \\Big/\\big/\\tfrac/\\le
  (use \\left..\\right, \\frac, \\leq).
- Extracted notebook figures are as-executed (~150 dpi) — flagged; can be
  re-rendered at 300 dpi on request by re-running the notebook code.
- PNGs not visually previewed (user previews, per convention).

## Not done / next (per user, upcoming in-conversation)

- Further plots, workflow diagrams, tables for the same meeting folder.
- Quantum-kick campaign completion (post-meeting): launch the 5 pending
  multi-k runs (0.0375, 0.100, 0.175, 0.250, 0.375; ~150 GPU-h) per campaign
  tasks 1–3, then rerun `quantum_kick/extract_omega_peaks.py` — it picks up
  new runs only via the MULTIK dict (add the new run dirs there).
- Possible confirmations/verifications of completed simulations.
- Pre-existing sibling folder `.../28-07-2026-meeting/energy_oscillation/`
  (empty at task start) — not this task's output.

## Conventions in force

- Python via `/local/data/public/skcb2/tddft/venv/bin/python3` (venv, VTK).
- Figures as .png; 2 s.f. rounding in human-facing numbers; σ always means
  σ_WP; no LaTeX in chat.

## 2026-07-28 — LIVE S(v) from the running replica_lz160_1cap CAP run (DONE)

- Running sim identified: `scripts/replica_lz160_1cap/run_jellium_replica.sh`
  → `./run` (GPU 0, since 03:16): LZ=160 one-sided-CAP replica, σ_WP=1,
  E=100 eV (k0=2.711), launch z=−20.5, slab |z|≤12.5, CAP z∈[60,80],
  dt=0.02, 8000 steps planned; at extraction step 2837 (t=56.7 a.u.).
- Method: Ehrenfest drag (#3) on the partial per-step record
  (wp_momentum_stats.csv every step + wp_real_space_stats.csv every 10) —
  the only refined method valid mid-run.
- Extractor (re-runnable as the run progresses):
  `/local/data/public/skcb2/tddft/docs/reports/28-07-2026-meeting/formulae/replica_lz160_1cap/live_sv_drag.py`
  → `live_sv_drag.png` (3 panels: identity+norm diagnostics; T_drift vs z_c
  with fits; sliding-window local S along path) + `live_sv_summary.json`.
- **Headline (2 s.f.): S_face-to-face = 0.20 eV/Bohr at v0 = 2.7** (pz
  2.742 at z=−12.5 → 2.675 at +12.5; image accel/decel cancels by slab
  symmetry); S_interior (|z_c|≤8) = 0.29 eV/Bohr at v̄ = 2.7.
- **CAP-selection onset measured at t = 17.2 a.u.** (norm < 0.999; norm now
  0.089, centroid stalled ~z=50): everything after is selection, excluded.
  The clean drag segment is ALREADY fully recorded — running further does
  NOT improve this extraction (run's remaining value = the energy plateau);
  the follow-on no-CAP run will give a selection-free record.
- Image acceleration confirmed again: pz 2.711 → peak 2.76 mid-slab
  (mirrors the v1p3 discovery; σ_pz²(0) = 0.500 = min-uncertainty ✓).

## 2026-07-28 — corrections explainer + comparison table (DONE)

- User asked for a clear explanation of the naive correction, orbital-free
  correction, and KS-orbital-based correction.
- `formulae/corrections_explained.md` — the explainer: one principle
  (S = −d(⟨p⟩²/2m)/ds; corrections differ only in WHO "the projectile" is);
  naive deposit S = [E(0)−E_plateau]/L over-counts (localisation 82 eV @
  σ=0.5 + capture binding; 59 eV "absorbed" vs 23 eV drift at k0=1.3);
  Scheme 1 = corridor excess density (methods #1, #2); Scheme 2 = injected
  KS orbital (methods #3, #4); cross-family agreement (0.30±0.42 vs
  0.38±0.27) validates both weak points simultaneously.
- `formulae/corrections_table.png` (+ gen_corrections_table.py) — 3-column
  slide table (projectile definition, observable, localisation counting,
  methods, assumptions, v1p3 results; naive column tinted red).
- FLAG for user: the qsp_phase5_momentum_stopping notebook PROSE quotes
  T_loc(0) = 1.5 Ha = 41 eV, but the measured gate (ekin 3.85 − drift 0.845
  = 3.0 Ha) and the analytic 3/(8σ_r²) with σ_r = 0.354 give 3.0 Ha = 82 eV
  — the notebook text looks like a σ_r-vs-σ_WP factor-2 slip; 82 eV used in
  all this conversation's resources.
- Bare equation PNGs (no captions, per user): `gen_correction_eqs.py` →
  `corr_eq0_principle.png` (S = −d(⟨p⟩²/2m)/ds), `corr_eq1_naive_deposit.png`,
  `corr_eq2_orbital_free.png` (T_drift = P_z²/2N with corridor N, P_z),
  `corr_eq3_ks_orbital.png` (⟨p⟩ from ψ_wp; E_kin,orb split).
