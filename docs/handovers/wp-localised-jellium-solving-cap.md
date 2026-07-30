# Handover: WP + localised-jellium CAP campaign (energy-plateau diagnostic)

Rolling handover. Task: build a campaign to test whether the plateauing
`energy_total` in high-density runs retains "too much" energy, by comparing a
localised-jellium WP run WITHOUT vs WITH a CAP. **Before the jellium runs**, the
user requested a quick vacuum single-WP warm-up experiment (in progress).

Session: `wp-localised-jellium-solving-cap`. Started 2026-07-22.

## ✅ COMPLETE (2026-07-23) — all four runs done, GS notebook added

The autonomous chain ran to completion unattended. CAP run finished 2026-07-23
19:55 (wall ≈ 13.6 h); comparison + final email done. All four RT runs report
`run_completed = true`.

**Headline:** no-CAP jellium plateau = −22462 eV, CAP = −22548 eV → **no-CAP sits
86 eV ABOVE CAP** = the trapped, should-have-radiated energy. Hypothesis holds.
(Vacuum warm-up gap was 36.6 eV; larger for jellium because the slab scatters more.)

**Notebooks (all exist, ordered for examination — see campaign doc "Results" table):**
1. GS: `.../scripts/wp_cap_energy_plateau/gs/report/gs_report.ipynb` **(added 2026-07-23**
   — was the one missing artefact; z-profile, 1D reconstructed potential, xz slice,
   WP-position validation, energy ledger. Builder: `gs/make_gs_report.py`.)
2–3. vacuum no-CAP / CAP: `systems/vacuum/scripts/wp_traversal_energy/results/{nocap,cap}/report/run_report.ipynb`
4–5. jellium no-CAP / CAP: `.../wp_cap_energy_plateau/wp/results/{nocap,cap}/report/run_report.ipynb`
6. comparison: `.../wp_cap_energy_plateau/wp/results/jellium_energy_compare.png`

**Logical-order guide with full clickable paths + the "why":** campaign doc
`docs/campaigns/localised_jellium/wp_cap_energy_plateau.md` → "Results — examine in
this order".

**Still open (not blocking):** physical decomposition of the 86 eV (reflected-WP KE
vs plasmon vs bound state — `scientific-panel`); E_ss/E_sp pairwise (twin decompose);
commit the machinery (two-commit split); catalogue the 4 runs.

## 2026-07-27 (RESOLUTION) — the CAP "energy artifact" is INQ reporting per-particle energy

**Root cause found in INQ source.** `inq/src/hamiltonian/energy.hpp:50-55` `occ_sum`
computes `sum occ*<psi|H|psi>/<psi|psi>` — it DIVIDES each orbital's energy by its
norm, i.e. reports the INTENSIVE (per-particle mean) energy. Normal calc: norm=1,
no-op. Under a CAP the WP norm decays -> reported `energy_total` stays ~E0 (even
climbs) as the packet is absorbed. See [[reference_inq_reports_normalized_energy]].

**Isolation experiments (vacuum, E=400 eV, all this session):**
- **Wrap is innocent:** longer no-CAP run (t=16, ~2 wraps) conserves energy to
  0.15 meV -> periodic re-entry is a clean unitary op; disproves "re-entry = new
  projectile". Plot `wrap_conservation.png` in the nocap notebook.
- **Not reflection:** cap_better (W=25, eta=-0.7), cap_fulllen (W=30 full +z half,
  eta=-1.0), and the earlier eta=-3.5 ALL give the same ~400 eV "residual" -> it is
  independent of CAP strength/width, so NOT leading-edge reflection.
- **Fix verified:** extensive energy `E_ext = E_reported * norm` decays cleanly with
  the norm. cap_fulllen: E_ext 402->1.65 eV, captured 400 eV = 99.6% (== norm
  absorbed); cap_better captured 386 eV (96%, 4% leaked). Energy IS captured
  correctly once un-normalized. There was never a physical artifact.

**Deliverables:** `energy_diagnostics.py` (decomposed-energy, compounded=extensive
vs reported energy, wrap-conservation; wrap-time marker t_wrap=(LZ/2-launch)/k0);
appended to nocap/cap/cap_better/cap_fulllen run notebooks. Runs: `results/nocap_long`
(t=16), `results/cap_better`, `results/cap_fulllen`. Dispatcher
`rerun_cap_experiments.sh`.

**Jellium correction — DONE (post-processing, no new run).**
`wp_kinetic_normalization_fix.py` computes <T_WP> and norm_WP from the saved complex
`wavefunction_wp` frames (501 each; validated t0 mean KE=120.4 eV=k0^2/2+3/4sigma^2),
and corrects E_total by `-<T_WP>*(1/norm-1)` (only the CAP run; nocap WP norm=1).
Result: **reported gap 93.5 eV -> corrected 115.9 eV** — the normalization ADDS ~22 eV
(peak 44 eV mid-absorption), it does NOT explain away the plateau. WHY only 22 eV: the
WP slows 120->22 eV in the bath BEFORE the CAP absorbs it (that ~98 eV kinetic is
deposited into the bath, correctly booked via the density-based terms). So the plateau
is LARGELY PHYSICAL; the WP-kinetic normalization is a real but sub-dominant systematic.
Plot `wp/results/comparison/corrected_plateau.png` appended to the jellium comparison
notebook; CSVs `{nocap,cap}_wp_kinetic.csv`.

## 2026-07-27 (latest) — low-spreading E=400 eV vacuum rerun (fixes self-interference)

- **User:** saw the WP self-interfere in the report GIFs; asked for higher energy /
  lower spreading. Diagnosis: transverse (x,y) is actually CLEAN (edge frac ~1e-9);
  the visible "interference" was the Z-WRAP (WP reaching the +z wall and re-entering
  at -z), amplified by the log-scale GIF panel.
- **New production run (dual-GPU, ~2-4 min each):** sigma0=3, **E=400 eV**
  (k0=5.421, k0*sigma0=16.3 -> **~5% transit spread**, down from 17%), grid **h=0.4**
  (cutoff guard PASS: k_max=7.85 > k0+4dk=6.36, E_cut=839 eV), dt=0.01. Box 30x30x45,
  one-sided +z CAP z in [7.5,22.5] (W=15). Launch z=-7.5 (5 sigma0 from CAP inner AND
  wrapped wall). Dispatcher `rerun_lowspread_dualgpu.sh` (build-once, nocap on GPU0 /
  cap on GPU1, thorough gates, notebooks, email).
- **CAP absorption tuning (validated):** survival = **exp(-|eta|W/v)** (the sin^2 CAP
  averages to 1/2, canceling the usual factor of 2 in exp(-2|eta|W/v)). Measured:
  eta=-1.0 -> 0.063 (6.3% leaked + wrapped), eta=-2.5 -> 0.0010, **eta=-3.5 -> 1e-4**
  (below the log-GIF floor -> invisible). One-sided CAP outer edge == box wall, so any
  leak wraps instantly -> must absorb hard. Reflection stayed **0.000** throughout
  (adiabatic W=15; neg-k/pos-k momentum weight check).
- **Final gates (all green):** N(t0)=1.0000; spread@transit=+5%; transverse edge
  ~1e-6 (no x/y wrap); CAP norm(tF)=1e-4 (99.99% absorbed); reflection 0.000; no-CAP
  control z_wrapped=False (350 steps stops it before the wall). The cap `z_wrapped`
  flag is a detector artifact on the 1e-4 ghost (negligible).
- **Notebooks (rebuilt):** `results/{nocap,cap}/report/run_report.ipynb` +
  `results/comparison/nocap_vs_cap_comparison.ipynb` + setup figure. run.cpp defaults
  now hold this design (k0=5.421, sigma=3, h=0.4, dt=0.01, LZ=45, CAP_L=15,
  launch=-7.5, eta=-3.5 via dispatcher). Supersedes sigma3/100eV (17%).

## 2026-07-27 (later) — compact non-dispersing vacuum rerun + WP-dispersion formula correction

- **User goal:** a COMPACT WP that does not expand appreciably over the sim.
  Physics: a free Gaussian disperses; the design must control it, not the box.
- **FORMULA CORRECTION (important):** the density width spreads as
  `sigma_dens(t)=sqrt(sigma0^2/2 + t^2/(2 sigma0^2))`, i.e. expansion factor
  `R=sqrt(1+(t/sigma0^2)^2)`, spreading time `tau=sigma0^2` — NOT `2 sigma0^2`
  (my first estimate had a 2x error). Verified against the vacuum runs to 3 dp
  (minimum-uncertainty, no bug). Now in [[reference_wp_dispersion_formula]] +
  `.claude/rules`-adjacent memory. Design rule: transit `R=sqrt(1+(5/(k0 sigma0))^2)`
  -> dispersion controlled ONLY by **k0*sigma0** (>=16 for ~5%, >=11 for ~10%).
- **Runs done (dual-GPU, GPU0=nocap GPU1=cap, ~5 min each):** sigma0=3, E=100 eV,
  box 30x30x40, one-sided +z CAP z in [10,20], launch z=-5 (5 sigma0 from CAP inner
  AND wrapped -z wall), 600 steps. Defaults baked into `run.cpp`. Dispatcher
  `rerun_compact_dualgpu.sh` (setsid, build-once then concurrent runs + notebooks +
  verify + email).
- **Verified:** N(t0)=1.0000 (true vacuum), sigma_wf(t0)=3.00 as designed, CAP
  norm(tF)=0.077 (92% absorbed), nocap norm conserved. The nocap sigma_z blow-up at
  t>7 is the CONTROL WP wrapping the periodic box (no absorber) — expected, not a bug.
- **Transit expansion ~17% (NOT the 5% I first promised)** because k0*sigma0=8.1, not
  16 — the formula error. Corrected options: sigma3/400eV (5%, same box), sigma6/100eV
  (5%, box 60x60x70), sigma4/100eV (10%, 40x40x50). **User chose to KEEP the current
  sigma3/100eV (~17%)** — no rerun. Notebooks stand.
- **Notebooks (all built + tabulated):** vacuum run+phase under
  `.../vacuum/scripts/wp_traversal_energy/results/{nocap,cap}/report/run_report.ipynb`
  and `.../results/comparison/nocap_vs_cap_comparison.ipynb`. Jellium notebooks (GS +
  2 runs + jellium comparison) unchanged from earlier this day, with the ΔE plots.
- **NEW builder:** `.../wp_cap_energy_plateau/compare_notebook.py` — the jellium
  cap-vs-nocap PHASE notebook (was missing; vacuum had one). Produces
  `.../wp/results/comparison/jellium_nocap_vs_cap_comparison.ipynb`.

## 2026-07-27 — true-vacuum + 30×30 transverse rerun; ΔE decomposition plots; autonomous pipeline

Two user-driven corrections + a new energy-decomposition deliverable. All done
autonomously on **GPU 1** (GPU 0 was occupied by another user's task; chosen via a
custom `/tmp/gpuprobe.cu` cudaMemGetInfo probe since `nvidia-smi`/NVML is broken —
cosmetic). Pipeline: `.../vacuum/scripts/wp_traversal_energy/autorun_pipeline.sh`
(setsid-detached, survives session; log `autorun_pipeline.log`; ran to COMPLETE
2026-07-27 01:49).

### Verified problems in the prior vacuum runs (user was right on both)
1. **Not vacuum — 2 background electrons.** `run.cpp` used
   `extra_states(1).extra_electrons(2.0)`: the WP was a 3rd electron on top of a
   uniform k=0 2-electron gas (measured total=3.0, WP=1.0, bg=2.0 uniform). In
   `non_interacting` this is kinetically inert but pollutes `density_total`.
2. **Transverse box 12 Bohr far too small.** A free Gaussian disperses to
   σ_dens≈11 Bohr by t=32 (>> half-box 6). Measured σ_x saturated ~3.5 (vs analytic
   11.3) and edge-probability hit 0.84 → the WP wrapped x/y and self-interfered
   (the "interference" the user saw).

### Fixes (run.cpp), user-chosen geometry (30×30, keep long z)
- **True vacuum:** `extra_states(0).extra_electrons(1.0)` — the WP REPLACES the one
  base electron (INQ needs ≥1 electron for GS AND validates the count in propagate;
  `extra_electrons(0)` throws "no electrons" from BOTH initial_guess/calculate AND
  propagate). GS relaxes the single electron to k=0, then
  `inject_into_last_extra_state` overwrites that (only) state → WP is the sole
  electron. **Verified N_total = 1.00000, norm_after=1.**
- **Transverse box:** `WP_LPERP` 12 → **30** (half-box 15). LZ=80, launch z=−30
  unchanged. **Verified:** transverse wrap edge-fraction **84% → 10.5%** (residual
  only in the final frames, after the WP is mostly CAP-absorbed — the compromise the
  user accepted vs a ~80-Bohr box). CAP run norm→0.229 (77% absorbed); nocap
  conserved.
- Old LZ=60/launch=−26 AND the 12-Bohr-transverse runs are both superseded.

### NEW: ΔE energy-decomposition plots (user request) — both jellium runs
- **`energy_decomposition.py`** (new, in `.../wp_cap_energy_plateau/`): per-run
  pairwise electrostatic decomposition reconstructed from the saved density VTIs via
  FFT-Poisson, mirroring `inqkit::jellium::interaction_energies.hpp`
  (`compute_coulomb_wp`). P=WP, S=slab e⁻ (n_total−n_wp), B=+background
  (`n0·½erfc((|z|−12.5)/0.5)`, n0=102/15625). Emits `interactions.csv`.
  **Closure-gated (validation):** reconstructed E_hartree matches INQ to
  **1.4e-08 eV** (Poisson convention exact) → E_ss/E_ps/E_pp trustworthy; E_external
  closes to ~14 eV (0.03%, analytic-n+/charged-cell G0), absorbed into the E_sb/E_pb
  split (their sum exact). E_pb absolute carries the charged-cell G=0 gauge (flagged
  on plot). E_ss+E_ps+E_pp≡E_hartree, E_sb+E_pb≡E_external enforced.
- **`analyse.py`** now also emits `energy_delta_components.png` (ΔE(t) per KS
  component) and `energy_delta_pairwise.png` (ΔE of E_ss/E_ps/E_pp/E_sb/E_pb) and
  embeds them in `run_report.ipynb`.
- Headline: nocap ΔE_total≈0 (conserved sanity check); **cap ΔE_total = −93 eV**
  (electronic energy the CAP removes).

### Caveat carried forward
- **Jellium transverse box is ALSO undersized** for the dispersing WP: 25 Bohr
  (half-box 12.5), 100 a.u. propagation → σ_dens→~35 Bohr, worse wrap than the old
  vacuum. NOT yet re-run. The jellium notebooks/ΔE plots are on the EXISTING runs; if
  the transverse wrap matters for the physics conclusion, the jellium runs need the
  same enlargement (expensive — 13.6 h each).
- Completion email must NOT attach the comparison notebook (~25 MB embedded GIFs >
  Gmail cap → bounces); pipeline patched to send text-only. Resent text-only OK.

## 2026-07-24 (later) — vacuum CAP is ONE-SIDED; 5σ-clearance rerun (80-Bohr box)

- **Correction (my earlier setup figure was wrong).** The vacuum CAP is
  **one-sided at the +z end only** — `perturbations::absorbing` (absorbing.hpp:44)
  makes a SINGLE band `mid−w/2 < z_frac < mid+w/2`; the vacuum run builds exactly
  one (run.cpp:100). Two-sided needs two summed bands (as the jellium runs do). My
  first `setup_vacuum_cap.png` drew crimson CAP bands on BOTH sides, which made the
  WP (launch z=−26) look like it started inside a −z CAP. It did NOT — norm=1.0 at
  t=0. The genuine issue: z=−26 was only 4σ from the −z periodic wall, and the wall
  ≡ +z CAP outer edge (periodicity), so the WP's wrapped tail grazed the CAP.
- **Rerun (user-directed, ≥5σ).** New corrected geometry, two identical runs
  (no-CAP + CAP), one-sided +z CAP:
  - box **12×12×80** (was 12×12×60); CAP one-sided **z∈[30,40]**; WP launch
    **z=−30** → **10σ** clear of both the −z wall/wrapped CAP edge and (60 Bohr
    from) the +z CAP inner edge; n_steps **1600** (t=32); h=0.5, σ=1, E=100 eV.
  - Verified: both norm=1.0000 at t=0 (no init absorption); no-CAP norm conserved
    1.0; CAP WP reaches +30.8 and is absorbed to **norm=0.23** (77% removed).
  - run.cpp defaults updated to this geometry (self-documenting, reproducible);
    old LZ=60/launch=−26 results removed (superseded).
- **Dispatch/regeneration scripts (skill-local, reusable):**
  `.../vacuum/scripts/wp_traversal_energy/rerun_5sigma.sh` (build inq-study + both
  runs) and `.../regen_notebooks.sh` (setup fig + both per-run notebooks +
  comparison).
- **One-sided CAP drawing fixed everywhere:** `make_density_gif_battery` gained a
  `cap_lines` override (density_gifs.py); `analyse.py --cap-lines 30,40`;
  `make_setup_figure.py` + `compare_notebook.py` draw a single +z band. No more
  spurious −z CAP line.

## 2026-07-24 — GIF scaling fix + vacuum comparison artifacts

- **Bug found & fixed (density GIF scaling).** Vacuum WP-CAP GIFs "showed no
  motion" — a colour-scale artefact: a σ=1 free WP disperses so its peak density
  collapses ~1/σ³ (~100×: 0.18→0.0017 a₀⁻³ by t≈6), and `_save_gif` locked ONE
  linear vmax for all frames → dispersed WP <2% of vmax → black. WP genuinely
  moves (peak-z −26→+20, verified from VTIs + `wp_real_space_stats.csv`).
  Fix (opt-in, backward-compatible): `make_density_gif_battery(...,
  per_frame_norm_wp=True)` + `analyse.py --per-frame-norm-wp` → WP linear panel
  per-frame normalised (n/nₘₐₓ(t)) + log panel widened to 4 decades. In
  `inq-stack/python/inqview/visualisation/density_gifs.py`. Memory:
  `reference_dispersing_wp_gif_scaling`.
- **Vacuum notebooks regenerated** with the fix (nocap + cap `run_report.ipynb`).
- **Setup figure (real density, per scientific-figures §4):**
  `.../vacuum/scripts/wp_traversal_energy/results/cap/report/setup_vacuum_cap.png`
  — t=0 total density xz slice, dashed CAP bands |z|∈[20,30] + WP launch z=−26.
  Builder: `.../wp_traversal_energy/make_setup_figure.py`.
- **no-CAP vs CAP comparison notebook (NEW, user-requested):**
  `.../vacuum/scripts/wp_traversal_energy/results/comparison/nocap_vs_cap_comparison.ipynb`
  — energy(t) overlay + side-by-side TOTAL density GIF with SHARED FIXED LOG scale
  (shows motion AND CAP absorption; a per-frame-norm twin is motion-only). Builder:
  `.../wp_traversal_energy/compare_notebook.py`.
- **OPEN (user message truncated at "Then, I …"):** the 3rd requested comparison
  item is unknown — ask the user before adding it.

## AUTONOMOUS RUN LAUNCHED (2026-07-22, detached, GPU 0)

Orchestrator: `ResearchProject/systems/localised_jellium/scripts/wp_cap_energy_plateau/orchestrate.sh`
launched `setsid`-detached (PPID=1, own session) → survives session end. Log:
`.../wp_cap_energy_plateau/orchestrate.log`. Emails to chiddukanna@gmail.com,
family `[wp-cap-energy-plateau]`, at each stage.

**Chain:** vacuum warm-up notebooks → jellium GS → self-validate WP pos →
WP smoke (40 steps) → WP no-CAP (5000 steps/100 a.u.) → WP CAP (η=−0.7) → compare.
Idempotent + resumable: re-launch the same command to resume (completed stages
skipped via `run_completed=true`; RT runs resume from `rt_ckpt`).

**Re-launch command (if killed):**
```
cd .../scripts/wp_cap_energy_plateau && setsid nohup bash orchestrate.sh >/dev/null 2>&1 </dev/null &
```
**Kill:** `kill <pgid of orchestrate.sh>` (find via `ps -eo pid,sess,cmd | grep orchestrate`).
GS checkpoint: `.../shared_gs/slab_n102_L25x25x140_w0p5_h0p5`.

### Live progress (2026-07-22, verified through every runtime gate)
- Vacuum warm-up DONE + emailed: no-CAP plateau 120.4 eV (conserved), CAP 83.8 eV,
  **gap 36.6 eV**. Notebooks in `systems/vacuum/scripts/wp_traversal_energy/results/*/report/`.
- Jellium **GS converged: E_GS = −830.0 Ha**, saved to shared_gs. (~3.5 min.)
- WP-position validation **PASS**: n(launch)=8e-7 (0.013% of centre), slab edge
  ≈−12.34, launch-to-face gap = 8.00 Bohr.
- WP **smoke** (40 steps) PASS: WP injected idx=74, norm=1.0, max_overlap=2e-4.
- **no-CAP full run (5000 steps) RUNNING** since 17:45. Then CAP, then compare.
  All subsequent stages emit stage emails; nothing else needs manual action.

### Verified before launch
- Engine settled empirically (stock inq can't compile CAP; inq-study absorbs).
- Vacuum WP no-CAP: E_total=4.42 Ha CONSERVED to 1e-5; CAP: drains 4.42→~0. Clean
  decomposition (total=kinetic, H=xc=0) after the `ground_state::calculate` bath-relax fix.
- GS run + WP run both COMPILE (compile-probe ok on inq-study).
- `analyse.py` verified on vacuum data: energy/momentum PNGs + 9 density GIFs +
  embedded-GIF `run_report.ipynb` (per notebook-density-gif rule).
- Email path (`notify.py`) verified working.

### Files (campaign machinery)
- Config: `shared/configs/slab_n102_L25x25x140_w0p5.hpp`
- GS run: `scripts/wp_cap_energy_plateau/gs/run.cpp`
- RT run (env CAP on/off): `scripts/wp_cap_energy_plateau/wp/run.cpp`
- Orchestrator: `orchestrate.sh`; validators/analysis: `validate_wp.py`,
  `analyse.py`, `compare.py`, `notify.py`
- Vacuum warm-up run: `systems/vacuum/scripts/wp_traversal_energy/run.cpp`
  (results/nocap, results/cap DONE + report/ notebooks built)

## Status

- **DONE — engine question settled empirically.** Stock `inq/` CANNOT compile a
  real-time CAP run (`double += complex` at `absorbing.hpp:45`; `vscalar=vion_` is
  REAL in `self_consistency.hpp:176`). Only `inq-study` compiles + absorbs
  (probe: WP norm 1.0→0.30, 70% absorbed). Probe:
  `ResearchProject/systems/localised_jellium/scripts/cap_engine_probe/run.cpp`.
  Both campaign runs build against inq-study. Memory:
  `reference_stock_inq_cannot_compile_cap`.
- **IN PROGRESS — vacuum warm-up experiment** (this is the current focus):
  single WP in vacuum, full traversal, WITH and WITHOUT CAP, non-interacting.
  Deliverables: run notebooks for both + total-energy-vs-time + all decomposed
  energies plotted. Location: `ResearchProject/systems/vacuum/scripts/wp_traversal_energy/`.

## Locked decisions — jellium campaign (deferred until after vacuum warm-up)

- **Geometry:** localised jellium = SLAB (fills periodic x,y face; localised in z).
  Box 25×25×140 Bohr. Slab 25 Bohr thick (half-width 12.5), centred z=0, faces ±12.5.
- **Density:** N=102 electrons in 25×25×25 = 15625 Bohr³ → n₀=0.00653, **r_s=3.32**,
  E_F=4.5 eV, ħω_p=7.8 eV. (User: keep N=102, moderate density, independent of the
  prior high-density run.)
- **Smoothing:** w = 0.5 Bohr → `edge_width=0.5` (erfc softening in
  `inqkit/jellium/localised_background.hpp`).
- **Grid spacing:** h = 0.5 Bohr (GS + both RT). Must pass `cutoff_guard.py`.
- **WP:** σ_WP=1, E=100 eV (k₀=2.71), mass=1 (electron), launched z=−20.5
  (8 Bohr from −12.5 slab face), moving +z.
- **Engine:** BOTH runs on inq-study (CAP-off vs CAP-on).
- **CAP (run 2):** two-sided sin² absorber, 10 Bohr/side at far ends
  (z∈[±60,±70], inner faces ±60, in vacuum → no bath over-drain), **η=−0.7 Ha**.
- **Diagnostic:** no-CAP (closed periodic) conserves energy_total → plateau =
  all deposited energy retained; CAP drains escaping flux → lower plateau. The
  **gap between plateaus = energy radiated to the boundaries** = the suspected
  "too much retained energy".
- **Still to resolve (jellium):** simulation length ("long"), observable cadence
  (user wants momentum distribution EVERY step, WP wavefunction every 10 steps),
  E_ss/E_sp decomposition (post-processing via twin decomposition), GS-first +
  GS notebook (z-density profile, potential+background vs z for WP-position
  validation), campaign doc + monitoring.

## Vacuum warm-up — design (current)

- **Theory:** non-interacting (user choice) → E_total = E_kinetic; Hartree/xc/
  external = 0 (bookkeeping check). ETRS propagator (correct for CAP; CN would
  renormalise and defeat absorption).
- **WP:** σ=1, E=100 eV (k₀=2.71), mass=1, matching the jellium WP.
- **Runs:** one env-driven binary; CAP_ETA=0 → no-CAP, CAP_ETA=−0.7 → CAP.
- **Record:** all energy components each step (total,kinetic,hartree,external,
  non_local,xc); momentum distribution; WP density frames (for density GIF).
- **Build:** `export INQ_SOURCE=.../inq-study` before `inq-run`.

## Key references

- Reuse: `ResearchProject/systems/vacuum/scripts/cap_sweep/run.cpp` (single WP +
  in-built CAP, inq-study, ETRS, energy/momentum writers).
- CAP mechanism: `inq-study` diff in `self_consistency.hpp` (complex vscalar +
  `real()` on external energy so CAP's imaginary term never inflates energetics).

## Milestone 2026-07-29 — IN-RUN extensive-kinetic observable (norm-division fix)

Prelude to the next localised-jellium run design: fix INQ's per-particle
(norm-divided) `energy_kinetic` (energy.hpp:55) at the WRITE layer, no engine
edit. Plan section: `docs/plans/norm-corrected-stopping-power.md` → "Extension
(2026-07-29)". Root-cause note: `docs/notes/inq-energy-normalization-error.md`.

### Done (verified)
- **`inqkit::observables::OrbitalKineticStats`**
  (`/local/data/public/skcb2/tddft/inq-stack/include/inqkit/observables/orbital_kinetic_stats.hpp`):
  per-orbital BARE T_i = ½(dV/N_grid)Σk²|ψ̃_i|² and norm_i (physical units;
  INQ's to_fourier is an unnormalized DFT — raw Parseval sum = (N_grid/dV)·∫|ψ|²dV,
  verified numerically: raw norm 13,183,593.75 = 843750/0.064 for norm-1 orbital).
  CSV: kin_bare_total (extensive), kin_normdiv_total (identity reconstruction of
  INQ's reported kinetic — matched 14.7769538333 Ha to all printed digits at t=0,
  = analytic ½k₀²+3/(4σ₀²)), norm_total, per-orbital columns, wall_ms self-timing.
- **Vacuum `run.cpp` extended**
  (`/local/data/public/skcb2/tddft/ResearchProject/systems/vacuum/scripts/wp_traversal_energy/run.cpp`):
  `WP_CAP2=1` double-sided CAP (−z band `absorbing(η,−0.5+w/2,w)` +z band via
  `perturbations::sum`; contravariant/fraction convention verified in source),
  `WP_EXTKIN`/`WP_EXTKIN_EVERY`, propagate wall-time in run_summary, final
  checkpoint (electrons.save + rt_state.txt).
- Cutoff guard PASS (σ₀=3, E=400 eV, h=0.4 → tail 0.00%, E_cut=839 eV).
- Comparison script:
  `/local/data/public/skcb2/tddft/ResearchProject/systems/vacuum/hypotheses/cap_norm_investigation/extensive_kinetic/compare_extkin.py`.

### ✅ VALIDATED (2026-07-29 02:11) — all acceptance checks PASS
- Pair completed (`dcap_extkin` + `dcap_baseline`, two-sided CAP, 700 steps).
  First TWO attempts crashed at step ~420: **/local/data was 100% FULL** (VTI
  writer ENOSPC — not the mid-run rebuild as first suspected). Fixed by
  WP_WF_EVERY=700 (t=0+final frames only; comparison is CSV-based). Disk later
  freed to 125G by someone else — still tight, USER SHOULD REVIEW.
- **Identity EXACT**: Σocc·T_i/norm_i == energies.csv:kinetic, 0.0 Ha at all
  701 steps. t=0 bare kinetic == analytic 14.777 Ha.
- **Fix works**: E_reported(final)=383 eV (pinned at remnant mean, norm 3.5e-6)
  vs E_corr → 0.00 eV == E0·norm (captured 100.0% of 402 eV). Post-hoc route
  (e_kin_ha·norm) agrees with in-run bare to 2.5e-9 eV.
- **Cost**: 0.42 ms/step self-timed (0.14% of the ~300 ms step, 1 orbital);
  run-level ON−OFF Δ = −15 ms/step (noise). Jellium-162: expect a few % at
  every-step cadence (one extra set-FFT vs ETRS's several per step); measure in
  pilot; WP_EXTKIN_EVERY available.
- Results in `docs/notes/inq-energy-normalization-error.md` (§ IN-RUN FIX
  VALIDATED), test-catalogue row added, figs+summary in
  `systems/vacuum/hypotheses/cap_norm_investigation/extensive_kinetic/`.

### Next
- Wire OrbitalKineticStats into the next localised-jellium run.cpp (all 163
  orbitals, WP_EXTKIN_EVERY from pilot timing); design the run per the
  2026-07-29 recap (extensive ledger + explicit CAP-sink bookkeeping).

## Milestone 2026-07-29 — replica_lz160_1cap CRASHED on full disk; ~130 GB freed (VTI prune, user-approved)

### What happened
- The LZ=160 one-sided-CAP replica run
  (`/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/scripts/replica_lz160_1cap/`,
  GPU 0, 8000 steps) **aborted 2026-07-28 20:54 at step 4382** (t = 87.6 a.u.,
  ~55%): `VTIImageDataWriter` failed writing
  `results/cap/raw/vti/density_delta/density_delta_t004382.vti` → std::runtime_error
  → abort. Cause: `/local/data` was **100% full (0 bytes free)**. Physics was clean
  (E ≈ −997.71 Ha, ~7 s/step, no instability). Failure email was sent by the
  dispatcher. The no-CAP leg never started.
- **Recoverable:** interior checkpoint at `results/cap/rt_ckpt/` with
  `last_step=3200`, t=64, 75 states. Resume: `WP_RESUME=1 WP_GS_DIR=... WP_CAP_ETA=-0.7
  WP_OUT=cap WP_N_STEPS=8000 WP_CKPT_EVERY=1600 inq-run` from the `wp/` dir
  (loses 1182 steps ≈ 2.3 h). Observables on disk up to the crash.

### Disk cleanup (DESTRUCTIVE, user-approved 2026-07-29: "Remove VTI dirs of superseded runs")
- Deleted **only `results/raw/vti/` directories** of 16 old
  `systems/jellium/run_*` runs, each strictly name-superseded by a completed
  `_v2`/`_v3`/`_attempt2` twin (twin verified to have non-empty observables) and
  with zero VTI references in hypotheses/, inqview, docs/. All CSVs, observables,
  run_summary, checkpoints kept. Freed ~130 GB → **125 GB now available (99%)**.
- Deleted (VTI only): classical L30 highdens E50/E300; classical L50 E50, E600;
  wp L30 highdens_sigma1 E50/E100/E200/E300; wp L50 E300, E50, E600,
  E20_sigma1, E25_sigma1, E50_sigma1, E200_sigma1, E300_sigma1.
- **SPARED deliberately:** `run_classical_n162_L50_E100` + `run_wp_n162_L50_E100`
  (E100 case-study pair — VTIs referenced by
  `docs/reports/classical-vs-wp-case-study.md` and
  `docs/reports/14-05-2026-meeting-emilio/case_study_E100eV.py`),
  `run_wp_n162_L50_E100_sigma1` (VTIs used by
  `docs/reports/2026-05-21-meeting-emilio/build_axial_gifs.py` / `build_density_diff.py`),
  `run_classical_n162_L30_E100_highdens` (its _v2 twin has EMPTY observables — old
  run may be the good one), all `_wf` σ-sweep pairs (old runs hold WP-included
  density, `_wf` hold bath-only — different fields, not redundant).

### Open / caution
- **125 GB may NOT cover both remaining legs**: CAP remainder (~4800 steps ≈
  40 GB VTI) + full no-CAP 8000-step leg (~85 GB) ≈ 125 GB — zero margin.
  Options before relaunch: prune Tier-2 (qsp_phase2/4 VTIs, ~80 GB, user approval
  needed) or reduce VTI cadence for the no-CAP leg.
- Ranked deletion-candidate list (Tiers 1–3) delivered in-session 2026-07-29;
  only Tier 1a executed.

### 2026-07-29 second prune pass (user-approved: "remove more of the wf, v1 superseded")
- Deleted VTI dirs of the three `_wf`-superseded σ runs (successors verified
  complete, zero VTI references): `run_wp_n162_L50_E100_sigma0p5` (15G),
  `_sigma3` (19G), `_sigma8` (20G) → **178 GB now free (98%)**.
- **`run_wp_n162_L50_E100_sigma1_v2` SPARED — do not treat v3 as its successor:**
  all report1 analysis scripts (`stopping_power_data.py`, `fig_matched_pair.py`,
  `fig_gs_decomposition*.py`) consume **v2**; `sigma1_v3` has ZERO analysis
  mentions (appears to be an abandoned later experiment, ~15G — candidate for a
  future pass if user confirms it is dead).
- Still spared (referenced by report/meeting figure scripts): E100 case-study
  pair + `run_wp_n162_L50_E100_sigma1` (v1) — ~47 GB more if user releases them.

### 2026-07-29 02:32 — replica RESUMED (launched, verified propagating)
- New script `.../replica_lz160_1cap/resume_jellium_replica.sh` (setsid-detached,
  survives session; log `jellium_replica_resume.log`): CAP leg `WP_RESUME=1` via
  the EXISTING `wp/run` binary (deliberately NOT inq-run — no rebuild of a
  possibly-drifted source; binary is the Jul-28 inq-study CAP build), then chains
  the fresh no-CAP 8000-step leg + COMPLETE email, mirroring the original script.
- Pre-flight: GPU 0 free (extkin test DONE marker), 178 GB disk free vs ~125 GB
  projected need.
- **Restore verified bit-consistent**: resumed step 3200 e = −997.341832913767
  == original log at step 3200 (all printed digits). ~10 s/step → CAP leg done
  ≈ +9–13 h, no-CAP ≈ +15–22 h after that. Segment CSVs: `*.from3200.csv`.
- Post-processing must CONCATENATE `observables*.csv`/segment files in step order
  (final-timestep-checkpoint rule).

### 2026-07-29 02:46 — resume attempt 1 crashed (VTI collision); FIXED, attempt 2 running
- Attempt 1 aborted at step 3210: `VTIImageDataWriter: file already exists and
  overwrite=false: .../wavefunction_wp/wavefunction_t003210.vti`. Cause: resume
  writers open with `overwrite=!RESUME`, but the ORIGINAL run died at 4382 —
  PAST the 3200 checkpoint — leaving orphan frames 3201–4382 that the resumed
  (recomputing) propagation must rewrite.
- Fix: deleted the 2599 orphan frames with step > 3200 across all five VTI dirs
  (density_delta 1182, density_delta_coarse 1181, density_total 59, density_wp 59,
  wavefunction_wp 118), incl. the truncated `density_delta_t004382.vti`.
- **General rule for any resume after a crash BEYOND the last checkpoint: prune
  every per-step output with step > last_step before relaunching.**
- Attempt 2 launched 02:46, verified past the collision point (step 3234,
  e drifting smoothly, ~10 s/step). Failure emails fired correctly on both aborts.
- 02:52 full health check PASSED: e(3250)/e(3270) bit-identical to original log;
  VTI frames rewriting past 3210 with no overwrite errors; segment CSVs
  (`*.from3200.csv`) start at 3201 with NO duplicate steps (attempt 2 truncated
  attempt 1's rows); ~7 s/step → CAP leg ETA ≈ 12:00 same day; 189 GB free.

## 2026-07-29 (later): extkin_plateau_E100 wired — first jellium run with the IN-RUN fix

User-interviewed design (full decision log in
`docs/plans/norm-corrected-stopping-power.md` "Run design (2026-07-29)"):
35×35×120 box h=0.5; slab N=92 / r_s=4.0 exact / thickness 20.13 Bohr
(faces ±10.07, edge 0.5); 46 occ + 16 extra states, T≈100 K; WP σ=1.5 /
100 eV, launch z=−17.5 (5σ rule met; ×1.58 free-dispersion at slab entry
ACCEPTED by user); two-sided CAP 15 Bohr/side η=−1.0 (inner edges ±45);
dt=0.04 × 1500 steps (t=60 a.u.). OrbitalKineticStats ALL states EVERY step.
USER SCOPE CUTS (recorded): no no-CAP twin, no dt=0.04 vacuum gate — first
CAP run at this dt, absorption quality unverified; E_plateau single-source.

Files (all new, untracked):
- `ResearchProject/systems/localised_jellium/shared/configs/slab_n92_L35x35x120_w0p5.hpp`
- `.../scripts/extkin_plateau_E100/{gs,wp}/run.cpp` (clones of
  wp_cap_energy_plateau + OrbitalKineticStats + η=−1.0 + per-step timing +
  final ckpt; resumable)
- `.../scripts/extkin_plateau_E100/run_extkin_plateau.sh` (autonomous:
  GS → CAP → notebook; emails; EXTKIN_GPU env; DONE marker)
- `.../hypotheses/extkin_plateau_E100/build_extkin_plateau_report.py`
  (partial-tolerant; density GIF + E_corr/plateau/identity/cost cells)

Gates: cutoff guard PASS (0.00% aliased, k_Nyq=6.28 vs p0+3σp=4.13).
Compile probes + launch: see below / next milestone.

GPU NOTE: the replica campaign RESUMED 02:39 on GPU 0 (CAP run from step
3200/8000, ~7–12 s/step → many hours + no-CAP after). extkin_plateau launches
on GPU 1 (EXTKIN_GPU=1). Revised wall estimate from replica timing: ~2–4 h,
not minutes. /local/data at 99% (125G free) — run writes ~2 GB.

### Launch (2026-07-29 ~03:20): extkin_plateau_E100 RUNNING on GPU 1

Compile probes PASS (gs + wp, inq-study). Dispatcher setsid-detached
(`run_extkin_plateau.sh`, pid 2845144, EXTKIN_GPU=1) — VERIFIED live: GS SCF
iterating all 62 states (46 occ + 16 extra), top eigenvalues ≈ −0.073 Ha.
Fully autonomous: GS → CAP (1500 steps, ckpt/200 + final) → notebook
auto-build; email milestones via notify.py. Expected ~2–4 h total.
Log: `.../scripts/extkin_plateau_E100/extkin_plateau.log`;
DONE marker: `EXTKIN_PLATEAU_DONE.txt`; results: `wp/results/cap/`.
To extend afterwards: WP_RESUME=1 + larger WP_N_STEPS (final ckpt present).
Reminder of open caveats: dt=0.04 CAP absorption ungated; no no-CAP twin.
GPU 0 still owned by the replica resume (step ~3300/8000 + no-CAP after).

### COMPLETE (2026-07-29 13:12): extkin_plateau_E100 results

Run finished (1500/1500 steps, t=60; ~10 h wall — ~25 s/step, ~2× slowed by
host/disk contention with the concurrent replica on GPU 0). Final ckpt at 1500
(extendable). Notebook built+executed 0 errors after fixing a builder API
mismatch (_nbreport md()/code() RETURN cells, no anchor kwarg; cells list ->
build(cells, path)): `hypotheses/extkin_plateau_E100/extkin_plateau_E100_study.ipynb`
(77 KB, GIF path-referenced, 76 frames).

Headline (2 s.f.):
- **E_plateau = 4.1 eV** (corrected E_corr − E_GS; windows 4.6/4.2/4.1 eV,
  final drift −0.006 eV/a.u. — converged; residual WP norm 9.7e-4 ⇒ ≤0.1 eV).
- Reported total would say 22 eV — the 17.8 eV norm-division artifact
  inflates deposition ~5×; removed in-run.
- Identity EXACT with 62 interacting states: 0.0 Ha, all 1501 steps.
- Absorption 99.9% (9.7e-4 vs ~8e-4 predicted) ⇒ dt=0.04 CAP behaves as the
  dt=0.01 calibration — ungated-dt caveat empirically retired.
- Bath norms 0.99999999999 (CAP does not drain the slab).
- OrbitalKineticStats cost 0.6%/step (143 ms of 23.7 s), every step, 62 states.

Open: no no-CAP twin (user cut). Replica still on GPU 0 (step ~5240/8000).

### 2026-07-29 (later): stopping-power sections added to the study notebook

Two new sections in extkin_plateau_E100_study.ipynb (builder updated, executed
0 errors): (1) "Orbital-free stopping power" — S = E_pl/L = 4.1/20.13 =
**0.20 eV/Bohr** vs Bethe bulk point-charge (ω_p²/v²)ln(2v²/ω_p) =
0.73 eV/Bohr (ratio 0.28; suppression attributed — inference — to
20-Bohr slab ≪ 79-Bohr wake wavelength, packet form factor, subbands).
(2) "Orbital-dependent stopping power" — position via ∫⟨p_z⟩dt (every step)
vs density centroid (VTI every 10 steps + z_mean every 100): agree to
0.25 Bohr while norm>0.995; orbital KE 109.9→108.4 eV across the slab →
S_orb = **0.07 eV/Bohr** (endpoint AND mid-slab gradient). S_orb ≈ ⅓·S_free
QUANTIFIES the standing caveat: the WP orbital's KE loss misses ~2/3 of the
deposition (Hartree/xc channels) — total-energy deposition remains the
primary quantum stopping measure. New figs: fig_wp_position/ke_time/ke_position.
