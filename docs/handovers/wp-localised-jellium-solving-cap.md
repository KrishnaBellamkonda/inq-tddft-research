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
