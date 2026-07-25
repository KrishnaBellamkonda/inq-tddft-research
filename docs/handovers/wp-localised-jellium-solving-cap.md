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
