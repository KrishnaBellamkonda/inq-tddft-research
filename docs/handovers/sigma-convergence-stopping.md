# Handover: σ-convergence of stopping power → point-charge Lindhard

Plan: `docs/plans/sigma-convergence-stopping.md`. Branch
`convergence-gaussian-electron-previous-scheme`.

---

## Milestone: 2026-06-15 23:28 BST — σ=3 large-width probe LAUNCHED + plot revisions

### Current status
**Done + running.** Plot revisions applied and verified; σ=3.0 Bohr run set built,
smoke-tested, and launched detached (PPID=1, survives session exit). 7-velocity
ladder running across both GPUs (~11–12 h); emails BOTH figures + table on
completion. Decided in a grill-with-docs session (plan §12).

### Plot revisions (`hypotheses/06_sigma_convergence/sigma_sweep_report.py`)
- Points now at **nominal launch v₀** (not window-mean); **vertical linregress
  stderr bar only**, no horizontal bar.
- ★ peak marker → **vertical dashed line annotated "Lindhard peak"**; existing
  dashed line annotated **"k_F"** (`_annotate_lines`).
- New **companion energy figure** `figures/sv_convergence_energy.png` (x = E=½v²
  in eV, **log x**; reference + both dashed lines mapped to energy). Velocity
  figure `figures/sv_convergence.png` retained as primary. Reference curve
  computed once (90 pts) and shared by both (`_reference_curve`).
- Regenerated clean on all 23 points; v1p0 + σ=3 rows correctly pending.

### σ=3.0 run set
- UPF `shared/pseudopotentials/electron_gaussian_sigma3p0.upf` (rel err 5e-11 to
  C·erf(r/(σ√2))/r, C=2 Ry, in-file V(0)=0.532 Ry ≈ 0.266 Ha — weak/diffuse).
- Cfg `SV_Ladder_L50_sigma3p0` in `shared/configs/sv_ladder_L50_sigma0p5.hpp`
  (launch_z=**−13** = standard 4σ rule; σ=3 UPF).
- Run dir `run_classical_n162_L50_sv_sigma3p0/` + run.cpp (= sigma0p15 run.cpp
  with the one Cfg-typedef line changed) + own `build/`. Compiles clean.
- **Smoke PASSED** (v=1.0, 30 steps, GPU1): run_completed=true, final_z=−12.40,
  final_vz=0.99975, 6 VTI frames ×2, observables.csv with cod_*+density_l2,
  track every step, overlap t=0/end, manifest jellium-classical. (24 s/step here
  is startup-inflated over 30 steps; asymptotic ≈17 s/step.)
- Ladder {3.0,2.0,1.3,**1.0** (new peak-refinement),0.8,0.6,0.2}; N_STEPS reused
  {3.0:300,2.0:450,1.3:700,1.0:700,0.8:700,0.6:700,0.2:1000}.
- Launcher: `hypotheses/06_sigma_convergence/orchestrate_sigma3.sh`, detached via
  `setsid nohup` (PID 4114726, PPID=1). GPU0=(v0p2,v1p3,v0p6)=2400 steps,
  GPU1=(v0p8,v1p0,v2p0,v3p0)=2150. Emails σ=3 result (both figures) at end.

### RESCHEDULED 23:35 BST — launch delayed 2 h (user request)
The immediate launch (PID 4114726) was **killed** (process-group kill) and the
two partial run dirs (results/v0p2, v0p8) removed. Re-armed via
`launch_sigma3_delayed.sh` (DELAY_SECONDS=7200), detached (setsid, **session
leader, PPID=1**), PID 4121311 sleeping until **01:35:43 BST 2026-06-16**, then
execs `orchestrate_sigma3.sh`. results/smoke retained as the validation artifact.

### Email content (satisfies the "energy-scale comparison" request)
`orchestrate_sigma3.sh` → `sigma_sweep_report.py --email "3.0"` attaches BOTH
`sv_convergence.png` (velocity) AND `sv_convergence_energy.png` (S vs projectile
KE in eV, log x) with σ=3 overlaid against σ={0.15,0.25,0.35,0.5} — i.e. the
stopping power at σ=3 across the energy scale vs the previously obtained results.

### Monitoring / kill
- Before 01:35: `pkill -f launch_sigma3_delayed` (or kill PID 4121311).
- During the run: `pkill -f sv_sigma3p0/run`.
- Progress: `hypotheses/06_sigma_convergence/{launch_sigma3.log, orchestrate_sigma3.log}`
  + each `run_classical_n162_L50_sv_sigma3p0/results/<vtag>/run.log`.

### Expectation
σ=3 should land **well below** point-charge Lindhard (strong e^{−q²σ²}
suppression) — the first σ in this study expected to deviate strongly. ~11–12 h
after the 01:35 launch ⇒ email ~midday 2026-06-16.

---

## Update: 2026-06-15 — sweep complete; σ=0.5 anchor added to convergence plot

Status: **done**. The autonomous sweep finished cleanly 2026-06-15 07:46Z (all 18
new runs, σ∈{0.15,0.25,0.35} × 6 v). Per-σ emails sent.

Changed: `hypotheses/06_sigma_convergence/sigma_sweep_report.py` now also plots the
pre-existing σ=0.5 anchor set (`run_sv_sigma0p5`). That set uses a *flat* layout
(`results/<vtag>/observables.csv`) and carries `run_completed  = false` (halted at
the boundary rule), so a `layout` tag was threaded through `_obs_path`/`_completed`:
"flat" gates on file presence, "nested" on the completion flag. Colour map inverted
so smaller σ reads darker (nearer the black point-charge curve). Figure
`figures/sv_convergence.png` regenerated → **23 points** (σ=0.5 has no v=0.2 run).

σ=0.5 numbers reproduce the original `build_sv_extraction_notebook.py` extraction
(sim/LR: v3.0→1.01, v2.0→1.05, v1.3→0.94, v0.8→0.59, v0.6→0.77).

**Observation for the user to interpret (no verdict rendered):** across
σ=0.5→0.15 the extracted S(v) is essentially σ-independent at low/mid v
(v≲1.3 points coincide to ~3-4 sig figs) and only separates at high v
(v=3: S_A=0.0089→0.0103 as σ 0.5→0.15). Physically consistent with the form
factor exp(−q²σ²) only biting at large q (high v); the low-v shortfall vs
point-charge Lindhard near the peak (sim/LR≈0.6–0.78) does **not** shrink with σ.
Whether that shortfall is physical (TDDFT vs linear response / finite-size) or a
reference artefact is the open question — flagged, not adjudicated.

Next: optional — a study `.ipynb` in this folder formalising S(v;σ) vs the point
reference + Method A/B parity; the deferred per-run `analyse.py` derived-observable
pipeline still runs on saved base data anytime.

---

## Milestone: 2026-06-14 02:26 — ARMED for autonomous 04:00 launch

### Current status
**Built, smoke-tested, and ARMED.** A detached launcher (PID under init) sleeps
until 04:00 BST then fires the orchestrator. Fully autonomous; survives session
exit. ARMED confirmation email sent (email path proven end-to-end).

### What is armed
- Launcher: `hypotheses/06_sigma_convergence/launch_at_4am.sh` (setsid+nohup,
  detached, parent=init) → at 04:00 execs `orchestrate_sigma_sweep.sh`.
- Orchestrator runs σ ∈ {0.15,0.25,0.35} **sequentially**, each splitting its 6
  velocities across GPU0/GPU1 (~10 h/σ, ~30 h total), then emails the cumulative
  S(v) plot (`sigma_sweep_report.py --email`) to chiddukanna@gmail.com.
- One shared binary: `run_classical_n162_L50_sv_sigma0p15/run`; outputs to each
  `run_classical_n162_L50_sv_sigma0p{15,25,35}/results/<vtag>/raw/...`.

### Smoke test (PASSED, 30 steps, GPU1)
build OK; ran normally; produced the full contract: manifest (jellium-classical),
**6 density VTI frames** (total+system, decoupled from dense energy), observables.csv
with cod_x/y/z_bohr + density_l2 (11 rows), electron_track (31), state_energies,
gs eigenvalues/occupations, overlap_full at t=0 & end only. Method A/B parser
verified on the real CSVs. ~15–18 s/step.

### Files added this milestone
- `run_classical_n162_L50_sv_sigma0p15/run.cpp` (+ build/, run binary; smoke in results/smoke)
- `run_classical_n162_L50_sv_sigma0p{25,35}/` (output dirs)
- `hypotheses/06_sigma_convergence/{orchestrate_sigma_sweep.sh, launch_at_4am.sh, sigma_sweep_report.py, figures/sv_convergence.png}`

### Deferred (NOT in overnight critical path)
Full per-run `analyse.py` derived-observable pipeline (animations, energy
decomposition, density_fourier — the last needs dense VTI we don't have).
Base observables ARE saved, so this runs on the saved data anytime post-sweep.

### Monitoring / kill
- Progress: `hypotheses/06_sigma_convergence/orchestrate.log` and each
  `results/<vtag>/run.log`.
- To stop: `pkill -f launch_at_4am.sh` (before 4 AM) or `pkill -f sv_sigma0p15/run`
  (during the sweep).

---

## Milestone: 2026-06-14 — design locked (grill-with-docs), reference + UPFs done; runs not yet built

### Current status
Design fully locked via a grill-with-docs session. The single point-charge
Lindhard reference is implemented, tested, and now overlaid on every plot; the
three new σ pseudopotentials are generated and verified. The 18-run overnight
sweep is **not yet built or launched** — awaiting go-ahead, then a smoke test.

### Locked decisions (see plan §4)
- New σ = {0.15, 0.25, 0.35} (reuse σ=0.5 anchor). Reference = **single
  point-charge Lindhard**, infinite system, no box correction — never σ=0.2.
- Velocities v₀ ∈ {3.0, 2.0, 1.3, 0.8, 0.6, 0.2}; **m_e for ALL runs** (firm).
- Full ADR-0006 jellium_classical observable set + COD + manifest, BUT:
  density VTI = **6 frames/run**, scalars every step, overlap t=0/end only,
  momentum N/A. ⇒ no spectral loss-function (S(v) unaffected).
- 3 run dirs (one per σ); analysis in `hypotheses/06_sigma_convergence/`.
  ~2.5 nights on 2 GPUs.

### What changed (committed code/figures)
- `inqview/analysis/lindhard_elf.py`: added `stopping_power_point(v, kF)` — the
  one canonical reference (natural kinematic qmax; no 1/σ blow-up).
- `test_lindhard_elf.py`: +`test_stopping_point_converged`,
  `test_stopping_point_above_finite_sigma` (16/16 pass, 29 s).
- `build_sv_extraction_notebook.py` rebuilt + executed → regenerated
  `figures/sv_preexisting_extraction.png`, `sv_method_crosscheck.png`,
  `sv_fit_example_v1p3.png` against the single point reference; Barkas removed.
- `run_sv_sigma0p5/analyse_sv.py` + `make_sv_comparison.py`: single point
  reference, Barkas annotations stripped, marked superseded.
- UPFs generated + verified: `shared/pseudopotentials/electron_gaussian_sigma0p{15,25,35}.upf`.

### Files touched
- `/local/data/public/skcb2/tddft/inq-stack/python/inqview/analysis/lindhard_elf.py`
- `/local/data/public/skcb2/tddft/inq-stack/tests/python/inqview/analysis/test_lindhard_elf.py`
- `/local/data/public/skcb2/tddft/docs/reports/overnight-gaussian-classical-jellium/build_sv_extraction_notebook.py` (+ executed .ipynb + 3 figures)
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/run_sv_sigma0p5/analyse_sv.py`, `make_sv_comparison.py`
- `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/pseudopotentials/electron_gaussian_sigma0p{15,25,35}.upf`
- `/local/data/public/skcb2/tddft/docs/plans/sigma-convergence-stopping.md` (new)

### Key numbers
- Reference S_LR^point (Ha/Bohr): v=0.20→0.0144, 0.62→0.0603 (peak), 1.00→0.0416,
  1.94→0.0170, 2.98→0.0088.
- σ=0.5 sim/LR: 1.01 (v=2.98) … 0.59 (v=0.62) — gap to close as σ↓.
- Per-step wall: full set 18.1 s/step (VTI+overlap +15 % over stripped 15–16).
- 2 GPUs (`/dev/nvidia0,1`); nvidia-smi NVML-broken (compute fine).

### Tests / validation
- DONE: lindhard_elf 16/16; UPF rel-err 5e-11; reference 0.00 % converged.
- PENDING: smoke test of the new full-observable run.cpp before full launch
  (manifest, 6 VTI frames, energy drift <1 mHa, Method A/B <10 %).
- TODO: add test-catalogue row for `stopping_power_point`.

### Exact next steps
1. Build config + 3 per-σ run dirs from the full classical template with
   decoupled VTI cadence (VTI_EVERY=N_STEPS/5) + manifest.
2. inq-run build per σ dir (GPU).
3. Smoke test (σ=0.25, v=1.0) → validate, then launch 6-velocity ladders on 2 GPUs.
4. Per-run analyse.py; then `hypotheses/06_sigma_convergence/` study notebook
   (S(v;σ) vs single point reference + Method A/B parity).
