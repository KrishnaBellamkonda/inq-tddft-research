# Handover — jellium slab→bulk L_slab sweep (`lz_bulk_sweep`)

**Rolling file. Latest milestone at top.**
Plan (the authoritative design): `/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/docs/plans/jellium-slab-extend-Lz.md`
Branch `quantum-stopping-power`. Machine CSD3, `ampere`,
account `mphil-nikiforakis-skcb2-sl2-gpu`.

---

## STATUS 2026-08-05 17:00 — PILOT PASSED AUTONOMOUSLY; PRODUCTION RUNNING; FIRST S(L) RESULT

**Chain state:** GS ×4 COMPLETED (3–9 min each, all gates passed); smoke
COMPLETED exit 0 (**first compile of all binaries clean, all 8 t=0 gates
passed**); pilot ×8 COMPLETED; **pilot gate PASSED 13:56 UTC (job 32886169)**
→ production auto-released. wp array: 3/12 done (all of s0p5_L15), 3 running;
cl array + vac PENDING on **AssocGrpGRESMinutes** — the account's GPU-minutes
quota is throttling, not an error; jobs start as budget frees (the sigma56
re-run 32880125 is also consuming it). Measured s/step 1.75–4.89 → production
projection 70 GPU·h (matches plan).

**First physics (v = 3.0, corrected deposit, eV/Bohr; PILOT_REPORT.md):**

| σ_WP | half | S(15) | S(25 anchor) | S(35) |
|---|---|---|---|---|
| 5 | classical | 0.131 | 0.133 | 0.133 |
| 5 | WP | 0.309 | 0.396 | 0.445 |
| 0.5 | classical | 0.373 | 0.381 | 0.387 |
| 0.5 | WP | 0.189 | 0.167 | 0.160 |

- σ=5 classical is L-INDEPENDENT (already bulk-like; surface term ~0).
- σ=5 WP RISES with L — the WP/classical gap GROWS with thickness (2.4× at
  L=15 → 3.3× at L=35); naive 1/L fit → S_bulk(WP) ≈ 0.55 vs classical 0.13,
  i.e. the excess is NOT a 1/L surface term (single-velocity, provisional;
  linearity residual at the L=25 point ~2 %).
- Bulk-likeness WARN fired exactly as predicted: L=15 boxes n(0) −3.2/−3.6 %
  off n0, interior Friedel peak-to-peak ~23 % of n0 (L=35: ~18–20 % — the
  oscillation is large at r_s=4.18 everywhere; the single-point n(0) check is
  luck-of-phase, peak-to-peak is the honest metric).
- Classical plateau-not-settled notes are the expected receding E_PS tail
  (endpoint is corrected); σ=5 WP plateaus are settled.

**Fix during this check (commit 708deca):** the E_PS(t_f) cut is now
CLASSICAL-ONLY in `lzb_stopping.measure`. The pilot exposed that subtracting
E_PS from the σ=0.5 WP runs — whose dispersed packet SURVIVES (norm ~3e-2)
with a large NEGATIVE E_PS — inflated S ~2.5× and made the trace spuriously
non-monotone vs its no-cut anchor. On the fixed convention all four traces are
monotone in L. At σ=5 the cut was a no-op (norm ~4e-10) so nothing else moved.

---

## STATUS 2026-08-05 (earlier) — MACHINERY BUILT + VALIDATED, PILOT-GATED CHAIN SUBMITTED

**User instructions this milestone:** "Submit this plan too and ensure
everything is autonomous. … do one or two runs of velocity for all the Lz.
Check if everything is alright before committing the massive number of GPU
hours." → implemented as a PILOT-GATED autonomous chain: v = 3.0 first at all
four boxes on both halves (~20 GPU·h), an automated gate job, and only on PASS
does SLURM release the remaining ~70 GPU·h (production arrays sit
`afterok:<gate>`; a gate failure leaves them DependencyNeverSatisfied with
`PILOT_REPORT.md` + email saying why).

### What was BUILT (all new; sigma56_sv clones, validated logic untouched)

- **Config:** `ResearchProject/systems/localised_jellium/shared/configs/lzb_boxes.hpp`
  — the four box presets selected at RUNTIME via env `LZB_CFG` (one binary per
  run type serves all boxes): s0p5_L15 (L_z 75, N 60, launch −19),
  s0p5_L35 (95, 140, −29), s5p0_L15 (95, 60, −22.5), s5p0_L35 (115, 140, −32.5).
  Per-family standoff 11.5 (σ=0.5) / 15 (σ=5) Bohr; face→CAP 17.5 / 27.5.
- **Binaries:** `scripts/lz_bulk_sweep/{gs,wp,classical,vac}/run.cpp`.
  gs+wp converted from compile-time Cfg to the runtime preset; classical was
  already env-driven (EXTRA_STATES made a knob — the sigma56 binary pinned 24,
  wrong for N = 60/140); vac is a verbatim copy (provenance strings only).
- **Dispatchers:** `shared/bin/lzb-params.sh` (bash mirror of the presets +
  LITERAL step tables), `run-lzb-{gs,smoke,wp,cl,vac,pilotgate,finalize}.slurm`,
  `submit-lz-bulk-sweep.sh` (the chain; `--exclude=gpu-q-2,gpu-q-25` by
  default, `LZB_EXCLUDE=none` to lift). wp/cl have three modes: `smoke`,
  velocity-index (pilot/repair), `--array=0-11%4` (production, cheap-first).
- **Analysis:** `hypotheses/lz_bulk_sweep/lzb_stopping.py` (adapter over the
  validated `e_absorbed.measure_dir`; per-box E_GS + L_slab divisor + E_PS
  monopole correction; `anchors()` reads all three L = 25 anchor sources),
  `pilot_gate.py` (HARD: completeness, ledger closure ≤ 1e-5 Ha, finite S;
  WARN: GS bulk-likeness from the density VTI, S(L) ordering, measured-cost
  projection), `finalize.py` (status → repair with per-run ~2× timeout cap →
  deliverables → CAMPAIGN_REPORT.md, email best-effort),
  `build_lzb_figures.py` (S vs 1/L + fits + S_bulk intercepts; every figure
  also written to `docs/reports/report2/drafts/draft1/figures/jellium_slab/`
  with `slab_` prefix at the house standard — per the standing feedback rule).

### VERIFIED this session (no GPU needed)

- `bash -n` clean on all 9 shell scripts; `py_compile` clean on all 4 python.
- `lzb_boxes.hpp` compiles standalone (g++ -fsyntax-only).
- `python lzb_stopping.py` self-test: r_s = 4.1815 for all four presets,
  LITERAL step tables == formula (all 16), dispersion table matches the plan
  (σ=5 L=35 slowest exit t = 25 a.u. < t_ov = 32.8 — every velocity is
  transversely clean, no velocity dropped).
- **All three anchor loaders read REAL data**: σ=5 from `s56_S_summary.csv`
  (e.g. classical v=2.0 → 0.2408), σ=0.5 WP from `sigma_sweep_S_deposit.csv`
  (0.239 at v=2.0), σ=0.5 classical from `S_of_v_cap.csv` with the monopole
  correction — reproduces the recorded **0.760** at v=2.0 exactly; its
  z_final/v_final are consistent with launch −24, which RESOLVES the plan's
  anchor-geometry verification item.
- `pilot_gate.py --report` end-to-end: fails all 8 pilot checks on missing
  runs (correct pre-run behaviour), writes PILOT_REPORT.md, email degrades to
  a log line (no Gmail creds), cost projection 67 GPU·h production.
- `finalize.py --status-only`: enumerates all 48 expected runs (32 production
  + 16 vacuum).
- `build_lzb_figures.py`: assembled 17 anchor rows, wrote `S_of_invL.png`
  locally AND `slab_S_of_invL.png` to the report-2 folder (anchor-only
  placeholder; overwritten with real data by the finalizer).

### NOT yet verified

- **The C++ has never compiled against INQ** — the chain's smoke stage
  (`run-lzb-smoke.slurm`) is deliberately the first compile and afterok-blocks
  everything if it fails (the sigma56 pattern). Most likely failure: a missed
  Cfg→B conversion in wp/run.cpp.
- Any physics in the new boxes; the GS SCF behaviour at N = 60/140.
- The measured s/step (pilot gate re-projects from real numbers).

### The submitted chain — jobs 32886156–32886177 (submitted 2026-08-05 07:25)

| stage | job(s) | dependency |
|---|---|---|
| gs s0p5_L15 (builds binary) | 32886156 | — |
| gs s0p5_L35 / s5p0_L15 / s5p0_L35 | 32886157/8/9 | afterok 156 |
| smoke (wp+cl builds, 8 t=0 gates) | 32886160 | afterok 156-159 |
| pilot v=3.0, wp+cl per box | 32886161–32886168 | afterok 160 |
| **pilot gate** | **32886169** | afterany 161-168 |
| wp production array 0-11%4 | 32886170 | **afterok 169** |
| cl production array 0-11%4 | 32886171 | **afterok 169** |
| vac ×4 boxes | 32886172–32886175 | afterok 169 |
| finalize ×2 | 32886176, 32886177 | afterany chain |

Dependency wiring VERIFIED in squeue after submission: production and vac hang
on `afterok:32886169`, so a pilot-gate failure leaves them
DependencyNeverSatisfied. Kill everything:
`scancel 32886156 … 32886177` (full list in the submit log / lzb-smoke output).
Node excludes active: gpu-q-2, gpu-q-25.

Inspect without changing anything:
    cd ResearchProject/systems/localised_jellium/hypotheses/lz_bulk_sweep
    python finalize.py --status-only        # run status
    python pilot_gate.py --report           # rebuild pilot report, no gating

Resume/extend a run:
    sbatch --export=ALL,LZB_CFG=<preset>,LJ_RESUME=1 shared/bin/run-lzb-wp.slurm <vidx>

### Open / follow-ups

- Run notebooks (density-GIF rule) for the 48 runs — builder not yet written;
  `finalize.py` already imports `build_run_notebooks` if it appears.
- Contingent L=50 point if L=15 proves non-bulk-like or off the 1/L line.
- Email: `python -m inqview.email setup` (one-time, interactive) would turn the
  gate/finalizer disk reports into emails.

---

## STATUS 2026-08-05 — DESIGN LOCKED VIA USER INTERVIEW

Locked: L_slab {15, 25, 35}; σ_WP {0.5, 5} twins; all 4 velocities; σ=0.5@25
NOT re-run (existing 85-box data anchors its family) ⇒ per-σ geometry;
VTI cadence ~3× coarser; fully autonomous + pilot-first. Full rationale,
geometry tables, caveat register and cost math: the plan file.
