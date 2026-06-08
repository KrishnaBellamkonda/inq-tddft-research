---
name: tddft-run-catalogue
description: >
  Maintain an automated, filterable CSV catalogue of every TDDFT run and the
  observables it produced. Use after a simulation completes to upsert its row,
  or to rebuild the whole catalogue, or to answer "which runs have observable X?"
  Complements the tddft-simulations skill (which produces the runs).
---

# TDDFT Run Catalogue Skill

A reproducible inventory of all runs under `ResearchProject/systems/` and which
observables each produced. Replaces hand-maintained run tables (e.g.
`docs/observables/catalogue.md` §3) with a script that reflects actual disk state.

## When to use
- A simulation just finished → **upsert** its catalogue row (post-run hook).
- You need to **filter runs by observable** ("which runs have a loss_function?",
  "which classical runs have electron_track + stopping_curve?").
- You want a **fresh full rebuild** of the catalogue.

## The artefact
- **CSV:** `docs/runs_catalogue.csv` (one row per run; key = `run_name`).
- **Scanner:** `scan_runs.py` (this skill dir). Pure stdlib — runs with any
  Python 3, but prefer the project venv for consistency:
  `/local/data/public/skcb2/tddft/venv/bin/python3`.

## Commands
```bash
PY=/local/data/public/skcb2/tddft/venv/bin/python3
SCAN=/local/data/public/skcb2/tddft/.claude/skills/tddft-run-catalogue/scan_runs.py

# rebuild the whole catalogue (run from repo root)
$PY $SCAN --all

# upsert one run after it completes (the post-run hook)
$PY $SCAN --run ResearchProject/systems/jellium/run_wp_n162_L50_E100_sigma1
```

## CSV schema
**Metadata** (parsed from `results/raw/run_summary.txt`): `system, run_name,
sim_type, run_completed, energy_ev, sigma_bohr, L_bohr, n_electrons, r_s, dt_au,
n_steps, total_time_au, write_every, norm_after, max_overlap, date_finished,
wall_time_s, run_path`. `r_s` is derived for cubic jellium as
`(3 L³ / 4πN)^{1/3}`.

**`n_observables`** then one 0/1 flag per observable:
- *raw:* `observables_csv, eigenvalues, wp_momentum_stats, wp_realspace_stats,
  state_energies, occupations, momentum_distribution, gamma_transitions,
  electron_track, density_total_vti, density_system_vti, density_wp_vti,
  density_delta_vti, wp_wavefunction_vti, overlap_wp, overlap_full,
  overlap_proxies, leed_screens`
- *post-processed:* `report_md, loss_function, energy_decomposition,
  gs_basis_decomposition, knudsen_ke, kl_divergence, bath_energy, stopping_curve,
  overlap_heatmap, momentum_before_after, planewave_decomposition,
  momentum_2d_map, density_fourier, density_gifs`

Detection: metadata from `run_summary.txt`; file-observables by basename under
`results/` (VTI frame contents are not walked — only series dir names recorded,
so it stays fast); dir-observables by characteristic `results/`-relative paths.
The scan reflects what is actually on disk — a flag is 0 if the artefact is
absent, even if it was "expected".

## Extending the observable set
When a new observable type is added to `tddft-simulations`, add one line to
`FILE_OBS` (basename test) or `DIR_OBS` (dir-suffix test) in `scan_runs.py` and
re-run `--all`. Keep this list in sync with `docs/observables/catalogue.md` §1–2.

## Filtering examples
```python
import csv
rows = list(csv.DictReader(open("docs/runs_catalogue.csv")))
# runs with a loss function:
[r["run_name"] for r in rows if r["loss_function"] == "1"]
# classical runs with a stopping curve at r_s≈5.69:
[r["run_name"] for r in rows
 if r["sim_type"]=="classical" and r["stopping_curve"]=="1" and r["r_s"]=="5.69"]
```

## Validation status
Validated 2026-05-31 on 80 runs: WP σ=5 → wp_momentum_stats=0, WP σ=1 → =1;
classical → electron_track/stopping_curve/bath_energy=1; coronene → leed_screens=1;
r_s: L50→5.69, L30→3.414. Detection matches simulation type in every spot check.
