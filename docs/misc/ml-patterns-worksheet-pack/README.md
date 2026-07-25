# ml-patterns worksheet pack

A self-contained collation of the `ml-patterns` campaign — interpretable / data-
driven pattern-finding on the **induced electron density** of a **classical vs
quantum-wavepacket projectile in bulk jellium** — assembled as **input for a
worksheet-building agent**. Its job is to turn this into a study worksheet
(foundational theory → tasks → results). This pack is a **neutral record**: it
states *what was done* and *the results*; it deliberately offers **no
interpretation** — that is left to the expert reader.

## What this pack contains

| Path | Contents |
|---|---|
| `01_algorithm_papers.md` | Seminal reference **per algorithm actually used** (POD, DMD, SINDy, PDE-FIND, form-factor/residual) + method→task map. Physics-context papers are intentionally out of scope. |
| `02_tasks.md` | Every task performed, in 4 chronological eras: input data, method, code path, output artifacts. |
| `03_results.md` | The recorded results per task: numeric tables + verdicts + figure/JSON pointers. Recorded flags (panel reviews, provenance measurements) are given as facts, no judgment. |
| `notebooks/` | 3 freshly-executed Jupyter notebooks with embedded outputs (figures + tables). |
| `figures/` | All result PNGs referenced by `03_results.md`. |
| `data/` | The `*.json` result files the notebooks and results doc are built from. |

## Notebooks (embedded outputs; snapshots)

- `notebooks/pdefind_recovery_demo.ipynb` — runs PDE-FIND **live** on synthetic
  fields with known governing PDEs (advection, diffusion, wave); the algorithm
  validation done in task T9.
- `notebooks/linear_response_residual.ipynb` — reconstructs the linear-response
  residual result (`|R(q)|` vs Gaussian form factor; temporal flatness) from the
  recorded per-σ JSON.
- `notebooks/pod_dmd_bath_structure.ipynb` — reconstructs the POD/DMD
  bath-structure sweep (POD rank, leading-mode energy, dominant DMD frequency;
  σ-sweep and velocity sweep) from the recorded summary JSON.

These are frozen snapshots. If re-execution is ever wanted, the live source lives
in `docs/campaigns/ml-patterns/` (kernels + runner scripts) and requires the
project venv (`/local/data/public/skcb2/tddft/venv/bin/python3`).

## The campaign in one paragraph (factual)

Across four eras the campaign applied interpretable/data-driven methods to the
GS-subtracted induced bath density `δn(r,t) = (n_total − n_wp) − n_bath^GS` of
matched classical and wavepacket projectiles in bulk jellium (`r_s ≈ 5.69`,
`L = 50`, grid `125³`, `ω_p ≈ 3.5 eV`): **(1)** signature gates T1–T7
(form-factor ratio, wake DMD, latent SINDy); **(2)** a governing-PDE discovery
redo T8–T14 (PDE-FIND/STRidge with three validation walls), which a scientific
panel then reviewed and a re-analysis re-measured (recorded facts: `density_system`
163 e WP-included vs `total − wp` 162 e bath-only); **(3)** a POD/DMD
bath-structure sweep on the blob-free bath (POD rank & DMD-frequency descriptors
across σ and velocity); **(4)** a time-domain linear-response residual test
(q-space ratio vs Gaussian form factor). See `02_tasks.md` and `03_results.md`
for the specifics.

## Source of record (for provenance / deeper reading)

- Campaign brief + task design: `docs/campaigns/ml-patterns/pattern-finding-in-wp-classical-runs.md`
- Deep-research method/physics survey: `docs/campaigns/ml-patterns/research/ml_induced_density_research.md`
- Methodology ADR: `docs/adr/0012-agnostic-pde-discovery-three-walls.md`
- Session handover (full narrative): `docs/handovers/ml-pattern-finding-wp-classical.md`
- Kernels (implementations): `docs/campaigns/ml-patterns/kernels/`
