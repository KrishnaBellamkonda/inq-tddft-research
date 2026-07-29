# Phase 1 — Vacuum exit test (classical-highdensity-sv campaign)

## Verdict: **PASS — CLEAN-EXIT, no wraparound**

A moving Gaussian charge, built by the real inqkit code path
(`inqkit::jellium::gaussian_density`, the exact function the moving-projectile
perturbation calls at `moving_gaussian_projectile_perturbation.hpp:48`) on a
z-open `periodicity(2)` grid, is **cleanly CLIPPED by the finite z-grid** as its
center leaves the box. Only the in-box portion survives; nothing wraps to the
opposite (-z) face. This de-risks the projectile-exit mechanism.

## Setup

- Cell: `orthorhombic(35, 35, 85)_b .periodicity(2)` (x,y periodic; **z open**), EMPTY
  (no jellium background, no projectile potential added to any Hamiltonian; no SCF,
  no propagation — a pure density-construction test).
- Basis: `dx = 0.5` → INQ grid `70 x 70 x 175`, z-spacing 0.4857 b, z-axis span
  physical `[-42.0, +42.5]` (far face `+Lz/2 = +42.5`).
- `sigma_pot = 0.35355` (= sigma_WP / sqrt(2) for sigma_WP = 0.5).
- Sweep `z_center = {-20, 0}` (deep-interior refs) then `30 → 50` in steps of 1.0
  (crossing the far face +42.5).
- Per z_center: `n_proj = gaussian_density(basis, {0,0,z}, sigma_pot)`,
  `integral = operations::integral(n_proj)`, and the "wrap witness" = max density in
  the near-face region `z < -38`.
- Run machinery: `scripts/classical_highdensity_sv/vac_exit/run.cpp`
  (`results/exit_scan.csv`, `results/nproj_z*.vti`, `results/run_summary.txt`).

## Results

**Integral(n_proj) vs z_center** (should be ~1 interior, half at +42.5, →0 beyond,
no secondary rise):

| region | z_center | ∫ n_proj |
|---|---|---|
| deep interior | −20 … +41 | ≈ 1.00 (0.9385 only once clipping starts at 42) |
| just below face | 42 | 0.94 |
| **at face** | **42.5** | **≈ 0.50** (between 42→0.94 and 43→0.062) |
| just past face | 43 | 0.062 |
| beyond | 44 | 2.9e−6 |
| far beyond | 45 … 50 | 4.7e−14 → 3.9e−105 (monotone decay to 0, **no rebound**) |

The integral holds at ≈1.0 while the Gaussian is inside, passes through ≈0.5 at the
open face `+42.5`, and decays monotonically to zero past it — exactly the signature
of a hard clip of a normalised Gaussian by the box boundary. (The tiny +0.0003
surplus interior is the discrete-grid quadrature of a Gaussian normalised in the
continuum; physically 1.)

**Wrap witness (max density in z < −38) vs z_center:**

- **= 0 (exactly) for EVERY z_center**, including all the far-exit points 43–50.
- VTI cross-check: for `z_center = 45` (well past the face), `max(z < −38) = 0`
  while the only surviving density (max 1.2e−13) sits at the **+z** boundary — a
  clipped remnant, not a wrap. Peak-z tracks the center (+40.1 for z_center=40,
  pinned at +42.5 once the center is outside).

Any non-zero near-face bump would have signalled wraparound = FAIL. There is none.

## Dashboard

`ResearchProject/systems/localised_jellium/hypotheses/classical_highdensity_sv/phase1_vac_exit/`

- `integral_vs_zcenter.png` — ∫ n_proj vs z_center, +42.5 marked, half-line at 0.5.
- `wrap_witness_vs_zcenter.png` — wrap witness flat at 0 across the whole sweep.
- `montage_xz_slices.png` — n_proj x-z slices (mid-y, `load_vti` physical order) at
  z_center = 40, 42.5, 45, 48: Gaussian clipping at the +z face, nothing at −z.
- `make_dashboard.py` — regenerates the above from `results/`.

## Key numbers (2 s.f.)

- interior ∫ n_proj = 1.0
- ∫ at the face (z=42.5) ≈ 0.5 (half clipped)
- wrap-witness peak over the entire sweep = **0** (no wraparound)

## Provenance

- Build: GPU 1 (`CUDA_VISIBLE_DEVICES=1`), `inq-run`, default `inq/` engine,
  `TMPDIR=/local/data/public/skcb2/tddft/.build_tmp`.
- One implementation note: the pure density-construction object uses
  `extra_electrons(2)` (a valid populated basis) rather than `extra_electrons(0)`;
  an empty (0-electron) object triggered a spurious `CUDA ERROR: invalid argument`
  at the first `gpu::run`. The electron count does NOT affect `gaussian_density`,
  which reads only the basis geometry — the tested clip/BC behaviour is INQ's own.
