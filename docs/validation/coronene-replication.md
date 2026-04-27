# Validation status: coronene replication framework

Per the plan (`docs/plans/coronene-replication.md` §5) and per
`.claude/skills/simulation-validation.md`, only Tier A is run in this round
on representative runs. Tiers B and C are documented but deferred.

---

## Tier A — Fast, always run on representative runs

Representative runs:

| What | Run(s) |
|---|---|
| Every new GS save (full Tier A) | `save_gs/gs_35x35x60_cut40`, `save_gs/gs_35x35x80_cut40`, `save_gs/gs_35x35x40_cut40` |
| Every propagation (GPU-exec + results-tree check only) | `run_base`, `run_E30`, ..., `run_35x35x40` |

### Checks (per the plan §5)

- [ ] **GPU execution** — run launched via `inq-run` (no `--cpu`); INQ log
      contains the GPU/CUDA initialisation banner; `nvidia-smi` shows the
      run's PID on the assigned GPU during execution.
- [ ] Coronene xyz parses to 36 atoms, all z = 0 within 1e-12.
- [ ] `cell.contains()` is true for every atom (defensive check is in each
      `save_gs/<sig>/run.cpp`).
- [ ] SCF converges to 1e-6 Ha.
- [ ] No NaN / Inf in printed energies.
- [ ] WP injection: `norm_after ∈ [0.97, 1.03]`, `max_overlap < 1e-3`
      against the occupied subspace.
- [ ] `RealField3DWriter` with `emit_raw=false, emit_vti=true` produces no
      zero-byte `.raw` files.
- [ ] `results/` tree passes the `find` checks at the bottom of
      `docs/results_folder_structure_spec.md` §22.
- [ ] `coronene_postprocess.py --results <dir>` completes every phase
      without error.

Status will be filled in below as runs complete.

### Per-run outcomes

| Run | GPU exec | xyz | cell | SCF | NaN | WP norm | WP overlap | spec tree | postprocess |
|---|---|---|---|---|---|---|---|---|---|
| save_gs/gs_35x35x60_cut40 | — | — | — | — | — | n/a | n/a | n/a | n/a |
| save_gs/gs_35x35x80_cut40 | — | — | — | — | — | n/a | n/a | n/a | n/a |
| save_gs/gs_35x35x40_cut40 | — | — | — | — | — | n/a | n/a | n/a | n/a |
| run_base                  | — | — | — | — | — | — | — | — | — |
| run_E30                   | — | — | — | — | — | — | — | — | — |
| run_E800                  | — | — | — | — | — | — | — | — | — |
| run_s0p33                 | — | — | — | — | — | — | — | — | — |
| run_s3                    | — | — | — | — | — | — | — | — | — |
| run_E800_s0p33            | — | — | — | — | — | — | — | — | — |
| run_E30_s3                | — | — | — | — | — | — | — | — | — |
| run_b18_35x35x80          | — | — | — | — | — | — | — | — | — |
| run_b6_35x35x80           | — | — | — | — | — | — | — | — | — |
| run_35x35x40              | — | — | — | — | — | — | — | — | — |

Legend: `—` = not yet run, `✓` = pass, `✗` = fail (with note).

---

## Tier B — Medium (deferred)

Tier B was deprioritised by the user for this round. To be revisited if
Tier A surfaces concerning behaviour, or before the framework is used for
paper-quality results.

- [ ] Restart check (load checkpoint, propagate a few steps, compare).
- [ ] Energy conservation: drift < 0.1 % over full propagation.
- [ ] Charge / norm conservation within 1e-3 over the run.
- [ ] CPU/GPU consistency on a short representative run.
- [ ] Coordinate-mapped vs raw-index LEED screen plots match (no
      four-corner split).

---

## Tier C — Expensive (deferred)

- [ ] Time-step convergence sweep (dt ∈ {0.020, 0.010, 0.005} a.u.).
- [ ] Cutoff convergence on the new geometry (30, 40, 50 Ha).
- [ ] Tsubonoya 2014 Fig. 2 quantitative reproduction at the paper window.
- [ ] Far-field b-scan trend monotonicity (would need additional b values
      between the existing `b6_35x35x80` and `b18_35x35x80`).
- [ ] Projectile vs paper comparison (`run_E800` vs `run_base`).
