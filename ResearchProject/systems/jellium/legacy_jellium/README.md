# legacy_jellium — archived runs and checkpoints

This folder holds jellium runs and ground-state checkpoints that were
either **flawed by design**, **never finished**, or **superseded** by a
later configuration. Nothing here is a current production run; everything
here is kept for reproducibility and post-mortem inspection only.

## Contents

### `runs/`
| Dir | Why archived |
|---|---|
| `run_base/` | L=60 cubic, **N = 128 (partial fill of |G|²=6 shell)**. Produces a non-uniform GS density (~20 % spatial modulation) and "blob"-like KS orbitals because the shell is not closed. The propagation completed (320 steps in ~83 min), so all observables exist, but the physical interpretation as "WP through uniform jellium" is compromised. |
| `run_base_n514/` | L=60 cubic, **N = 514 (closed-shell |G|²=16, 4× density)**. Propagator never completed even one step in 78 minutes wall time on a free A30 — INQ's per-step cost (likely the orthonormalisation Cholesky and overlap matmul) scales super-linearly with `n_states`, and at 257 occupied states it is prohibitively slow. Replaced by the lower-`n_states` `run_base_n138` path. |

### `save_gs/`
| Dir | Why archived |
|---|---|
| `gs_L40_cubic_N38/` | Legacy L=40 cubic, N=38 closed shell. Pre-redesign; superseded by L=60 runs. |
| `gs_L40_cubic_N40/` | Legacy L=40 cubic, N=40 open shell. Pre-redesign. |
| `gs_L60_cubic_N128/` | GS for the partial-shell N=128 run — see `runs/run_base/`. |
| `gs_L60_cubic_N135/` | Open-shell N=135 sketch; never used in any propagation. |
| `gs_L60_cubic_N514/` | GS for the high-density N=514 run — see `runs/run_base_n514/`. |

### `checkpoints/`
Saved electron states corresponding to the archived `save_gs/` entries.
- `gs_L60_cubic_N128/`
- `gs_L60_cubic_N514/`

(Legacy `gs_L40_*` checkpoints were never present on disk — only the
`save_gs/` runner sources.)

## What is *not* archived

- **`save_gs/gs_L60_cubic_N138_dx0p55/`** and **`checkpoints/gs_L60_cubic_N138_dx0p55/`**
  — current production GS (closed-shell, low-energy WP, coarsened grid).
- The six `run_E*_*/` variant skeletons (`run_E50_s0p53/`, `run_E200_s0p265/`,
  `run_E200_s0p53_N40/`, `run_E200_s0p53_tilt45/`, `run_E200_s2p0/`,
  `run_E400_s0p53/`) — **awaiting a base-run sign-off** and may need their
  configs updated to track the new Base_N138 conventions before launch.
- The shared configs (`shared/configs/base*.hpp`) are kept in place even
  though some of them (e.g. `base_highN.hpp` for the N=514 run) drive only
  archived runs. They would need re-targeting if the wave-packet redesign
  picks them up again.

## Recovery

To resurrect any archived run, move the relevant dir back to
`ResearchProject/systems/jellium/` (and any matching checkpoint back to
`checkpoints/`). All paths inside the run.cpp files are absolute so the
binary will still find its checkpoint after the move.
