# Rule: Every run/analysis notebook ships a density-matrix GIF

Apply to: every run-notebook and twin/analysis notebook of a TDDFT run that has
saved density frames — `run-notebook` + `twin-run-analysis` skills, their builders
(`run_notebook_builder.py`, `twin_notebook_builder.py`), and any `.ipynb` under
`ResearchProject/systems/**/hypotheses/`. Always on.

## The one rule

**Every notebook describing a run (or a classical-vs-WP pair) MUST include an
animated density-matrix GIF — the real-space density n(r,t) evolving over time —
not only static density carpets/snapshots.** (User decision, 2026-07-15.)

- **The density GIF is the 2D field in the propagation x–z plane** (mid-y slice),
  NOT a 1D z-lineout. It animates the **matrix of three kinds** —
  `density` n(x,z,t), `induced` Δn=n(t)−n(0), `instantaneous` Δn=n(t)−n(t−Δt) —
  rendered by the canonical inqview functions
  `inqview.visualisation.make_density_gif_battery` (per-run: kinds × {total, wp,
  bath}) and `make_twin_density_matrix` (classical vs WP vs WP−classical for a
  twin pair). Mid-y xz slice, physical-order VTIs (never fftshift'd — the
  vti-coordinate-mapping rule), LINEAR | LOG panels, slab faces dashed. Written to
  the notebook's own directory. (A 1D n(z,t) lineout is an OPTIONAL companion, not
  the requirement — user clarified 2026-07-15 that "density matrix" means the full
  xz field, all three kinds, for classical, WP, and WP−classical.)
- **It MUST be DISPLAYED inline, near the TOP of the notebook** (a "visual
  intuition" section), not merely written to disk. Saving the file and printing
  "wrote …gif" is NOT sufficient — the reader must SEE the animation when they
  open the `.ipynb`. Emit it as cell output via
  `IPython.display.Image(filename="density_evolution.gif")`, which base64-embeds
  the bytes as `image/gif` into the stored outputs so it animates on reopen
  without the sidecar file. (User decision, 2026-07-15: GIFs go at the top and
  serve as the first-read visual intuition for the run.)
- Fixed axis limits across frames (global max for n, symmetric for Δn) so the
  animation is comparable frame-to-frame; mark the slab faces.
- It is IN ADDITION to (not a replacement for) the static density carpets/maps and
  the pairwise-energy GIF.
- If a run saved NO density frames, the notebook must still emit the GIF cell (it
  prints "no density frames"), and the run should be re-run with frame saving
  (`LJ_SAVE_EVERY>0`) enabled — a WP/classical run intended for a notebook SHOULD
  save frames.

## Why

The density evolution is the most direct picture of the quantum effect (dispersion,
reflection, tunnelling, capture) — a static carpet compresses time onto one axis and
hides the moment-to-moment behaviour. An animated n(r,t) is what the user reads first
to understand how the wavepacket departs from the classical projectile.

## How to apply

- `twin_notebook_builder.py` and `run_notebook_builder.py` emit the
  `density_evolution.gif` cell as a MANDATORY section (see the "Density-matrix GIF"
  cell in `twin_notebook_builder.py`).
- When adding a new notebook builder or hand-writing a notebook, include the density
  GIF cell; a notebook without it is incomplete.
- Ensure the producing run.cpp saves density frames (`LJ_SAVE_EVERY`); the classical
  `proj_dyn/run.cpp` and WP `phase5_wp/run.cpp` already do.
