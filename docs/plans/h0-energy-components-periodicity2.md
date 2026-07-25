# Plan: stream all energy components + re-run H0 at periodicity 2

**Requested (2026-07-07):** make inqkit stream *every* INQ energy component each
step; re-run the H0 insertion-energy experiment (WP + classical vs radius r) at
**periodicity 2** (open-z slab) using the new columns; update the H0 notebook.

## Why
The campaign only streamed `total/kinetic/hartree/xc`, so `E_ext = ∫ n·v_ext`
(the term the classical ghost changes) was invisible — folded into `total`.
INQ's `energy` object already holds all components with public accessors
(`inq/src/hamiltonian/energy.hpp`); the gap is purely in the inqkit wrapper, so
**no `inq/` or `inq-study/` edit is needed** (immutability preserved).

Periodicity 2 removes the full-PBC `G=0`/uniform-background ambiguity that makes
absolute charged-cell energies box-dependent (E_GS p3 −108.5 Ha vs p2 +60.4 Ha).
Slab → periodicity 2 (open only in z), NOT 0 (which would add lateral edges).

## Components streamed (all public accessors)
Sum to `total()`:
`kinetic + external + non_local + hartree + xc + exact_exchange + ion + ion_kinetic`.
Diagnostics (not in total): `nvxc`, `eigenvalues`.

## Steps
1. **inqkit (additive, header-only):**
   - `real_time/step_context.hpp` — add `energy_{external,nonlocal,ion,ion_kinetic,exact_exchange,nvxc,eigenvalues}` fields.
   - `real_time/real_time_session.hpp` — copy each from `data.energy().<accessor>()`.
   - `io/observables_writer.hpp` — matching `ObservableSelection` flags (default
     false → other runs unchanged) + `col()`/`val()` lines + doc block.
2. **run drivers:** `campaign_autorun/{classical,wp}/run.cpp` — enable all new
   `sel.energy_*` flags. Periodicity stays env-driven (`LJ_PERIODICITY`).
3. **Rebuild** both binaries via `inq-run` (INQ_SOURCE=inq-study, GPU). Editing
   run.cpp forces recompile → header changes picked up.
4. **Re-run H0 at periodicity 2:** new dispatcher `rerun_h0_p2.py`; radii
   (4,12,20,28,36,40) × {wp, cl}; `LJ_PERIODICITY=2`, GS =
   `runs/h2/gs_p2_lz120/checkpoint` (the p2 GS); output `runs/h0_p2/{tag}_r{r}_p2`.
   Single-point 3-step runs → minutes on GPU.
5. **Validate:** per-row `total ≈ kinetic+external+non_local+hartree+xc+exact_exchange+ion+ion_kinetic`
   to solver precision (record in `docs/validation/test-catalogue.md`).
6. **Notebooks:** regenerate H0 study + run-evidence + representative notebooks to
   show the p2 energy-component decomposition (E_ext(r), etc.). Neutral — no
   interpretation; provisional box only. Reference GS = p2.

## Status — COMPLETE (2026-07-07)
- [x] 1 inqkit  [x] 2 run.cpp  [x] 3 rebuild  [x] 4 re-run  [x] 5 validate  [x] 6 notebooks

Results: 12 periodicity-2 runs in `runs/h0_p2/{wp,cl}_r{4..40}_p2`; every row's
`total` == Σ(8 components) to ~1e-13 Ha. WP p2 excess 81.2→79.6 eV (flat, near the
81.6 eV zero-point); classical p2 excess 185→12 eV (decays). E_ext now measured:
classical r-dependence sits in E_ext (153.6→147.2 Ha); WP shows E_ext↑ vs U_H↓.
Notebook additions: H0 study "periodicity-2 measured decomposition" section +
`runs/H0_p2_runs.ipynb` evidence table (both executed, 0 errors).

### Follow-up (2026-07-07, later) — interpretation aids + extended r
- `H0_p2_interpretation.ipynb` (builder `build_h0_p2_interpretation.py`): waterfall
  decomposition (sums to total), GS charge-distribution (n_-, n_+, diff xz/yz/xy +
  profile), extended-r excess plot. 4 figs, 0 errors.
- Extended-r sweep `rerun_h0_p2_far.py` (Lz=200 box, own p2 GS E_GS=60.25 Ha),
  r={4..76}, wp+cl, sum-checked 1e-13. Classical excess bottoms at 0.7 eV near
  r≈52 (not r=40, still 12 eV) then edges up to 1.8 eV by r=76; WP flat 75.5–77.4 eV.
