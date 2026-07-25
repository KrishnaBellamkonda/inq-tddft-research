# Rule: `inq/` is immutable

Apply to: `inq/` (the upstream INQ source tree) — always on, every session.

## The rule

**Never modify, add, delete, or move any file inside `/local/data/public/skcb2/tddft/inq/`.**
It is the unmodified upstream INQ library (header-only C++17 GPU engine,
gitignored but present on disk). Treat it as strictly read-only.

This includes — but is not limited to — editing source/headers, adding new
files, applying patches, "quick fixes", build-config tweaks, or reverting the
existing CUB workaround in `inq/external_libs/gpurun/include/gpu/reduce.hpp`
(that change is pre-existing and must stay; do not touch it either way).

## Why

`inq/` must remain a clean, reproducible mirror of upstream INQ so that builds,
benchmarks, and physics results are attributable to *our* code, not to local
engine edits. Divergence from upstream silently invalidates every comparison and
makes the engine impossible to update.

## What to do instead

- Engine-level analysis or experimentation that needs a *modified* INQ goes in
  **`inq-study/`** — a byte-identical replica created exactly for this purpose.
  Modify `inq-study/`, never `inq/`. (Verify they are still byte-identical before
  relying on `inq/` as the reference, e.g. `diff -q inq/src/... inq-study/src/...`.)
- New library code that builds *on top of* INQ goes in **`inq-stack/`**
  (`inq-stack/include/inqkit/` for C++, `inq-stack/python/inqview/` for Python).
- Per-step hooks, observables, wavepacket/absorber logic, etc. can almost always
  be done in the wrapper (`inqkit`) by owning the `electrons` object and acting
  in the `real_time::propagate` callback — no engine edit required. See
  [[reference_inq_propagator_mask_absorber]] for a worked example (the mask
  absorber: applied entirely in the callback, `inq/` and `inq-study` untouched).

## On exception requests

If a task appears to *require* editing `inq/`, stop and surface it to the user
with the specific reason and the `inq-study/` alternative — do not edit `inq/`
on your own initiative. Only an explicit, unambiguous user instruction to modify
the upstream tree overrides this rule, and even then prefer `inq-study/` first.
