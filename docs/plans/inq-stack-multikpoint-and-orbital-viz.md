# Plan: inq-stack multi-kpoint compatibility + per-kpoint orbital visualisation

Branch: `features/inqkit-multikpoint` cut from `features/jellium-ks-energy-observables` (the
current main-repo branch). Independent of `features/li-extensive` in the QKE
repo, but uses its GS checkpoint (`QuantumKickExtension/inq-codebase/Li/checkpoints/li_54_2x2x2_T200/`) as the heaviest smoke test.

## Context

The inqkit C++ library (`inq-stack/include/inqkit/`) was developed against
single-kpoint use cases — coronene (Γ-only finite cell) and jellium. Several
modules silently assume `electrons.kpin().size() == 1`:

| Module | Assumption | Behaviour with multi-kpoint |
|---|---|---|
| `fields/density.hpp::total` | k-sum is automatic via `electrons.density()` | ✅ Already correct (verified Phase 1). |
| `fields/density.hpp::orbital` | Takes `kpoint_index` arg, default 0 | ⚠️ Most callers pass nothing → only kpoint 0 is read. |
| `fields/orbital.hpp` | Same | Same |
| `observables/orbital_overlap.hpp` | Single-kpoint tracker | ❌ Multi-kpoint silently broken. |
| `observables/density_delta.hpp` | Operates on `RealField3D` (k-summed) | ✅ Already correct. |
| `wavepacket/wavepacket.hpp` | `inject_into_last_extra_state` writes one orbital | ❌ Only injects at kpoint 0. |
| `screens/leed_pattern_accumulator.hpp` | Reads `electrons.density()` | ✅ Should be correct. |
| `io/observables_writer.hpp` | Reads `StepContext` scalars | ✅ Energies and currents already k-summed by INQ. |
| `real_time/real_time_session.hpp` | Generic step callback | ✅ Generic. |

The `features/li-extensive` Phase 3 propagation works for energies and
density, but per-kpoint orbital information is not exposed by the current
inqkit, and there is no diagnostic helper for visualising
ψ_{n,k}(r) variation across the BZ.

The user has asked for two deliverables:

1. **Audit and harden inqkit for multi-kpoint use** so future runs (Li,
   any periodic crystal) work without surprises.
2. **A per-kpoint orbital visualiser** to build intuition for what
   ψ_{n,k}(r) looks like at different k-points and bands — a teaching tool,
   essential for interpreting the FFT energy spectra in Phase 5.

Both must be smoke-tested on systems simpler than 54-atom Li before they are
trusted on the production checkpoint.

## Smoke-test ladder (in order, lightest first)

| Tier | System | Cells / atoms | k-grid | Cost target | Validates |
|---|---|---|---|---|---|
| T0 | Atomic-H in cubic box | 1 H | Γ-only (1 kpoint) | < 30 s | Round-trip baseline (orbital writer, slice plotter) on the simplest possible case. |
| T1 | Atomic-H in cubic box | 1 H | 2×2×2 shifted MP | < 60 s | First multi-kpoint case. Single electron → analytical-ish ground state. |
| T2 | He in cubic box | 1 He | 2×2×2 shifted MP | < 60 s | Two-electron closed shell. |
| T3 | Li BCC primitive | 2 Li | 4×4×4 shifted MP | ~3 min | Metallic, dense BZ, but small cell. Already a tutorial (`Tutorial/li-bcc/li_bcc.cpp`). |
| T4 | Li BCC 3×3×3 supercell | 54 Li | 2×2×2 shifted MP | use existing GS checkpoint | Production-scale test. No new GS run needed. |

All smoke tests live in `Tutorial/inqkit-multikpoint-smoke/` (one subdir per
tier), follow the standard `inq-run` layout, and produce a small pass/fail
summary in their `results/run_summary.txt`.

Pass criteria for **each** tier:
- Density round-trip: ∫ρ d³r matches expected electron count to 1e-3.
- Orbital round-trip: ∫|ψ_{n,k}|² d³r ≈ 1 for every (n, k) we write
  (this is the explicit numerical check the user's `development-feedback-loop`
  rule mandates).
- Orbital symmetry: for k=0 (Γ), ψ_{n,0}(r) is real (or purely imaginary up
  to a global phase) — its imaginary-part L² norm should be < 1e-3.
- For k ≠ 0, ψ_{n,k}(r) shows non-trivial spatial phase modulation (visible
  in a real/imag side-by-side slice).

---

## Phase A — multi-kpoint audit and inqkit hardening

### A0. Audit (read-only first)

Walk every header under `inq-stack/include/inqkit/` and tag each public
function with one of: `kpoint-summed`, `requires-kpoint-index`,
`kpoint-broadcasts`, or `single-kpoint-assumed`. Output:
`docs/handovers/inqkit_multikpoint_audit.md` listing every function with its
classification and a one-line justification (the pattern in the source that
proves it).

Scope: `fields/`, `observables/`, `wavepacket/`, `screens/`, `io/`,
`real_time/`. Skip `core/` and `config/` (no k-point logic).

### A1. New helper — `fields/orbital_kpoint.hpp`

Add a small helper that returns ψ_{n,k}(r) **as a complex field** (not
|ψ|²) so the visualiser can show real and imaginary parts. Signature:

```cpp
namespace inqkit::fields {
  ComplexField3D orbital_complex(
      inq::systems::electrons const& electrons,
      int orbital_index,
      int kpoint_index);  // no default — caller must be explicit
}
```

The existing `density::orbital(electrons, n, kpoint_index=0)` keeps its
default for backward compat but gets a deprecation comment recommending
explicit `kpoint_index`. Both helpers loop the same FFT-shift pattern as
`density::total`.

Smoke test (T0): write ψ_{0,0} for atomic H, check it's a real Gaussian-like
1s orbital. Tier T1 (k=Γ shifted): same orbital should now be complex with
a smooth spatial phase.

### A2. Multi-kpoint orbital-overlap (`observables/orbital_overlap.hpp`)

The current implementation tracks a single orbital index against the entire
ground-state set at one k-point. Generalise to:

```cpp
class OrbitalOverlapMatrix {
  // Accept either a single (kpoint, orbital) pair, or "all kpoints".
  // For each kpoint independently, compute |<psi_GS_{n,k} | psi_t_{m,k}>|.
  // Cross-kpoint overlaps are zero by Bloch's theorem; we don't compute them.
};
```

Smoke test (T3): kick a Li bcc 2-atom cell, propagate 200 steps, verify the
diagonal of the per-kpoint overlap matrix is ≈ 1 in the first few steps and
deviates smoothly after.

### A3. Wavepacket injection (`wavepacket/wavepacket.hpp`)

`inject_into_last_extra_state` currently writes into one slot (state index
`n_states - 1`) at kpoint 0. We don't need this in the production Li runs
(no wavepacket), but for future multi-kpoint use the function should:
- Accept an optional `kpoint_index` argument with no default.
- Throw a clear error if called on a multi-kpoint electrons object without
  an explicit kpoint_index.

Smoke test: keep the existing single-kpoint test green; add a new T2 He test
that injects a wavepacket at kpoint 1 of a 2×2×2 MP grid and verifies the
norm via the already-available `injection_report`.

### A4. Documentation

Update `inq-stack/README.md` (if absent, create) with a short "k-point
support" section listing per-module behaviour and the explicit
`kpoint_index` requirement for orbital-level helpers.

---

## Phase B — per-kpoint orbital visualiser

### B1. C++ side: `inq-codebase/Li/run_save_gs_2x2x2_T200/dump_per_kpoint_orbitals.cpp`

VTI is the primary output. Folder layout is **band-major** so a user opening
ParaView on a single band directory immediately sees all 8 k-points side by
side as a series.

A small driver that:
1. Loads the Li 54-atom GS checkpoint.
2. Iterates over a user-selected set of band indices (default: {1, 40, 81, 100}
   — one deep, one mid-occupied, one Fermi-level, one empty), and for each
   band over all 8 MP k-points.
3. For each `(band, k)` pair, writes:
   - `Re ψ_{n,k}(r)` as a binary VTI (real-valued, ParaView-native)
   - `Im ψ_{n,k}(r)` as a binary VTI
   - `|ψ_{n,k}(r)|²` as a binary VTI
   plus appends a row to `orbital_index.csv`:
   `band | kpoint_index | kx | ky | kz | occ | evalue_ha | re_l2 | im_l2 | density_l2`.

Output (band-major — one folder per band, k as the time-series index):
```
results/analysis/ground_state/orbitals_per_kpoint/
├── orbital_index.csv
├── band_001/
│   ├── re_psi_k000.vti, re_psi_k001.vti, … re_psi_k007.vti
│   ├── im_psi_k000.vti, … im_psi_k007.vti
│   └── density_k000.vti, … density_k007.vti
├── band_040/
│   └── (same 24 VTIs)
├── band_081/
│   └── (same 24 VTIs)
└── band_100/
    └── (same 24 VTIs)
```
Total: 4 bands × 8 k-points × 3 components = 96 VTIs (~few MB each, ~300 MB).

Why band-major: opening `band_001/` in ParaView and animating "step" gives
you the k-sweep for that band — the natural comparison. Cross-band
comparison stays available by opening multiple folders. The CSV index lets
the Python helper join (band, k) ↔ eigenvalue ↔ kpoint-vector cleanly.

Implementation note: `inqkit::io::RealField3DWriter` is used for all three
components; the complex field is split into Re/Im before writing rather than
adding a new writer. This avoids depending on
`inq-stack/include/inqkit/io/complex_field_3d_writer.hpp` (which may not yet
support binary VTI). The new helper from A1, `orbital_complex(electrons, n,
k)`, is the source of the complex field.

Smoke-test driver versions for each tier (T0–T3) live in
`Tutorial/inqkit-multikpoint-smoke/<tier>/dump_orbital.cpp`, configured
inline. The smoke tests use a smaller band set (e.g. {0, 1} for H/He) but
the same band-major folder layout so the visualisation script is identical.

### B2. Python visualiser: `scripts/plot_per_kpoint_orbitals.py`

The VTI files in B1 are the primary scientific deliverable — open them in
ParaView for full 3D inspection. The Python helper produces *complementary*
2D summary sheets that are easier to put in a report than a 3D screenshot.

Reads `orbital_index.csv` plus the per-band VTI folders, produces (one
figure per band, plus a summary):

- `band_NNN_re.png` — for one band, all 8 k-points as a 2×4 grid of mid-cell
  slices of Re ψ. Fixed colour scale across the 8 panels (per band) so the
  Bloch-phase variation is visually comparable.
- `band_NNN_im.png` — same for Im ψ.
- `band_NNN_density.png` — same for |ψ|². Shows where ψ has support; should
  be more k-invariant than Re/Im, illustrating that the Bloch phase carries
  the k-dependence.
- `band_fragments_vs_kpoint.png` — single figure: evalue vs k-point index
  for each chosen band (toy band-structure-style summary).
- `parview_recipe.md` — short note explaining how to open one
  `band_NNN/re_psi_k*.vti` series in ParaView, set "step" as the time index,
  and animate to sweep across the BZ (this is the headline visualisation).

### B3. Smoke-test harness

`Tutorial/inqkit-multikpoint-smoke/run_all.sh` — sequential dispatcher:
1. Build + run T0 → assert pass.
2. T1 → T2 → T3 → T4 (T4 only if all earlier tiers green).
3. For each tier, exercise A1 (orbital_complex), A2 (overlap), B1+B2
   (visualiser).
4. A single `summary.md` is appended with per-tier pass/fail.

Failure handling: stop the harness on first tier that fails; print the
offending output to stderr; do not run subsequent tiers. The user sees a
clear failure point rather than a wall of confusing output from a downstream
test that depended on the broken upstream.

---

## Critical files

### To read

| File | Why |
|---|---|
| `/local/data/public/skcb2/tddft/inq-stack/include/inqkit/fields/density.hpp` | Reference k-sum + FFT-shift pattern (lines 43-102). |
| `/local/data/public/skcb2/tddft/inq-stack/include/inqkit/observables/orbital_overlap.hpp` | Current single-kpoint overlap; needs generalising. |
| `/local/data/public/skcb2/tddft/inq-stack/include/inqkit/wavepacket/wavepacket.hpp` | Injection logic. |
| `/local/data/public/skcb2/tddft/inq/src/systems/electrons.hpp` (lines 600-720) | `kpin()`, `kpoint_index()`, save/load semantics with k-points. |
| `/local/data/public/skcb2/tddft/Tutorial/li-bcc/li_bcc.cpp` | Reference primitive-cell Li for tier T3. |

### To create

| File | Phase |
|---|---|
| `inq-stack/include/inqkit/fields/orbital_kpoint.hpp` | A1 |
| Updates to `inq-stack/include/inqkit/observables/orbital_overlap.hpp` | A2 |
| Updates to `inq-stack/include/inqkit/wavepacket/wavepacket.hpp` | A3 |
| `inq-stack/README.md` (or section update) | A4 |
| `Tutorial/inqkit-multikpoint-smoke/T{0..3}/{dump_orbital,smoke}.cpp` | A1–B3 |
| `Tutorial/inqkit-multikpoint-smoke/run_all.sh` | B3 |
| `QuantumKickExtension/inq-codebase/Li/run_save_gs_2x2x2_T200/dump_per_kpoint_orbitals.cpp` | B1 |
| `QuantumKickExtension/inq-codebase/Li/scripts/plot_per_kpoint_orbitals.py` | B2 |
| `docs/handovers/inqkit_multikpoint_audit.md` | A0 (rolling) |
| `docs/handovers/inqkit_multikpoint_changes.md` | A1–B3 (rolling) |

---

## Verification

Per-phase, the **write → known-case-test → fix → confirm** loop:

| Step | Test | Pass criterion |
|---|---|---|
| A0 | Audit completeness | Every public function in scope appears in the audit doc with a tag and one-line justification. |
| A1, T0 | H atom Γ-only `orbital_complex` | ψ_{0,0} is real (Im L² < 1e-6). \|ψ\|² ≈ 1. |
| A1, T1 | H atom 2×2×2 MP | For k≠Γ, Im ψ has visible structure (not numerical noise). \|ψ\|² ≈ 1 at every (n,k). |
| A2, T3 | Li primitive bcc overlap | Diagonal ≈ 1 at t=0; off-diagonal < 1e-3. After 100 steps with no kick, diagonal still > 0.99 (energy-conserving evolution). |
| A3 | Wavepacket inject at k≠0 | T2 He inject at kpoint 1: norm_after ∈ [0.97, 1.03]. Single-kpoint case still passes. |
| B1, T0–T3 | VTI round-trip | For each (k,n), read back the VTI and recompute \|ψ\|² → matches density::orbital(electrons, n, k) to 1e-6. |
| B2 | Visual grid for T3 | Grid figure shows real-imag structure that is symmetric under k → −k (time-reversal); confirm by eye. |
| T4 | Production Li smoke | All A2/B1/B2 tests green on the existing `li_54_2x2x2_T200` checkpoint. No SCF re-run. |

End-to-end:

```bash
# in a worktree on features/inqkit-multikpoint
bash Tutorial/inqkit-multikpoint-smoke/run_all.sh
# expect "ALL GREEN" final line on stdout

# production check (uses existing GS checkpoint, no new SCF)
cd QuantumKickExtension/inq-codebase/Li/run_save_gs_2x2x2_T200
inq-run dump_per_kpoint_orbitals.cpp
python ../scripts/plot_per_kpoint_orbitals.py .
```

---

## Sequencing and isolation

- This work happens on `features/inqkit-multikpoint` in the main repo, in a
  separate **worktree**. The Phase 3 propagation on
  `features/li-extensive` (in QKE) keeps running undisturbed; it doesn't
  share code with the new branch.
- Build artefacts: each smoke-test tier has its own `inq-run` build dir, all
  gitignored.
- The Phase 3 propagation reaches step 1000 around 22:08 BST and finishes
  around 16:50 BST tomorrow. Phase A and B can complete well within that
  window.

## Risks and fallbacks

1. **`electrons.kpin()` distribution.** When inq is launched with k-point
   parallelism (multi-rank), `phi.set_part().contains(orbital_index)`
   filters out orbitals not on the local rank. Our visualiser must either
   (a) gather across ranks, or (b) fail loudly when run with > 1 rank. We'll
   start with (b) and document; (a) is a future extension.
2. **Complex-field VTI**. INQ's VTI writers in inqkit are real-valued
   (`RealField3DWriter`). We need `ComplexField3DWriter` —
   already exists at `inq-stack/include/inqkit/io/complex_field_3d_writer.hpp`
   (verified earlier in this session). If it's a stub, fall back to writing
   two real VTIs (re_psi.vti, im_psi.vti) per orbital.
3. **Per-orbital VTI count for the 54-atom run.** 8 k-points × 4 bands × 1
   complex-field-as-2-real-files = 64 VTIs ~ a few MB each → ~200 MB. Fine.
4. **GPU contention.** Phase 3 owns GPU 1. Smoke tests T0–T2 are cheap
   enough to run on CPU (`inq-run --cpu`). T3 (Li primitive bcc 4×4×4) is
   borderline; run on CPU if GPU 1 is busy. T4 reuses an existing
   checkpoint and only does I/O — no propagation, GPU not strictly needed.

## Out of scope

- Modifying INQ itself. All changes live in `inq-stack/include/inqkit/` and
  in user `.cpp` files.
- Multi-rank gather/scatter for the orbital visualiser (future work; flagged
  in Risks).
- Spin-polarised electrons (no use case in current research).
- Re-running the Li 54-atom GS at different k-grids (a separate question
  about k-grid convergence; can be added as a follow-up plan if needed).

---

## Approvals required before execution

1. Tier ladder T0 → T4 — the user may want a different order (e.g. start at
   T3 if they trust the simpler tests).
2. Number and choice of bands for the Li 54-atom visualiser (default:
   1, 40, 81, 100 — adjustable).
3. Whether to put the smoke tests under `Tutorial/` (the current proposal)
   or under a new `inq-stack/tests/` directory.
4. Confirm primary deliverable is the VTI series under `band_NNN/` (one
   folder per band, k as the series index in ParaView). PNG sheets are the
   secondary, paper-ready summary; cross-band comparison is left to the
   user opening multiple folders in ParaView.
