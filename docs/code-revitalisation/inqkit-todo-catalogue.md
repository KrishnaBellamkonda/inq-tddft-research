# inqkit TODO catalogue (review extraction)

Verbatim extraction of every `TODO` review comment in `inqkit`
(`inq-stack/include/inqkit/`), the backbone for the grilling that produces
**tests + a proposed architecture** (`docs/code-revitalisation/
inqkit-review-and-next-steps.md`). 30 substantive comments across 11 headers
(+ 5 "write this file" placeholder stubs, listed last).

**Type** key: `Q` question to answer via a test · `BUG?` suspected defect ·
`REFAC` restructuring · `EXPT` one-off experiment (visualise → decide) ·
`DOC` documentation/understanding · `NAME` naming/convention.
**Test class** (per the plan): `unit` · `integration` · `char` (characterization
/ golden-master baseline) · `oneoff` (answers a hypothesis, not in CI).

## Cross-cutting themes (recur across files)
- **Θ-coord** — real-space index ↔ physical-coordinate convention (FFT-natural
  vs −L/2-origin) and the `fft_shift` helper. → T01, T05, T20, overall #5, #12.
- **Θ-parallel** — GPU `gpu::run`/`reduce` + MPI `all_reduce` correctness,
  esp. **missing Allreduce in `plane_screen`**. → T17, T18, T22, T24, T28.
- **Θ-vector** — three scalars `*_x/_y/_z` should be a `vector3`-style unit.
  → T07, T09, overall #7.
- **Θ-naming** — `f→field`, `w→weight`, `m?`, the `function_` suffix
  convention. → T10, T13, T19, overall #1, #3.
- **Θ-density-semantics** — does INQ `electrons.density()` give total or
  target(system, WP-excluded)? → T02, T04.

## detail/ & fields/
| ID | File:line | Type | Verbatim essence | Implied test / change | Class |
|---|---|---|---|---|---|
| T01 | fields/density.hpp:34 | REFAC | `fft_shift`/index helper "seems like an important utility function" — move to a shared file. | Extract `fft_shift_index` to a `detail/` util; re-point density+orbital. | char→unit |
| T02 | fields/density.hpp:71 | Q/BUG? | **IMPORTANT:** confirm INQ `electrons.density()` returns *total* vs *target/system* (WP-excluded) density. | Validation test on a tiny system w/ + w/o WP occupation → compare integral. | integration |
| T03 | fields/density.hpp:313 | DOC/REFAC | `total_excluding_orbital` (rho_bath = total − occ·orbital): find where used, document, maybe remove. | Usage audit; decide keep/delete. | — |
| T04 | fields/orbital.hpp:14 | Q | Is the *full complex* wavefunction used when building observables (momentum dist. etc.)? | Test observables consume complex ψ, not |ψ|. | unit |
| T05 | fields/orbital.hpp:20 | REFAC | orbital.hpp imports density.hpp *only* for `fft_shift` → move it. | Same extraction as T01. | char→unit |

## io/
| ID | File:line | Type | Verbatim essence | Implied test / change | Class |
|---|---|---|---|---|---|
| T06 | io/vti_image_data_writer.hpp:56 | Q/BUG? | Check the array indexing convention — right coords map to right indices (fft-indexing ↔ iz-fastest). | VTI coord round-trip test (links overall #12). | char + integration |
| T07 | io/observables_writer.hpp:153 | REFAC | current & dipole must be tracked as a **vector unit**, not x/y/z. | Introduce `vector3`; refactor selection. | char→unit |
| T08 | io/real_field_3d_writer.hpp:367 | DOC | Meta sidecar: future expansion to full sim config (ions, electrons, GS/RT info). | Define richer meta schema (later). | — |

## observables/center_of_density.hpp
| ID | File:line | Type | Verbatim essence | Implied test / change | Class |
|---|---|---|---|---|---|
| T09 | :33 | Q/REFAC | Is `CenterOfDensityResult` (vector struct) used for COD in all sims? | Vector-type consolidation (Θ-vector). | — |
| T10 | :43 | NAME | Rename `f→field`, `w→weight`; what is `m` in `mx/my/mz` (moment?). | Rename pass; behaviour unchanged. | char |
| T11 | :416 | Q | What unit is the `_L` suffix? Shouldn't it be Bohr? | Pin unit; assert in test. | unit |
| T12 | :433 | DOC | Explain the `total_weight` condition; write a succinct comment. | Document; no behaviour change. | — |

## observables/density_delta.hpp
| ID | File:line | Type | Verbatim essence | Implied test / change | Class |
|---|---|---|---|---|---|
| T13 | :670 | NAME | Is `function_` suffix consistent across inq-stack C++/Python? | Convention audit → coding-standard skill. | — |
| T14 | :489 | EXPT | Verify first snapshot is t=0 vs t=dt; does t=dt remove the "deep hole behind"? Run both, visualise, decide. | One-off A/B run + visual. | oneoff |
| T15 | :554 | Q | Does δn at t+1 use t as base, or is t=0 the base for all steps? | Test base-frame semantics. | unit |
| T16 | :555 | EXPT | View step-by-step changes (δn within each individual timestep). | New observable variant (later). | oneoff |

## observables/eigenvalue_dump.hpp & wp_real_space_stats.hpp
| ID | File:line | Type | Verbatim essence | Implied test / change | Class |
|---|---|---|---|---|---|
| T17 | eigenvalue_dump.hpp:696 | DOC | Understand the data-structure types of `electrons` properties. | Document the INQ accessors used. | — |
| T18 | wp_real_space_stats.hpp:778 | Q | Is `write_every` passed from the real-time session? | Trace config plumbing; test. | unit |
| T19 | wp_real_space_stats.hpp:779 | NAME | Avoid `w`; use ix/iy/iz, x/y/z, wx/wy/wz — worth fixing throughout? | Naming convention decision (Θ-naming). | char |
| T20 | wp_real_space_stats.hpp:892 | DOC | "What does this syntax really do?" (destructor closing file). | Explain RAII destructor. | — |
| T21 | wp_real_space_stats.hpp:1032 | Q/BUG? | "Need to test that the parallelisation is working as expected." (state-partition local check) | Parallelism test case (Θ-parallel). | integration |
| T22 | wp_real_space_stats.hpp:1085 | Q | MPI reduction — is this GPU comms? How are GPU+MPI compatible? Is GPU parallelised? | Understand + test the 7-partial-sum reduce. | integration |

## screens/plane_screen.hpp
| ID | File:line | Type | Verbatim essence | Implied test / change | Class |
|---|---|---|---|---|---|
| T23 | :4 | REFAC | Generalise the class beyond the z-axis to all dims. | API generalisation + tests. | char→unit |
| T24 | :5 | EXPT | Implement time-averaged total density. | New feature (later). | — |
| T25 | :1203 | **BUG?** | **No Allreduce on the slice** after the loop — multi-rank state-parallel runs accumulate only local states. Add Allreduce if multi-rank needed. | Multi-rank correctness test (Θ-parallel). | integration |

## wavepacket/wavepacket.hpp
| ID | File:line | Type | Verbatim essence | Implied test / change | Class |
|---|---|---|---|---|---|
| T26 | :1276 | EXPT/REFAC | Make a **momentum-space** Gram-Schmidt; compare observables vs real-space; decide. (links overall #10) | Extract GS (overall #2) + A/B compare. | oneoff |
| T27 | :1349 | Q/BUG? | Check the parallelisation code; write tests proving `gpu::run`+`reduce` understanding; same for MPI. | Parallelism proof tests (Θ-parallel). | integration |
| T28 | :120 | Q | Unknown units of `k` (guess 1/Bohr) — sanity check via free-WP / jellium-WP propagation. | k-unit test: inject k₀ → momentum peak at k₀. | integration |
| T29 | :1396 | Q/BUG? | Test orthogonalisation rigorously. | Post-ortho overlaps ≈ 0; norm ≈ 1. | integration |
| T30 | :1419 | REFAC | Generalise injection to any k-point config (currently gamma-only throws). | Multi-kpoint support (later). | — |
| T31 | :1528 | Q/DOC | Limitations of the GS protocol? Should overlap emerge during subtraction (KS orbitals orthonormal → no)? | Test residual overlaps during GS. | unit/integration |

> (IDs run T01–T31 with T17/T18 sharing a row group; 30 distinct comments —
> the count differs from a naive line count because two `wavepacket` blocks and
> two `wp_real_space_stats` blocks each carry stacked TODOs.)

## Placeholder stubs (`// TODO: Write this file`) — out of scope, restructure phase
`detail/validation.hpp`, `detail/filesystem.hpp`, `detail/text_io.hpp`,
`io/manifest_writer.hpp`, `io/text_summary_writer.hpp`. See
`docs/notes/inqkit-rejuvenation-ideas.md`.

## The grid_layout coordinate test (overall comment #12 — flagship)
Run a tiny sim with an ion at an **off-centre, arbitrary** position → use the
`grid_layout` path to get its 3D position field → take xy/xz/yz slices → assert
the peak sits at the specified position (checks the FFT-shift mapping). Then
confirm the **VTI** output read in Python places the same coordinate correctly.
This single test exercises Θ-coord end-to-end (T01, T05, T06) across C++→VTI→
Python.
