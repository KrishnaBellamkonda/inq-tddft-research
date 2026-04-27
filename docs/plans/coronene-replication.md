# Plan — Coronene TDDFT replication framework



## 1. Context

A large number of coronene wave-packet TDDFT runs were scattered across the workspace
during the effort to reproduce Tsubonoya, Hu & Watanabe, *Phys. Rev. B* **90**, 035416 (2014).
Many of the older runs used a buggy `coronene.xyz` whose carbon plane was placed at
`z = L/2` rather than `z = 0` — they implicitly assumed a `[0, L]` cell, but INQ uses
`[-L/2, +L/2]` (cf. `inq/src/systems/cell.hpp:212`, `inqkit/config/tsubonoya_2014_coronene.hpp:15-17`).

This plan does three things:

1. Catalogues every coronene TDDFT configuration ever attempted (Tables A and B), so the
   user can pick which to replicate with the corrected geometry.
2. Designs a single, shared infrastructure under
   `ResearchProject/systems/coronene/shared/` — one canonical `coronene.xyz`
   (already at z = 0), a reusable C++ baseline config, a corrected `run.cpp` template,
   and one centralised Python post-processor that consumes a `results/` directory.
3. Defines the mapping from each replicated run to a `results/` tree that complies with
   `docs/results_folder_structure_spec.md` and `docs/visualisation-instructions-v1.md`.

### 1.1 Units and folder-structure precedence

- **All lengths and filenames in this plan are in Bohr.** INQ stores cells and
  positions in Bohr internally (`inq-stack/include/inqkit/config/tsubonoya_2014_coronene.hpp:32-34`
  defines `LX_BOHR = 18.4 * ANG_TO_BOHR`; the `_b` literal in
  `run_propagate_paper_replica/run.cpp:78` is INQ's Bohr unit). Earlier drafts
  showed Å values pulled from each run's source comments; those have been
  converted in Tables A and B below. Conversion factor: 1 Å = 1.8897259886 Bohr.
- **The folder structure laid out in §4 supersedes
  `.claude/rules/file-placement.md`.** The latter is auxiliary / outdated; this
  plan is authoritative for coronene replication work and will be mirrored to
  `docs/plans/coronene-replication.md`.

### 1.2 Tsubonoya 2014 base configuration (the reference for everything)

| Quantity | Value (Bohr / a.u. / eV) | Comes from |
|---|---|---|
| Cell L<sub>x</sub> × L<sub>y</sub> × L<sub>z</sub> | **35 × 35 × 60 Bohr** (precise: 34.77 × 34.77 × 59.90) | paper, via `tsubonoya_2014_coronene.hpp:32-34` |
| Cutoff | **40 Ha** | INQ-converged value (locked in §7) |
| Extra states | 8 | paper / `tsubonoya_2014_coronene.hpp:39` |
| WP σ | **1.0 Bohr** (precise: 1.0015 = 0.53 Å) | paper Eq. 1, `tsubonoya_2014_coronene.hpp:49` |
| WP offset b | **12 Bohr** (precise: 12.00 = 6.35 Å) | paper Eq. 1, `tsubonoya_2014_coronene.hpp:50` |
| WP energy E | **200 eV** | paper, `tsubonoya_2014_coronene.hpp:51` |
| WP direction k<sub>0</sub> | along −z | paper, `tsubonoya_2014_coronene.hpp:63` |
| dt | 0.020 a.u. | paper Eq. 5, `tsubonoya_2014_coronene.hpp:68` |
| N<sub>steps</sub> | 600 | `tsubonoya_2014_coronene.hpp:72` |
| LEED window [t<sub>1</sub>, t<sub>2</sub>] | [0.077, 0.25] fs = [3.18, 10.34] a.u. | paper Eq. 5, `tsubonoya_2014_coronene.hpp:79-82` |
| Screens | 20, in [−L<sub>z</sub>/2, +L<sub>z</sub>/2] | `tsubonoya_2014_coronene.hpp:75` |

Throughout the rest of this plan, "the base run" or "Tsubonoya base" means this
exact configuration. Run names omit any quantity that matches it (see §4.4).

---

## 2. Table A — Parameter inventory (all known coronene TDDFT configurations)

**All lengths in Bohr**, rounded to 2 sig figs in this table for readability;
internal storage in INQ uses higher precision. Geometry tag:
`✓ centred` = molecule at z = 0 (correct for INQ);
`✗ z=L/2` = old-bug shifted geometry (atoms placed at +L<sub>z</sub>/2 assuming `[0, L]`);
`GS-only` = no real-time, no wave packet.

| # | Run path | Type | Cell L<sub>x</sub>×L<sub>y</sub>×L<sub>z</sub> (Bohr) | Cutoff (Ha) | Extra states | xyz file | Geom. tag | XC | dt (a.u.) | N<sub>steps</sub> | t<sub>tot</sub> (fs) | WP σ (Bohr) | WP offset b (Bohr) | E<sub>kin</sub> (eV) | k<sub>0</sub> dir | Screens | Snap every | Window [t<sub>1</sub>,t<sub>2</sub>] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | RP/01_geometry | GS-only | 35×35×60 | 30 | 3 | 01_geometry/coronene.xyz | ✓ centred | LDA | — | — | — | — | — | — | — | — | — | — |
| 2 | RP/02_ground_state_analysis | GS-only | 35×35×60 | 30 | 3 | 01_geometry/coronene.xyz | ✓ centred | LDA | — | — | — | — | — | — | — | — | — | — |
| 3 | RP/03_ecut_convergence | GS sweep | 35×35×60 | **20–60 (9 pts)** | 3 | 01_geometry/coronene.xyz | ✓ centred | LDA | — | — | — | — | — | — | — | — | — | — |
| 4 | Tut/run_01_tight_scf | GS-only | 35×35×60 | ≈54 | 8 | coronene_leed.xyz | ✗ z=L/2 | PBE | — | — | — | — | — | — | — | — | — | — |
| 5 | Tut/run_02_benzene | GS-only | 35×35×60 | ≈54 | 8 | benzene.xyz | ✗ z=L/2 | PBE | — | — | — | — | — | — | — | — | — | — |
| 6 | Tut/run_03_coronene_hardcoded | GS-only | 35×35×60 | ≈54 | 8 | hard-coded coords | ✗ z=L/2 | PBE | — | — | — | — | — | — | — | — | — | — |
| 7 | Tut/run_04_graphene | GS-only | 35×35×60 | ≈54 | 12 | graphene_nanoflake.xyz | ✗ z=L/2 | PBE | — | — | — | — | — | — | — | — | — | — |
| 8 | Tut/run_05_quarter_coronene | GS-only | 35×35×60 | ≈54 | 8 | quarter_coronene.xyz | ✗ z=L/2 | PBE | — | — | — | — | — | — | — | — | — | — |
| 9 | Tut/run_06_centred_writer_check | GS+writer probe | 35×35×60 | 54 | 8 | coronene_centred.xyz | ✓ centred | PBE | — | — | — | — | — | — | — | — | — | — |
| 10 | Tut/run_07_paper_replica | **Full RT-LEED** | 35×35×60 | 54 | 8 | coronene_centred.xyz | ✓ centred | ALDA | 0.020 | 600 | 0.290 | 1.0 | 12 | 200 | −z | 20 | 30 | [0.077, 0.25] fs |
| 11 | Tut/run_08_gs_only_wp_check | GS+WP norm | 35×35×60 | 54 | 8 | coronene_centred.xyz | ✓ centred | ALDA | — | — | — | 1.0 | 12 | 200 | −z | — | — | — |
| 12 | Tut/run_09_gs_vti_writer | GS VTI test | 35×35×60 | 54 | 8 | coronene_centred.xyz | ✓ centred | ALDA | — | — | — | — | — | — | — | — | — | — |
| 13 | RP/04_leed_simulation/run_001 | Full RT-LEED | 35×35×60 | 40 | 3 | coronene.xyz (uncentred) | ✗ uncentred | LDA | 0.020 | ≈517 | ≈0.250 | 1.0 | 12 | 200 | −z | n/a | n/a | n/a |
| 14 | RP/04_leed_simulation/run_002 | Full RT-LEED | 35×35×60 | 54 | 3 | coronene_centered.xyz | ✗ z=L/2 | LDA | 0.020 | 516 | 0.250 | **2.6** | 12 | 200 | −z | n/a | n/a | n/a |
| 15 | RP/04_leed_simulation/run_003 | Full RT-LEED (validated) | 35×35×60 | 40 | 3 | coronene_centered.xyz | ✗ z=L/2 | LDA | 0.020 | 516 | 0.250 | **2.6** | 12 | 200 | −z | n/a | n/a | n/a |
| 16 | RP/04_leed_simulation/run_004 | Full RT-LEED | 35×35×**90** | 40 | 3 | coronene_centered.xyz | ✗ z=L/2 | LDA | 0.020 | 1561 | 0.756 | 1.0 | **40** | 200 | −z | n/a | n/a | n/a |
| 17 | RP/04_leed_simulation/run_005 | Setup only | 35×35×60 | n/a | n/a | coronene_centered.xyz | ✗ z=L/2 | n/a | n/a | n/a | n/a | 1.0 | n/a | 200 | −z | n/a | n/a | n/a |
| 18 | RP/coronene-wp-rt/run_01_d635_base | Full RT-LEED | 35×35×90 | 40 | 3 | coronene_centered.xyz | ✗ z=L/2 | LDA | 0.020 | 10000 | 4.84 | 1.0 | 12 | 200 | −z | 4 | n/a | n/a |
| 19 | RP/coronene-wp-rt/run_02_d3 | Full RT-LEED | 35×35×90 | 40 | 3 | coronene_centered.xyz | ✗ z=L/2 | LDA | 0.020 | 10000 | 4.84 | 1.0 | **5.7** | 200 | −z | 4 | n/a | n/a |
| 20 | RP/coronene-wp-rt/run_03_d10 | Full RT-LEED | 35×35×90 | 40 | 3 | coronene_centered.xyz | ✗ z=L/2 | LDA | 0.020 | 10000 | 4.84 | 1.0 | **19** | 200 | −z | 4 | n/a | n/a |
| 21 | RP/coronene-wp-rt/run_04_d15 | Full RT-LEED | 35×35×90 | 40 | 3 | coronene_centered.xyz | ✗ z=L/2 | LDA | 0.020 | 10000 | 4.84 | 1.0 | **28** | 200 | −z | 4 | n/a | n/a |
| 22 | RP/coronene-wp-rt/run_05_d20 | Full RT-LEED | 35×35×90 | 40 | 3 | coronene_centered.xyz | ✗ z=L/2 | LDA | 0.020 | 10000 | 4.84 | 1.0 | **38** | 200 | −z | 4 | n/a | n/a |
| 23 | RP/coronene-wp-rt/run_06_projectile | Full RT-LEED | 35×35×90 | 40 | 3 | coronene_centered.xyz | ✗ z=L/2 | LDA | 0.020 | 10000 | 4.84 | **0.50** | 12 | **800** | −z | 4 | n/a | n/a |
| 24 | RP/run_save_gs_paper_replica | GS save (checkpoint) | 35×35×60 | 40 | 3 | configurations/coronene_centred.xyz | ✓ centred | LDA | — | — | — | — | — | — | — | — | — | — |
| 25 | RP/run_propagate_paper_replica | RT load + propagate | 35×35×60 | 40 | 3 | (from checkpoint) | ✓ centred | LDA | 0.020 | 600 | 0.290 | 1.0 | 12 | 200 | −z | 20 | 30 | [0.077, 0.25] fs |

`RP` = `ResearchProject/systems/coronene`; `Tut` = `Tutorial/coronene-leed/run_diagnoses`.
Numbers shown for runs 1–12 come from `inqkit/config/tsubonoya_2014_coronene.hpp` and the per-run `run.cpp`. Numbers for runs 13–25 come from each run's local `run.cpp`/`config.hpp`/`run_summary.txt` (read by Explore agents — please flag if any cell-size or cutoff disagrees with what is on disk so I can re-verify).

### Geometry-bug summary

- **Correct (atoms at z = 0)**: rows 1–3, 9–12, 24, 25.
- **Buggy (atoms at z = L/2)**: rows 4–8 (Tutorial diagnostics 01–05) and rows 13–23 (every `04_leed_simulation` and `coronene-wp-rt` run that produced LEED data).

The buggy runs are **the entire set of completed full RT-LEED scattering simulations** with the exception of `run_07_paper_replica`. That is exactly why the user wants to re-run.

---

## 3. Table B — Hypothesis / motivation per run

| # | Run | Hypothesis being tested |
|---|---|---|
| 1 | RP/01_geometry | Verify the bare geometry: pseudopotentials load, single-point forces are sane. Used as the source-of-truth xyz for the early runs. |
| 2 | RP/02_ground_state_analysis | Baseline GS energy / forces / KS eigenvalues at the paper geometry, against which all later runs are compared. |
| 3 | RP/03_ecut_convergence | E<sub>cut</sub> sweep 20–60 Ha to find the energy-minimising cutoff for INQ + pseudo-dojo norm-conserving PSPs (settled on 40 Ha). |
| 4 | Tut/run_01_tight_scf | Does tightening SCF to 1e-6 Ha (vs 1e-4) restore D6h symmetry of nearly-degenerate orbitals near HOMO? |
| 5 | Tut/run_02_benzene | Benzene reference: D6h with only 15 occupied → if the orbital-cross artefact persists, it isn't degeneracy-specific. |
| 6 | Tut/run_03_coronene_hardcoded | Bypass the xyz parser: hard-code coronene coords. If the artefact persists → parser is innocent. |
| 7 | Tut/run_04_graphene | Honeycomb fragment without H. If the cross artefact appears → it's not aromaticity-specific. |
| 8 | Tut/run_05_quarter_coronene | Symmetry-breaking 7C+3H fragment with odd electrons. If cross disappears → it was a symmetry artefact, not a code bug. |
| 9 | Tut/run_06_centred_writer_check | Diagnose the four-corner artefact: is it the writer (FFT-natural index 0 ↔ origin) plus a centred geometry, vs the geometry itself? Confirmed: writer-side issue when `origin = -L/2` is set explicitly. |
| 10 | Tut/run_07_paper_replica | **Closest existing thing to a full Tsubonoya 2014 reproduction**: ALDA, σ = 1.0 Bohr WP, b = 12 Bohr, 600 steps, 20 screens, paper window. |
| 11 | Tut/run_08_gs_only_wp_check | Quantify how Modified-Gram-Schmidt + renormalisation of the injected WP deforms the Gaussian (norm before/after). |
| 12 | Tut/run_09_gs_vti_writer | Validate the new C++-native VTI writer (binary + ASCII) against the previous Python `.raw → .vti` converter. |
| 13 | RP/04_leed_simulation/run_001 | First end-to-end paper attempt; **buggy** uncentred geometry + paper-narrow WP — produced a four-copy LEED ghost. |
| 14 | RP/04_leed_simulation/run_002 | Geometry recentred to `z = L/2` (bug, but better than uncentred) + WP widened to σ = 2.6 Bohr + cutoff bumped to paper's 54 Ha → SCF instability. |
| 15 | RP/04_leed_simulation/run_003 | Same as 14 but cutoff reverted to 40 Ha (the converged minimum). Marked "validated" in old handover; gives a 6-fold LEED. |
| 16 | RP/04_leed_simulation/run_004 | Long-Z cell (L<sub>z</sub> = 90 Bohr) + paper WP narrowness + maximised offset (b = 40 Bohr) → far-field LEED + background-subtracted accumulator. |
| 17 | RP/04_leed_simulation/run_005 | Variant template; never executed. |
| 18 | RP/coronene-wp-rt/run_01_d635_base | Parameter-sweep baseline: paper b = 12 Bohr on the long-Z cell, very long propagation (10 000 steps). |
| 19 | RP/coronene-wp-rt/run_02_d3 | Close approach (b = 5.7 Bohr): does the WP visibly hybridise with the molecule? |
| 20 | RP/coronene-wp-rt/run_03_d10 | Intermediate b = 19 Bohr for trend interpolation. |
| 21 | RP/coronene-wp-rt/run_04_d15 | Intermediate b = 28 Bohr. |
| 22 | RP/coronene-wp-rt/run_05_d20 | Far-field b = 38 Bohr: pure diffraction regime. |
| 23 | RP/coronene-wp-rt/run_06_projectile | High-energy (800 eV) + narrow WP (σ = 0.50 Bohr) → high-resolution scattering test. |
| 24 | RP/run_save_gs_paper_replica | Compute the paper-replica GS once, save to disk (decouples GS cost from many propagations). |
| 25 | RP/run_propagate_paper_replica | Load checkpoint + inject WP + propagate. **Cleanest existing template** — but its `results/` layout does not yet match the new spec. |

The user reviews Tables A and B and tells me which numbered rows to replicate.

---

## 4. Implementation design — shared/ folder

```
ResearchProject/systems/coronene/
├── shared/
│   ├── geometry/
│   │   └── coronene.xyz                  # canonical, z = 0 centred
│   ├── configs/
│   │   ├── tsubonoya_2014_base.hpp       # paper baseline (cell, dt, N_steps, σ, b, k0…)
│   │   └── (per-variant headers as needed; see §4.2)
│   ├── cpp/
│   │   ├── run_template.hpp              # helpers used by every run.cpp
│   │   ├── results_paths.hpp             # compile-time results subpaths per the spec
│   │   └── leed_screen_layout.hpp        # screen z-positions (parametrised by Lz)
│   └── python/
│       ├── postprocess.py                # CLI entry point (argparse)
│       └── coronene_pipeline/            # internal modules; thin wrapper over inqview
│           ├── __init__.py
│           ├── ground_state.py
│           ├── observables.py
│           ├── density.py
│           ├── screens.py
│           ├── overlap.py
│           ├── orbitals.py
│           ├── vti.py                    # only if inqview's VTI helpers fall short
│           └── run_summary.py            # builds run_summary.txt section by section
```

### 4.1 `shared/geometry/coronene.xyz`

Already verified (read above): `01_geometry/coronene.xyz` is centred — atoms symmetric about
x = y = 0, all z = 0. **Action**: copy this file verbatim to `shared/geometry/coronene.xyz`,
keep the existing one as the legacy source, and make every new `run.cpp` point at the shared
copy. No re-centring needed.

> Reusability check (Section B of the audit): `inq/src/parse/xyz.hpp:57` reads coordinates verbatim, `inq/src/systems/cell.hpp` builds an orthorhombic cell `[-L/2, +L/2]`. A single z = 0 xyz therefore works for any orthorhombic cell whose half-extent is large enough to contain the molecule **and** the WP centre at b = 6.35 Å (so L<sub>z</sub> ≥ 2(b + a few σ) ≈ 25–30 Å minimum).

### 4.2 `shared/configs/`

One baseline header (`tsubonoya_2014_base.hpp`) — essentially the existing
`inqkit/config/tsubonoya_2014_coronene.hpp`, but **moved into the coronene
tree** and **fixed at cutoff 40 Ha** (locked in §7) so the runs are
self-contained:

- Cell: 35 × 35 × 60 Bohr (Tsubonoya base; precise 34.77 × 34.77 × 59.90).
- DFT: ALDA, cutoff 40 Ha, extra_states 8.
- WP: σ = 1.0 Bohr, b = 12 Bohr, E = 200 eV, k<sub>0</sub> along −z.
- RT: dt = 0.020 a.u., N<sub>steps</sub> = 600, write_every = 10, screen_snap_every = 30.
- LEED: 20 screens, window [0.077, 0.25] fs.

Per-variant overrides (one tiny header per actual run; named to mirror the
folder name from §4.4):

| Variant header | Differs from base in |
|---|---|
| `tsubonoya_2014_base.hpp` | (the baseline) |
| `E30.hpp` | E = 30 eV |
| `E800.hpp` | E = 800 eV |
| `s0p33.hpp` | σ = 0.33 Bohr (= base/3) |
| `s3.hpp` | σ = 3.0 Bohr (= 3 × base) |
| `E800_s0p33.hpp` | E = 800 eV, σ = 0.33 Bohr |
| `E30_s3.hpp` | E = 30 eV, σ = 3.0 Bohr |
| `b18_35x35x80.hpp` | L<sub>z</sub> = 80 Bohr, b = 18 Bohr |
| `b6_35x35x80.hpp` | L<sub>z</sub> = 80 Bohr, b = 6 Bohr |
| `35x35x40.hpp` | L<sub>z</sub> = 40 Bohr (= 2/3 × base) |

Only the headers actually used by the §10 investigations are created. Any other
variant header is added on demand.

### 4.3 `shared/cpp/results_paths.hpp`

A header that fixes the per-spec paths once, so each `run.cpp` just calls
`results::density_rt_total_dir()`, `results::vti_density_rt_total_dir()`, etc.
This is what makes "every run produces the same results tree".

### 4.4 Each run lives in its own folder + run-naming convention

```
ResearchProject/systems/coronene/<run_name>/         # flat sibling under coronene/
├── run.cpp                              # the propagation
├── paths.hpp                            # paths to GS checkpoint, results dir
├── README.md                            # 5-line description: hypothesis, parameter delta vs base
└── results/                             # populated at runtime; structure per spec
```

(Cross-run analysis goes elsewhere — see §10's `hypotheses/` folder.)

**Naming convention** (per user instruction):

- Prefix: `run_`.
- Symbols, in canonical order: `b<value>` (WP offset, Bohr), `s<value>` (WP σ, Bohr), `E<value>` (WP energy, eV), `<Lx>x<Ly>x<Lz>` (cell dims, Bohr, 2 sig figs).
- Each token separated by `_`.
- **Lowercase `p`** = decimal point in a numerical value (e.g. σ = 0.33 Bohr → `s0p33`, σ = 1.5 Bohr → `s1p5`).
- **Uppercase `P`** = any literal period (punctuation) needed in a filename. None are needed in this round, but the convention is reserved for future use.
- **Only quantities that differ from the Tsubonoya base config (§1.2) appear in the name.** A run that changes only σ omits b, E and the cell.
- The base config itself is `run_base/` (no parameters in the name — they all match the base).

Examples derived from the base values (b<sub>0</sub> = 12, σ<sub>0</sub> = 1.0,
E<sub>0</sub> = 200, cell<sub>0</sub> = 35×35×60, all Bohr / eV):

| Folder name | Means |
|---|---|
| `run_base` | The Tsubonoya base run. |
| `run_E30` | E = 30 eV; everything else at base. |
| `run_E800` | E = 800 eV. |
| `run_s0p33` | σ = 0.33 Bohr. |
| `run_s3` | σ = 3.0 Bohr. |
| `run_E800_s0p33` | E = 800 eV and σ = 0.33 Bohr (fast projectile). |
| `run_E30_s3` | E = 30 eV and σ = 3.0 Bohr (electron capture). |
| `run_b18_35x35x80` | b = 18 Bohr, L<sub>z</sub> = 80 Bohr. |
| `run_b6_35x35x80` | b = 6 Bohr, L<sub>z</sub> = 80 Bohr. |
| `run_35x35x40` | L<sub>z</sub> = 40 Bohr (2/3 × base). |

### 4.5 Ground-state checkpoints

A separate subtree for save runs:

```
ResearchProject/systems/coronene/save_gs/<gs_signature>/run.cpp
ResearchProject/systems/coronene/checkpoints/<gs_signature>/   # produced by save_gs run
```

`<gs_signature>` encodes `(cell, cutoff)` for the coronene molecule
(atoms + extra_states are constant across all selected runs). Three checkpoints
are needed for the §10 investigations:

| GS signature | Used by |
|---|---|
| `gs_35x35x60_cut40/` | base + every run that keeps the base cell (E30, E800, s0p33, s3, E800_s0p33, E30_s3) |
| `gs_35x35x80_cut40/` | b18_35x35x80, b6_35x35x80 |
| `gs_35x35x40_cut40/` | 35x35x40 |

Three GS save runs total, each saving once and reused by every propagation that
matches the cell.

**Per the user's instruction**, every `run.cpp` pastes the corresponding GS construction code
in a top-of-file `/* ... */` block as documentation:

```cpp
/* Ground state used (compiled and saved by save_gs/<gs_signature>/run.cpp):
 *
 *   auto cell = systems::cell::orthorhombic(LX_BOHR*1.0_b, LY_BOHR*1.0_b, LZ_BOHR*1.0_b).finite();
 *   auto ions = systems::ions::parse(SHARED_GEOMETRY_XYZ, cell);
 *   auto electrons = systems::electrons(ions,
 *       options::electrons{}.cutoff(54.0_Ha).extra_states(8));
 *   ground_state::calculate(ions, electrons,
 *       options::theory{}.lda(),
 *       inq::options::ground_state{}.energy_tolerance(1e-6_Ha).max_steps(1000));
 *   electrons.save("checkpoints/<gs_signature>/");
 */
```

### 4.6 Corrections to the `run_propagate_paper_replica/run.cpp` template

The user flagged this file as the closest reference. To bring it in line with
`docs/results_folder_structure_spec.md`, the new template differs in:

1. **All output paths prefixed with `results/raw/...` or `results/analysis/...`** (writers currently default to `results/density_rt_target/`).
2. **Three density categories**, not two:
   - `density_rt_system` ← `density::total(electrons)` (occupied orbitals only — INQ's `density::total` excludes the WP extra state, per `docs/observables_reference.md:27`).
   - `density_rt_wp` ← `density::orbital(electrons, wp_idx)`.
   - `density_rt_total` ← `system + wp` (computed at write-time by adding the two `RealField3D`s).
3. **Flat instantaneous-screen filenames** (`screen_NN_tXXXXXX.dat`), not nested per-step folders.
4. **Wave-packet artefacts**: emit `results/raw/wavepacket/wavepacket_config.txt`, `injection_report.txt`, `orthogonality_report.csv`, `density_wp_initial/`, and `wavefunction_wp_initial/`.
5. **Ground-state artefacts**: also write the GS density and the per-orbital GS densities into `results/raw/ground_state/density_system/`, `density_gs_orbitals/`, plus `results/raw/vti/density_gs_system/` and `density_gs_orbitals/` (flat filenames).
6. **Time-windowed screens**: in addition to the single paper window `[T1, T2]`, emit a small set of contiguous windows (e.g. five equal windows over the WP-transit time), each saved as `screen_NN_tAAAAAA_to_tBBBBBB.dat`. The paper-window result is one of those, kept under a documented subfolder.
7. **`results/run_summary.txt`** filled with the full template from `docs/results_folder_structure_spec.md` §3.2, written at the end of every run **and** also written stub-style at the start so a crashed run still has a partial summary.
8. **VTI directly**: every density writer is configured `emit_raw=false, emit_vti=true, vti_format=binary` and the path is `results/raw/vti/<series>/`. Confirmed in §C of the audit that `RealField3DWriter` skips the `.raw` write cleanly when `emit_raw=false` and produces flat sibling filenames.
9. **Overlap optimisation (per user instruction)**: do **not** compute the full
   `O_ij(t)` matrix between every evolved KS orbital and every GS KS orbital.
   Instead, at each timestep compute only the row corresponding to the WP state:
   `O_i,wp(t) = |⟨ψᵢ<sup>GS</sup> | ψ<sub>wp</sub>(t)⟩|²` for `i ∈ [0, n_occupied)`.
   This costs `O(n_occupied)` inner products per step instead of
   `O(n_states × n_occupied)`. A new helper
   `inqkit::observables::WPOverlapVector` (a stripped-down version of
   `OrbitalOverlapMatrix`) is added under
   `inq-stack/include/inqkit/observables/` and called from the RT step lambda
   in place of the full matrix snapshot. Output goes to
   `results/raw/overlap/wp_overlap_with_gs_orbitals.csv` (one row per saved step,
   columns = GS orbital index 0…n_occupied−1) — exactly the source the spec's
   `wp_overlap_with_gs_orbitals.gif` consumes.

### 4.7 Centralised Python post-processor (modularised)

Per the user's decision: **generalisable core logic lives in inqview**; the
coronene-specific entry point lives at
`ResearchProject/systems/coronene/scripts/coronene_postprocess.py` (note: under
`systems/coronene/scripts/`, **not** `shared/python/`).

Core (generalisable) modules added to `inq-stack/python/inqview/postprocess/`:

```
inqview/postprocess/
├── __init__.py
├── pipeline.py             # phase dispatch, --rebuild, --skip, --only
├── run_summary.py          # writes results/run_summary.txt per the spec template
├── ground_state.py         # SCF history, eigenvalues, occupations, GS density figs
├── observables.py          # time-domain plots + FFT (energy, current, dipole)
├── density.py              # 2D slice GIFs (xy/xz/yz), 3D ParaView volume GIFs
├── screens.py              # total / instantaneous / time-windowed LEED + coordinate checks
├── overlap.py              # overlap matrix heatmap GIF + WP-overlap bar GIF
└── orbitals.py             # orbital galleries + WP orbital density GIF
```

These wrap the existing `inqview` primitives (`fields`, `vti`, `paraview`,
`overlap`, `fourier`, `screens`, `plots`) and codify the visualisation rules
from `docs/visualisation-instructions-v1.md` (3 sig figs, fixed colour scales,
`step k/N` + `t = X.XX fs` in every frame title, run-name in titles, axis labels
with units).

The system-specific wrapper is a thin CLI:

```
ResearchProject/systems/coronene/scripts/coronene_postprocess.py
    --results <results_dir>
    [--phases gs,observables,density,screens,overlap,orbitals,paraview,summary]
    [--rebuild] [--skip-paraview]
```

It supplies coronene-specific defaults (run-name parsing, colour-scale percentiles
appropriate for diffraction-pattern intensities, screen-z layout for the
coordinate-check helper) and otherwise delegates to
`inqview.postprocess.pipeline.run(results_dir, ...)`. Any other system gets its
own thin wrapper in `systems/<name>/scripts/`.

Outputs produced (per the spec):

- `results/run_summary.txt` (full)
- `results/analysis/ground_state/{scf_convergence,eigenvalue_spectrum,occupations,density_gs_system,gs_orbital_gallery,ground_state_summary}.png`
- `results/analysis/observables/*.png` — all the time-domain + FFT plots listed in spec §15
- `results/analysis/density/{total,system,wp}_{xy,xz,yz}.gif` + volume GIFs (via ParaView)
- `results/analysis/screens/{total,instantaneous,time_windowed,filtered,spectra,coordinate_checks}/...`
- `results/analysis/overlap/wp_overlap_with_gs_orbitals.gif` (animated bar chart per `docs/visualisation-instructions-v1.md` §5; the only overlap visualisation, since the full O<sub>ij</sub> matrix is no longer computed — see §4.6 item 9)
- `results/analysis/orbitals/{homo_lumo_density.png, selected_orbital_gallery.png, ...}` (only if orbital VTIs exist)

LEED coordinate-mapping safeguard (per `docs/results_folder_structure_spec.md` §9.4 and §17.6): each screen is plotted both as a raw-index image *and* a coordinate-mapped image, with explicit handling of the writer's index-0-at-cell-centre convention surfaced by `run_06_centred_writer_check`.

---

## 5. Validation — Tier A only, on representative runs

Per the user's instruction, **only Tier A is run** in this round; Tiers B and C
are deferred. Tier A is run only on **representative** runs — not on every run.
Specifically:

- Each new GS save run (one per signature in §4.5 — three total) gets the full
  Tier-A pass, since loading a *new* electron configuration is the most likely
  point of failure.
- Each propagation run gets the **GPU-execution check** (below) and the
  results-tree `find` check, but the heavy injection / postprocess checks are
  not repeated for every variant in a parameter scan.

**Tier A checks**

- **GPU execution** (new, per user instruction): the run was launched with
  `inq-run` (not `inq-run --cpu`), and its log contains the INQ "GPU device" /
  CUDA initialisation banner. `nvidia-smi` shows the run's PID on the assigned
  GPU during execution. If either fails, the run is treated as invalid.
- Coronene xyz parses to 36 atoms, all z = 0 within 1e-12.
- `cell.contains()` is true for every atom (molecule fits the chosen cell).
- SCF converges to 1e-6 Ha; no NaN / Inf in printed energies.
- WP injection: `norm_after ∈ [0.97, 1.03]`, `max_overlap < 1e-3` against the occupied subspace.
- `RealField3DWriter` with `emit_raw=false, emit_vti=true` produces no zero-byte `.raw` files.
- `results/` tree passes the `find` checks at the bottom of `docs/results_folder_structure_spec.md` §22.
- `coronene_postprocess.py --results <dir>` completes every phase without error.

Tier B (energy-conservation, restart, CPU/GPU consistency) and Tier C (dt /
cutoff convergence, paper-figure comparison) are documented but not executed in
this pass; they are left as future work entries in
`docs/validation/coronene-replication.md`.

---

## 6. Critical files to read / write during implementation

To create or update:

- `ResearchProject/systems/coronene/shared/geometry/coronene.xyz` (copy of `01_geometry/coronene.xyz`)
- `ResearchProject/systems/coronene/shared/configs/tsubonoya_2014_base.hpp` (copy and adapt `inqkit/config/tsubonoya_2014_coronene.hpp`)
- `ResearchProject/systems/coronene/shared/configs/<variant>.hpp` (one per chosen run delta)
- `ResearchProject/systems/coronene/shared/cpp/results_paths.hpp` (new)
- `ResearchProject/systems/coronene/scripts/coronene_postprocess.py` (new thin CLI; per user instruction, **not** under `shared/`)
- `inq-stack/python/inqview/postprocess/*.py` (new generalisable modules — see §4.7)
- `ResearchProject/systems/coronene/save_gs/<gs_sig>/run.cpp` (one per unique GS tuple)
- `ResearchProject/systems/coronene/<run_name>/run.cpp` (one per chosen replication; flat sibling under `coronene/`, no `runs/` parent)
- `ResearchProject/systems/coronene/<run_name>/paths.hpp`
- `ResearchProject/systems/coronene/hypotheses/<NN>_*/README.md` (one per investigation; comparison artefacts produced alongside)
- `docs/sources/tsubonoya-2014-coronene-leed.md` (source note per `.claude/skills/literature-review.md`; every replica grounds in this paper)
- `docs/plans/coronene-replication.md` (mirror of this plan, per project rule `.claude/rules/file-placement.md`)
- `docs/handovers/coronene-replication.md` (rolling handover during execution, per `.claude/skills/handover-update.md`)
- `docs/validation/coronene-replication.md` (Tier-A/B/C validation log per the skill)
- `docs/todo_later.md` (append the legacy-cleanup item recorded in §7)

To reuse (read-only):

- `inq-stack/include/inqkit/io/real_field_3d_writer.hpp` (VTI writer config)
- `inq-stack/include/inqkit/io/complex_field_3d_writer.hpp` (for WP wavefunction)
- `inq-stack/include/inqkit/io/observables_writer.hpp`
- `inq-stack/include/inqkit/wavepacket/wavepacket.hpp`
- `inq-stack/include/inqkit/observables/orbital_overlap.hpp`
- `inq-stack/include/inqkit/screens/leed_pattern_accumulator.hpp`
- `inq-stack/python/inqview/*` (all modules listed in audit §A)
- `Tutorial/coronene-leed/run_diagnoses/run_07_paper_replica/{run.cpp, analysis.py}` (current closest reference)
- `ResearchProject/systems/coronene/run_propagate_paper_replica/run.cpp` (current load-from-checkpoint reference)
- `docs/results_folder_structure_spec.md` and `docs/visualisation-instructions-v1.md` (the binding specifications)

---

## 7. Decisions locked in

| Topic | Decision |
|---|---|
| Postprocess location | Generalisable core in `inq-stack/python/inqview/postprocess/`; coronene-specific wrapper at `ResearchProject/systems/coronene/scripts/coronene_postprocess.py`. |
| Cutoff for new replicas | **40 Ha** for everything (INQ-converged minimum on dojo PSPs). |
| GS checkpoint strategy | One GS save run per unique `(cell, cutoff)` tuple — three checkpoints total (§4.5). |
| Old buggy run directories | Leave untouched on disk for now. Wholesale legacy reorganisation deferred — entry written to `docs/todo_later.md` during implementation: *"Lump every legacy buggy coronene run (`04_leed_simulation/`, `coronene-wp-rt/`, the buggy `run_diagnoses` rows 4–8) under a single `legacy/` subtree once the new framework has reproduced the important ones."* |
| Units / filenames | Bohr throughout (§1.1). |
| Folder-structure precedence | This plan supersedes `.claude/rules/file-placement.md` for coronene work (§1.1). |
| Validation scope | Tier A only, on representative runs, with explicit GPU-execution check (§5). |
| Overlap computation | WP-only (§4.6 item 9); full O<sub>ij</sub> matrix dropped. |
| GPU dispatch | Two simultaneous runs across the two GPUs via `CUDA_VISIBLE_DEVICES`; queue refills as runs finish (§9). |
| Run selection | The five investigations of §10, *not* a re-run of any specific Table-A row. |

## 8. (Removed — superseded by §10)

---

## 9. Multi-GPU dispatch

**Hardware**: two GPUs available on the host (identical for practical
purposes). At any moment **one or both GPUs may already be busy with someone
else's work**, so the dispatcher must:

1. Detect which GPUs are *free* before queuing anything.
2. Hold at most one job per *free* GPU.
3. Re-poll for free GPUs as jobs finish, re-acquiring a GPU if it becomes
   available later.
4. If both GPUs are busy on entry, sleep-poll until at least one frees up.

**Free-GPU detection**: a GPU is considered free when, over a 3-poll window
(2 s apart), `nvidia-smi --query-compute-apps=pid --format=csv,noheader` shows
no PIDs on that GPU **and** memory utilisation is below 200 MiB. Any other
state means busy.

**Per-run launch**: each run's process is pinned via `CUDA_VISIBLE_DEVICES`:

```bash
CUDA_VISIBLE_DEVICES=<gpu_id> inq-run run.cpp >run.log 2>&1
```

The Tier-A "GPU execution" check (§5) confirms the run actually used the
intended GPU — `nvidia-smi` shows the run's PID on the assigned GPU and the
INQ log prints the GPU initialisation banner.

**Driver**: a Python script at
`ResearchProject/systems/coronene/scripts/dispatch_runs.py`:

- Reads a queue file (one run-directory path per line, e.g. `run_queue.txt`).
- Maintains an in-memory map `gpu_id → child_pid_or_None`.
- Main loop: every poll interval (default 30 s), for each GPU not currently
  running one of *our* children, check whether it is free; if so, pop the next
  run from the queue and launch it on that GPU.
- Logs `[gpu N] launched <run_dir>`, `[gpu N] finished <run_dir> exit=K
  walltime=...`, and `[gpu N] still busy (other process)` lines to a single
  `dispatch.log` so progress is visible from one tail.
- On Ctrl-C: cleanly waits for in-flight children before exiting; does not
  kill them.
- Includes a `--dry-run` mode that prints the launch commands without
  executing them — used for the user's review pass.
- Per-run argv: `dispatch_runs.py run_queue.txt [--poll-seconds 30]
  [--gpus 0,1] [--dry-run]`.

**Queue ordering**: the three GS save runs (§4.5) go first and are serialised
into the queue at the top, followed by the ten propagation runs (§10) in a
fixed order. The dispatcher itself is order-preserving: it pops the queue's
head whenever a free GPU appears.

**Operational note**: the user reports one GPU is currently occupied, so the
dispatcher is expected to start with one free GPU and ramp to two as the
other clears. No special-casing is needed beyond the free-GPU poll above.

---

## 10. Run set + hypotheses layout

Per the user's clarification:

- **Each propagation run is a flat sibling under `systems/coronene/`** —
  e.g. `systems/coronene/run_E30/`, `systems/coronene/run_b18_35x35x80/`. There
  is no `investigations/` or `runs/` parent.
- **Cross-run analysis lives in a new sibling folder
  `systems/coronene/hypotheses/`**, with one subfolder per investigation. Each
  hypotheses subfolder *only collates and compares* the existing runs'
  `results/` trees — no `run.cpp`, no GS, no propagation happens here.

### 10.1 Top-level layout

```
ResearchProject/systems/coronene/
├── shared/                                # geometry, configs, C++ helpers (§4.1–4.3)
├── scripts/                               # coronene_postprocess.py, dispatch_runs.sh (§4.7, §9)
├── save_gs/                               # one save_gs/<sig>/run.cpp per checkpoint (§4.5)
├── checkpoints/                           # produced by save_gs runs
│
├── run_base/                              # Tsubonoya base (b=12, σ=1.0, E=200, 35×35×60)
├── run_E30/                               # low energy
├── run_E800/                              # high energy
├── run_s0p33/                             # narrow WP
├── run_s3/                                # wide WP
├── run_E800_s0p33/                        # fast projectile (E=800, σ=0.33)
├── run_E30_s3/                            # electron-capture probe (E=30, σ=3)
├── run_b18_35x35x80/                      # large box, far b
├── run_b6_35x35x80/                       # large box, near b
├── run_35x35x40/                          # smaller box, base b
│
└── hypotheses/                            # comparison-only; no simulations
    ├── 00_base/                           # diagnostic summary of run_base
    ├── 01_wp_energy_spread/               # collates run_base, run_E30, run_E800
    ├── 02_wp_sigma_spread/                # collates run_base, run_s0p33, run_s3
    ├── 03_fast_projectile_classical/      # collates run_E800_s0p33 (vs run_base)
    ├── 04_electron_capture/               # collates run_E30_s3 (vs run_base)
    └── 05_box_length_and_distance/        # collates run_b18_35x35x80, run_b6_35x35x80, run_35x35x40 (vs run_base)
```

### 10.2 Run / GS-checkpoint mapping

| Run | GS checkpoint | σ (Bohr) | b (Bohr) | E (eV) | Cell (Bohr) | Hypothesis being tested |
|---|---|---|---|---|---|---|
| `run_base` | `gs_35x35x60_cut40` | 1.0 | 12 | 200 | 35×35×60 | Tsubonoya replica with corrected geometry. |
| `run_E30` | `gs_35x35x60_cut40` | 1.0 | 12 | **30** | 35×35×60 | Low-energy WP probe. |
| `run_E800` | `gs_35x35x60_cut40` | 1.0 | 12 | **800** | 35×35×60 | High-energy WP probe. |
| `run_s0p33` | `gs_35x35x60_cut40` | **0.33** | 12 | 200 | 35×35×60 | Narrow WP. |
| `run_s3` | `gs_35x35x60_cut40` | **3.0** | 12 | 200 | 35×35×60 | Wide WP. |
| `run_E800_s0p33` | `gs_35x35x60_cut40` | **0.33** | 12 | **800** | 35×35×60 | Classical fast-electron limit. |
| `run_E30_s3` | `gs_35x35x60_cut40` | **3.0** | 12 | **30** | 35×35×60 | Electron-capture probe. |
| `run_b18_35x35x80` | `gs_35x35x80_cut40` | 1.0 | **18** | 200 | **35×35×80** | Large box, far b. |
| `run_b6_35x35x80` | `gs_35x35x80_cut40` | 1.0 | **6** | 200 | **35×35×80** | Large box, near b. |
| `run_35x35x40` | `gs_35x35x40_cut40` | 1.0 | 12 | 200 | **35×35×40** | Smaller box, base b. |

Total: **9 propagation runs + the base = 10 propagations**, plus **3 GS save
runs** → 13 launches across the two-GPU queue.

### 10.3 What goes inside each `hypotheses/<NN>_*/` folder

Each hypothesis folder contains comparison artefacts only:

```
hypotheses/01_wp_energy_spread/
├── README.md                                       # which runs are compared, hypothesis statement, key takeaways
├── leed_screen_comparison_E30_E200_E800.png        # 3-panel LEED at the brightest screen
├── peak_intensity_vs_energy.png                    # log-log
├── energy_spectrum_comparison.png                  # FFT(total_energy) for the three runs
├── overlap_evolution_E30.gif / E200.gif / E800.gif # WP-overlap-with-GS-orbitals bar charts
└── ...                                             # any custom plot the hypothesis needs
```

The comparison plots are produced by a dedicated submodule
`inqview/postprocess/compare.py`, called from a new sibling subcommand:

```
coronene_postprocess.py hypothesis --hypothesis-dir <hypotheses/NN_*> \
    [--runs <abs_path_1> <abs_path_2> ...]
```

If `--runs` is omitted, the hypothesis dir's `README.md` (or a `runs.txt` next
to it) lists which `run_*/results/` trees feed into the comparison. The
submodule:

1. Reads each listed `results/` per the spec.
2. Aligns runs on a common time axis where applicable.
3. Renders side-by-side LEED panels, peak-intensity-vs-parameter plots,
   overlap-vs-time line plots, and any investigation-specific custom plots
   (e.g. for `04_electron_capture`, the asymptotic
   `O_i,wp(t → t_final)` for each occupied GS orbital — non-trivial residual
   population there is the electron-capture signature).

Hypotheses folders are flat (just PNGs / GIFs / README), not nested under
`raw/` + `analysis/` — those substructures only apply to actual simulation
`results/` trees.
