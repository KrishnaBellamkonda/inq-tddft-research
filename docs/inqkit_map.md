# inqkit library map (subtask 1 deliverable)

Canonical component map of the `inqkit` C++ header library
(`inq-stack/include/inqkit/`), produced for the unit-testing rejuvenation
(`docs/plans/inqkit-rejuvenation.md`). Built from the `understand-anything`
knowledge graph (`inq-stack/include/inqkit/.understand-anything/
knowledge-graph.json`, 37 files · 88 nodes · 112 edges · 10 layers · 14-step
tour — explore live with `/understand-dashboard`) and augmented with the
testing columns the plan requires (**Tier** = pure/engine, **Formula** =
formula-bearing, **Test approach**).

## Legend

- **Tier** — `pure` = no INQ engine, builds with a C++17 compiler alone
  (cheap CI lane, `ctest -L pure`); `engine` = links INQ (`inq::systems::
  electrons` / GPU / MPI), runs only where INQ is built.
- **Formula** — ✚ = formula-bearing (needs the two-agent formula
  verification before a test is locked); — = no nontrivial formula.
- **Status** — `target` = in the ~24-header testing surface; `stub` =
  0-byte / `// TODO` placeholder, **deferred** (handled in restructuring).

> **Refinement to the plan (logged):** the plan's bottom-up order assumed
> "fields/observables/io = engine". The graph + signature checks show the
> **I/O writers and two observables are actually `pure`** — they operate on
> POD `RealField3D`/`ComplexField3D`, not on INQ objects. Verified:
> `center_of_density(RealField3D const&)`, `RealField3DWriter::write(
> RealField3D const&)`, `density_delta::snapshot(RealField3D const&)` (host-
> only) vs `density::total(inq::systems::electrons const&)` (engine). This
> widens the pure lane — good for CI. Confirm the writer trio
> (`observables_writer`/`occupations_writer`/`state_energy_writer`) at their
> subtask-2 turn.

## Testing surface (bottom-up)

### 1 · Foundational schema & grid (`detail/`)
| Header | Tier | Formula | Role · API · test approach |
|---|---|---|---|
| `detail/grid_layout.hpp` | pure | ✚ | Raw on-disk field schemas + the flatten convention `flat=((ix·ny)+iy)·nz+iz` and `step_suffix`. API: `ComplexField3DRawSchema`, `RealField3DRawSchema`, `flatten_index`, `step_suffix`. **Test:** assert `flatten_index` against hand-computed indices on a tiny grid; schema string fields fixed. |

### 2 · Field data model & extraction (`fields/`)
| Header | Tier | Formula | Role · API · test approach |
|---|---|---|---|
| `fields/real_field_3d.hpp` | pure | — | POD 3D real field (dims, origin, spacing, flat buffer). API: `RealField3D`. **Test:** construction + flat-index accessors on a toy grid. |
| `fields/complex_field_3d.hpp` | pure | — | POD 3D complex field. API: `ComplexField3D`. **Test:** as above. |
| `fields/density.hpp` | engine | ✚ | Extracts ρ(r) from INQ `electrons` into `RealField3D` (total / per-orbital / total-excluding-orbital) with MPI gather + FFT-shift reindex. API: `total`, `orbital`, `total_excluding_orbital`. **Test:** tiny INQ system, constant/known density → integral and FFT-shift placement. |
| `fields/orbital.hpp` | engine | ✚ | Extracts one complex KS orbital from INQ into centred-box `ComplexField3D`. API: `wavefunction`. **Test:** known orbital on small grid → reindex + norm. |

### 3 · Serialization & writers (`io/`)
| Header | Tier | Formula | Role · API · test approach |
|---|---|---|---|
| `io/real_field_3d_writer.hpp` | pure | — | Writes `RealField3D` → raw binary + `.meta.txt` (+ optional VTI), with per-step series overloads. API: `RealField3DWriter`. **Test:** write toy field, read back (Python `inqview`), assert round-trip + layout. |
| `io/complex_field_3d_writer.hpp` | pure | — | Writes `ComplexField3D` → `_real.raw`/`_imag.raw` + meta (+ optional VTI). API: `ComplexField3DWriter`. **Test:** round-trip real/imag parts. |
| `io/vti_image_data_writer.hpp` | pure | — | Serialises real/complex fields to VTK ImageData `.vti` (ASCII + base64 appended). API: `VTIImageDataWriter`, `write_real`, `write_complex`. **Test:** toy field → valid `.vti` parsed by VTK/ParaView; base64 decode matches. |
| `io/observables_writer.hpp` | pure | — | Streams a scalar-observable selection (energies, current, dipole, CoD, density-L2) to CSV. API: `ObservableSelection`, `ObservablesWriter`. **Test:** known scalars → exact CSV header + precision. |

### 4 · Observables & analysis (`observables/`)
| Header | Tier | Formula | Role · API · test approach |
|---|---|---|---|
| `observables/center_of_density.hpp` | pure | ✚ | Density-weighted centroid `∫r·n/∫n` of a `RealField3D` → Bohr coords + total weight. API: `center_of_density`. **Test:** offset Gaussian → known centroid; constant field → box centre. |
| `observables/density_delta.hpp` | pure | ✚ | `n − n_ref` raw + coarse-grained VTI snapshots + optional L2. API: `DensityDelta`. **Test:** known ref/current pair → delta values + L2 norm. |
| `observables/momentum_distribution.hpp` | engine | ✚ | Per-step FFT of WP orbital → radial `|ψ̃(k)|²` histogram → CSV. API: `MomentumDistributionConfig`, `MomentumDistribution`. **Test:** plane-wave orbital → peak at known `k₀`. |
| `observables/wp_momentum_stats.hpp` | engine | ✚ | GPU moments `⟨k⟩,⟨k²⟩,σ_k` over Fourier orbital → CSV. API: `WPMomentumStatsConfig`, `WPMomentumStats`. **Test:** Gaussian WP → analytic `⟨k⟩=k₀`, `σ_k=1/σ_r`. |
| `observables/wp_real_space_stats.hpp` | engine | ✚ | GPU real-space moments (CoD, spread) over orbital amplitude → CSV. API: `WPRealSpaceStatsConfig`, `WPRealSpaceStats`. **Test:** Gaussian WP → `⟨r⟩`, `σ_r`. |
| `observables/orbital_overlap.hpp` | engine | ✚ | `|⟨φ_refᵢ|ψ_evolⱼ⟩|²` matrices vs cached GS refs → per-step CSV. API: `OrbitalOverlapMatrix`, `snapshot{,_wp_only,_proxies}`. **Test:** identity at t=0 (diag≈1, off≈0). |
| `observables/state_energy_writer.hpp` | engine | ✚ | Per-state `⟨ψ|H|ψ⟩` (+ optional variance) each step → CSV. API: `StateEnergyWriter`. **Test:** eigenstate → energy = eigenvalue, variance≈0. |
| `observables/eigenvalue_dump.hpp` | engine | — | Per-k KS eigenvalues + occupations from INQ → CSV. API: `dump_eigenvalues`. **Test:** known small system → values match SCF. |
| `observables/occupations_writer.hpp` | engine | — | Per-state occupations (iter, time cols) → CSV during RT. API: `OccupationsWriter`. **Test:** known occupations → exact CSV. |

### 5 · Wavepacket injection (`wavepacket/`)
| Header | Tier | Formula | Role · API · test approach |
|---|---|---|---|
| `wavepacket/wavepacket.hpp` | engine | ✚ | Gaussian × plane-wave WP injected into last extra KS state; GPU norm check, optional Gram-Schmidt, returns report. API: `WavePacket`, `inject_into_last_extra_state`. **Test:** inject → norm≈1, momentum peak at `k₀`. |
| `wavepacket/injection_report.hpp` | pure | — | POD result (target k/state, norms, max overlap, ortho/tol flags). API: `InjectionReport`. **Test:** field defaults + flag logic. |

### 6 · Real-time orchestration (`real_time/`)
| Header | Tier | Formula | Role · API · test approach |
|---|---|---|---|
| `real_time/step_context.hpp` | pure | — | POD per-step bundle (step, time, ions, electrons, energies, current, dipole, WP centre, density-L2). API: `StepContext`. **Test:** construction + field wiring. |
| `real_time/real_time_session.hpp` | engine | — | Drives RT TDDFT: builds `StepContext` per step, dispatches to registered observable tasks. API: `RealTimeSession`. **Test (integration):** few-step run dispatches tasks in order. |

### 7 · LEED screens (`screens/`, in-progress — test completed paths only)
| Header | Tier | Formula | Role · API · test approach |
|---|---|---|---|
| `screens/plane_screen.hpp` | engine | ✚ | 2D density slice at grid plane nearest target z (sum of occupied amplitudes). API: `PlaneScreen`, `extract`, `save`. **Note:** carries internal `// TODO`. **Test:** plane-pick + slice integral on small grid. |
| `screens/leed_pattern_accumulator.hpp` | engine | ✚ | Time-integrates plane density into a LEED intensity map → disk. API: `LeedPatternAccumulator`, `accumulate`, `save`. **Test:** constant slice → linear-in-time accumulation. |

### Jellium analytics (`jellium/`)
| Header | Tier | Formula | Role · API · test approach |
|---|---|---|---|
| `jellium/shells.hpp` | pure | ✚ | Enumerates free-electron `|G|²` plane-wave shells + degeneracies; picks unoccupied proxy states; writes shell CSV. API: `ShellInfo`, `default_shell_table`, `enumerate_for_n_states`, `write_shells_csv`. **Test:** first shells `|G|²=0,1,2,3…` with degeneracies vs `docs/sources/free-electron-gas-magic-numbers.md`. |

### Configuration (`config/`)
| Header | Tier | Formula | Role · API · test approach |
|---|---|---|---|
| `config/tsubonoya_2014_coronene.hpp` | pure | — | Compile-time constants reproducing the coronene LEED setup of Tsubonoya et al. (2014). **Test:** constants equal the cited paper values. |

## Deferred placeholders (0-byte / `// TODO` — NOT tested this pass)

`core/pipeline.hpp`, `core/session_context.hpp`, `core/task.hpp`,
`jellium/analytics.hpp`, `config/simulation_config.hpp`,
`ground_state/ground_state_tasks.hpp`, `detail/validation.hpp`,
`detail/filesystem.hpp`, `detail/text_io.hpp`, `io/manifest_writer.hpp`,
`io/text_summary_writer.hpp`. Confirmed never implemented on any branch; see
`docs/notes/inqkit-rejuvenation-ideas.md`. **CLAUDE.md drift:** advertises
populated `jellium/analytics` + `core/*` that are empty.

## Count reconciliation

37 headers = 13 pure targets + 11 engine targets + (13 deferred stubs).
The "~24 testing surface" in the plan = 13 pure + 11 engine. Note the pure
count rose vs the plan (writers + `center_of_density` + `density_delta`
reclassified pure).
