# Plan: inqkit real-time density writing — N2 tutorial

*Created: 2026-04-19 | Branch: features/real-time-orbital-density-writing*

---

## Context

The inqkit library currently supports ground-state field writing (GS snapshot → raw file).
The next library milestone is writing fields during real-time TDDFT propagation: the total
electronic density must be saved at a regular stride throughout the simulation, producing a
time-indexed file series that inqview can load as a `FieldSeries` and render as a movie.

This task implements the minimal real-time writing infrastructure and validates it on N2.
It does NOT implement orbital or wavepacket density writing (those come later).

---

## Physics system

| Parameter | Value | Rationale |
|---|---|---|
| System | N₂ molecule, slightly stretched | Simple diatomic; 5 occupied orbitals; well-studied |
| Bond length | 1.15 Å (stretch from equil. 1.098 Å) | Creates non-equilibrium geometry for richer dynamics |
| Cell | finite cubic, L = 20 bohr | Adequate vacuum (~9 bohr each side); smaller than n2-with-inqkit (30 bohr) |
| E_cut | 20 Ha (40 Ry) | Coarse; allows dt ≈ 0.1 au while keeping grid manageable |
| dt | 0.1 au | Near ETRS stability limit for 20 Ha cutoff |
| N_steps | 6100 | ≈ 610 au total time; one N₂ vibrational period (T ≈ 583 au at ~2359 cm⁻¹) |
| Perturbation | impulsive kick (`.impulsive()`) | Excites the system; generates density oscillations |
| Theory | LDA | Adequate for I/O validation |

**Expected oscillations**: the kick excites the N₂ electronic system. The total electron
density will oscillate at electronic excitation frequencies. Over 610 au of total time,
multiple oscillation cycles should be visible in the written density frames.

**Atom placement** (centered at L/2 = 10 bohr, not at corner):
```cpp
auto L = 20.0_bohr;
auto half_bond = 1.15_angstrom / 2;   // 1.087 bohr
ions.insert("N", {L/2, L/2, L/2 - half_bond});   // (10, 10, 8.913) bohr
ions.insert("N", {L/2, L/2, L/2 + half_bond});   // (10, 10, 11.087) bohr
```

---

## File size estimate

| Parameter | Value |
|---|---|
| Grid dimensions | ~40 × 40 × 40 (from N = L√(2 E_cut)/π = 40.3) |
| Grid points | ~64,000 |
| Bytes per frame (.raw, float64) | ~512 KB |
| Write stride | every 100 steps |
| Frames written | 6100 / 100 = 61 |
| Total density data | 61 × 512 KB ≈ 31 MB |
| Observables CSV | 6100 rows ≈ 0.6 MB |
| **Total** | **~32 MB** |

**Permission requested**: user must approve before the simulation is launched (it runs
for ~10–20 minutes on GPU).

---

## inqkit C++ work required

### 1. `real_time/step_context.hpp` (implement stub)

```cpp
#include <inq/inq.hpp>

namespace inqkit {

struct StepContext {
    int    step   = 0;
    double time_au = 0.0;
    inq::systems::ions     const* ions      = nullptr;
    inq::systems::electrons const* electrons = nullptr;
};

}
```

### 2. `real_time/real_time_session.hpp` (implement stub)

Minimal implementation: owns a stride-aware callback that:
- Increments the step counter
- Builds `StepContext`
- Calls each registered task's `on_step(ctx)` method

```cpp
class RealTimeSession {
public:
    RealTimeSession(inq::systems::ions& ions,
                    inq::systems::electrons& electrons,
                    int write_every = 1);

    // Add a task that will be called every write_every steps.
    // Task must have void operator()(StepContext const&).
    template <typename TTask>
    void add(TTask task);

    // Call this inside the propagate lambda.
    template <typename RTData>
    void step(RTData const& data);
};
```

### 3. `io/real_field_3d_writer.hpp` — add real-time write overload

Add an overload that generates step-indexed filenames:

```cpp
// Real-time write: basename = field_name + "_t" + zero-padded step
void write(RealField3D const& field, double time_au, int step) const;
```

Filename convention: `total_density_t000000.raw` / `total_density_t000000.meta.txt`

The meta.txt gains a `time_au = ...` field to record the simulation time.

---

## Tutorial: `Tutorial/n2-real-time-with-inqkit/`

### Files to create
- `run.cpp`
- `analysis.py`

### `run.cpp` sketch

```cpp
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/real_time/real_time_session.hpp>

using namespace inq;
using namespace inq::magnitude;

int main() {
    auto L = 20.0_bohr;
    auto cell = systems::cell::cubic(L).finite();
    auto half_bond = 1.15_angstrom / 2;

    systems::ions ions(cell);
    ions.insert("N", {L/2, L/2, L/2 - half_bond});
    ions.insert("N", {L/2, L/2, L/2 + half_bond});

    systems::electrons electrons(ions,
        options::electrons{}.cutoff(20.0_Ha));

    // Ground state
    ground_state::initial_guess(ions, electrons);
    ground_state::calculate(ions, electrons,
        options::theory{}.lda(),
        options::ground_state{}.energy_tolerance(1e-6_Ha));

    // Real-time session
    inqkit::RealTimeSession rt(ions, electrons, /*write_every=*/100);

    rt.add(inqkit::io::RealField3DWriter(
        "results/real_time/density",
        {.field_name = "total_density", .include_meta = true}
    ));

    real_time::propagate(
        ions, electrons,
        [&](auto const& data) { rt.step(data); },
        options::theory{}.lda(),
        options::real_time{}
            .num_steps(6100)
            .dt(0.1_atomictime)
            .impulsive()
    );
}
```

### `analysis.py` sketch

```python
sim = SimulationData("results")
density_series = sim.field_series("real_time/density")
# density_series.files: list of 61 meta.txt paths, time-ordered

# Validate norms at several timesteps
for i, meta_path in enumerate(density_series.files[::10]):
    rho = load_real_field(meta_path=meta_path)
    n_elec = rho.array.sum() * rho.meta.voxel_volume_bohr3
    print(f"  step {i*100*10}: N_elec = {n_elec:.3f}")

# Density slice at z-midplane, several time points
# (shows how density evolves over time)
# ...

# VTI series → ParaView movie
pv.render_density_from_meta_series(
    density_series,
    vti_output_dir=...,
    render=VolumeRenderSpec(array_name="total_density"),
    animation=AnimationSpec(output_frames_dir=..., image_size=(1600, 1200)),
    atoms=AtomSpec(
        positions=[[10, 10, 8.913], [10, 10, 11.087]],
        symbols=["N", "N"],
    ),
)
pv.build_gif(..., fps=12)
```

---

## Output schema

```
results/
  real_time/
    density/
      total_density_t000000.raw
      total_density_t000000.meta.txt
      total_density_t000100.raw
      total_density_t000100.meta.txt
      ... (61 files total)
```

The meta.txt for each frame gains a `time_au` field matching the INQ simulation time.

---

## Validation criteria

| Check | Pass condition |
|---|---|
| N₂ electrons | N_elec = 10.0 ± 0.1 at every written frame |
| Density norm stability | N_elec varies by < 0.01 across all 61 frames (energy conservation) |
| 61 frames written | All `total_density_t{step:06d}.raw` files present |
| meta.txt time field | `time_au` field present and monotonically increasing |
| ParaView movie | 61-frame GIF shows density evolution; N atoms visible as CPK spheres |

---

## Execution order

1. Implement `StepContext`, `RealTimeSession` (minimal), `write(field, time_au, step)` overload
2. Write `Tutorial/n2-real-time-with-inqkit/run.cpp`
3. Build + run (`inq-run`, GPU, ~10–20 min) — **requires user permission**
4. Write `analysis.py`
5. Run analysis, capture output, validate norms
6. Build ParaView GIF, visually inspect
7. Record results in handover

---

## Out of scope (this task)

- Orbital density writing during real-time
- Complex wavefunction writing during real-time
- Observables CSV (`ObservablesWriter`) — add in a follow-up if time allows
- Ehrenfest (ionic) dynamics
- TDDFT with laser field
