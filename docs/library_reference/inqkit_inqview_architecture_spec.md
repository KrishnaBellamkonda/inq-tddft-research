# inqkit and inqview Architecture Specification

This document fixes the current agreed structure of the C++ **inqkit** library and the Python **inqview** library.
It is an API and architecture plan, not a compile-ready implementation. Names, method signatures, and example code are intended to stabilise the design before coding begins.

---

# 1. Library summaries

## inqkit (C++)
inqkit is a helper library layered on top of INQ that abstracts recurring simulation tasks without hiding the underlying INQ workflow.
Its main roles are: typed simulation configuration, structured input/output, ground-state and real-time task orchestration, wavepacket injection, detector/screen accumulation, and reusable analytics helpers.

## inqview (Python)
inqview is a lightweight post-processing and visualisation library that reads the standardised outputs produced by inqkit.
Its main roles are: loading and validating simulation data, plotting observables, performing Fourier/post-processing analysis, generating LEED and other diagnostic plots, and automating ParaView-based volume-rendering workflows.

---

# 2. Design principles

## 2.1 Common principles
- Keep the **native INQ workflow explicit**. The user still writes ordinary C++/INQ code.
- Standardise **output schema**, **writer APIs**, and **analysis entry points**.
- Separate **human-readable outputs** from **performance-mode outputs** and from **INQ restart outputs**.
- Prefer **identity names for classes** and **action names for functions**.
- Prefer **typed structs** over magic strings wherever possible.

## 2.2 Writer API rule
All writer classes should follow the same broad pattern:
1. **First constructor argument**: output path.
2. **Second constructor argument**: a typed layout/spec/selection object.
3. **Third constructor argument**: options.
4. Main action method: `.write(...)`.
5. Optional streaming method: `.append(...)`.
6. Optional end-of-run method: `.finish()`.

## 2.3 Output format rule
### Human-readable outputs
Use:
- `.txt` for key-value summaries and metadata.
- `.csv` for time-series/tabular outputs.
- `.txt` grid-based field dumps where readability is required.

### Performance-mode outputs
Use:
- raw binary arrays such as `.raw` for field data.
- small sidecar `.meta.txt` files describing dimensions, spacing, type, layout, units, and time.

### Restart outputs
Use:
- native INQ restart data via `electrons.save(...)`.

---

# 3. Output schema

## 3.1 Human-readable schema

```text
results/
  manifest.txt
  simulation_summary.txt

  config/
    simulation_config.txt
    derived_values.txt

  ground_state/
    summary.txt
    observables.csv
    density/
    orbitals/
    checkpoint/

  real_time/
    observables.csv
    density/
    wavepacket_density/
    orbitals/
    screens/
    leed/
    diagnostics/

  visualisation/
    vti/
    frames/
    movies/

  analysis/
    plots/
    spectra/
```

## 3.2 Performance schema

### Real scalar field
```text
total_density_t000250.raw
total_density_t000250.meta.txt
```

### Complex field
```text
homo_t000260_real.raw
homo_t000260_imag.raw
homo_t000260.meta.txt
```

### Example sidecar metadata
```text
type = real_field_3d
dtype = float64
nx = 160
ny = 160
nz = 160
origin_bohr = 0 0 0
spacing_bohr = 0.25 0.25 0.25
layout = x_slowest_z_fastest
time_au = 5.0
field_name = total_density
units = bohr^-3
```

## 3.3 VTI schema

```text
results/visualisation/vti/
  total_density/
    total_density_t000000.vti
    total_density_t000010.vti
    ...

  wavepacket_density/
    wavepacket_density_t000000.vti
    wavepacket_density_t000010.vti
    ...

  orbitals/
    homo_t000000_real.vti
    homo_t000000_imag.vti
    lumo_t000000_real.vti
    lumo_t000000_imag.vti
```

---

# 4. inqkit (C++) structure

## 4.1 Directory tree

```text
inqkit/
  include/inqkit/
    config/
      simulation_config.hpp

    core/
      task.hpp
      pipeline.hpp
      session_context.hpp

    io/
      manifest_writer.hpp
      text_summary_writer.hpp
      observables_writer.hpp
      real_field_3d_writer.hpp
      complex_field_3d_writer.hpp
      vti_image_data_writer.hpp

    ground_state/
      ground_state_tasks.hpp

    real_time/
      step_context.hpp
      real_time_session.hpp

    wavepacket/
      wavepacket.hpp
      injection_report.hpp

    screens/
      plane_screen.hpp
      leed_pattern_accumulator.hpp

    jellium/
      analytics.hpp

    detail/
      filesystem.hpp
      validation.hpp
      text_io.hpp
      grid_layout.hpp
```

---

## 4.2 File-by-file specification

## `config/simulation_config.hpp`
**Function:** Defines the typed simulation configuration object used by the simulation and written to disk.

### Main types
- `SimulationConfig`
- `SystemConfig`
- `GroundStateConfig`
- `RealTimeConfig`
- `WavePacketConfig`
- `OutputConfig`
- `ObservableSelection`

### Responsibilities
- Hold the parameters for the whole run.
- Provide derived quantities such as timestep count, output cadence, packet momentum magnitude, and selected orbital indices.
- Provide a single source of truth for the simulation and its output writers.

### API sketch
```cpp
struct SystemConfig {
    double lx_bohr = 0.0;
    double ly_bohr = 0.0;
    double lz_bohr = 0.0;
    double cutoff_ha = 0.0;
    int extra_states = 0;
    double temperature_ev = 0.0;
    bool gamma_only = true;
};

struct GroundStateConfig {
    double energy_tolerance_ha = 1e-8;
    double mixing_alpha = 0.1;
    int max_steps = 300;
    bool use_broyden = true;
};

struct RealTimeConfig {
    int num_steps = 0;
    double dt_au = 0.0;
    int write_every_density = 10;
    int write_every_orbitals = 20;
    bool observe_current = true;
    bool observe_dipole = false;
};

struct WavePacketConfig {
    bool enabled = false;
    double center_x_bohr = 0.0;
    double center_y_bohr = 0.0;
    double center_z_bohr = 0.0;
    double sigma_bohr = 1.0;
    double k0_x_bohr_inv = 0.0;
    double k0_y_bohr_inv = 0.0;
    double k0_z_bohr_inv = 0.0;
    double occupation = 1.0;
};

struct OutputConfig {
    std::string root = "results";
    bool human_readable = true;
    bool performance_mode = false;
    bool write_vti = true;
};

struct ObservableSelection {
    bool step = true;
    bool time_au = true;
    bool energy_total = true;
    bool energy_kinetic = true;
    bool energy_hartree = false;
    bool energy_xc = false;
    bool current_x = true;
    bool current_y = true;
    bool current_z = true;
    bool dipole_x = false;
    bool dipole_y = false;
    bool dipole_z = false;
};

struct SimulationConfig {
    SystemConfig system;
    GroundStateConfig gs;
    RealTimeConfig rt;
    WavePacketConfig wp;
    OutputConfig output;
    ObservableSelection observables;

    int homo_index() const;
    int lumo_index() const;
    std::string root_path() const;
};
```

---

## `core/task.hpp`
**Function:** Defines the minimal task interface shared by ground-state and real-time task objects.

### Main types
- `Task`
- `TaskStage`

### Responsibilities
- Provide a common interface for objects that participate in a session.
- Allow optional lifecycle hooks without forcing every task to implement all of them.

### API sketch
```cpp
enum class TaskStage {
    run_start,
    ground_state_complete,
    real_time_step,
    run_end
};

class Task {
public:
    virtual ~Task() = default;

    virtual void on_run_start() {}
    virtual void on_ground_state_complete() {}
    virtual void on_real_time_step() {}
    virtual void on_run_end() {}
};
```

---

## `core/pipeline.hpp`
**Function:** Owns an ordered collection of tasks and dispatches lifecycle events to them.

### Main types
- `Pipeline`

### Responsibilities
- Store task objects.
- Ensure tasks are executed in a deterministic order.
- Provide a shared composition mechanism used by GS and RT orchestration.

### API sketch
```cpp
class Pipeline {
public:
    template <class TTask>
    void add(TTask task);

    void run_start();
    void ground_state_complete();
    void real_time_step();
    void run_end();
};
```

---

## `core/session_context.hpp`
**Function:** Defines shared contextual state available to session-level code.

### Main types
- `SessionContext`

### Responsibilities
- Hold references to config, ions, electrons, and optional shared payloads.
- Provide a clean way to extend future tasks without changing every task signature.

### API sketch
```cpp
struct SessionContext {
    SimulationConfig const* config = nullptr;
    inq::systems::ions* ions = nullptr;
    inq::systems::electrons* electrons = nullptr;
    void* user_context = nullptr;
};
```

---

## `io/manifest_writer.hpp`
**Function:** Writes the minimal machine-readable run manifest and schema/provenance information.

### Main types
- `ManifestWriter`
- `ManifestOptions`

### Responsibilities
- Write `manifest.txt`.
- Record schema version, run type, file layout version, and provenance.
- Point Python code to the canonical root of the run.

### API sketch
```cpp
struct ManifestOptions {
    std::string schema_version = "0.1";
    std::string code_name = "inqkit";
    std::string run_type;
};

class ManifestWriter {
public:
    ManifestWriter(std::string path, ManifestOptions options = {});
    void write(SimulationConfig const& config);
};
```

---

## `io/text_summary_writer.hpp`
**Function:** Writes human-readable key-value summaries and diagnostic reports.

### Main types
- `TextSummaryWriter`
- `SummaryOptions`

### Responsibilities
- Write `simulation_summary.txt`, GS summaries, and wavepacket diagnostics.
- Write simple text summaries that can be checked quickly by a human.

### API sketch
```cpp
struct SummaryOptions {
    bool overwrite = true;
};

class TextSummaryWriter {
public:
    TextSummaryWriter(std::string path, SummaryOptions options = {});

    void write_simulation_summary(SimulationConfig const& config);
    void write_ground_state_summary(/* gs result */);
    void write_wavepacket_report(/* InjectionReport const& report */);
    void write_key_value(std::string key, std::string value);
};
```

---

## `io/observables_writer.hpp`
**Function:** Streams selected scalar and vector observables to a standard CSV file.

### Main types
- `ObservablesWriter`
- `ObservableSelection`
- `ObservablesWriterOptions`

### Responsibilities
- Write `observables.csv` with a consistent column order.
- Use typed selection rather than string column names.
- Support appending one row per timestep.

### API sketch
```cpp
struct ObservablesWriterOptions {
    bool overwrite = true;
    char separator = ',';
};

class ObservablesWriter {
public:
    ObservablesWriter(
        std::string path,
        ObservableSelection selection,
        ObservablesWriterOptions options = {}
    );

    void write_header();
    void append(/* StepContext const& ctx */);
    void finish();
};
```

---

## `io/real_field_3d_writer.hpp`
**Function:** Writes one real-valued 3D scalar field to the standard text, raw, or both output formats.

### Main types
- `RealField3DWriter`
- `RealField3DLayout`
- `RealField3DWriteOptions`

### Responsibilities
- Write total density, wavepacket density, charge-density slices, or any other real field.
- Support human-readable text output, raw binary output, or both.
- Standardise metadata and naming.

### API sketch
```cpp
struct RealField3DLayout {
    std::string field_name;
    bool write_text = true;
    bool write_raw = false;
    bool include_meta = true;
};

struct RealField3DWriteOptions {
    bool overwrite = true;
    int stride = 1;
};

class RealField3DWriter {
public:
    RealField3DWriter(
        std::string path,
        RealField3DLayout layout,
        RealField3DWriteOptions options = {}
    );

    void write(/* field object */, double time_au, int step);
};
```

---

## `io/complex_field_3d_writer.hpp`
**Function:** Writes one complex-valued 3D field into separate real and imaginary outputs, plus metadata.

### Main types
- `ComplexField3DWriter`
- `ComplexField3DLayout`
- `ComplexField3DWriteOptions`

### Responsibilities
- Write HOMO/LUMO orbitals or any other complex orbital field.
- Write separate real and imaginary outputs in text or raw form.
- Standardise naming and metadata for orbital outputs.

### API sketch
```cpp
struct ComplexField3DLayout {
    std::string field_name;
    bool write_text = true;
    bool write_raw = false;
    bool include_meta = true;
};

struct ComplexField3DWriteOptions {
    bool overwrite = true;
    int stride = 1;
};

class ComplexField3DWriter {
public:
    ComplexField3DWriter(
        std::string path,
        ComplexField3DLayout layout,
        ComplexField3DWriteOptions options = {}
    );

    void write(/* complex field object */, double time_au, int step);
};
```

---

## `io/vti_image_data_writer.hpp`
**Function:** Exports real scalar fields to VTI image-data files for ParaView workflows.

### Main types
- `VTIImageDataWriter`
- `VTIImageDataLayout`
- `VTIWriteOptions`

### Responsibilities
- Write `.vti` files for volume rendering and image-data time series.
- Standardise array naming and time-index naming.
- Prepare outputs for ParaView automation in inqview.

### API sketch
```cpp
struct VTIImageDataLayout {
    std::string array_name;
};

struct VTIWriteOptions {
    bool overwrite = true;
};

class VTIImageDataWriter {
public:
    VTIImageDataWriter(
        std::string path,
        VTIImageDataLayout layout,
        VTIWriteOptions options = {}
    );

    void write(/* real field object */, double time_au, int step);
};
```

---

## `ground_state/ground_state_tasks.hpp`
**Function:** Orchestrates tasks that run once the ground state has been found.

### Main types
- `GroundStateTasks`

### Responsibilities
- Own and execute GS-related writers and post-SCF tasks.
- Keep GS orchestration separate from the writers themselves.
- Save the checkpoint and any selected GS outputs.

### API sketch
```cpp
class GroundStateTasks {
public:
    GroundStateTasks(SimulationConfig const& config);

    template <class TTask>
    void add(TTask task);

    void run(/* gs result */, inq::systems::ions& ions, inq::systems::electrons& electrons);
};
```

---

## `real_time/step_context.hpp`
**Function:** Defines the typed per-step context passed to real-time tasks.

### Main types
- `StepContext`

### Responsibilities
- Carry the current step data plus references to config, ions, electrons, and optional GS references.
- Avoid fragile callback signatures that have to be rewritten whenever new information is needed.

### API sketch
```cpp
struct StepContext {
    // INQ time-step data
    /* real-time data view */ const* data = nullptr;

    // Shared session state
    SimulationConfig const* config = nullptr;
    inq::systems::ions const* ions = nullptr;
    inq::systems::electrons const* electrons = nullptr;
    inq::systems::electrons const* gs_reference = nullptr;
    void const* user_context = nullptr;
};
```

---

## `real_time/real_time_session.hpp`
**Function:** Owns the real-time task pipeline, persistent writer state, and end-of-run accumulation logic.

### Main types
- `RealTimeSession`
- `RealTimeSessionOptions`

### Responsibilities
- Hold references to config, ions, electrons, screens, and optional extra data.
- Convert the raw INQ callback into a stable task-driven session.
- Provide `.step(...)` and `.finish()` as the public API.

### API sketch
```cpp
struct RealTimeSessionOptions {
    inq::systems::electrons const* gs_reference = nullptr;
    void const* user_context = nullptr;
};

class RealTimeSession {
public:
    RealTimeSession(
        SimulationConfig const& config,
        inq::systems::ions& ions,
        inq::systems::electrons& electrons,
        RealTimeSessionOptions options = {}
    );

    template <class TTask>
    void add(TTask task);

    void start();
    void step(/* real-time data view */ const& data);
    void finish();
};
```

---

## `wavepacket/wavepacket.hpp`
**Function:** Defines the wavepacket object, its builder-style parameter API, orthogonalisation methods, and injection methods.

### Main types
- `WavePacket`

### Responsibilities
- Store packet parameters in bohr/bohr^-1 conventions.
- Support builder-style definition.
- Orthogonalise against occupied or all states.
- Inject into the chosen target state and return an insertion report.

### API sketch
```cpp
class WavePacket {
public:
    WavePacket& center(double x_bohr, double y_bohr, double z_bohr);
    WavePacket& sigma(double sigma_bohr);
    WavePacket& k0(double kx_bohr_inv, double ky_bohr_inv, double kz_bohr_inv);

    WavePacket& orthogonalise_against_occupied(
        inq::systems::electrons const& electrons,
        double tolerance
    );

    WavePacket& orthogonalise_against_all_states(
        inq::systems::electrons const& electrons,
        double tolerance
    );

    /* InjectionReport */ auto inject_into_last_extra_state(
        inq::systems::electrons& electrons,
        double occupation
    ) const;

    /* InjectionReport */ auto inject_into_state(
        inq::systems::electrons& electrons,
        int kpoint_index,
        int state_index,
        double occupation
    ) const;
};
```

---

## `wavepacket/injection_report.hpp`
**Function:** Stores the diagnostics of a wavepacket insertion and orthogonalisation operation.

### Main types
- `InjectionReport`

### Responsibilities
- Record the target state.
- Record norms and overlaps before and after orthogonalisation.
- Record pass/fail information relative to tolerance.

### API sketch
```cpp
struct InjectionReport {
    int kpoint_index = 0;
    int state_index = -1;
    double norm_before = 0.0;
    double norm_after = 0.0;
    double max_overlap_before = 0.0;
    double max_overlap_after = 0.0;
    bool orthogonalised = false;
    bool passed_tolerance = false;
};
```

---

## `screens/plane_screen.hpp`
**Function:** Defines a detector plane in Cartesian space.

### Main types
- `PlaneScreen`

### Responsibilities
- Represent measurement surfaces such as detector planes.
- Start with `at_z(...)` and support future extensions such as `at_x(...)` and `at_y(...)`.

### API sketch
```cpp
class PlaneScreen {
public:
    static PlaneScreen at_z(double z_bohr);
    static PlaneScreen at_x(double x_bohr);
    static PlaneScreen at_y(double y_bohr);

    double position_bohr() const;
};
```

---

## `screens/leed_pattern_accumulator.hpp`
**Function:** Accumulates time-integrated density on a plane screen and writes the final LEED pattern.

### Main types
- `LEEDPatternAccumulator`
- `LEEDAccumulatorOptions`

### Responsibilities
- Evaluate the density on the detector plane during the requested time window.
- Accumulate the integrated LEED pattern.
- Write the final pattern and metadata at the end of the run.

### API sketch
```cpp
struct LEEDAccumulatorOptions {
    double start_time_au = 0.0;
    double end_time_au = 0.0;
};

class LEEDPatternAccumulator {
public:
    LEEDPatternAccumulator(
        std::string path,
        PlaneScreen screen,
        LEEDAccumulatorOptions options = {}
    );

    void append(StepContext const& ctx);
    void finish();
};
```

---

## `jellium/analytics.hpp`
**Function:** Provides analytic reference functions for jellium and related validation tasks.

### Main types
- Prefer namespace functions, not a heavy class, in the first implementation.

### Responsibilities
- Compute `r_s`, mean density, `k_F`, `E_F`, `omega_p`, shell fillings, and related quantities.
- Support validation and interpretation of simulations.

### API sketch
```cpp
namespace inqkit::jellium {
    double wigner_seitz_radius(int n_electrons, double l_bohr);
    double mean_density(int n_electrons, double l_bohr);
    double fermi_wavevector(int n_electrons, double l_bohr);
    double fermi_energy(int n_electrons, double l_bohr);
    double plasmon_frequency(int n_electrons, double l_bohr);
}
```

---

## `detail/filesystem.hpp`
**Function:** Holds shared internal filesystem utilities.

### Main types
- namespace-only helpers

### Responsibilities
- Create directories safely.
- Check existence.
- Enforce overwrite policies.

### API sketch
```cpp
namespace inqkit::detail::filesystem {
    void ensure_directory(std::string const& path);
    bool exists(std::string const& path);
    void ensure_parent_directory(std::string const& path);
}
```

---

## `detail/validation.hpp`
**Function:** Holds shared internal validation routines.

### Responsibilities
- Validate path existence, grid sizes, field sizes, and selection consistency.
- Centralise precondition checks used by the public API.

### API sketch
```cpp
namespace inqkit::detail::validation {
    void require(bool condition, std::string const& message);
    void require_file_exists(std::string const& path);
    void require_same_grid_shape(/* a */, /* b */);
}
```

---

## `detail/text_io.hpp`
**Function:** Holds shared helpers for writing standard text formats.

### Responsibilities
- Write key-value files, metadata sidecars, and standard headers.
- Avoid repeated formatting code across writers.

### API sketch
```cpp
namespace inqkit::detail::text_io {
    void write_key_value(std::string const& path, std::string const& key, std::string const& value);
    void write_lines(std::string const& path, std::vector<std::string> const& lines);
}
```

---

## `detail/grid_layout.hpp`
**Function:** Holds shared helpers for field layout, indexing, and metadata conventions.

### Responsibilities
- Standardise axis ordering, spacing/origin metadata, and filename conventions.
- Ensure C++ writers and Python readers use the same assumptions.

### API sketch
```cpp
namespace inqkit::detail::grid_layout {
    std::string default_layout_name();
    std::string step_suffix(int step);
}
```

---

# 5. inqview (Python) structure

## 5.1 Directory tree

```text
inqview/
  inqview/
    data.py
    plots.py
    fourier.py
    paraview.py
    defaults.py
    config.py
```

---

## 5.2 File-by-file specification

## `data.py`
**Function:** Loads simulation data from the standard schema and validates that required files are present.

### Main classes
- `SimulationData`
- `FieldSeries`
- `DataError`

### Responsibilities
- Represent one simulation directory.
- Load observables, LEED patterns, text summaries, field metadata, and VTI series.
- Raise clear errors if files or fields are missing.

### API sketch
```python
from dataclasses import dataclass
from pathlib import Path

class DataError(RuntimeError):
    pass

@dataclass
class FieldSeries:
    root: Path
    files: list[Path]
    field_name: str

class SimulationData:
    def __init__(self, root: str | Path):
        ...

    def require(self, relative_path: str) -> Path:
        ...

    def load_observables(self, relative_path: str = "real_time/observables.csv"):
        ...

    def load_leed_pattern(self, relative_path: str = "real_time/leed/leed_pattern.txt"):
        ...

    def field_series(self, relative_dir: str) -> FieldSeries:
        ...

    def summary(self, relative_path: str = "simulation_summary.txt") -> dict:
        ...
```

---

## `plots.py`
**Function:** Provides reusable scientific plotting functions for observables, LEED patterns, and field-derived diagnostics.

### Main classes
- Mostly functions.
- Optional `PlotStyle` dataclass if needed.

### Responsibilities
- Plot time-series observables.
- Plot LEED patterns in linear and log scale.
- Plot orbital norms, overlaps, and other standard diagnostics.
- Use the defaults defined in `config.py`.

### API sketch
```python
from dataclasses import dataclass

@dataclass
class PlotStyle:
    dpi: int = 150
    figsize: tuple[float, float] = (6.0, 4.0)


def plot_observables(observables, columns, output_path, style: PlotStyle | None = None):
    ...


def plot_leed_pattern(leed, output_path, log_scale: bool = True):
    ...


def plot_leed_fft(leed_fft, output_path):
    ...


def plot_field_slice(field, output_path):
    ...
```

---

## `fourier.py`
**Function:** Provides FFT-based post-processing for energy, current, dipole, LEED, and any other sampled observable.

### Main classes
- `FourierResult`
- mostly functions

### Responsibilities
- Compute frequency-domain transforms from time-series data.
- Provide light preprocessing such as detrending and windowing.
- Return data objects that can be plotted with `plots.py`.

### API sketch
```python
from dataclasses import dataclass
import numpy as np

@dataclass
class FourierResult:
    frequency_axis: np.ndarray
    amplitude: np.ndarray
    units: str


def fft_observable(time_au, values, window=None, detrend=True) -> FourierResult:
    ...


def fft_energy(observables, column="E_total") -> FourierResult:
    ...


def fft_current(observables, column="Jx") -> FourierResult:
    ...


def fft_leed_pattern(leed_pattern):
    ...
```

---

## `paraview.py`
**Function:** Automates ParaView volume-rendering workflows and scripted animation export for VTI field series.

### Main classes
- `VolumeRenderSpec`
- `AnimationSpec`
- `ParaViewPipeline`

### Responsibilities
- Load a VTI field series.
- Apply a repeatable volume-render preset.
- Export PNG frames or an animation using a scripted ParaView pipeline.
- Hide ParaView string-valued details behind typed Python dataclasses.

### API sketch
```python
from dataclasses import dataclass
from pathlib import Path

@dataclass
class VolumeRenderSpec:
    array_name: str
    opacity_preset: str = "default_density"
    color_preset: str = "default_density"

@dataclass
class AnimationSpec:
    output_frames_dir: Path
    image_size: tuple[int, int] = (1600, 1200)
    fps: int = 12
    frame_stride: int = 1

class ParaViewPipeline:
    def __init__(self, data: "SimulationData"):
        ...

    def render_volume_series(
        self,
        series: "FieldSeries",
        render: VolumeRenderSpec,
        animation: AnimationSpec,
    ) -> list[Path]:
        ...

    def build_gif(self, frames_dir: str | Path, output_path: str | Path, fps: int = 12):
        ...
```

---

## `defaults.py`
**Function:** Provides standard high-level analysis and visualisation pipelines so most runs can be processed with a few function calls.

### Main functions
- `default_density_movie`
- `default_wavepacket_movie`
- `default_leed_plots`
- `default_observables_plots`
- `default_quantum_kick_analysis`

### Responsibilities
- Bundle together common workflows.
- Reduce boilerplate in per-run `analysis.py` files.
- Provide readable, conventional defaults while remaining overridable.

### API sketch
```python
def default_density_movie(data: "SimulationData"):
    ...


def default_wavepacket_movie(data: "SimulationData"):
    ...


def default_leed_plots(data: "SimulationData"):
    ...


def default_observables_plots(data: "SimulationData"):
    ...


def default_quantum_kick_analysis(data: "SimulationData"):
    ...
```

---

## `config.py`
**Function:** Stores shared plotting, colour, rendering, and style defaults used throughout inqview.

### Main classes
- `PlotDefaults`
- `RenderDefaults`

### Responsibilities
- Centralise colour maps, line styles, figure sizes, and render presets.
- Keep a consistent visual identity across all plots and movies.

### API sketch
```python
from dataclasses import dataclass

@dataclass
class PlotDefaults:
    dpi: int = 150
    figsize: tuple[float, float] = (6.0, 4.0)
    leed_cmap: str = "hot"
    observable_linewidth: float = 1.6

@dataclass
class RenderDefaults:
    image_size: tuple[int, int] = (1600, 1200)
    density_opacity_preset: str = "default_density"
    density_color_preset: str = "default_density"
```

---

# 6. How the APIs integrate together

The following examples are intentionally rough. They are architecture tests showing how the public API fits together.

---

# 7. Example 1: LEED simulation workflow

## 7.1 C++ simulation sketch

**Task requirements:**
- Ground state saved.
- Wavepacket injected into the last extra state.
- LEED screen at fixed `z`.
- Observables recorded as a function of time.
- Wavepacket density written every 10th timestep.
- Total density written every 10th timestep.
- HOMO and LUMO orbitals written every 20th timestep.
- LEED pattern written at end.
- VTI files written for volume rendering.

```cpp
#include <inq/inq.hpp>

#include <inqkit/config/simulation_config.hpp>
#include <inqkit/io/manifest_writer.hpp>
#include <inqkit/io/text_summary_writer.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/io/complex_field_3d_writer.hpp>
#include <inqkit/io/vti_image_data_writer.hpp>
#include <inqkit/ground_state/ground_state_tasks.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>
#include <inqkit/screens/plane_screen.hpp>
#include <inqkit/screens/leed_pattern_accumulator.hpp>

using namespace inq;
using namespace inq::magnitude;

int main() {
    // 1. Configuration
    inqkit::SimulationConfig cfg;
    cfg.system.lx_bohr = 60.0;
    cfg.system.ly_bohr = 60.0;
    cfg.system.lz_bohr = 120.0;
    cfg.system.cutoff_ha = 40.0;
    cfg.system.extra_states = 4;
    cfg.system.temperature_ev = 0.0862;

    cfg.gs.energy_tolerance_ha = 1e-8;
    cfg.gs.mixing_alpha = 0.1;
    cfg.gs.max_steps = 300;
    cfg.gs.use_broyden = true;

    cfg.rt.num_steps = 600;
    cfg.rt.dt_au = 0.02;
    cfg.rt.write_every_density = 10;
    cfg.rt.write_every_orbitals = 20;
    cfg.rt.observe_current = true;

    cfg.wp.enabled = true;
    cfg.wp.center_x_bohr = 30.0;
    cfg.wp.center_y_bohr = 30.0;
    cfg.wp.center_z_bohr = 100.0;
    cfg.wp.sigma_bohr = 1.0;
    cfg.wp.k0_x_bohr_inv = 0.0;
    cfg.wp.k0_y_bohr_inv = 0.0;
    cfg.wp.k0_z_bohr_inv = -3.83;
    cfg.wp.occupation = 1.0;

    cfg.output.root = "results";
    cfg.output.human_readable = true;
    cfg.output.performance_mode = false;
    cfg.output.write_vti = true;

    cfg.observables.step = true;
    cfg.observables.time_au = true;
    cfg.observables.energy_total = true;
    cfg.observables.energy_kinetic = true;
    cfg.observables.current_x = true;
    cfg.observables.current_y = true;
    cfg.observables.current_z = true;

    // 2. INQ system construction
    auto cell = systems::cell::orthorhombic(
        cfg.system.lx_bohr * 1.0_bohr,
        cfg.system.ly_bohr * 1.0_bohr,
        cfg.system.lz_bohr * 1.0_bohr
    ).finite();

    auto ions = systems::ions::parse("geometry/system.xyz", cell);

    systems::electrons electrons(
        ions,
        options::electrons{}
            .cutoff(cfg.system.cutoff_ha * 1.0_Ha)
            .extra_states(cfg.system.extra_states)
            .temperature(cfg.system.temperature_ev * 1.0_eV),
        input::kpoints::gamma()
    );

    // 3. Top-level run metadata
    inqkit::ManifestWriter manifest_writer(
        cfg.root_path() + "/manifest.txt",
        {.schema_version = "0.1", .code_name = "inqkit", .run_type = "leed"}
    );
    manifest_writer.write(cfg);

    inqkit::TextSummaryWriter summary_writer(cfg.root_path() + "/simulation_summary.txt");
    summary_writer.write_simulation_summary(cfg);

    // 4. Ground state
    ground_state::initial_guess(ions, electrons);
    auto gs = ground_state::calculate(
        ions,
        electrons,
        options::theory{}.pbe(),
        options::ground_state{}
            .energy_tolerance(cfg.gs.energy_tolerance_ha * 1.0_Ha)
            .mixing(cfg.gs.mixing_alpha)
            .max_steps(cfg.gs.max_steps)
            .broyden_mixing()
    );

    inqkit::GroundStateTasks gs_tasks(cfg);
    gs_tasks.add(inqkit::TextSummaryWriter(cfg.root_path() + "/ground_state/summary.txt"));
    gs_tasks.add(inqkit::RealField3DWriter(
        cfg.root_path() + "/ground_state/density",
        {.field_name = "total_density", .write_text = true, .write_raw = false, .include_meta = true}
    ));
    gs_tasks.run(gs, ions, electrons);

    electrons.save(cfg.root_path() + "/ground_state/checkpoint");

    // 5. Wavepacket
    auto wp = inqkit::WavePacket{}
        .center(cfg.wp.center_x_bohr, cfg.wp.center_y_bohr, cfg.wp.center_z_bohr)
        .sigma(cfg.wp.sigma_bohr)
        .k0(cfg.wp.k0_x_bohr_inv, cfg.wp.k0_y_bohr_inv, cfg.wp.k0_z_bohr_inv)
        .orthogonalise_against_occupied(electrons, 1e-8);

    auto wp_report = wp.inject_into_last_extra_state(electrons, cfg.wp.occupation);
    summary_writer.write_wavepacket_report(wp_report);

    // 6. LEED detector
    auto detector = inqkit::PlaneScreen::at_z(100.0);

    // 7. Real-time session
    inqkit::RealTimeSession rt(cfg, ions, electrons);

    rt.add(inqkit::ObservablesWriter(
        cfg.root_path() + "/real_time/observables.csv",
        cfg.observables
    ));

    rt.add(inqkit::RealField3DWriter(
        cfg.root_path() + "/real_time/density",
        {.field_name = "total_density", .write_text = true, .write_raw = false, .include_meta = true},
        {.stride = cfg.rt.write_every_density}
    ));

    rt.add(inqkit::RealField3DWriter(
        cfg.root_path() + "/real_time/wavepacket_density",
        {.field_name = "wavepacket_density", .write_text = true, .write_raw = false, .include_meta = true},
        {.stride = cfg.rt.write_every_density}
    ));

    rt.add(inqkit::ComplexField3DWriter(
        cfg.root_path() + "/real_time/orbitals/homo",
        {.field_name = "homo", .write_text = true, .write_raw = false, .include_meta = true},
        {.stride = cfg.rt.write_every_orbitals}
    ));

    rt.add(inqkit::ComplexField3DWriter(
        cfg.root_path() + "/real_time/orbitals/lumo",
        {.field_name = "lumo", .write_text = true, .write_raw = false, .include_meta = true},
        {.stride = cfg.rt.write_every_orbitals}
    ));

    rt.add(inqkit::VTIImageDataWriter(
        cfg.root_path() + "/visualisation/vti/total_density",
        {.array_name = "density"}
    ));

    rt.add(inqkit::VTIImageDataWriter(
        cfg.root_path() + "/visualisation/vti/wavepacket_density",
        {.array_name = "wavepacket_density"}
    ));

    rt.add(inqkit::LEEDPatternAccumulator(
        cfg.root_path() + "/real_time/leed/leed_pattern.txt",
        detector,
        {.start_time_au = 3.0, .end_time_au = 10.0}
    ));

    rt.start();

    real_time::propagate(
        ions,
        electrons,
        [&](auto const& data) {
            rt.step(data);
        },
        options::theory{}.pbe(),
        options::real_time{}
            .num_steps(cfg.rt.num_steps)
            .dt(cfg.rt.dt_au * 1.0_atomictime)
            .impulsive()
            .observables_current()
    );

    rt.finish();
}
```

---

## 7.2 Python analysis and visualisation sketch

**Task requirements:**
- Load the standardised outputs.
- Make a ParaView-driven density movie from the VTI series.
- Plot the LEED pattern and its FFT.
- Plot observables over time.
- Plot or otherwise inspect wavepacket/HOMO/LUMO outputs as needed.

```python
from pathlib import Path

from inqview.data import SimulationData
from inqview.paraview import ParaViewPipeline, VolumeRenderSpec, AnimationSpec
from inqview.plots import plot_observables, plot_leed_pattern, plot_leed_fft
from inqview.fourier import fft_leed_pattern

# 1. Load the run
sim = SimulationData("results")

# 2. Observables plot
obs = sim.load_observables("real_time/observables.csv")
plot_observables(
    obs,
    columns=["E_total", "E_kinetic", "Jx", "Jy", "Jz"],
    output_path="results/analysis/plots/observables.png"
)

# 3. LEED pattern plot
leed = sim.load_leed_pattern("real_time/leed/leed_pattern.txt")
plot_leed_pattern(
    leed,
    output_path="results/analysis/plots/leed_pattern.png",
    log_scale=True
)

leed_fft = fft_leed_pattern(leed)
plot_leed_fft(
    leed_fft,
    output_path="results/analysis/plots/leed_fft.png"
)

# 4. ParaView volume-render movie for total density
pv = ParaViewPipeline(sim)

total_density_series = sim.field_series("visualisation/vti/total_density")
wavepacket_series = sim.field_series("visualisation/vti/wavepacket_density")

pv.render_volume_series(
    total_density_series,
    render=VolumeRenderSpec(array_name="density"),
    animation=AnimationSpec(
        output_frames_dir=Path("results/visualisation/frames/total_density"),
        image_size=(1600, 1200),
        fps=12,
        frame_stride=1,
    )
)

pv.build_gif(
    frames_dir="results/visualisation/frames/total_density",
    output_path="results/visualisation/movies/total_density.gif",
    fps=12,
)

# 5. ParaView volume-render movie for the wavepacket density
pv.render_volume_series(
    wavepacket_series,
    render=VolumeRenderSpec(array_name="wavepacket_density"),
    animation=AnimationSpec(
        output_frames_dir=Path("results/visualisation/frames/wavepacket_density"),
        image_size=(1600, 1200),
        fps=12,
        frame_stride=1,
    )
)

pv.build_gif(
    frames_dir="results/visualisation/frames/wavepacket_density",
    output_path="results/visualisation/movies/wavepacket_density.gif",
    fps=12,
)
```

---

# 8. Example 2: Quantum kick workflow

## 8.1 C++ simulation sketch

**Task requirements:**
- Prepare a system in the ground state.
- Apply a quantum kick: ions begin moving impulsively at constant velocity at `t = 0`.
- Record observables such as current and energy versus time.
- Save outputs for later Fourier analysis.

```cpp
#include <inq/inq.hpp>

#include <inqkit/config/simulation_config.hpp>
#include <inqkit/io/manifest_writer.hpp>
#include <inqkit/io/text_summary_writer.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/ground_state/ground_state_tasks.hpp>
#include <inqkit/real_time/real_time_session.hpp>

using namespace inq;
using namespace inq::magnitude;

int main() {
    // 1. Configuration
    inqkit::SimulationConfig cfg;
    cfg.system.lx_bohr = 20.0;
    cfg.system.ly_bohr = 20.0;
    cfg.system.lz_bohr = 20.0;
    cfg.system.cutoff_ha = 54.0;
    cfg.system.extra_states = 0;
    cfg.system.temperature_ev = 0.0;

    cfg.gs.energy_tolerance_ha = 1e-8;

    cfg.rt.num_steps = 4000;
    cfg.rt.dt_au = 0.04;
    cfg.rt.write_every_density = 100;
    cfg.rt.observe_current = true;

    cfg.output.root = "results";

    cfg.observables.step = true;
    cfg.observables.time_au = true;
    cfg.observables.energy_total = true;
    cfg.observables.energy_kinetic = true;
    cfg.observables.current_x = true;
    cfg.observables.current_y = true;
    cfg.observables.current_z = true;

    // 2. Example system
    auto ions = systems::ions(
        systems::cell::cubic(cfg.system.lx_bohr * 1.0_bohr).periodic()
    );
    ions.insert_fractional("Li", {0.0, 0.0, 0.0});
    ions.insert_fractional("Li", {0.5, 0.5, 0.5});

    systems::electrons electrons(
        ions,
        options::electrons{}
            .cutoff(cfg.system.cutoff_ha * 1.0_Ha)
            .temperature(cfg.system.temperature_ev * 1.0_eV),
        input::kpoints::gamma()
    );

    // 3. Metadata
    inqkit::ManifestWriter manifest_writer(
        cfg.root_path() + "/manifest.txt",
        {.schema_version = "0.1", .code_name = "inqkit", .run_type = "quantum_kick"}
    );
    manifest_writer.write(cfg);

    inqkit::TextSummaryWriter summary_writer(cfg.root_path() + "/simulation_summary.txt");
    summary_writer.write_simulation_summary(cfg);

    // 4. Ground state
    ground_state::initial_guess(ions, electrons);
    auto gs = ground_state::calculate(
        ions,
        electrons,
        options::theory{}.pbe(),
        options::ground_state{}.energy_tolerance(cfg.gs.energy_tolerance_ha * 1.0_Ha)
    );

    // 5. Apply impulsive motion to ions at t=0
    //    Exact low-level implementation will depend on the chosen INQ setup for ionic velocities.
    //    This block is intentionally schematic.
    //    Example intent: all ions move at constant +x velocity after t=0.
    // ions.set_velocity_for_all({vx_bohr_per_au, 0.0, 0.0});

    // 6. Real-time session
    inqkit::RealTimeSession rt(cfg, ions, electrons);

    rt.add(inqkit::ObservablesWriter(
        cfg.root_path() + "/real_time/observables.csv",
        cfg.observables
    ));

    rt.add(inqkit::RealField3DWriter(
        cfg.root_path() + "/real_time/density",
        {.field_name = "total_density", .write_text = true, .write_raw = false, .include_meta = true},
        {.stride = cfg.rt.write_every_density}
    ));

    rt.start();

    real_time::propagate(
        ions,
        electrons,
        [&](auto const& data) {
            rt.step(data);
        },
        options::theory{}.pbe(),
        options::real_time{}
            .num_steps(cfg.rt.num_steps)
            .dt(cfg.rt.dt_au * 1.0_atomictime)
            .impulsive()
            .observables_current()
    );

    rt.finish();
}
```

---

## 8.2 Python post-processing sketch

**Task requirements:**
- Load the observables.
- Plot current and energy versus time.
- Fourier transform the energy response and current response.
- Save the frequency-domain plots.

```python
from inqview.data import SimulationData
from inqview.plots import plot_observables
from inqview.fourier import fft_energy, fft_current
from inqview.plots import plot_field_slice

sim = SimulationData("results")
obs = sim.load_observables("real_time/observables.csv")

# 1. Time-domain observables
plot_observables(
    obs,
    columns=["E_total", "E_kinetic"],
    output_path="results/analysis/plots/energy_vs_time.png"
)

plot_observables(
    obs,
    columns=["Jx", "Jy", "Jz"],
    output_path="results/analysis/plots/current_vs_time.png"
)

# 2. Frequency-domain analysis
energy_fft = fft_energy(obs, column="E_total")
current_fft = fft_current(obs, column="Jx")

# These would be plotted by dedicated plotting helpers in plots.py.
# Example names are shown here for clarity.
from inqview.plots import plot_observables as plot_fft_like  # placeholder

plot_fft_like(
    {"frequency": energy_fft.frequency_axis, "amplitude": energy_fft.amplitude},
    columns=["amplitude"],
    output_path="results/analysis/spectra/energy_fft.png"
)

plot_fft_like(
    {"frequency": current_fft.frequency_axis, "amplitude": current_fft.amplitude},
    columns=["amplitude"],
    output_path="results/analysis/spectra/current_fft.png"
)
```

---

# 9. Recommended first implementation order

## inqkit
1. `simulation_config.hpp`
2. `manifest_writer.hpp`
3. `text_summary_writer.hpp`
4. `observables_writer.hpp`
5. `real_field_3d_writer.hpp`
6. `complex_field_3d_writer.hpp`
7. `ground_state_tasks.hpp`
8. `wavepacket.hpp` + `injection_report.hpp`
9. `plane_screen.hpp`
10. `leed_pattern_accumulator.hpp`
11. `step_context.hpp`
12. `real_time_session.hpp`
13. `vti_image_data_writer.hpp`
14. `jellium/analytics.hpp`

## inqview
1. `data.py`
2. `config.py`
3. `plots.py`
4. `fourier.py`
5. `paraview.py`
6. `defaults.py`

---

# 10. Notes on scope
- The public API names in this document are intended to be memorable and stable.
- Writer classes live only in `io/`.
- Orbitals are treated as complex 3D fields, so no separate orbital writer class is required.
- `WavePacket` is one public class with builder-style definition, orthogonalisation methods, and injection methods.
- `SimulationData` is the central Python data handle.
- `fourier.py` is the dedicated home for FFT-based post-processing.
- ParaView automation is kept behind `ParaViewPipeline`, so ParaView-specific details do not leak into the user-facing API.

