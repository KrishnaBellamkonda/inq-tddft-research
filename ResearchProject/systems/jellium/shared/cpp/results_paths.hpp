// ============================================================================
// shared/cpp/results_paths.hpp  (jellium)
//
// Compile-time results-tree paths. Direct port of
// ResearchProject/systems/coronene/shared/cpp/results_paths.hpp into namespace
// `jellium::results`. The on-disk schema is identical to coronene's so a
// single post-processor structure can serve both systems.
//
// All paths are relative to the run's working directory (the run_*/ folder
// when launched by inq-run). Helpers create directories on first access so
// writers can write straight away.
// ============================================================================
#pragma once

#include <filesystem>
#include <string>

namespace jellium::results {

inline std::string ensure_dir(std::string const &p) {
    std::filesystem::create_directories(p);
    return p;
}

inline std::string ensure_parent(std::string const &file_path) {
    auto parent = std::filesystem::path(file_path).parent_path();
    if (!parent.empty()) std::filesystem::create_directories(parent);
    return file_path;
}

// ---- Top-level -----------------------------------------------------------
inline std::string root()                  { return ensure_dir("results"); }
inline std::string raw()                   { return ensure_dir("results/raw"); }
inline std::string analysis()              { return ensure_dir("results/analysis"); }
inline std::string run_summary_path()      { return ensure_parent("results/run_summary.txt"); }

// ---- raw/ground_state/ ---------------------------------------------------
inline std::string gs_dir()                       { return ensure_dir("results/raw/ground_state"); }
inline std::string gs_density_system_dir()        { return ensure_dir("results/raw/ground_state/density_system"); }
inline std::string gs_density_orbitals_dir()      { return ensure_dir("results/raw/ground_state/density_gs_orbitals"); }
inline std::string gs_summary_path()              { return ensure_parent("results/raw/ground_state/summary.txt"); }
inline std::string gs_eigenvalues_path()          { return ensure_parent("results/raw/ground_state/eigenvalues.csv"); }
inline std::string gs_occupations_path()          { return ensure_parent("results/raw/ground_state/occupations.csv"); }

// ---- raw/wavepacket/ -----------------------------------------------------
inline std::string wp_dir()                       { return ensure_dir("results/raw/wavepacket"); }
inline std::string wp_density_initial_dir()       { return ensure_dir("results/raw/wavepacket/density_wp_initial"); }
inline std::string wp_wavefunction_initial_dir()  { return ensure_dir("results/raw/wavepacket/wavefunction_wp_initial"); }
inline std::string wp_config_path()               { return ensure_parent("results/raw/wavepacket/wavepacket_config.txt"); }
inline std::string wp_injection_report_path()     { return ensure_parent("results/raw/wavepacket/injection_report.txt"); }
inline std::string wp_orthogonality_report_path() { return ensure_parent("results/raw/wavepacket/orthogonality_report.csv"); }

// ---- raw/density/ --------------------------------------------------------
// Three categories: total, system (jellium background), wp.
inline std::string density_total_dir()  { return ensure_dir("results/raw/density/density_rt_total"); }
inline std::string density_system_dir() { return ensure_dir("results/raw/density/density_rt_system"); }
inline std::string density_wp_dir()     { return ensure_dir("results/raw/density/density_rt_wp"); }

// ---- raw/screens/ --------------------------------------------------------
inline std::string screens_total_dir()         { return ensure_dir("results/raw/screens/total"); }
inline std::string screens_instantaneous_dir() { return ensure_dir("results/raw/screens/instantaneous"); }
inline std::string screens_time_windowed_dir() { return ensure_dir("results/raw/screens/time_windowed"); }
inline std::string screens_config_path()       { return ensure_parent("results/raw/screens/screen_config.csv"); }
inline std::string screens_window_ranges_path(){ return ensure_parent("results/raw/screens/window_ranges.csv"); }

// ---- raw/overlap/ --------------------------------------------------------
inline std::string overlap_dir() { return ensure_dir("results/raw/overlap"); }

// ---- raw/observables/ ----------------------------------------------------
inline std::string observables_dir()      { return ensure_dir("results/raw/observables"); }
inline std::string observables_csv_path() { return ensure_parent("results/raw/observables/observables.csv"); }
inline std::string state_energies_csv_path() { return ensure_parent("results/raw/observables/state_energies.csv"); }
inline std::string occupations_csv_path()    { return ensure_parent("results/raw/observables/occupations_vs_time.csv"); }
inline std::string momentum_distribution_csv_path() { return ensure_parent("results/raw/observables/momentum_distribution.csv"); }

// ---- raw/vti (delta) -----------------------------------------------------
// Co-located with the other density VTI series so the existing density
// postprocessor can pick them up via the same _CATEGORIES sweep.
inline std::string vti_density_delta_dir()        { return ensure_dir("results/raw/vti/density_rt_delta"); }
inline std::string vti_density_delta_coarse_dir() { return ensure_dir("results/raw/vti/density_rt_delta_coarse"); }

// ---- raw/vti/ ------------------------------------------------------------
inline std::string vti_density_gs_system_dir()   { return ensure_dir("results/raw/vti/density_gs_system"); }
inline std::string vti_density_gs_orbitals_dir() { return ensure_dir("results/raw/vti/density_gs_orbitals"); }
inline std::string vti_density_total_dir()       { return ensure_dir("results/raw/vti/density_rt_total"); }
inline std::string vti_density_system_dir()      { return ensure_dir("results/raw/vti/density_rt_system"); }
inline std::string vti_density_wp_dir()          { return ensure_dir("results/raw/vti/density_rt_wp"); }
inline std::string vti_orbitals_dir()            { return ensure_dir("results/raw/vti/orbitals"); }

}  // namespace jellium::results
