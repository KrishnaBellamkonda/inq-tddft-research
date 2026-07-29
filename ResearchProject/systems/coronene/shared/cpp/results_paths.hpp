
// ============================================================================
// shared/cpp/results_paths.hpp
//
// Compile-time results-tree paths per docs/results_folder_structure_spec.md.
// Every coronene run.cpp uses these helpers so the output layout is identical
// across runs, and so the post-processor in scripts/coronene_postprocess.py
// can rely on a fixed structure.
//
// All paths are relative to the run's working directory (which is the run_*/
// folder itself when launched by inq-run). Helpers create directories on first
// access so writers can write straight away.
// ============================================================================
#pragma once

#include <filesystem>
#include <string>

namespace coronene::results {

inline std::string ensure_dir(std::string const &p) {
    std::filesystem::create_directories(p);
    return p;
}

// Helper for *_path() returns: makes sure the file's parent directory exists
// before the run.cpp opens an ofstream. Without this, writing to a file in a
// subtree whose dir was not yet created (e.g. raw/wavepacket/...) silently
// drops the data.
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
// Three categories: total, system, wp. See docs/results_folder_structure_spec.md §7.
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

// ---- raw/vti/ ------------------------------------------------------------
inline std::string vti_density_gs_system_dir()   { return ensure_dir("results/raw/vti/density_gs_system"); }
inline std::string vti_density_gs_orbitals_dir() { return ensure_dir("results/raw/vti/density_gs_orbitals"); }
inline std::string vti_density_total_dir()       { return ensure_dir("results/raw/vti/density_rt_total"); }
inline std::string vti_density_system_dir()      { return ensure_dir("results/raw/vti/density_rt_system"); }
inline std::string vti_density_wp_dir()          { return ensure_dir("results/raw/vti/density_rt_wp"); }
inline std::string vti_orbitals_dir()            { return ensure_dir("results/raw/vti/orbitals"); }

}  // namespace coronene::results
