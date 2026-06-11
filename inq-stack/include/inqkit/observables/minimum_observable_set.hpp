/*
 * inqkit::observables — the minimum observable set + manifest writer (ADR 0006).
 *
 * Single source of truth mapping a run-type to the primary observables it MUST
 * (required) and MAY (optional) produce. A run consults `minimum_set()` and
 * writes `results/observables_manifest.json` at startup via `write_manifest()`,
 * committing to the contract the inqview validator later checks (4 tiers).
 *
 * Pure: depends only on the standard library (string/vector/fstream). No INQ.
 * Tested by tests/cpp/test_minimum_observable_set.cpp.
 */
#pragma once

#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace inqkit::observables {

enum class RunType { coronene, jellium_wp, jellium_classical, free_wp };

inline std::string run_type_name(RunType t) {
  switch (t) {
    case RunType::coronene:          return "coronene";
    case RunType::jellium_wp:        return "jellium-wp";
    case RunType::jellium_classical: return "jellium-classical";
    case RunType::free_wp:           return "free-wp";
  }
  return "unknown";
}

// A declared physical invariant (tier 4). kind=="" means none. params are
// (name,value) pairs serialised into the manifest (e.g. drift_max value_mHa=1).
struct Invariant {
  std::string kind;
  std::string col;
  std::vector<std::pair<std::string, double>> params;
};

struct ObservableSpec {
  std::string name;
  bool required = true;
  std::string format = "csv";              // csv | vti | text
  std::string file;                        // for csv/text (relative to results/)
  std::string path;                        // for vti glob (relative to results/)
  std::string column;                      // optional CSV column
  std::string cadence = "write_every";     // step | write_every
  std::vector<std::string> schema;         // optional declared columns
  Invariant invariant;                     // optional
};

// ── helpers to build specs tersely ─────────────────────────────────────────
inline ObservableSpec csv(std::string name, std::string file, std::string column,
                          std::string cadence, Invariant inv = {}, bool required = true) {
  ObservableSpec s;
  s.name = std::move(name); s.format = "csv"; s.file = std::move(file);
  s.column = std::move(column); s.cadence = std::move(cadence);
  s.invariant = std::move(inv); s.required = required;
  return s;
}
inline ObservableSpec vti(std::string name, std::string path, bool required = true) {
  ObservableSpec s;
  s.name = std::move(name); s.format = "vti"; s.path = std::move(path);
  s.cadence = "write_every"; s.required = required;
  return s;
}
inline ObservableSpec text(std::string name, std::string file, bool required = true) {
  ObservableSpec s;
  s.name = std::move(name); s.format = "text"; s.file = std::move(file);
  s.cadence = "once"; s.required = required;
  return s;
}

// ── UNIVERSAL CORE (every run) ──────────────────────────────────────────────
inline std::vector<ObservableSpec> universal_core() {
  return {
    csv("energy_total", "raw/observables/observables.csv", "energy_total", "step",
        {"drift_max", "energy_total", {{"value_mHa", 1.0}}}),
    csv("energy_kinetic", "raw/observables/observables.csv", "energy_kinetic", "step"),
    csv("energy_hartree", "raw/observables/observables.csv", "energy_hartree", "step"),
    csv("energy_xc", "raw/observables/observables.csv", "energy_xc", "step"),
    csv("density_l2", "raw/observables/observables.csv", "density_l2", "write_every",
        {"zero_at_t0", "density_l2", {{"atol", 1e-9}}}),
    csv("gs_eigenvalues", "raw/observables/eigenvalues/eigenvalues.csv", "", "once"),
    csv("gs_occupations", "raw/observables/eigenvalues/occupations.csv", "", "once"),
    vti("gs_system_density", "raw/vti/density_gs_system/*.vti"),
    text("run_summary", "run_summary.txt"),
  };
}

// ── per run-type required∪optional set (core ∪ type) ────────────────────────
inline std::vector<ObservableSpec> minimum_set(RunType type) {
  auto set = universal_core();
  Invariant norm_band{"norm_band", "norm_check", {{"lo", 0.97}, {"hi", 1.03}}};

  auto add = [&set](std::vector<ObservableSpec> v) {
    for (auto &s : v) set.push_back(std::move(s));
  };

  if (type == RunType::jellium_wp || type == RunType::coronene ||
      type == RunType::free_wp) {
    add({
      // wp_momentum_stats.norm_check is a k-space integral (NOT normalised to 1),
      // so no norm-band invariant; existence/schema/finite suffice.
      csv("wp_momentum_stats", "raw/observables/wp_momentum_stats.csv", "", "write_every"),
      // wp_real_space_stats.norm_check is the real-space ∫|ψ|²dV ≈ 1.
      csv("wp_real_space_stats", "raw/observables/wp_real_space_stats.csv", "", "write_every", norm_band),
      vti("density_wp_rt", "raw/vti/density_wp/*.vti", /*required=*/false),
    });
  }
  if (type == RunType::jellium_wp) {
    add({
      csv("momentum_distribution", "raw/observables/momentum_distribution.csv", "", "write_every"),
      csv("state_energies", "raw/observables/state_energies.csv", "", "write_every"),
      csv("occupations_vs_time", "raw/observables/occupations_vs_time.csv", "", "write_every"),
      vti("density_system_rt", "raw/vti/density_system/*.vti"),
      vti("density_total_rt", "raw/vti/density_total/*.vti"),
    });
  }
  if (type == RunType::jellium_classical) {
    add({
      csv("electron_track", "raw/observables/electron_track.csv", "", "step"),
      csv("state_energies", "raw/observables/state_energies.csv", "", "write_every"),
      vti("density_system_rt", "raw/vti/density_system/*.vti"),
    });
  }
  if (type == RunType::coronene) {
    add({
      text("leed_screen_config", "raw/screens/screen_config.csv"),
      vti("density_total_rt", "raw/vti/density_rt_total/*.vti"),
    });
  }
  return set;
}

// ── JSON serialisation ──────────────────────────────────────────────────────
namespace detail {
inline void json_str(std::ostream &o, std::string const &s) {
  o << '"';
  for (char c : s) {
    if (c == '"' || c == '\\') o << '\\' << c;
    else o << c;
  }
  o << '"';
}
inline void json_invariant(std::ostream &o, Invariant const &inv) {
  o << "\"invariant\":{\"kind\":"; json_str(o, inv.kind);
  if (!inv.col.empty()) { o << ",\"col\":"; json_str(o, inv.col); }
  for (auto const &p : inv.params) { o << ",\""; o << p.first << "\":" << p.second; }
  o << "}";
}
inline void json_observable(std::ostream &o, ObservableSpec const &s) {
  o << "{\"name\":"; json_str(o, s.name);
  o << ",\"required\":" << (s.required ? "true" : "false");
  o << ",\"format\":"; json_str(o, s.format);
  if (!s.file.empty()) { o << ",\"file\":"; json_str(o, s.file); }
  if (!s.path.empty()) { o << ",\"path\":"; json_str(o, s.path); }
  if (!s.column.empty()) { o << ",\"column\":"; json_str(o, s.column); }
  o << ",\"cadence\":"; json_str(o, s.cadence);
  if (!s.schema.empty()) {
    o << ",\"schema\":[";
    for (std::size_t i = 0; i < s.schema.size(); ++i) {
      if (i) o << ',';
      json_str(o, s.schema[i]);
    }
    o << "]";
  }
  if (!s.invariant.kind.empty()) { o << ','; json_invariant(o, s.invariant); }
  o << "}";
}
}  // namespace detail

// Serialise the manifest to a JSON string.
inline std::string manifest_json(RunType type, int write_every, int n_steps) {
  std::ostringstream o;
  o << "{\"run_type\":"; detail::json_str(o, run_type_name(type));
  o << ",\"schema_version\":1";
  o << ",\"write_every\":" << write_every;
  o << ",\"n_steps\":" << n_steps;
  o << ",\"observables\":[";
  auto set = minimum_set(type);
  for (std::size_t i = 0; i < set.size(); ++i) {
    if (i) o << ',';
    detail::json_observable(o, set[i]);
  }
  o << "]}";
  return o.str();
}

// Write results/observables_manifest.json. Returns true on success.
inline bool write_manifest(std::string const &path, RunType type,
                           int write_every, int n_steps) {
  std::ofstream f(path);
  if (!f) return false;
  f << manifest_json(type, write_every, n_steps) << "\n";
  return static_cast<bool>(f);
}

}  // namespace inqkit::observables
