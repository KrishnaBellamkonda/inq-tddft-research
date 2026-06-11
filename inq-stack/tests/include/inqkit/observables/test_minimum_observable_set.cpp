// Pure-tier test of inqkit::observables::MinimumObservableSet + manifest writer
// (ADR 0006). No INQ — pure std-library data + JSON serialisation.

#include <catch2/catch_test_macros.hpp>

#include <inqkit/observables/minimum_observable_set.hpp>

#include <algorithm>
#include <string>

using namespace inqkit::observables;

static bool has(std::vector<ObservableSpec> const &s, std::string const &name) {
  return std::any_of(s.begin(), s.end(),
                     [&](ObservableSpec const &o) { return o.name == name; });
}
static ObservableSpec const &get(std::vector<ObservableSpec> const &s,
                                 std::string const &name) {
  return *std::find_if(s.begin(), s.end(),
                       [&](ObservableSpec const &o) { return o.name == name; });
}

TEST_CASE("MinimumObservableSet: universal core in every run-type", "[obs_set][pure]") {
  for (auto t : {RunType::coronene, RunType::jellium_wp,
                 RunType::jellium_classical, RunType::free_wp}) {
    auto s = minimum_set(t);
    CHECK(has(s, "energy_total"));
    CHECK(has(s, "energy_hartree"));   // promoted to core (was default-off)
    CHECK(has(s, "energy_xc"));
    CHECK(has(s, "density_l2"));
    CHECK(has(s, "gs_occupations"));
    CHECK(has(s, "run_summary"));
  }
}

TEST_CASE("MinimumObservableSet: per-type required observables", "[obs_set][pure]") {
  auto wp = minimum_set(RunType::jellium_wp);
  CHECK(has(wp, "wp_momentum_stats"));
  CHECK(has(wp, "momentum_distribution"));
  CHECK(get(wp, "energy_total").required);

  auto cl = minimum_set(RunType::jellium_classical);
  CHECK(has(cl, "electron_track"));
  CHECK_FALSE(has(cl, "wp_momentum_stats"));   // classical has no WP

  auto co = minimum_set(RunType::coronene);
  CHECK(has(co, "leed_screen_config"));
  CHECK(has(co, "wp_momentum_stats"));         // coronene has a WP
}

TEST_CASE("MinimumObservableSet: energy_total carries the drift invariant", "[obs_set][pure]") {
  auto wp = minimum_set(RunType::jellium_wp);   // bind: get() returns a ref into it
  auto const &e = get(wp, "energy_total");
  CHECK(e.invariant.kind == "drift_max");
  REQUIRE(e.invariant.params.size() == 1);
  CHECK(e.invariant.params[0].first == "value_mHa");
  CHECK(e.invariant.params[0].second == 1.0);
}

TEST_CASE("MinimumObservableSet: manifest JSON is well-formed + declares the set",
          "[obs_set][pure]") {
  std::string j = manifest_json(RunType::jellium_wp, 2, 190);
  CHECK(j.find("\"run_type\":\"jellium-wp\"") != std::string::npos);
  CHECK(j.find("\"write_every\":2") != std::string::npos);
  CHECK(j.find("\"n_steps\":190") != std::string::npos);
  CHECK(j.find("\"name\":\"energy_total\"") != std::string::npos);
  CHECK(j.find("\"name\":\"wp_momentum_stats\"") != std::string::npos);
  CHECK(j.find("\"kind\":\"drift_max\"") != std::string::npos);
  CHECK(j.find("\"kind\":\"norm_band\"") != std::string::npos);
  CHECK(j.find("\"value_mHa\":1") != std::string::npos);
  // balanced braces / brackets => well-formed structure
  int br = 0, sq = 0;
  for (char c : j) {
    if (c == '{') br++; else if (c == '}') br--;
    else if (c == '[') sq++; else if (c == ']') sq--;
  }
  CHECK(br == 0);
  CHECK(sq == 0);
}
