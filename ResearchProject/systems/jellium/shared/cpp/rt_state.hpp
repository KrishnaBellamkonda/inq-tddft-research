// ============================================================================
// shared/cpp/rt_state.hpp  (jellium)
//
// Minimal reader/writer for results/rt_state.txt — the sidecar file that carries
// the per-run dynamical state which electrons.save() does NOT hold, so a run can
// be resumed and EXTENDED instead of recomputed
// (.claude/rules/final-timestep-checkpoint.md).
//
// electrons.save() stores the collective RT state (all orbitals, including an
// injected wavepacket). It does NOT store: which state index the WP occupies,
// the step/time reached, or a classical projectile's position and velocity.
// Those live here.
//
// Format: one `key = value` per line, so it stays greppable from the shell and
// parseable from post-processing without a schema.
// ============================================================================
#pragma once

#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

namespace jellium::rt_state {

struct State {
    int    last_step = 0;
    double time_au   = 0.0;
    double dt        = 0.0;
    int    wp_idx    = -1;      // wavepacket state index (-1 = no WP)
    // Classical projectile state (unused, left at 0, when there is no ion).
    double proj_z    = 0.0;
    double proj_vz   = 0.0;
    double proj_mass = 0.0;
    double proj_charge = 0.0;
};

inline void write(std::string const& path, State const& s) {
    std::ofstream f(path);
    f << std::setprecision(16);
    f << "last_step = "   << s.last_step   << "\n"
      << "time_au = "     << s.time_au     << "\n"
      << "dt = "          << s.dt          << "\n"
      << "wp_idx = "      << s.wp_idx      << "\n"
      << "proj_z = "      << s.proj_z      << "\n"
      << "proj_vz = "     << s.proj_vz     << "\n"
      << "proj_mass = "   << s.proj_mass   << "\n"
      << "proj_charge = " << s.proj_charge << "\n";
}

inline State read(std::string const& path) {
    State s;
    std::ifstream f(path);
    std::string line;
    while (std::getline(f, line)) {
        auto eq = line.find('=');
        if (eq == std::string::npos) continue;
        std::string key = line.substr(0, eq);
        std::string val = line.substr(eq + 1);
        // trim
        auto trim = [](std::string& t) {
            while (!t.empty() && (t.front() == ' ' || t.front() == '\t')) t.erase(t.begin());
            while (!t.empty() && (t.back()  == ' ' || t.back()  == '\t' ||
                                  t.back()  == '\r' || t.back() == '\n')) t.pop_back();
        };
        trim(key); trim(val);
        std::istringstream is(val);
        if      (key == "last_step")   is >> s.last_step;
        else if (key == "time_au")     is >> s.time_au;
        else if (key == "dt")          is >> s.dt;
        else if (key == "wp_idx")      is >> s.wp_idx;
        else if (key == "proj_z")      is >> s.proj_z;
        else if (key == "proj_vz")     is >> s.proj_vz;
        else if (key == "proj_mass")   is >> s.proj_mass;
        else if (key == "proj_charge") is >> s.proj_charge;
    }
    return s;
}

}  // namespace jellium::rt_state
