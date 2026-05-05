// ============================================================================
// inqkit/observables/occupations_writer.hpp
//
// Per-step dump of the KS-orbital occupations f_i to a long-format CSV.
//
// In TDDFT propagation INQ holds the occupations f_i FROZEN — they're
// the GS values, set once at electrons.load() and never updated.  This
// observable therefore records what is, in principle, a flat trace.
// Two reasons to record it anyway:
//
//   1. Audit: a non-flat trace would be a numerics red flag (e.g.
//      something inadvertently rewriting electrons.occupations()).
//      Cheap to verify.
//   2. Companion to state_energies.csv: knowing f_i at every recorded
//      time-step lets the energy-balance postprocess compute the
//      occupation-weighted bath energy change ΣΔE_i · f_i without
//      having to join with a separate ground-state CSV.
//
// Cost: trivial — n_states scalar writes per snapshot.
// ============================================================================
#pragma once

#include <inq/inq.hpp>

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <utility>

namespace inqkit::observables {

class OccupationsWriter {
public:
    OccupationsWriter(std::string path)
        : path_(std::move(path)) {
        namespace fs = std::filesystem;
        if (auto parent = fs::path(path_).parent_path(); !parent.empty())
            fs::create_directories(parent);
        file_.open(path_);
        if (!file_)
            throw std::runtime_error(
                "OccupationsWriter: cannot open '" + path_ + "'");
        file_ << "step,time_au,kpoint_index,state_index,occupation\n";
        file_.flush();
    }

    template <typename Viewables>
    void snapshot(Viewables const& data) {
        using namespace inq;
        auto const& electrons = data.electrons();
        const int step    = data.iter();
        const double time = data.time();

        const int nk_local = electrons.kpin_size();
        const long kp_offset = electrons.kpin_part().start();

        file_ << std::setprecision(15);
        for (int ik = 0; ik < nk_local; ++ik) {
            auto const& phi = electrons.kpin()[ik];
            auto const& occs = electrons.occupations()[ik];
            const long state_start = phi.set_part().start();
            const auto nst = static_cast<long>(occs.size());
            for (long ist = 0; ist < nst; ++ist) {
                const double f_i = static_cast<double>(occs[ist]);
                const long state_index = state_start + ist;
                file_ << step << ',' << time << ','
                      << (kp_offset + ik) << ',' << state_index << ','
                      << f_i << '\n';
            }
        }
        file_.flush();
    }

    ~OccupationsWriter() {
        if (file_.is_open()) file_.close();
    }

private:
    std::string path_;
    std::ofstream file_;
};

} // namespace inqkit::observables
