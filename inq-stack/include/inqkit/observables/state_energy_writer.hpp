// ============================================================================
// inqkit/observables/state_energy_writer.hpp
//
// Writes per-state energy expectation values <phi_i|H|phi_i>(t) and the
// energy variance <phi_i|H^2|phi_i> - <phi_i|H|phi_i>^2 over the course
// of a real-time propagation, in long format CSV
//
//     step, time_au, kpoint_index, state_index, weight, occupation,
//     E_expect_ha, E_variance_ha2
//
// Implementation adapted from the professor-supplied reference at
//     ResearchProject/literature/misc/viewables.hpp
// (methods state_energy_expectations() and state_energy_variance()),
// reorganised so all caching/IO lives outside the propagator's
// `viewables` class. The Hamiltonian is accessed via the local extension
// `viewables::ham()` we added to inq/src/real_time/viewables.hpp.
//
// Cost: 2 ham() applications + 2 overlap_diagonal calls per snapshot
// (one for <H>, one extra for <H^2>). ~O(n_states) per step.
// ============================================================================
#pragma once

#include <inq/inq.hpp>
#include <operations/overlap_diagonal.hpp>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace inqkit::observables {

class StateEnergyWriter {
public:
    StateEnergyWriter(std::string path, bool emit_variance = true)
        : path_(std::move(path)), emit_variance_(emit_variance) {
        namespace fs = std::filesystem;
        if (auto parent = fs::path(path_).parent_path(); !parent.empty())
            fs::create_directories(parent);
        file_.open(path_);
        if (!file_)
            throw std::runtime_error(
                "StateEnergyWriter: cannot open '" + path_ + "'");
        file_ << "step,time_au,kpoint_index,state_index,weight,occupation,"
                 "E_expect_ha";
        if (emit_variance_) file_ << ",E_variance_ha2";
        file_ << '\n';
    }

    // Snapshot is templated so it can take any concrete viewables<...> type
    // produced by INQ's propagator.
    template <typename Viewables>
    void snapshot(Viewables const& data) {
        using namespace inq;

        auto const& electrons = data.electrons();
        auto const& ham       = data.ham();

        const int step    = data.iter();
        const double time = data.time();

        // Local-rank kpin loop. Each rank writes its own slice; if needed
        // a downstream MPI gather can be added (jellium runs gamma-only,
        // single rank, so this is fine for the current target.)
        const int nk_local = electrons.kpin_size();
        const long kp_offset = electrons.kpin_part().start();

        file_ << std::setprecision(15);
        for (int ik = 0; ik < nk_local; ++ik) {
            auto const& phi = electrons.kpin()[ik];

            auto hphi    = ham(phi);
            auto h_diag  = operations::overlap_diagonal(phi, hphi);

            // Variance requires <H^2> = <phi|H H phi>
            std::vector<double> var_local;
            if (emit_variance_) {
                auto hhphi   = ham(hphi);
                auto h2_diag = operations::overlap_diagonal(phi, hhphi);
                var_local.resize(h2_diag.size());
                for (std::size_t i = 0; i < h2_diag.size(); ++i) {
                    const double e1 = real(h_diag[i]);
                    const double e2 = real(h2_diag[i]);
                    var_local[i] = e2 - e1 * e1;
                }
            }

            const auto nstates = h_diag.size();
            const double weight = electrons.kpin_weights()[ik];
            auto const& occs    = electrons.occupations()[ik];

            const long state_start = phi.set_part().start();

            for (std::size_t ist = 0; ist < nstates; ++ist) {
                const double E_exp = real(h_diag[ist]);
                const long state_index = state_start + ist;
                const double occ = (ist < occs.size())
                    ? static_cast<double>(occs[ist]) : 0.0;

                file_ << step << ',' << time << ','
                      << (kp_offset + ik) << ',' << state_index << ','
                      << weight << ',' << occ << ',' << E_exp;
                if (emit_variance_) {
                    const double v = std::max(0.0, var_local[ist]);
                    file_ << ',' << v;
                }
                file_ << '\n';
            }
        }
        file_.flush();
    }

    ~StateEnergyWriter() {
        if (file_.is_open()) file_.close();
    }

private:
    std::string path_;
    bool emit_variance_;
    std::ofstream file_;
};

} // namespace inqkit::observables
