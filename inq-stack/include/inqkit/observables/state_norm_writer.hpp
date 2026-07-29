// inqkit::observables::StateNormWriter — per-orbital norm diagnostic.
//
// Writes ∫|ψ_i(r,t)|² dV for EVERY evolved orbital (occupied + extra, incl. the
// WP) at a chosen cadence. All KS orbitals are normalised, so each norm should
// stay ≈ 1; a drift beyond ~1e-3 flags numerical leakage during propagation
// (the diagnostic's purpose — todo_later.md "norm-conservation diagnostic").
//
// CSV: step,time_au,state_index,norm
//
// Mirrors the compute()/accumulate() split of the WP stats observables so the
// per-state norms are unit-testable directly (no CSV, no RT Viewables).
// Gamma-only (single k-point); single-rank set/domain partition (MPI reduction
// is a TODO, matching the sibling observables).

#ifndef INQKIT_OBSERVABLES_STATE_NORM_WRITER_HPP
#define INQKIT_OBSERVABLES_STATE_NORM_WRITER_HPP

#include <inq/inq.hpp>

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace inqkit {
namespace observables {

struct StateNorm {
    int    state_index = -1;
    double norm        = 0.0;   // ∫|ψ_i|² dV  (≈ 1 for a normalised orbital)
};

struct StateNormWriterConfig {
    int write_every = 1;        // accumulate every Nth iteration; <=0 disables
};

class StateNormWriter {
public:
    StateNormWriter(std::string csv_path, StateNormWriterConfig cfg = {})
        : csv_path_(std::move(csv_path)), cfg_(cfg) {
        namespace fs = std::filesystem;
        if (auto parent = fs::path(csv_path_).parent_path(); !parent.empty())
            fs::create_directories(parent);
        file_.open(csv_path_);
        if (!file_)
            throw std::runtime_error("StateNormWriter: cannot open '" + csv_path_ + "'");
        file_ << "# write_every=" << cfg_.write_every << '\n';
        file_ << "step,time_au,state_index,norm\n";
    }

    ~StateNormWriter() { if (file_.is_open()) file_.close(); }

    // Per-state norms from the current electrons state. Split out so it is
    // directly unit-testable.
    std::vector<StateNorm> compute(inq::systems::electrons const& electrons) const {
        using namespace inq;

        if (electrons.kpin_size() != 1)
            throw std::runtime_error(
                "StateNormWriter: only single-kpoint (gamma-only) runs supported.");

        auto const& phi   = electrons.kpin()[0];
        auto const& basis = phi.basis();
        const double dV   = basis.volume_element();
        const long st_start = phi.set_part().start();
        const long st_size  = phi.set_part().local_size();
        auto const sizes    = basis.local_sizes();
        auto phic           = begin(phi.hypercubic());

        std::vector<StateNorm> out;
        out.reserve(static_cast<std::size_t>(st_size));
        for (int il = 0; il < static_cast<int>(st_size); ++il) {
            const double norm = gpu::run(
                gpu::reduce(sizes[2]), gpu::reduce(sizes[1]), gpu::reduce(sizes[0]),
                0.0,
                [phic, il, dV] GPU_LAMBDA (auto iz, auto iy, auto ix) {
                    auto v = phic[ix][iy][iz][il];
                    return dV * (inq::real(v) * inq::real(v)
                               + inq::imag(v) * inq::imag(v));
                });
            out.push_back({static_cast<int>(st_start + il), norm});
        }
        return out;
    }

    template <typename Viewables>
    void maybe_accumulate(Viewables const& data) {
        if (cfg_.write_every <= 0) return;
        if (data.iter() % cfg_.write_every != 0) return;
        accumulate(data);
    }

    template <typename Viewables>
    void accumulate(Viewables const& data) {
        auto const norms = compute(data.electrons());
        const int    step = data.iter();
        const double t_au = data.time();
        for (auto const& s : norms)
            file_ << step << ',' << t_au << ',' << s.state_index << ','
                  << s.norm << '\n';
        file_.flush();
    }

private:
    std::string csv_path_;
    StateNormWriterConfig cfg_;
    std::ofstream file_;
};

}  // namespace observables
}  // namespace inqkit

#endif  // INQKIT_OBSERVABLES_STATE_NORM_WRITER_HPP
