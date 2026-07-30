/*
 * OrbitalKineticStats — per-orbital BARE (extensive) kinetic energy + norm.
 *
 * Motivation (docs/notes/inq-energy-normalization-error.md)
 * ---------------------------------------------------------
 * INQ's reported `energy_kinetic` is INTENSIVE: energy.hpp:83 reduces the
 * per-orbital kinetic expectation values as Σ_i occ_i ⟨ψ_i|T|ψ_i⟩ / ⟨ψ_i|ψ_i⟩
 * (occ_sum with the norms argument, energy.hpp:55). For unitary dynamics the
 * division is a no-op, but under a norm-losing absorber (CAP) it pins the
 * reported kinetic at the surviving remnant's MEAN energy — the "energy shoots
 * up" artifact. The physically-extensive kinetic content is the BARE sum
 *
 *   E_kin_bare = Σ_i occ_i · T_i ,   T_i = ½ Σ_k |k|² |ψ̃_i(k)|²   [Ha]
 *
 * with no norm division anywhere. This observable computes, per KS orbital,
 *   norm_i = (dV/N_grid) Σ_k |ψ̃_i(k)|²        (= ∫|ψ_i|² dV; 1 at t=0)
 *   T_i    = ½ (dV/N_grid) Σ_k |k|² |ψ̃_i(k)|² (bare kinetic, occ NOT applied)
 * (dV/N_grid because INQ's to_fourier is an unnormalized DFT — Parseval then
 * reads Σ_k|ψ̃|² = N_grid Σ_r|ψ|² = (N_grid/dV)∫|ψ|²dV)
 * from ONE to_fourier of the whole orbital set per invocation, and logs
 *   kin_bare_total_ha    = Σ_i occ_i T_i             (the extensive kinetic)
 *   kin_normdiv_total_ha = Σ_i occ_i T_i / norm_i    (reconstruction of INQ's
 *                          reported energy_kinetic — per-step identity check;
 *                          division deliberately unguarded to match energy.hpp)
 *   norm_total           = Σ_i occ_i norm_i          (surviving electron count)
 * plus per-orbital norm_i / T_i columns and the wall-clock cost of the
 * evaluation (wall_ms) so the observable's overhead is measured in-run.
 *
 * Same Fourier-space math as ham.kinetic_expectation_value (=
 * operations::laplacian_expectation_value with factor −½) at gamma with zero
 * vector potential; validated by the t=0 identity kin_bare == kin_normdiv ==
 * energies.csv:kinetic. Gamma-only (single k-point), like WPMomentumStats.
 *
 * No engine edit: inq/ and inq-study/ untouched; the SCF Rayleigh-quotient
 * convention in occ_sum stays as upstream intends.
 *
 * Usage
 * -----
 *   OrbitalKineticStats obs("output/orbital_kinetic_stats.csv", {.write_every=1});
 *   // inside the real-time callback:
 *   obs.maybe_accumulate(data);    // no-op on skipped steps
 */
#pragma once

#include <inq/inq.hpp>
#include <operations/transform.hpp>
#include <basis/fourier_space.hpp>
#include <math/vector3.hpp>

#include <chrono>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace inqkit::observables {

struct OrbitalKineticStatsConfig {
    int write_every = 1;   // accumulate every Nth iteration; <=0 disables
};

// One evaluation — per-GLOBAL-orbital arrays plus the three totals.
struct OrbitalKineticSnapshot {
    std::vector<double> occ;        // occupations (fixed during RT)
    std::vector<double> norm;       // Parseval norm per orbital
    std::vector<double> tkin;       // bare ½Σk²|ψ̃|² per orbital (occ NOT applied)
    double kin_bare_total    = 0.0; // Σ occ_i·T_i        (extensive kinetic, Ha)
    double kin_normdiv_total = 0.0; // Σ occ_i·T_i/norm_i (== INQ reported kinetic)
    double norm_total        = 0.0; // Σ occ_i·norm_i     (electrons remaining)
};

class OrbitalKineticStats {
public:
    OrbitalKineticStats(std::string csv_path, OrbitalKineticStatsConfig cfg = {})
        : csv_path_(std::move(csv_path)), cfg_(cfg) {
        namespace fs = std::filesystem;
        if (auto parent = fs::path(csv_path_).parent_path(); !parent.empty())
            fs::create_directories(parent);
        file_.open(csv_path_);
        if (!file_)
            throw std::runtime_error(
                "OrbitalKineticStats: cannot open '" + csv_path_ + "'");
    }

    ~OrbitalKineticStats() {
        if (file_.is_open()) file_.close();
    }

    template <typename Viewables>
    void maybe_accumulate(Viewables const& data) {
        if (cfg_.write_every <= 0) return;
        if (data.iter() % cfg_.write_every != 0) return;
        accumulate(data);
    }

    // compute() split from accumulate() so it is testable without RT Viewables.
    OrbitalKineticSnapshot compute(inq::systems::electrons const& electrons) const {
        using namespace inq;

        if (electrons.kpin_size() != 1)
            throw std::runtime_error(
                "OrbitalKineticStats: only single-kpoint (gamma-only) "
                "runs are supported.");

        // INQ's to_fourier is an UNNORMALIZED DFT: Σ_k|ψ̃|² = N_grid·Σ_r|ψ|²
        // (verified numerically 2026-07-29: raw norm = N_grid/dV for a
        // unit-normalized orbital). Physical scale: ∫|ψ|²dV = (dV/N_grid)·Σ_k|ψ̃|².
        auto const& rbasis = electrons.kpin()[0].basis();
        const double scale = rbasis.volume_element() / rbasis.size();

        auto fphi = inq::operations::transform::to_fourier(electrons.kpin()[0]);
        auto const& fbasis = fphi.basis();

        const long n_global = fphi.set_part().size();
        const long st_start = fphi.set_part().start();
        const long st_size  = fphi.set_part().local_size();

        // [norm_0..norm_{n-1}, k2sum_0..k2sum_{n-1}] — zeros for non-local states,
        // filled by all_reduce (same pattern as WPMomentumStats).
        std::vector<double> buf(2 * n_global, 0.0);

        auto const sizes  = fbasis.local_sizes();
        auto fhc          = begin(fphi.hypercubic());
        auto point_op     = fbasis.point_op();

        for (long ist_l = 0; ist_l < st_size; ++ist_l) {
            const long ig = st_start + ist_l;

            buf[ig] = gpu::run(
                gpu::reduce(sizes[2]), gpu::reduce(sizes[1]), gpu::reduce(sizes[0]),
                0.0,
                [fhc, ist_l] GPU_LAMBDA (auto iz, auto iy, auto ix) {
                    auto v = fhc[ix][iy][iz][ist_l];
                    return inq::real(v)*inq::real(v) + inq::imag(v)*inq::imag(v);
                });

            buf[n_global + ig] = gpu::run(
                gpu::reduce(sizes[2]), gpu::reduce(sizes[1]), gpu::reduce(sizes[0]),
                0.0,
                [fhc, ist_l, point_op] GPU_LAMBDA (auto iz, auto iy, auto ix) {
                    auto v = fhc[ix][iy][iz][ist_l];
                    double w = inq::real(v)*inq::real(v) + inq::imag(v)*inq::imag(v);
                    auto k = point_op.gvector_cartesian(ix, iy, iz);
                    return (k[0]*k[0] + k[1]*k[1] + k[2]*k[2]) * w;
                });
        }

        if (fbasis.comm().size() > 1)
            fbasis.comm().all_reduce_in_place_n(buf.data(), buf.size(), std::plus<>{});
        if (fphi.set_comm().size() > 1)
            fphi.set_comm().all_reduce_in_place_n(buf.data(), buf.size(), std::plus<>{});

        // Occupations live on the state partition, replicated across the basis
        // (domain) ranks — gather over the set comm ONLY.
        std::vector<double> occ(n_global, 0.0);
        for (long ist_l = 0; ist_l < st_size; ++ist_l)
            occ[st_start + ist_l] = electrons.occupations()[0][ist_l];
        if (fphi.set_comm().size() > 1)
            fphi.set_comm().all_reduce_in_place_n(occ.data(), occ.size(), std::plus<>{});

        OrbitalKineticSnapshot s;
        s.occ  = std::move(occ);
        s.norm.resize(n_global);
        s.tkin.resize(n_global);
        for (long ig = 0; ig < n_global; ++ig) {
            s.norm[ig] = scale * buf[ig];
            s.tkin[ig] = 0.5 * scale * buf[n_global + ig];
            s.kin_bare_total    += s.occ[ig] * s.tkin[ig];
            s.kin_normdiv_total += s.occ[ig] * s.tkin[ig] / s.norm[ig];
            s.norm_total        += s.occ[ig] * s.norm[ig];
        }
        return s;
    }

    // accumulate() = timed compute() + one CSV row.
    template <typename Viewables>
    void accumulate(Viewables const& data) {
        const auto t0 = std::chrono::steady_clock::now();
        auto const s = compute(data.electrons());
        const double wall_ms =
            std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - t0).count();

        if (!header_written_) {
            file_ << "# n_states=" << s.norm.size()
                  << "  write_every=" << cfg_.write_every << '\n';
            file_ << "# occ =";
            for (auto o : s.occ) file_ << ' ' << o;
            file_ << '\n';
            file_ << "step,time_au,wall_ms,"
                     "kin_bare_total_ha,kin_normdiv_total_ha,norm_total";
            for (std::size_t i = 0; i < s.norm.size(); ++i)
                file_ << ",norm_" << i << ",tkin_" << i;
            file_ << '\n';
            header_written_ = true;
        }

        file_ << std::setprecision(12);
        file_ << data.iter() << ',' << data.time() << ',' << wall_ms << ','
              << s.kin_bare_total << ',' << s.kin_normdiv_total << ','
              << s.norm_total;
        for (std::size_t i = 0; i < s.norm.size(); ++i)
            file_ << ',' << s.norm[i] << ',' << s.tkin[i];
        file_ << '\n';
        file_.flush();
    }

private:
    std::string csv_path_;
    OrbitalKineticStatsConfig cfg_;
    std::ofstream file_;
    bool header_written_ = false;
};

} // namespace inqkit::observables
